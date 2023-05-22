Diagnostic_Plot <- function(fit_num, fit_Model_list, IFRvals,
                            Return_likelihoods_only = F,
                            Return_data_only_plots = F,
                            get_data = F
){
  gc();print(fit_num)
  fit_Model <- fit_Model_list[[fit_num]]
  if(is.null(fit_Model)){return(NA)}

  ################################################
  ################################################
  ### First, load inputs
  pars_obs <- fit_Model$pmcmc_results$inputs$pars_obs
  pcr_det <- pars_obs$pcr_det
  index <- squire:::odin_index(fit_Model$model)
  combined_data <- cbind(expand.grid(Week_gr = 1:16, Age_gr = 1:17, Replicate = 1:100), pars_obs$combined_data)
  dfj_mcmc_data <- pars_obs$drj_mcmc
  frac_reg <- pars_obs$frac_reg

  ################################################
  ################################################
  ## Get Susceptibles, Infections, Cumulative infections, pcr_positive and pcr_perc
  Mod_Age_Deaths_Lus_pcr <- fit_Model$output[,c(index$D,index$S),] %>%
    melt(varnames = c("date","var","Replicate"), value.name = "value") %>%
    mutate(Age_gr = as.numeric(gsub(".*?([0-9]+).*", '\\1', var)),
           var = substr(var, 1, 1),
           date = as.Date(date)) %>%
    dcast(... ~ var, value.var="value") %>%
    dplyr::rename(Mod_cd_Lus = "D", Sus = "S") %>%
    group_by(Age_gr, Replicate) %>%
    tidyr::replace_na(list(Mod_cd_Lus = 0)) %>%
    mutate(Sus = as.integer(ifelse(is.na(Sus), max(Sus, na.rm=T), Sus))) %>%
    mutate(cum_infs = max(Sus, na.rm = T)-Sus,
           infs = c(0, diff(max(Sus, na.rm = T)-Sus)),
           pcr_pos = cma:::roll_func_10(infs, pcr_det),
           pcr_perc = pcr_pos/max(Sus,na.rm = T))

  Mod_Age_Deaths_Lus <- Mod_Age_Deaths_Lus_pcr %>%
    filter(date > "2020-06-10") %>%
    mutate(Week_gr = as.numeric(cut.Date(date, breaks = "1 week"))) %>%
    group_by(Age_gr, Week_gr, Replicate) %>%
    summarise(Mod_cd_Lus = max(Mod_cd_Lus),
              pcr_perc = pcr_perc[4], # Thursday - middle of the study week
              infs = sum(infs),
              date = date[1]) %>%
    ungroup() %>% group_by(Age_gr, Replicate) %>%
    mutate(Mod_cd_Lus = c(0,diff(Mod_cd_Lus))) %>% #,
    filter(Week_gr != 1) %>%
    mutate(Week_gr = Week_gr-1) %>%
    filter(Week_gr <= 16) %>%
    merge(combined_data) %>%
    arrange(as.numeric(Replicate), Age_gr, Week_gr) %>%
    mutate(Mod_cd_morgue = Mod_cd_Lus*frac_reg)


  ################################################
  ################################################
  ### Calculate Likelihoods
  if(Return_likelihoods_only){
    ## Integrate over mcmc samples

    LL_Integrated <- lapply(unique(Mod_Age_Deaths_Lus$Replicate), function(x){

      tmp_data <- filter(Mod_Age_Deaths_Lus, Replicate == x)

      if(pars_obs$lld =="no scaling no weeks 4-5"){
        # Likelihood function for sensitivity analysis D (no weekly scaling)
        dpois_vec <- colSums(dpois(x = tmp_data$BurRegs,
                                   lambda = tmp_data$Mod_cd_morgue + pars_obs$drj_mcmc$drj_mcmc_data_baseline,
                                   log = T))
        dbinom_vec <- colSums(dbinom(x = tmp_data$PosTests,
                                     size = tmp_data$Samples,
                                     prob = (tmp_data$Mod_cd_morgue + pars_obs$drj_mcmc$drj_mcmc_data_baseline*tmp_data$pcr_perc)/
                                       (tmp_data$Mod_cd_morgue + pars_obs$drj_mcmc$drj_mcmc_data_baseline), log = T))
      } else {
        # General likelihood function
        dpois_vec <- colSums(dpois(x = tmp_data$BurRegs,
                                   lambda = tmp_data$Mod_cd_morgue * pars_obs$drj_mcmc$drj_mcmc_data_baseline_agstd + pars_obs$drj_mcmc$drj_mcmc_data_baseline,
                                   log = T))
        dbinom_vec <- colSums(dbinom(x = tmp_data$PosTests,
                                     size = tmp_data$Samples,
                                     prob = (tmp_data$Mod_cd_morgue*pars_obs$drj_mcmc$drj_mcmc_data_baseline_agstd + pars_obs$drj_mcmc$drj_mcmc_data_baseline*tmp_data$pcr_perc)/
                                       (tmp_data$Mod_cd_morgue*pars_obs$drj_mcmc$drj_mcmc_data_baseline_agstd + pars_obs$drj_mcmc$drj_mcmc_data_baseline), log = T))
      }

      LL_pois_binom <- dpois_vec + dbinom_vec
      LL_pois_rep <- log(Brobdingnag::sum(exp(Brobdingnag::as.brob(dpois_vec)))/length(dpois_vec))
      LL_binom_rep <- log(Brobdingnag::sum(exp(Brobdingnag::as.brob(dbinom_vec)))/length(dbinom_vec))
      LL_pois_binom_rep <- log(Brobdingnag::sum(exp(Brobdingnag::as.brob(LL_pois_binom)))/length(LL_pois_binom))

      return(data.frame(Replicate = x, LL_pois_binom_rep = LL_pois_binom_rep, LL_pois_rep = LL_pois_rep, LL_binom_rep = LL_binom_rep))

    })

    LL_Integrated <- data.table::rbindlist(LL_Integrated)


    ################################################
    ################################################
    ## PCR prevalence  and seroprevalence likelihoods
    pcr_vals <- lapply(1:dim(fit_Model$output)[3], function(x){
      pcr_tmp <- ll_prev_func(df = pars_obs$pcr_df, det = pars_obs$pcr_det, Incidence = "Infections", out = fit_Model$output[,,x], index = index)
    })

    sero_vals <- lapply(1:dim(fit_Model$output)[3], function(x){
      sero_tmp <- ll_prev_func(df = pars_obs$sero_df, det = pars_obs$sero_det, Incidence = "Symptoms", out = fit_Model$output[,,x], index = index, model_params = fit_Model$parameters)
    })

    sero_pcr_ll <- data.frame(Replicate = 1:100,
                              ll_pcr = unlist(pcr_vals),
                              ll_sero = unlist(sero_vals))

    ################################################
    ################################################
    ## Total likelihood
    ll_all <- merge(ll_ov_pois_bin, sero_pcr_ll) %>% ungroup() %>%
      mutate(ll_total = ll_pois_mean_overall + ll_bin_mean_overall + ll_pcr + ll_sero)

    LL_overall <- data.frame(Pois = log(Brobdingnag::sum(exp(Brobdingnag::as.brob(ll_all$LL_pois_rep)))/nrow(ll_all)),
                             Binom = log(Brobdingnag::sum(exp(Brobdingnag::as.brob(ll_all$LL_binom_rep)))/nrow(ll_all)),
                             Pois_Binom = log(Brobdingnag::sum(exp(Brobdingnag::as.brob(ll_all$LL_pois_binom_rep)))/nrow(ll_all)),
                             PCR = log(Brobdingnag::sum(exp(Brobdingnag::as.brob(ll_all$ll_pcr)))/nrow(ll_all)),
                             Sero = log(Brobdingnag::sum(exp(Brobdingnag::as.brob(ll_all$ll_sero)))/nrow(ll_all)),
                             Total = log(Brobdingnag::sum(exp(Brobdingnag::as.brob(ll_all$ll_total)))/nrow(ll_all)))

    lls <- list(
      ll_all = ll_all,
      ll_summary = LL_overall)

    return(lls)}


  ################################################
  ################################################
  ## Wrangle dataframes for plotting
  ################################################
  ################################################

  ## Get data by week
  Mod_Age_Deaths_Lus_Av_Week <- cbind(expand.grid(Week_gr = 1:16, Age_gr = 1:17), pars_obs$combined_data) %>% group_by(Week_gr) %>%
    summarise(across(c(Samples, PosTests, BurRegs), sum)) %>%
    merge(Mod_Age_Deaths_Lus %>% select(Week_gr, date) %>% unique())

  ## Get data by age
  Mod_Age_Deaths_Lus_Av_Age <- cbind(expand.grid(Week_gr = 1:16, Age_gr = 1:17), pars_obs$combined_data) %>% group_by(Age_gr) %>%
    summarise(across(c(Samples, PosTests, BurRegs), sum)) %>%
    mutate(Age_gr_label = ifelse(Age_gr %in% 1:16, paste0(Age_gr*5-5, "-",Age_gr*5-1),"80+"),
           Age_gr_label <- factor(Age_gr_label, levels = Age_gr_label))

  ## Get modelled PCR and seroprevalence
  pcr_sero_data <- fit_Model$output[,c(index$S, index$E2),] %>%
    melt(varnames = c("date","var","Replicate"), value.name = "value") %>%
    mutate(Age_gr = as.numeric(gsub(".*?([0-9]+).", '\\1', var)),
           var = substr(var, 1, 1),
           date = as.Date(date)) %>%
    dcast(... ~ var, value.var="value") %>%
    replace_na(list(E = 0)) %>%
    group_by(Age_gr) %>% mutate(S = as.integer(ifelse(is.na(S), max(S, na.rm = T),S))) %>% ungroup() %>%
    dplyr::rename(Sus = "S", Exp = "E") %>%
    group_by(date, Replicate) %>% summarise(Exp = sum(Exp), Sus = sum(Sus)) %>% ungroup() %>%
    group_by(Replicate) %>%
    mutate(cum_infs = max(Sus, na.rm = T)-Sus,
           attack_rate = cum_infs/sum(max(Sus)),

           infs = c(0, diff(max(Sus, na.rm = T)-Sus)),
           pcr_pos = cma:::roll_func_10(infs, pcr_det),
           pcr_perc = pcr_pos/max(Sus,na.rm = T),

           Symps = as.integer(Exp) * fit_Model$parameters$gamma_E,
           sero_pos = cma:::roll_func_10(Symps, pars_obs$sero_det),
           sero_perc = sero_pos/max(Sus,na.rm = T)) %>%
    group_by(date) %>%
    summarise_at(c("attack_rate", "pcr_perc","sero_perc"), list(median = median, ci = bayestestR::ci), na.rm = TRUE)


  #################################
  #################################
  ## Make data only plots first  ##
  #################################
  #################################

  # Figure 5A: weekly burial registrations
  Poisson_Figure_weeks <- ggplot(Mod_Age_Deaths_Lus_Av_Week, aes(x = date)) +
    geom_point(aes(y = BurRegs, color = "Burial registrations")) +
    xlab("Date") + ylab("Deaths") +
    theme_minimal(base_size = 7, base_family = "Helvetica") +
    theme(legend.key = element_rect(fill = "white", linetype = 0),
          plot.title = element_text(size = 7, face = "bold"))

  # Figure 5B: weekly burial registrations
  Poisson_Figure_age <- ggplot(data = Mod_Age_Deaths_Lus_Av_Age, aes(x = Age_gr_label)) +
    geom_point(aes(y = BurRegs, color = "Burial registrations")) +
    xlab("Age") + ylab("Deaths") +
    theme_minimal(base_size = 7, base_family = "Helvetica") +
    theme(legend.key = element_rect(fill = "white", linetype = 0),
          axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
          plot.title = element_text(size = 7, face = "bold"))

  # Figure 5C: PCR prevalence and seroprevalence
  PCR_sero_prev_plot <- ggplot(pcr_sero_data, aes(x = date, y = pcr_perc_median)) +
    geom_point(aes(x= as.Date(pars_obs$pcr_df$date_start) + 0.5*(as.Date(pars_obs$pcr_df$date_end)-as.Date(pars_obs$pcr_df$date_start)),
                   y=pars_obs$pcr_df$pos_tests/pars_obs$pcr_df$samples, color = "PCR %"),  size = 2) +
    geom_errorbar(aes(ymin=Hmisc::binconf(pars_obs$pcr_df$pos_tests,pars_obs$pcr_df$samples)[,"Lower"],
                      ymax=Hmisc::binconf(pars_obs$pcr_df$pos_tests,pars_obs$pcr_df$samples)[,"Upper"],
                      x=as.Date(pars_obs$pcr_df$date_start) + 0.5*(as.Date(pars_obs$pcr_df$date_end)-as.Date(pars_obs$pcr_df$date_start)),
                      width=10,color = "PCR %")) +
    geom_errorbarh(aes(xmin=as.Date(pars_obs$pcr_df$date_start),xmax=as.Date(pars_obs$pcr_df$date_end),y=pars_obs$pcr_df$pos_tests/pars_obs$pcr_df$samples, height=0,color = "PCR %")) +

    geom_point(aes(x= as.Date(pars_obs$sero_df$date_start) + 0.5*(as.Date(pars_obs$sero_df$date_end)-as.Date(pars_obs$sero_df$date_start)),
                   y=pars_obs$sero_df$pos_tests/pars_obs$sero_df$samples, color = "Sero %"), size = 2) +
    geom_errorbar(aes(ymin=Hmisc::binconf(pars_obs$sero_df$pos_tests,pars_obs$sero_df$samples)[,"Lower"],
                      ymax=Hmisc::binconf(pars_obs$sero_df$pos_tests,pars_obs$sero_df$samples)[,"Upper"],
                      x=as.Date(pars_obs$sero_df$date_start) + 0.5*(as.Date(pars_obs$sero_df$date_end)-as.Date(pars_obs$sero_df$date_start)),
                      width=10,color = "Sero %")) +
    geom_errorbarh(aes(xmin=as.Date(pars_obs$pcr_df$date_start),xmax=as.Date(pars_obs$pcr_df$date_end),y=pars_obs$sero_df$pos_tests/pars_obs$sero_df$samples, height=0,color = "Sero %")) +
    ylab(paste0("Population prevalence")) +
    coord_cartesian(xlim = c(as.Date("2020-04-15"), as.Date("2020-10-01")),
                    ylim = c(0, 0.28)) +
    xlab("Date")+
    scale_color_manual(name=NULL,
                       breaks = c("Modelled attack rate", "Modelled PCR %", "Modelled Sero %","PCR %","Sero %"),
                       values = c("black","darkgoldenrod2","chartreuse4","darkgoldenrod2","chartreuse4")) +
    scale_y_continuous(labels = scales::percent) +
    theme_minimal(base_size = 7, base_family = "Helvetica") +
    theme(plot.title = element_text(size = 7, face = "bold"))



  # Figure 5D: weekly post-mortem prevalence
  Week_prev_plot <- ggplot(data = Mod_Age_Deaths_Lus_Av_Week, aes(x = date)) +
    # Prevalence of positive in sample
    geom_point(aes(y = PosTests/Samples, color = "Post-mortem prevalence")) +
    geom_errorbar(aes(ymin = Hmisc::binconf(PosTests,Samples)[,"Lower"],
                      ymax = Hmisc::binconf(PosTests,Samples)[,"Upper"])) +
    coord_cartesian(ylim = c(0,1)) +
    xlab("Date") +
    ylab("Post-mortem prevalence") +
    scale_y_continuous(labels = scales::percent) +
    theme_minimal(base_size = 7, base_family = "Helvetica") +
    theme(plot.title = element_text(size = 7, face = "bold"))


  ## Plot 5E: Age group post-mortem prevalence fits
  Age_prev_plot <- ggplot(data = Mod_Age_Deaths_Lus_Av_Age, aes(x = Age_gr_label)) +
    # Prevalence of positive in sample
    geom_point(aes(y = PosTests/Samples, color = "Post-mortem prevalence")) +
    geom_errorbar(aes(ymin = Hmisc::binconf(PosTests,Samples)[,"Lower"],
                      ymax = Hmisc::binconf(PosTests,Samples)[,"Upper"])) +
    xlab("Age") +
    scale_x_discrete(limits = c(paste0(1:16*5-5, "-",1:16*5-1),"80+")) +
    scale_y_continuous(labels = scales::percent) +
    ylab("Post-mortem prevalence") +
    theme_minimal(base_size = 7, base_family = "Helvetica") +
    theme(plot.title = element_text(size = 7, face = "bold"),
          axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))


  if(Return_data_only_plots){

    return(list(Poisson_Figure_weeks = Poisson_Figure_weeks + scale_color_manual(name=NULL, breaks = c("Burial registrations"),values = c("black"),guide = guide_legend(override.aes = list(linetype = c(0), shape = c(16), alpha = c(1)))),
                Poisson_Figure_age = Poisson_Figure_age + scale_color_manual(name=NULL, breaks = c("Burial registrations"),values = c("black"),guide = guide_legend(override.aes = list(linetype = c(0), shape = c(16), alpha = c(1)))),
                Week_prev_plot = Week_prev_plot + scale_color_manual(name = NULL, breaks = c("Total +ve deaths"),values = c("black")),
                Age_prev_plot = Age_prev_plot + scale_color_manual(name = NULL, breaks = c("Total +ve deaths"),values = c("black")),
                PCR_sero_prev_plot = PCR_sero_prev_plot))
  }


  ####################################
  ####################################
  # Get summary statistics for plots #
  ####################################
  ####################################

  ## Get integrated summary statistics for each week
  Plotting_data_integrated_summarised_Weeks <- lapply(unique(Mod_Age_Deaths_Lus$Week_gr), function(x){

    # For each week group:
    # Get index vectors of week
    Week_vector <- Mod_Age_Deaths_Lus$Week_gr==x
    Week_vector_272 <- filter(Mod_Age_Deaths_Lus, Replicate ==1)$Week_gr==x

    # Filter data to week group and reshape
    Baseline <- as.data.frame(pars_obs$drj_mcmc$drj_mcmc_data_baseline)[Week_vector_272,]
    colnames(Baseline) <- 1:4000
    Baseline$Age_gr <- 1:17
    Baseline <- pivot_longer(Baseline, cols = 1:4000, names_to = "Sample", values_to = "Baseline")

    Age_standardisation <- as.data.frame(pars_obs$drj_mcmc$drj_mcmc_data_baseline_agstd)[Week_vector_272,]
    colnames(Age_standardisation) <- 1:4000
    Age_standardisation$Age_gr <- 1:17
    Age_standardisation <- pivot_longer(Age_standardisation, cols = 1:4000, names_to = "Sample", values_to = "Age_std")

    # Merge together
    Baselines <- merge(Baseline, Age_standardisation)
    Merged_Baselines <- merge(Mod_Age_Deaths_Lus[Week_vector,], Baselines)

    # Modelled burial registrations
    Merged_Baselines$Mod_tot_ds_morgue <- Merged_Baselines$Baseline + # Baseline
      Merged_Baselines$Age_std * Merged_Baselines$Mod_cd_morgue # Scaled COVID

    # Modelled coincidental PCR positive morgue deaths
    Merged_Baselines$Mod_coin_cpos_ds_morgue <- Merged_Baselines$Baseline * Merged_Baselines$pcr_perc # Baseline multiplied by PCR prevalence

    # Modelled total PCR morgue deaths
    Merged_Baselines$Mod_pos_ds_morgue <- Merged_Baselines$Mod_coin_cpos_ds_morgue + # Coincidental covid deaths that got registered
      Merged_Baselines$Age_std * Merged_Baselines$Mod_cd_morgue # Actual scaled COVID

    # Sum for each samples and replicate
    Merged_Baselines <- Merged_Baselines %>% ungroup() %>% group_by(Replicate, Sample) %>%
      summarise(across(c(Baseline, Mod_tot_ds_morgue, Mod_coin_cpos_ds_morgue, Mod_pos_ds_morgue), sum)) %>%
      ungroup()

    # Generate morgue prevalences from deaths
    Merged_Baselines$Pos_prev <- Merged_Baselines$Mod_pos_ds_morgue/Merged_Baselines$Mod_tot_ds_morgue
    Merged_Baselines$Pos_prev_coin <- Merged_Baselines$Mod_coin_cpos_ds_morgue/Merged_Baselines$Mod_tot_ds_morgue
    Merged_Baselines$Pos_prev_causal <- Merged_Baselines$Pos_prev - Merged_Baselines$Pos_prev_coin


    # Get summary statistics
    Median <- Merged_Baselines %>% summarise(across(c(Baseline, Mod_tot_ds_morgue, Mod_coin_cpos_ds_morgue, Mod_pos_ds_morgue, Pos_prev, Pos_prev_causal, Pos_prev_coin), median))
    Low_CI <- Merged_Baselines %>% summarise(across(c(Baseline, Mod_tot_ds_morgue, Mod_coin_cpos_ds_morgue, Mod_pos_ds_morgue, Pos_prev, Pos_prev_causal, Pos_prev_coin), ~bayestestR::ci(.)$CI_low))
    High_CI<- Merged_Baselines %>% summarise(across(c(Baseline, Mod_tot_ds_morgue, Mod_coin_cpos_ds_morgue, Mod_pos_ds_morgue, Pos_prev, Pos_prev_causal, Pos_prev_coin), ~bayestestR::ci(.)$CI_high))

    # Return results
    return(data.frame(Week_gr = x, Measure = c("Baseline", "Mod_tot_ds_morgue", "Mod_coin_cpos_ds_morgue", "Mod_pos_ds_morgue", "Pos_prev", "Pos_prev_causal", "Pos_prev_coin"),
                      Median = t(Median), Low_CI = t(Low_CI), High_CI = t(High_CI)))

  })

  ## Get integrated summary statistics for each age group
  Plotting_data_integrated_summarised_Ages <- lapply(unique(Mod_Age_Deaths_Lus$Age_gr), function(x){

    # For each age group:
    # Get index vectors of age
    Age_vector <- Mod_Age_Deaths_Lus$Age_gr==x
    Age_vector_272 <- filter(Mod_Age_Deaths_Lus, Replicate ==1)$Age_gr==x

    # Filter data to age group and reshape
    Baseline <- as.data.frame(pars_obs$drj_mcmc$drj_mcmc_data_baseline)[Age_vector_272,]
    colnames(Baseline) <- 1:4000
    Baseline$Week_gr <- 1:16
    Baseline <- pivot_longer(Baseline, cols = 1:4000, names_to = "Sample", values_to = "Baseline")

    Age_standardisation <- as.data.frame(pars_obs$drj_mcmc$drj_mcmc_data_baseline_agstd)[Age_vector_272,]
    colnames(Age_standardisation) <- 1:4000
    Age_standardisation$Week_gr <- 1:16
    Age_standardisation <- pivot_longer(Age_standardisation, cols = 1:4000, names_to = "Sample", values_to = "Age_std")

    # Merge together
    Baselines <- merge(Baseline, Age_standardisation)
    Merged_Baselines <- merge(Mod_Age_Deaths_Lus[Age_vector,], Baselines)

    # Modelled burial registrations
    Merged_Baselines$Mod_tot_ds_morgue <- Merged_Baselines$Baseline + # Baseline
      Merged_Baselines$Age_std * Merged_Baselines$Mod_cd_morgue # Scaled COVID

    # Modelled coincidental PCR positive morgue deaths
    Merged_Baselines$Mod_coin_cpos_ds_morgue <- Merged_Baselines$Baseline * Merged_Baselines$pcr_perc # Baseline multiplied by PCR prevalence

    # Modelled total PCR morgue deaths
    Merged_Baselines$Mod_pos_ds_morgue <- Merged_Baselines$Mod_coin_cpos_ds_morgue + # Coincidental covid deaths that got registered
      Merged_Baselines$Age_std * Merged_Baselines$Mod_cd_morgue # Actual scaled COVID

    # Sum for each samples and replicate
    Merged_Baselines <- Merged_Baselines %>% ungroup() %>% group_by(Replicate, Sample) %>%
      summarise(across(c(Baseline, Mod_tot_ds_morgue, Mod_coin_cpos_ds_morgue, Mod_pos_ds_morgue), sum)) %>%
      ungroup()

    # Generate morgue prevalences from deaths
    Merged_Baselines$Pos_prev <- Merged_Baselines$Mod_pos_ds_morgue/Merged_Baselines$Mod_tot_ds_morgue
    Merged_Baselines$Pos_prev_coin <- Merged_Baselines$Mod_coin_cpos_ds_morgue/Merged_Baselines$Mod_tot_ds_morgue
    Merged_Baselines$Pos_prev_causal <- Merged_Baselines$Pos_prev - Merged_Baselines$Pos_prev_coin

    # Get summary statistics
    Median <- Merged_Baselines %>% summarise(across(c(Baseline, Mod_tot_ds_morgue, Mod_coin_cpos_ds_morgue, Mod_pos_ds_morgue, Pos_prev, Pos_prev_causal, Pos_prev_coin), median))
    Low_CI <- Merged_Baselines %>% summarise(across(c(Baseline, Mod_tot_ds_morgue, Mod_coin_cpos_ds_morgue, Mod_pos_ds_morgue, Pos_prev, Pos_prev_causal, Pos_prev_coin), ~bayestestR::ci(.)$CI_low))
    High_CI<- Merged_Baselines %>% summarise(across(c(Baseline, Mod_tot_ds_morgue, Mod_coin_cpos_ds_morgue, Mod_pos_ds_morgue, Pos_prev, Pos_prev_causal, Pos_prev_coin), ~bayestestR::ci(.)$CI_high))

    # Return results
    return(data.frame(Age_gr = x, Measure = c("Baseline", "Mod_tot_ds_morgue", "Mod_coin_cpos_ds_morgue", "Mod_pos_ds_morgue", "Pos_prev", "Pos_prev_causal", "Pos_prev_coin"),
                      Median = t(Median), Low_CI = t(Low_CI), High_CI = t(High_CI)))

  })


  Plotting_data_integrated_summarised_Weeks <- data.table::rbindlist(Plotting_data_integrated_summarised_Weeks) %>% merge(Mod_Age_Deaths_Lus %>% select(Week_gr, date) %>% unique())
  Plotting_data_integrated_summarised_Ages <- data.table::rbindlist(Plotting_data_integrated_summarised_Ages) %>% mutate(Age_gr_label = ifelse(Age_gr %in% 1:16, paste0(Age_gr*5-5, "-",Age_gr*5-1),"80+"))


  if(get_data){
    return(list(fit_Model = fit_Model,
                pcr_sero_data = pcr_sero_data,
                Mod_Age_Deaths_Lus_Av_Age = Mod_Age_Deaths_Lus_Av_Age,
                Mod_Age_Deaths_Lus_Av_Week = Mod_Age_Deaths_Lus_Av_Week,
                Plotting_data_integrated_summarised_Weeks = Plotting_data_integrated_summarised_Weeks,
                Plotting_data_integrated_summarised_Ages = Plotting_data_integrated_summarised_Ages

    ))
  }

  #################################
  #################################
  ######### Add Model fit #########
  #################################
  #################################

  ## Plot 5A: Burial registrations by week
  Poisson_Figure_weeks_2 <- Poisson_Figure_weeks +
    geom_line(data = Plotting_data_integrated_summarised_Weeks %>% filter(Measure == "Mod_tot_ds_morgue"), aes(y=Median, color = "Model fit")) +
    geom_ribbon(data = Plotting_data_integrated_summarised_Weeks %>% filter(Measure == "Mod_tot_ds_morgue"), aes(ymin=Low_CI, ymax=High_CI), alpha=0.3) +
    scale_color_manual(name=NULL,
                       breaks = c("Burial registrations", "Model fit"),
                       values = c("black","black"),
                       guide = guide_legend(override.aes = list(linetype = c(0,1),
                                                                shape = c(16,NA),
                                                                alpha = c(1,1)),
                                            nrow = 1))

  ## Plot 5B: Burial registrations by age
  Poisson_Figure_age_2 <- Poisson_Figure_age +
    geom_line(data = Plotting_data_integrated_summarised_Ages %>% filter(Measure == "Mod_tot_ds_morgue"), aes(y=Median, color = "Model fit", group = 1)) +
    geom_ribbon(data = Plotting_data_integrated_summarised_Ages %>% filter(Measure == "Mod_tot_ds_morgue"), aes(ymin=Low_CI, ymax=High_CI, group = 1), alpha=0.3) +
    scale_color_manual(name=NULL,
                       breaks = c("Burial registrations", "Model fit"),
                       values = c("black","black"),
                       guide = guide_legend(override.aes = list(linetype = c(NA,1),
                                                                shape = c(16,NA),
                                                                alpha = c(1,1)),
                                            nrow = 1))

  ## Plot 5C: Population PCR prevalence and seroprevalence
  PCR_sero_prev_plot_2 <- PCR_sero_prev_plot +
    geom_line(aes(x=date, y=pcr_perc_median, color = "Modelled PCR %"), linetype=2) +
    geom_ribbon(aes(x=date,ymin=pcr_perc_ci$CI_low, ymax=pcr_perc_ci$CI_high), alpha=0.3, fill = "darkgoldenrod1")+
    geom_line(aes(x=date, y=sero_perc_median, color = "Modelled Sero %"), linetype=2) +
    geom_ribbon(aes(x=date,ymin=sero_perc_ci$CI_low, ymax=sero_perc_ci$CI_high), alpha=0.3, fill = "chartreuse4") +
    geom_line(aes(x = date, y = attack_rate_median, color = "Modelled attack rate"), inherit.aes = F) +
    geom_ribbon(aes(x = date, ymin = attack_rate_ci$CI_low, ymax = attack_rate_ci$CI_high), inherit.aes = F, alpha = 0.3) +
    guides(colour = guide_legend(override.aes = list(shape = c(NA,NA,NA,19,19),
                                                     linetype = c(1,2,2,NA,NA))))

  ## Plot 5D: Post-mortem PCR prevalence by week
  Week_prev_plot_2 <- Week_prev_plot +
    geom_line(data = Plotting_data_integrated_summarised_Weeks %>% filter(Measure == "Pos_prev_coin"), aes(y = Median, col = "Coincidental deaths"), linetype="dashed") +
    geom_ribbon(data = Plotting_data_integrated_summarised_Weeks %>% filter(Measure == "Pos_prev_coin"), aes(ymin = Low_CI, ymax= High_CI), alpha=0.3, fill = "darkblue") +
    geom_line(data = Plotting_data_integrated_summarised_Weeks %>% filter(Measure == "Pos_prev_causal"), aes(y = Median, col = "Causal deaths"),linetype="dashed") +
    geom_ribbon(data = Plotting_data_integrated_summarised_Weeks %>% filter(Measure == "Pos_prev_causal"), aes(ymin=Low_CI, ymax= High_CI), alpha=0.3, fill = "darkred") +
    geom_line(data = Plotting_data_integrated_summarised_Weeks %>% filter(Measure == "Pos_prev"), aes(y = Median, col = "Total +ve deaths"),linetype="dashed") +
    geom_ribbon(data = Plotting_data_integrated_summarised_Weeks %>% filter(Measure == "Pos_prev"), aes(ymin=Low_CI, ymax= High_CI), alpha=0.3) +
    scale_color_manual(name = NULL, breaks = c("Total +ve deaths","Causal deaths","Coincidental deaths","Post-mortem prevalence"),
                       values = c("black","darkred","darkblue","black"),
                       guide = guide_legend(override.aes = list(linetype = c(2,2,2,NA),
                                                                shape = c(NA,NA,NA,19)),
                                            nrow = 4, byrow =T))

  ## Plot 5E: Post-mortem PCR prevalence by age
  Age_prev_plot_2 <- Age_prev_plot +
    geom_line(data = Plotting_data_integrated_summarised_Ages %>% filter(Measure == "Pos_prev_coin"), aes(y = Median, col = "Coincidental deaths", group = 1), linetype="dashed") +
    geom_ribbon(data = Plotting_data_integrated_summarised_Ages %>% filter(Measure == "Pos_prev_coin"), aes(ymin = Low_CI, ymax= High_CI, group = 1), alpha=0.3, fill = "darkblue") +
    geom_line(data = Plotting_data_integrated_summarised_Ages %>% filter(Measure == "Pos_prev_causal"), aes(y = Median, col = "Causal deaths", group = 1),linetype="dashed") +
    geom_ribbon(data = Plotting_data_integrated_summarised_Ages %>% filter(Measure == "Pos_prev_causal"), aes(ymin=Low_CI, ymax= High_CI, group = 1), alpha=0.3, fill = "darkred") +
    geom_line(data = Plotting_data_integrated_summarised_Ages %>% filter(Measure == "Pos_prev"), aes(y = Median, col = "Total +ve deaths", group = 1),linetype="dashed") +
    geom_ribbon(data = Plotting_data_integrated_summarised_Ages %>% filter(Measure == "Pos_prev"), aes(ymin=Low_CI, ymax= High_CI, group = 1), alpha=0.3) +
    scale_color_manual(name = NULL, breaks = c("Total +ve deaths","Causal deaths","Coincidental deaths","Post-mortem prevalence"),
                       values = c("black","darkred","darkblue","black"),
                       guide = guide_legend(override.aes = list(linetype = c(2,2,2,NA),
                                                                shape = c(NA,NA,NA,19)),
                                            nrow = 4, byrow =T))

  ## Plot 5F: Rt trends
  p_Rt_Reff_a <- rt_plot_immunity(fit_Model)$plot +
    scale_fill_manual("",values = c("Reff" = "#48996b","Rt" ="#3f8da7"), labels = c(expression("R"[eff]),expression("R"[t]))) +
    guides(fill = guide_legend(nrow = 1))


  # Clean up
  rm(list = ls()[!ls() %in% c("lls","Poisson_Figure_weeks_2","Poisson_Figure_age_2",
                              "Week_prev_plot_2","Age_prev_plot_2","PCR_sero_prev_plot_2","p_Rt_Reff_a","pars_obs")])
  gc()

  # Return plots
  return(list(
    Poisson_Figure_weeks = Poisson_Figure_weeks_2,
    Poisson_Figure_age = Poisson_Figure_age_2,
    Week_prev_plot = Week_prev_plot_2,
    Age_prev_plot = Age_prev_plot_2,
    PCR_sero_prev_plot = PCR_sero_prev_plot_2,
    p_Rt_Reff_a = p_Rt_Reff_a))
}
