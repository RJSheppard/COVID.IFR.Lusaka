Diagnostic_Plot <- function(fit_num, fit_Model_list, IFRvals,
                               Return_likelihoods_only = F,
                               Return_data_only_plots = F
                               ){
  gc();print(fit_num)
  fit_Model <- fit_Model_list[[fit_num]]
  if(is.null(fit_Model)){return(NA)}

  ################################################
  ################################################
  ### First, load inputs
  pars_obs <- fit_Model$pmcmc_results$inputs$pars_obs
  pcr_det <- pars_obs$pcr_det
  pcr_det_PM <- pars_obs$pcr_det_PM
  if(is.null(pcr_det_PM)){pcr_det_PM <- pars_obs$pcr_det}
  index <- squire:::odin_index(fit_Model$model)
  FullData <- pars_obs$combined_data$Comb_data
  dfj_mcmc_data <- pars_obs$combined_data$dfj_mcmc_data

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
           pcr_pos = cma:::roll_func_10(infs, pcr_det_PM),
           pcr_perc = pcr_pos/max(Sus,na.rm = T))

  Mod_Age_Deaths_Lus <- Mod_Age_Deaths_Lus_pcr %>%
    filter(date > "2020-06-10") %>%
    mutate(Week_gr = as.numeric(cut.Date(date, breaks = "1 week"))) %>%
    group_by(Age_gr, Week_gr, Replicate) %>%
    summarise(Mod_cd_Lus = max(Mod_cd_Lus),
              pcr_perc = pcr_perc[4], # Thursday - middle of the study week
              infs = sum(infs)) %>%
    ungroup() %>% group_by(Age_gr, Replicate) %>%
    mutate(Mod_cd_Lus = c(0,diff(Mod_cd_Lus))) %>% #,
    filter(Week_gr != 1) %>%
    mutate(Week_gr = Week_gr-1) %>%
    merge(FullData)

  ################################################
  ################################################
  ### Calculate Likelihoods
  ## Integrate over mcmc samples
  if(pars_obs$lld =="no scaling no weeks 4-5"){
    LL_Distributions_diag <- lapply(1:length(pars_obs$combined_data$dfj_mcmc_data), function(y){
      tmp_df <- merge(Mod_Age_Deaths_Lus, pars_obs$combined_data$dfj_mcmc_data[[y]]) %>%
        arrange(Week_gr, Age_gr)
      # Calculate the poisson likelihood of burial registrations from Modelled Lusaka Covid deaths (factored by the proportional weekly changes in <5)
      # Added to the mcmc estimates for Lusaka level baseline deaths (these are already factored by the weekly changes based on <5s)
      tmp_df$Mod_cd_morgue <- tmp_df$Mod_cd_Lus
      tmp_df$Mod_tot_ds_morgue <-  tmp_df$Mod_cd_morgue + tmp_df$Mort_ncd_mcmc
      tmp_df$Mod_coin_cpos_ds_morgue <- tmp_df$Mort_ncd_mcmc*tmp_df$pcr_perc
      tmp_df$Mod_pos_ds_morgue <- tmp_df$Mod_cd_morgue + tmp_df$Mod_coin_cpos_ds_morgue

      tmp_df$l_pois <- dpois(x = tmp_df$Bur_regs, lambda = tmp_df$Mod_tot_ds_morgue, log = F)
      tmp_df$l_binom <- dbinom(tmp_df$PosTests, tmp_df$Samples, tmp_df$Mod_pos_ds_morgue/tmp_df$Mod_tot_ds_morgue, log = F)

      tmp_df$Sample_no <- y
      tmp_df <- tmp_df[c("Sample_no", "Week_gr", "Age_gr","Replicate", "l_pois", "l_binom",
                         "Mod_cd_morgue","Mort_ncd_mcmc",#"infs","pcr_perc",
                         "Mod_tot_ds_morgue","Mod_coin_cpos_ds_morgue","Mod_pos_ds_morgue",
                         "Bur_regs")]
    })

  } else if (pars_obs$lld =="no post mortem") {

    LL_Distributions_diag <- lapply(1:length(pars_obs$combined_data$dfj_mcmc_data), function(y){
      tmp_df <- merge(Mod_Age_Deaths_Lus, pars_obs$combined_data$dfj_mcmc_data[[y]]) %>%
        arrange(Week_gr, Age_gr)
      # Calculate the poisson likelihood of burial registrations from Modelled Lusaka Covid deaths (factored by the proportional weekly changes in <5)
      # Added to the mcmc estimates for Lusaka level baseline deaths (these are already factored by the weekly changes based on <5s)
      tmp_df$Mod_cd_morgue <- tmp_df$Mod_cd_Lus*tmp_df$ag1std
      tmp_df$Mod_tot_ds_morgue <-  tmp_df$Mod_cd_morgue + tmp_df$Mort_ncd_mcmc
      tmp_df$Mod_coin_cpos_ds_morgue <- tmp_df$Mort_ncd_mcmc*tmp_df$pcr_perc
      tmp_df$Mod_pos_ds_morgue <- tmp_df$Mod_cd_morgue + tmp_df$Mod_coin_cpos_ds_morgue

      tmp_df$l_pois <- dpois(x = tmp_df$Bur_regs, lambda = tmp_df$Mod_tot_ds_morgue, log = F)
      tmp_df$l_binom <- dbinom(tmp_df$PosTests, tmp_df$Samples, tmp_df$Mod_pos_ds_morgue/tmp_df$Mod_tot_ds_morgue, log = F)

      tmp_df$Sample_no <- y
      tmp_df <- tmp_df[c("Sample_no", "Week_gr", "Age_gr","Replicate", "l_pois", "l_binom",
                         "Mod_cd_morgue","Mort_ncd_mcmc",
                         "Mod_tot_ds_morgue","Mod_coin_cpos_ds_morgue","Mod_pos_ds_morgue",
                         "Bur_regs")]
    })

  } else {
    LL_Distributions_diag <- lapply(1:length(pars_obs$combined_data$dfj_mcmc_data), function(y){
      tmp_df <- merge(Mod_Age_Deaths_Lus, pars_obs$combined_data$dfj_mcmc_data[[y]], all = F) %>%
        arrange(Week_gr, Age_gr)
      # Calculate the poisson likelihood of burial registrations from Modelled Lusaka Covid deaths (factored by the proportional weekly changes in <5)
      # Added to the mcmc estimates for Lusaka level baseline deaths (these are already factored by the weekly changes based on <5s)
      tmp_df$Mod_cd_morgue <- tmp_df$Mod_cd_Lus*tmp_df$ag1std
      tmp_df$Mod_tot_ds_morgue <-  tmp_df$Mod_cd_morgue + tmp_df$Mort_ncd_mcmc
      tmp_df$Mod_coin_cpos_ds_morgue <- tmp_df$Mort_ncd_mcmc*tmp_df$pcr_perc
      tmp_df$Mod_pos_ds_morgue <- tmp_df$Mod_cd_morgue + tmp_df$Mod_coin_cpos_ds_morgue

      tmp_df$l_pois <- dpois(x = tmp_df$Bur_regs, lambda = tmp_df$Mod_tot_ds_morgue, log = F)
      tmp_df$l_binom <- dbinom(tmp_df$PosTests, tmp_df$Samples, tmp_df$Mod_pos_ds_morgue/tmp_df$Mod_tot_ds_morgue, log = F)

      tmp_df$Sample_no <- y
      tmp_df <- tmp_df[c("Sample_no", "Week_gr", "Age_gr","Replicate", "l_pois", "l_binom",
                         "Mod_cd_morgue","Mort_ncd_mcmc",#"infs","pcr_perc",
                         "Mod_tot_ds_morgue","Mod_coin_cpos_ds_morgue","Mod_pos_ds_morgue",
                         "Bur_regs")]
    })
  }

  Res_2d <- data.table::rbindlist(LL_Distributions_diag)
  rm(LL_Distributions_diag)
  gc()
  # Get overall ll for poisson and bin
  if(pars_obs$lld =="no scaling no weeks 4-5"){
    ll_ov_pois <- Res_2d %>% filter(!Week_gr %in% c(4,5)) %>% group_by(Replicate, Sample_no) %>%
      summarise(ll_pois_prod = log(Brobdingnag::prod(Brobdingnag::as.brob(l_pois))))

    ll_ov_bin <- Res_2d %>% group_by(Replicate, Sample_no) %>%
      summarise(l_binom_prod = prod(l_binom))

    ll_ov_pois_bin <- merge(ll_ov_pois, ll_ov_bin) %>%
      ungroup() %>% group_by(Replicate) %>%
      summarise(ll_pois_mean_overall = log(Brobdingnag::sum(exp(Brobdingnag::as.brob(ll_pois_prod)))/length(ll_pois_prod)),
                ll_bin_mean_overall = log(mean(l_binom_prod)))
  } else{
    ll_ov_pois_bin <- Res_2d %>% group_by(Replicate, Sample_no) %>%
      summarise(ll_pois_prod = log(Brobdingnag::prod(Brobdingnag::as.brob(l_pois))),
                l_binom_prod = prod(l_binom)) %>%
      ## Then take the mean of everything
      ungroup() %>% group_by(Replicate) %>%
      summarise(ll_pois_mean_overall = log(Brobdingnag::sum(exp(Brobdingnag::as.brob(ll_pois_prod)))/length(ll_pois_prod)),
                ll_bin_mean_overall = log(mean(l_binom_prod)))

  }

  ################################################
  ################################################
  ## pcr/sero prevalence likelihoods
  pcr_vals <- lapply(1:dim(fit_Model$output)[3], function(x){
    pcr_tmp <- ll_prev_func(df = pars_obs$pcr_df, det = pars_obs$pcr_det, Incidence = "Infections", out = fit_Model$output[,,x], index = index)
  })

  sero_vals <- lapply(1:dim(fit_Model$output)[3], function(x){
    sero_tmp <- ll_prev_func(df = pars_obs$sero_df, det = pars_obs$sero_det, Incidence = "Symptoms", out = fit_Model$output[,,x], index = index, model_params = fit_Model$parameters)
  })

  sero_pcr_ll <- data.frame(Replicate = 1:100,
                            ll_pcr = unlist(pcr_vals),
                            ll_sero = unlist(sero_vals))

  if (pars_obs$lld =="no post mortem") {
    ll_all <- merge(ll_ov_pois_bin, sero_pcr_ll) %>% ungroup() %>%
      mutate(ll_total = ll_pois_mean_overall + ll_pcr + ll_sero)
  } else{
    ll_all <- merge(ll_ov_pois_bin, sero_pcr_ll) %>% ungroup() %>%
      mutate(ll_total = ll_pois_mean_overall + ll_bin_mean_overall + ll_pcr + ll_sero)
  }

  lls <- list(
    ll_all = ll_all,
    ll_bin = log(mean(exp(ll_all$ll_bin_mean_overall))),
    ll_pois = log(Brobdingnag::sum(exp(Brobdingnag::as.brob(ll_all$ll_pois_mean_overall)))/length(ll_all$ll_pois_mean_overall)),
    ll_pcr = log(mean(exp(ll_all$ll_pcr))),
    ll_sero = log(mean(exp(ll_all$ll_sero))),
    ll_tot = log(Brobdingnag::sum(exp(Brobdingnag::as.brob(ll_all$ll_total)))/length(ll_all$ll_total)))

  ################################################
  ################################################
  if(Return_likelihoods_only){return(lls)}
  ################################################
  ################################################

  ll_week <- Res_2d %>% select(Sample_no, Week_gr, Age_gr, Replicate, Mod_cd_morgue, Mort_ncd_mcmc, Mod_tot_ds_morgue, Mod_coin_cpos_ds_morgue, Mod_pos_ds_morgue) %>%
    group_by(Sample_no, Week_gr, Replicate) %>%
    summarise(Mod_cd_morgue = sum(Mod_cd_morgue),
              Mort_ncd_mcmc = sum(Mort_ncd_mcmc),
              Mod_tot_ds_morgue = sum(Mod_tot_ds_morgue),
              Mod_coin_cpos_ds_morgue = sum(Mod_coin_cpos_ds_morgue),
              Mod_pos_ds_morgue = sum(Mod_pos_ds_morgue)
    ) %>% ungroup() %>%
    mutate(Pos_prev = Mod_pos_ds_morgue/Mod_tot_ds_morgue,
           Pos_prev_causal = Mod_cd_morgue/Mod_tot_ds_morgue,
           Pos_prev_coin = Mod_coin_cpos_ds_morgue/Mod_tot_ds_morgue
    ) %>%
    group_by(Week_gr) %>%
    summarise_at(c("Pos_prev", "Pos_prev_causal", "Pos_prev_coin", "Mod_tot_ds_morgue","Mort_ncd_mcmc"), list(median = median, ci = bayestestR::ci), na.rm = TRUE) %>%
    merge(unique(Mod_Age_Deaths_Lus %>% select(Week_gr, date)))

  # Group by week
  Mod_Age_Deaths_Lus_Av_Week <- ll_week %>%
    merge(FullData %>%
            group_by(Week_gr) %>%
            select(-date, -Age_gr) %>%
            summarise_all(sum, na.rm=T), all = T)


  ll_age <- Res_2d %>% select(Sample_no, Week_gr, Age_gr, Replicate, Mod_cd_morgue, Mort_ncd_mcmc, Mod_tot_ds_morgue, Mod_coin_cpos_ds_morgue, Mod_pos_ds_morgue) %>%
    group_by(Sample_no, Age_gr, Replicate) %>%
    summarise(Mod_cd_morgue = sum(Mod_cd_morgue),
              Mort_ncd_mcmc = sum(Mort_ncd_mcmc),
              Mod_tot_ds_morgue = sum(Mod_tot_ds_morgue),
              Mod_coin_cpos_ds_morgue = sum(Mod_coin_cpos_ds_morgue),
              Mod_pos_ds_morgue = sum(Mod_pos_ds_morgue)
    ) %>% ungroup() %>%
    mutate(Pos_prev = Mod_pos_ds_morgue/Mod_tot_ds_morgue,
      Pos_prev_causal = Mod_cd_morgue/Mod_tot_ds_morgue,
      Pos_prev_coin = Mod_coin_cpos_ds_morgue/Mod_tot_ds_morgue
    ) %>%
    group_by(Age_gr) %>%
    summarise_at(c("Pos_prev", "Pos_prev_causal", "Pos_prev_coin", "Mod_tot_ds_morgue","Mort_ncd_mcmc"), list(median = median, ci = bayestestR::ci), na.rm = TRUE) %>%
    merge(unique(Mod_Age_Deaths_Lus %>% select(Age_gr)))


  Mod_Age_Deaths_Lus_Av_Age <- ll_age %>%
    merge(FullData %>%
            group_by(Age_gr) %>%
            select(-date, -Week_gr) %>%
            summarise_all(sum, na.rm=T), all = T) %>%
    mutate(Age_gr_label = ifelse(Age_gr %in% 1:16,
                                 paste0(Age_gr*5-5, "-",Age_gr*5-1),
                                 "80+"))




  ################################################
  ################################################
  ## Plot data first  ##
  ################################################
  ################################################
  # Facet labels: Simple
  Age_groups.labs_Simple <-  Mod_Age_Deaths_Lus_pcr %>% group_by(Age_gr) %>% select(Age_gr) %>% unique() %>%
    mutate(Age_gr_lab = ifelse(Age_gr %in% 1:16, paste0("Age: ",Age_gr*5-5,"-",Age_gr*5-1),
                               "Age: 80+")) %>%
    ungroup() %>% arrange(Age_gr) %>% select(Age_gr_lab) %>%
    unlist()
  names(Age_groups.labs_Simple) <- 1:17


  ## Summarise by age and week
  Poisson_Figure_weeks <- ggplot(Mod_Age_Deaths_Lus_Av_Week, aes(x = date)) +
    geom_point(aes(y = Bur_regs, color = "Burial registrations")) +
    xlab("Date") + ylab("Deaths") +
    theme_minimal() +
    theme(legend.key = element_rect(fill = "white", linetype = 0))

  Mod_Age_Deaths_Lus_Av_Age$Age_gr_label <- factor(Mod_Age_Deaths_Lus_Av_Age$Age_gr_label, levels = Mod_Age_Deaths_Lus_Av_Age$Age_gr_label)

  Poisson_Figure_age <- ggplot(data = Mod_Age_Deaths_Lus_Av_Age, aes(x = Age_gr_label)) +
    geom_point(aes(y = Bur_regs, color = "Burial registrations")) +
    xlab("Age") + ylab("Deaths") +
    theme_minimal()
    theme(legend.key = element_rect(fill = "white", linetype = 0),
          axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))




  ################################################
  ## Plot 6: Weekly prevalence fit
  Week_prev_plot <- ggplot(data = Mod_Age_Deaths_Lus_Av_Week, aes(x = date)) +
    # Prevalence on positive in sample
    geom_point(aes(y = PosTests/Samples, color = "Total +ve deaths")) +
    geom_errorbar(aes(ymin = Hmisc::binconf(PosTests,Samples)[,"Lower"],
                      ymax = Hmisc::binconf(PosTests,Samples)[,"Upper"])) +
    coord_cartesian(ylim = c(0,1)) +
    xlab("Date") +
    theme(plot.title = element_text(size = 10)) +
    theme_minimal() +
    ylab("Post-mortem prevalence") +
    scale_y_continuous(labels = scales::percent)


  ## Plot 7: Age group prevalence fits
  Age_prev_plot <- ggplot(data = Mod_Age_Deaths_Lus_Av_Age, aes(x = Age_gr_label)) +
    # Prevalence on positive in sample
    geom_point(aes(y = PosTests/Samples, color = "Total +ve deaths")) +
    geom_errorbar(aes(ymin = Hmisc::binconf(PosTests,Samples)[,"Lower"],
                      ymax = Hmisc::binconf(PosTests,Samples)[,"Upper"])) +
    xlab("Age") +
    theme(plot.title = element_text(size = 10)) +
    theme_minimal() +
    scale_x_discrete(limits = c(paste0(1:16*5-5, "-",1:16*5-1),"80+")) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) +
    scale_y_continuous(labels = scales::percent) +
    ylab("Post-mortem prevalence")

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
    theme_minimal() +
    scale_color_manual(name=NULL,
                       breaks = c("Attack rate", "PCR %", "Sero %"),
                       values = c("black","darkgoldenrod2","chartreuse4")) +
    scale_y_continuous(labels = scales::percent)


  if(Return_data_only_plots){
    return(list(Poisson_Figure_weeks = Poisson_Figure_weeks + scale_color_manual(name=NULL, breaks = c("Burial registrations"),values = c("black"),guide = guide_legend(override.aes = list(linetype = c(0), shape = c(16), alpha = c(1)))),
                Poisson_Figure_age = Poisson_Figure_age + scale_color_manual(name=NULL, breaks = c("Burial registrations"),values = c("black"),guide = guide_legend(override.aes = list(linetype = c(0), shape = c(16), alpha = c(1)))),
                Week_prev_plot = Week_prev_plot + scale_color_manual(name = NULL, breaks = c("Total +ve deaths"),values = c("black")),
                Age_prev_plot = Age_prev_plot + scale_color_manual(name = NULL, breaks = c("Total +ve deaths"),values = c("black")),
                PCR_sero_prev_plot = PCR_sero_prev_plot))
  }

  if(get_data){
    return(list(lls = lls,
                fit_Model = fit_Model,
                pcr_sero_data = pcr_sero_data,
                Mod_Age_Deaths_Lus_Av_Age = Mod_Age_Deaths_Lus_Av_Age,
                Mod_Age_Deaths_Lus_Av_Week = Mod_Age_Deaths_Lus_Av_Week,
                Mod_Age_Deaths_Lus_Av_Age = Mod_Age_Deaths_Lus_Av_Age,
                Mod_Age_Deaths_Lus_Av_Week = Mod_Age_Deaths_Lus_Av_Week
    ))
  }


  Poisson_Figure_weeks <- Poisson_Figure_weeks +
    geom_line(aes(y=Mod_tot_ds_morgue_median, color = "Model fit")) +
    geom_ribbon(aes(ymin=Mod_tot_ds_morgue_ci$CI_low, ymax=Mod_tot_ds_morgue_ci$CI_high), alpha=0.3) +
    scale_color_manual(name=NULL,
                       breaks = c("Burial registrations", "Model fit"),
                       values = c("black","black"),
                       guide = guide_legend(override.aes = list(linetype = c(0,1),
                                                                shape = c(16,NA),
                                                                alpha = c(1,1)),
                                            nrow = 1))

  Poisson_Figure_age <- Poisson_Figure_age +
    geom_line(aes(y=Mod_tot_ds_morgue_median, color = "Model fit", group = 1)) +
    geom_ribbon(aes(ymin=Mod_tot_ds_morgue_ci$CI_low, ymax=Mod_tot_ds_morgue_ci$CI_high, group = 1), alpha=0.3) +
    scale_color_manual(name=NULL,
                       breaks = c("Burial registrations", "Model fit"),
                       values = c("black","black"),
                       guide = guide_legend(override.aes = list(linetype = c(0,1),
                                                                shape = c(16,NA),
                                                                alpha = c(1,1)),
                                            nrow = 1)) +
    guides(colour = guide_legend(nrow = 1))


  Week_prev_plot <- Week_prev_plot +
    geom_line(aes(y = Pos_prev_coin_median, col = "Coincidental deaths"),linetype="dashed") +
    geom_ribbon(aes(ymin=Pos_prev_coin_ci$CI_low, ymax=Pos_prev_coin_ci$CI_high), alpha=0.3, fill = "darkblue") +
    geom_line(aes(y = Pos_prev_causal_median,col = "Causal deaths"),linetype="dashed") +
    geom_ribbon(aes(ymin=Pos_prev_causal_ci$CI_low, ymax=Pos_prev_causal_ci$CI_high), alpha=0.3, fill = "darkred") +
    geom_line(aes(y = Pos_prev_median, col = "Total +ve deaths"),linetype="dashed") +
    geom_ribbon(aes(ymin=Pos_prev_ci$CI_low, ymax=Pos_prev_ci$CI_high), alpha=0.3) +
    scale_color_manual(name = NULL, breaks = c("Total +ve deaths","Causal deaths","Coincidental deaths"),
                       values = c("black","darkred","darkblue")) +
    guides(colour = guide_legend(nrow = 1))

  Age_prev_plot <- Age_prev_plot +
    geom_line(aes(y = Pos_prev_coin_median, col = "Coincidental deaths", group = 1),linetype="dashed") +
    geom_ribbon(aes(ymin=Pos_prev_coin_ci$CI_low, ymax=Pos_prev_coin_ci$CI_high, group = 1), alpha=0.3, fill = "darkblue") +
    geom_line(aes(y = Pos_prev_causal_median,col = "Causal deaths", group = 1),linetype="dashed") +
    geom_ribbon(aes(ymin=Pos_prev_causal_ci$CI_low, ymax=Pos_prev_causal_ci$CI_high, group = 1), alpha=0.3, fill = "darkred") +
    geom_line(aes(y = Pos_prev_median, col = "Total +ve deaths", group = 1),linetype="dashed") +
    geom_ribbon(aes(ymin=Pos_prev_ci$CI_low, ymax=Pos_prev_ci$CI_high, group = 1), alpha=0.3) +
    scale_color_manual(name = NULL, breaks = c("Total +ve deaths","Causal deaths","Coincidental deaths"),
                       values = c("black","darkred","darkblue")) +
    guides(colour = guide_legend(nrow = 1))

  PCR_sero_prev_plot <- PCR_sero_prev_plot +
    geom_line(aes(x=date, y=pcr_perc_median, color = "PCR %"),linetype="dashed") +
    geom_ribbon(aes(x=date,ymin=pcr_perc_ci$CI_low, ymax=pcr_perc_ci$CI_high), alpha=0.3, fill = "darkgoldenrod1")+
    geom_line(aes(x=date, y=sero_perc_median, color = "Sero %"),linetype="dashed") +
    geom_ribbon(aes(x=date,ymin=sero_perc_ci$CI_low, ymax=sero_perc_ci$CI_high), alpha=0.3, fill = "chartreuse4") +
    geom_line(aes(x = date, y = attack_rate_median, color = "Attack rate"), inherit.aes = F) +
    geom_ribbon(aes(x = date, ymin = attack_rate_ci$CI_low, ymax = attack_rate_ci$CI_high), inherit.aes = F, alpha = 0.3) +
    guides(colour = guide_legend(nrow = 1, override.aes = list(shape = c(NA,20,20))))


  ###############################################################
  ###############################################################
  ## Plot 2: Rt and Reff
  p_Rt_Reff_a <- rt_plot_immunity(fit_Model)$plot +
    scale_fill_manual("",values = c("Reff" = "#48996b","Rt" ="#3f8da7"), labels = c(expression("R"[eff]),expression("R"[t]))) +
    guides(fill = guide_legend(nrow = 1))

  rm(list = ls()[!ls() %in% c("lls","Poisson_Figure_weeks","Poisson_Figure_age",
                              "Week_prev_plot","Age_prev_plot","PCR_sero_prev_plot","p_Rt_Reff_a","pars_obs")])
  gc()

  return(list(lls = lls,
              Poisson_Figure_weeks = Poisson_Figure_weeks,
              Poisson_Figure_age = Poisson_Figure_age,
              Week_prev_plot = Week_prev_plot,
              Age_prev_plot = Age_prev_plot,
              PCR_sero_prev_plot = PCR_sero_prev_plot,
              p_Rt_Reff_a = p_Rt_Reff_a))
}

# Roll function
roll_func_10 <- function(x, det) {
  ret <- rep(0, length(x))
  for(i in seq_along(ret)) {
    to_sum <- x[seq_len(i)][max(1,i-length(det)+1):i]
    ret[i] <- sum(rev(to_sum)*det[seq_along(to_sum)], na.rm = T)
  }
  return(ret)
}


## Prev function
ll_prev_func <- function(df, det, Incidence, out, index, model_params = model_params){

  Sus <- rowSums(out[,index$S])
  Infs <- c(0,as.integer(diff(max(Sus, na.rm = T)-Sus)))

  if(Incidence == "Infections"){
    positives <- roll_func_10(Infs, det)
  }
  if(Incidence == "Symptoms"){
    Symps <- as.integer(rowSums(out[,index$E2]) * model_params$gamma_E)
    positives <- roll_func_10(Symps, det)
  }

  percs <- positives/max(Sus, na.rm = T)

  test_dates <- seq.Date(from = df$date_start, to = df$date_end, by = 1)
  perc_dates <- percs[as.Date(rownames(out)) %in% test_dates]

  ll <- dbinom(df$pos_tests, df$samples, mean(perc_dates), log = TRUE)

  return(ll)
}
