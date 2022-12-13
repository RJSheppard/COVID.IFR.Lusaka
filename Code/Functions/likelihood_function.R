# Run odin model to calculate log-likelihood
# return: Set to 'll' to return the log-likelihood (for MCMC) or to
#
######################
calc_loglikelihood_10_pois_bin_bin_ag1std_agRR <- function(pars, data, squire_model, model_params,
                                                           pars_obs, n_particles,
                                                           forecast_days = 0, return = "ll",
                                                           Rt_args,
                                                           interventions,...) {
  #----------------..
  # specify particle setup
  #----------------..
  switch(return,
         "full" = {
           save_particles <- TRUE
           full_output <- TRUE
           pf_return <- "sample"
         },
         "ll" = {
           save_particles <- FALSE
           forecast_days <- 0
           full_output <- FALSE
           pf_return <- "single"
         },
         {
           stop("Unknown return type to calc_loglikelihood")
         }
  )

  #----------------..
  # (potentially redundant) assertion
  #----------------..
  squire:::assert_in(c("R0", "start_date"), names(pars),
                     message = "Must specify R0, start date to infer")

  #----------------..
  # unpack current params
  #----------------..

  R0 <- pars[["R0"]]
  start_date <- pars[["start_date"]]

  # reporting fraction par if in pars list
  if("rf" %in% names(pars)) {
    assert_numeric(pars[["rf"]])
    pars_obs$phi_death <- pars[["rf"]]
  }


  #----------------..
  # more assertions
  #----------------..
  squire:::assert_pos(R0)
  squire:::assert_date(start_date)

  #----------------..
  # setup model based on inputs and interventions
  #----------------..
  R0_change <- interventions$R0_change
  date_R0_change <- interventions$date_R0_change
  date_contact_matrix_set_change <- interventions$date_contact_matrix_set_change
  date_ICU_bed_capacity_change <- interventions$date_ICU_bed_capacity_change
  date_hosp_bed_capacity_change <- interventions$date_hosp_bed_capacity_change
  date_vaccine_change <- interventions$date_vaccine_change
  date_vaccine_efficacy_infection_change <- interventions$date_vaccine_efficacy_infection_change
  date_vaccine_efficacy_disease_change <- interventions$date_vaccine_efficacy_disease_change

  # change betas
  if (is.null(date_R0_change)) {
    tt_beta <- 0
  } else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_R0_change,
                                                    change = R0_change,
                                                    start_date = start_date,
                                                    steps_per_day = round(1/model_params$dt),
                                                    starting_change = 1)
    model_params$tt_beta <- tt_list$tt
    R0_change <- tt_list$change
    date_R0_change <- tt_list$dates
  }

  # and contact matrixes
  if (is.null(date_contact_matrix_set_change)) {
    tt_contact_matrix <- 0
  } else {

    # here just provide positions for change and then use these to index mix_mat_set
    tt_list <- squire:::intervention_dates_for_odin(dates = date_contact_matrix_set_change,
                                                    change = seq_along(interventions$contact_matrix_set)[-1],
                                                    start_date = start_date,
                                                    steps_per_day = round(1/model_params$dt),
                                                    starting_change = 1)


    model_params$tt_matrix <- tt_list$tt
    model_params$mix_mat_set <- model_params$mix_mat_set[tt_list$change,,]
  }

  # and icu beds
  if (is.null(date_ICU_bed_capacity_change)) {
    tt_ICU_beds <- 0
  } else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_ICU_bed_capacity_change,
                                                    change = interventions$ICU_bed_capacity[-1],
                                                    start_date = start_date,
                                                    steps_per_day = round(1/model_params$dt),
                                                    starting_change = interventions$ICU_bed_capacity[1])
    model_params$tt_ICU_beds <- tt_list$tt
    model_params$ICU_beds <- tt_list$change
  }

  # and hosp beds
  if (is.null(date_hosp_bed_capacity_change)) {
    tt_hosp_beds <- 0
  } else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_hosp_bed_capacity_change,
                                                    change = interventions$hosp_bed_capacity[-1],
                                                    start_date = start_date,
                                                    steps_per_day = round(1/model_params$dt),
                                                    starting_change = interventions$hosp_bed_capacity[1])
    model_params$tt_hosp_beds <- tt_list$tt
    model_params$hosp_beds <- tt_list$change
  }

  # and vaccine coverage
  if (is.null(date_vaccine_change)) {
    tt_vaccine <- 0
  } else {
    tt_list <- squire:::intervention_dates_for_odin(dates = date_vaccine_change,
                                                    change = interventions$max_vaccine[-1],
                                                    start_date = start_date,
                                                    steps_per_day = round(1/model_params$dt),
                                                    starting_change = interventions$max_vaccine[1])
    model_params$tt_vaccine <- tt_list$tt
    model_params$max_vaccine <- tt_list$change
  }

  # and vaccine efficacy infection
  if (is.null(date_vaccine_efficacy_infection_change)) {
    tt_vaccine_efficacy_infection <- 0
  } else {

    # here we just pass the change as a position vector as we need to then
    # index the array of vaccine efficacies
    tt_list <- squire:::intervention_dates_for_odin(dates = date_vaccine_efficacy_infection_change,
                                                    change = seq_along(interventions$vaccine_efficacy_infection)[-1],
                                                    start_date = start_date,
                                                    steps_per_day = round(1/model_params$dt),
                                                    starting_change = 1)

    model_params$tt_vaccine_efficacy_infection <- tt_list$tt

    # here we have to not index the array by the postion vectors that are reutrned by intervention_dates_for_odin
    model_params$vaccine_efficacy_infection <- model_params$vaccine_efficacy_infection[tt_list$change,,]
  }

  # and vaccine efficacy disease
  if (is.null(date_vaccine_efficacy_disease_change)) {
    tt_vaccine_efficacy_disease <- 0
  } else {

    # here we just pass the change as a position vector as we need to then
    # index the array of vaccine efficacies
    tt_list <- squire:::intervention_dates_for_odin(dates = date_vaccine_efficacy_disease_change,
                                                    change = seq_along(interventions$vaccine_efficacy_disease)[-1],
                                                    start_date = start_date,
                                                    steps_per_day = round(1/model_params$dt),
                                                    starting_change = 1)

    model_params$tt_vaccine_efficacy_disease <- tt_list$tt

    # here we have to not index the array by the position vectors that are returned by intervention_dates_for_odin
    model_params$prob_hosp <- model_params$prob_hosp[tt_list$change,,]
  }

  #--------------------..
  # update new R0s based on R0_change and R0_date_change, and Meff_date_change
  #--------------------..
  # and now get new R0s for the R0
  R0 <- squire:::evaluate_Rt_pmcmc(R0_change = R0_change,
                                   R0 = R0,
                                   date_R0_change = date_R0_change,
                                   pars = pars,
                                   Rt_args = Rt_args)

  # which allow us to work out our beta
  beta_set <- squire:::beta_est(squire_model = squire_model,
                                model_params = model_params,
                                R0 = R0)

  #----------------..
  # update the model params accordingly from new inputs
  #----------------..
  model_params$beta_set <- beta_set


  #####################################################
  #####################################################

  #----------------..
  # run the particle filter
  #----------------..
  if (inherits(squire_model, "stochastic")) {

    pf_result <- squire:::run_particle_filter(data = data,
                                              squire_model = squire_model,
                                              model_params = model_params,
                                              model_start_date = start_date,
                                              obs_params = pars_obs,
                                              n_particles = n_particles,
                                              forecast_days = forecast_days,
                                              save_particles = save_particles,
                                              full_output = full_output,
                                              return = pf_return)

  } else if (inherits(squire_model, "deterministic")) {

    pf_result <- run_deterministic_10_pois_bin_bin_ag1std_agRR(data = data,
                                                               squire_model = squire_model,
                                                               model_params = model_params,
                                                               model_start_date = start_date,
                                                               obs_params = pars_obs,
                                                               forecast_days = forecast_days,
                                                               save_history = save_particles,
                                                               return = pf_return)

  }

  # out
  pf_result
}


#############################################################################
#############################################################################
#############################################################################
run_deterministic_10_pois_bin_bin_ag1std_agRR <- function(data,
                                                          squire_model,
                                                          model_params,
                                                          model_start_date = "2020-02-02",
                                                          obs_params = list(phi_cases = 0.1,
                                                                            k_cases = 2,
                                                                            phi_death = 1,
                                                                            k_death = 2,
                                                                            exp_noise = 1e6),
                                                          forecast_days = 0,
                                                          save_history = FALSE,
                                                          return = "ll") {


  # parameter checks
  if (!(return %in% c("full", "ll", "sample", "single"))) {
    stop("return argument must be full, ll, sample", "single")
  }
  if (as.Date(data$date[data$deaths > 0][1], "%Y-%m-%d") < as.Date(model_start_date, "%Y-%m-%d")) {
    stop("Model start date is later than data start date")
  }

  # convert data into particle-filter form
  data <- squire:::particle_filter_data(data = data,
                                        start_date = model_start_date,
                                        steps_per_day = round(1 / model_params$dt))

  model_params$tt_beta <- round(model_params$tt_beta*model_params$dt)
  model_params$tt_contact_matrix <- round(model_params$tt_contact_matrix*model_params$dt)
  model_params$tt_hosp_beds <- round(model_params$tt_hosp_beds*model_params$dt)
  model_params$tt_ICU_beds <- round(model_params$tt_ICU_beds*model_params$dt)

  #set up model
  model_func <- squire_model$odin_model(user = model_params,
                                        unused_user_action = "ignore")

  # steps for the deterministic
  steps <- c(0, data$day_end)
  fore_steps <- seq(data$day_end[nrow(data)], length.out = forecast_days + 1L)
  steps <- unique(c(steps,fore_steps))

  # model run
  if("atol" %in% names(obs_params) && "rtol" %in% names(obs_params)) {
    squire:::assert_numeric(obs_params$atol)
    atol <- obs_params$atol
    squire:::assert_numeric(obs_params$rtol)
    rtol <- obs_params$rtol
  } else {
    atol <- 1e-6
    rtol <- 1e-6
  }

  ##############################################################################
  ##############################################################################
  ######################## Likelihood calculations #############################
  ##############################################################################
  ##############################################################################
  ### Get model output

  if(!(obs_params$lld %in% c("","no post mortem","no scaling no weeks 4-5"))){
  stop("likelihood not found")
  }
  out <- model_func$run(t = seq(0, tail(steps,1), 1), atol = atol, rtol = rtol)
  index <- squire:::odin_index(model_func)

  ### Assign dates to the output
  date <- data$date[[1]] + seq_len(nrow(out)) - 1L
  rownames(out) <- as.character(date)
  attr(out, "date") <- date

  ### Data for comparison: age and weeks
  pcr_det <- obs_params$pcr_det
  pcr_det_PM <- obs_params$pcr_det_PM
  frac_reg <- obs_params$frac_reg # Mortuary captures 80% of deaths in Lusaka

  ### Get model data for comparison
  Days_for_comparison <-c(seq.Date(from = as.Date("2020-06-14"),
                                   to = as.Date("2020-10-02"),
                                   by = "week"),as.Date("2020-10-04"))


  Mod_Deaths_Age <- apply(out[as.Date(rownames(out)) %in% Days_for_comparison,index$D], 2, diff)
  colnames(Mod_Deaths_Age) <- 1:17
  rownames(Mod_Deaths_Age) <- 1:16
  Mod_Deaths_Age <- reshape2::melt(data = Mod_Deaths_Age, value.name = "Mod_cd_Lus", varnames= c("Week_gr","Age_gr")) # 16 time periods, 17 age groups.

  # Comb_data <- obs_params$combined_data$Comb_data
  Mod_Deaths_Age <- merge(Mod_Deaths_Age, obs_params$combined_data$Comb_data)

  # Mod_Deaths_Age %>% filter(Age_gr ==2, Week_gr ==1)

  ##############################################################################
  ##############################################################################

  # Get modelled pcr prevalence by age group
  # system.time(
  pcr_pos_age <- apply(out[,index$S],2,function(x){
    x <- as.integer(x)
    infs <- c(0,diff(max(x, na.rm = T)-x))
    pcr_positive <- roll_func_10(infs, pcr_det_PM)
    pcr_perc <- pcr_positive/max(x, na.rm = T)
  })
  # )

  # Get modelled pcr prevalence by week and shape for comparison
  pcr_perc_age <- pcr_pos_age[as.Date(rownames(out)) %in% (Days_for_comparison[-1]-3),] # This could be -3 or -4 to get Weds or Thurs.
  colnames(pcr_perc_age) <- 1:17
  pcr_perc_age <- reshape2::melt(data = pcr_perc_age, value.name = "pcr_perc", varnames= c("Week_gr","Age_gr")) # 16 time periods, 17 age groups.

  Mod_Deaths_Age <- merge(Mod_Deaths_Age, pcr_perc_age)
  Mod_Deaths_Age$Mod_cd_frac <- Mod_Deaths_Age$Mod_cd_Lus*frac_reg

  ##############################################################################
  ##############################################################################

  # if(obs_params$lld =="remove_weeks_4_and_5"){
    # Mod_Deaths_Age <- Mod_Deaths_Age %>% dplyr::filter(Week_gr != 4, Week_gr != 5)
    # obs_params$combined_data$dfj_mcmc_data <- lapply(obs_params$combined_data$dfj_mcmc_data, function(x){return(x %>% dplyr::filter(Week_gr != 4, Week_gr != 5))})}
  # if(obs_params$lld =="remove_under_5s"){
    # Mod_Deaths_Age <- Mod_Deaths_Age %>% dplyr::filter(Age_gr != 1)
    # obs_params$combined_data$dfj_mcmc_data <- lapply(obs_params$combined_data$dfj_mcmc_data, function(x){return(x %>% dplyr::filter(Age_gr != 1))})}
  # if(obs_params$lld =="remove_weeks_4_and_5_and_age_under_5s"){
    # Mod_Deaths_Age <- Mod_Deaths_Age %>% dplyr::filter(Age_gr != 1,Week_gr != 4,Week_gr != 5)}
    # obs_params$combined_data$dfj_mcmc_data <- lapply(obs_params$combined_data$dfj_mcmc_data, function(x){return(x %>% dplyr::filter(Age_gr != 1,Week_gr != 4,Week_gr != 5))})}

  ll_aw_pois <- 0
  ll_aw_bin <- 0

  if(obs_params$lld == "no post mortem"){
    LL_Distributions <- lapply(1:length(obs_params$combined_data$dfj_mcmc_data), function(y){
      tmp_df <- merge(Mod_Deaths_Age, obs_params$combined_data$dfj_mcmc_data[[y]], all = F)
      ll_pois_tmp <- sum(dpois(x = tmp_df$Bur_regs, lambda = (tmp_df$Mod_cd_frac * tmp_df$ag1std + tmp_df$Mort_ncd_mcmc), log = T))
      return(data.frame(ll_pois = ll_pois_tmp))
    })

    ll_aw_pois <- as.numeric(log(Brobdingnag::sum(exp(Brobdingnag::as.brob(sapply(LL_Distributions,"[[",1))))/length(LL_Distributions)))

  } else if(obs_params$lld == "no scaling no weeks 4-5"){
    LL_Distributions <- lapply(1:length(obs_params$combined_data$dfj_mcmc_data), function(y){
      tmp_df <- merge(Mod_Deaths_Age, obs_params$combined_data$dfj_mcmc_data[[y]], all = F)
      ll_binom_tmp <- sum(dbinom(x = tmp_df$PosTests, size =  tmp_df$Samples, prob =  (tmp_df$Mod_cd_frac + tmp_df$Mort_ncd_mcmc*tmp_df$pcr_perc)/(tmp_df$Mod_cd_frac + tmp_df$Mort_ncd_mcmc), log = T))
      tmp_df <- dplyr::filter(tmp_df, !Week_gr %in% c(4,5))
      ll_pois_tmp <- sum(dpois(x = tmp_df$Bur_regs, lambda = (tmp_df$Mod_cd_frac + tmp_df$Mort_ncd_mcmc), log = T))
      return(data.frame(ll_pois = ll_pois_tmp, ll_binom = ll_binom_tmp))
    })
    ll_aw_pois <- as.numeric(log(Brobdingnag::sum(exp(Brobdingnag::as.brob(sapply(LL_Distributions,"[[",1))))/length(LL_Distributions)))
    ll_aw_bin <- as.numeric(log(Brobdingnag::sum(exp(Brobdingnag::as.brob(sapply(LL_Distributions,"[[",2))))/length(LL_Distributions)))

  } else {
    LL_Distributions <- lapply(1:length(obs_params$combined_data$dfj_mcmc_data), function(y){
      tmp_df <- merge(Mod_Deaths_Age, obs_params$combined_data$dfj_mcmc_data[[y]], all = F)
      ll_pois_tmp <- sum(dpois(x = tmp_df$Bur_regs, lambda = (tmp_df$Mod_cd_frac*tmp_df$ag1std + tmp_df$Mort_ncd_mcmc), log = T))
      ll_binom_tmp <- sum(dbinom(x = tmp_df$PosTests, size =  tmp_df$Samples, prob =  (tmp_df$Mod_cd_frac * tmp_df$ag1std + tmp_df$Mort_ncd_mcmc * tmp_df$pcr_perc)/(tmp_df$Mod_cd_frac*tmp_df$ag1std + tmp_df$Mort_ncd_mcmc), log = T))
      # tmp_df$Mod_ncd_Lus <- tmp_df$Bur_regs/tmp_df$ag1std - tmp_df$Mod_cd_Lus
      # tmp_df$Mod_ncd_Lus <- ifelse(tmp_df$Mod_ncd_Lus<0,0,tmp_df$Mod_ncd_Lus)
      # tmp_df$Mod_pos_ds_Lus <- tmp_df$Mod_cd_Lus+tmp_df$Mod_ncd_Lus*tmp_df$pcr_perc # total positive deaths (covid and coincidental)
      # tmp_df$Pos_prev <- ifelse(tmp_df$Mod_pos_ds_Lus<=tmp_df$Bur_regs/tmp_df$ag1std,tmp_df$Mod_pos_ds_Lus/(tmp_df$Bur_regs/tmp_df$ag1std), 0.999)
      # ll_binom_tmp <- sum(dbinom(tmp_df$PosTests, tmp_df$Samples, prob = tmp_df$Pos_prev, log = T))
      return(data.frame(ll_pois = ll_pois_tmp, ll_binom = ll_binom_tmp))
    })

    ll_aw_pois <- as.numeric(log(Brobdingnag::sum(exp(Brobdingnag::as.brob(sapply(LL_Distributions,"[[",1))))/length(LL_Distributions)))
    ll_aw_bin <- as.numeric(log(Brobdingnag::sum(exp(Brobdingnag::as.brob(sapply(LL_Distributions,"[[",2))))/length(LL_Distributions)))

  }

  ##############################################################################
  ##############################################################################

  ################################################
  ################################################
  ## And the ll for the combined prevalence.
  # llc <- 0
  # if(all(unlist(lapply(list(obs_params$comb_df, obs_params$comb_det),function(x){!is.null(x)})))){
  #   comb_df <- obs_params$comb_df
  #   comb_det <- obs_params$comb_det
  #
  #   ### Model PCR: Number of Infections/combined prevalence
  #   Sus <- rowSums(out[,index$S])
  #   infs <- c(0,as.integer(diff(max(Sus)-Sus)))
  #   comb_positive <- roll_func_10(infs, comb_det)
  #   comb_perc <- comb_positive/max(Sus)
  #
  #   # comb_dates <- c(comb_df$date_start, comb_df$date_end, comb_df$date_start + as.integer((comb_df$date_end - comb_df$date_start)/2))#list(comb_df$date_end, comb_df$date_start, comb_df$date_start + as.integer((comb_df$date_end - comb_df$date_start)/2))
  #   comb_dates <- seq.Date(from = comb_df$date_start, to = comb_df$date_end, by = 1)
  #   comb_perc <- comb_perc[as.Date(rownames(out)) %in% comb_dates]
  #
  #   # likelihood of model obvs
  #   llc <- dbinom(comb_df$comb_pos, comb_df$samples, mean(comb_perc), log = TRUE)
  # }
  ################################################
  ################################################

  ################################################
  ################################################
  ## ll for pcr prevalence.
  ll_pcr <- 0
  if(all(unlist(lapply(list(obs_params$pcr_df, obs_params$pcr_det),function(x){!is.null(x)})))){
    ll_pcr <- ll_prev_func(df = obs_params$pcr_df, det = obs_params$pcr_det, Incidence = "Infections", out = out, index = index)
  }

  ## ll for pcr prevalence.
  ll_sero <- 0
  if(all(unlist(lapply(list(obs_params$sero_df, obs_params$sero_det),function(x){!is.null(x)})))){
    ll_sero <- ll_prev_func(df = obs_params$sero_df, det = obs_params$sero_det, Incidence = "Symptoms", out = out, index = index, model_params = model_params)
  }

  ################################################
  ################################################

  ##############################################################################
  ##############################################################################
  ll <- sum(ll_aw_pois, ll_aw_bin, ll_pcr, ll_sero)
  ##############################################################################
  ##############################################################################

  # format similar to particle_filter nomenclature
  pf_results <- list()
  pf_results$log_likelihood <- ll

  # pf_results$ll_pois <- LL_int %>% dplyr::filter(LL =="ll_pois") %>% dplyr::pull(ll)
  # pf_results$ll_bin <- LL_int %>% dplyr::filter(LL =="ll_bin") %>% dplyr::pull(ll)
  # pf_results$ll_comb <- llc

  # single returns final state
  if (save_history) {
    pf_results$states <- out
  } else if (return == "single") {
    pf_results$sample_state <- out[nrow(out), ]
  }

  # create returned object
  if (return == "ll") {
    ret <- pf_results$log_likelihood
  } else if (return == "sample") {
    ret <- pf_results$states
  } else if (return == "single" || return == "full") {
    ret <- pf_results
  }

  ret
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
    Symps <- rowSums(out[,index$E2]) * model_params$gamma_E
    positives <- roll_func_10(Symps, det)
  }

  percs <- positives/max(Sus, na.rm = T)

  test_dates <- seq.Date(from = df$date_start, to = df$date_end, by = 1)
  perc_dates <- percs[as.Date(rownames(out)) %in% test_dates]

  ll <- dbinom(df$pos_tests, df$samples, mean(perc_dates), log = TRUE)

  return(ll)
}
