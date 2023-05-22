
#' Fit squire model using rw splines
#'
#' @inheritParams squire::pmcmc
#' @inheritParams squire::parameters_explicit_SEEIR
#' @param hosp_beds General Hospital Beds
#' @param icu_beds ICU Beds
#' @param rw_duration Random Walk/Spline Duration. Default = 14 days
#' @param reporting_fraction_bounds Region to vary reporting fraction
#' @param combined_data Data to fit to
#'
#' @return Model fit from [squire:::pmcmc]
#' @export
#'
fit_spline_rt <- function(data,
                          country,
                          population,
                          baseline_contact_matrix = baseline_contact_matrix,
                          reporting_fraction=1,
                          reporting_fraction_bounds=NULL,
                          n_mcmc = 10000,
                          replicates = 100,
                          rw_duration = 14,
                          hosp_beds = 10000000000,
                          icu_beds = 10000000000,
                          sero_sens = 0.9,
                          pcr_sens = 0.95,
                          comb_df = NULL,
                          pcr_df = NULL,
                          pcr_det = NULL,
                          sero_df = NULL,
                          IFR_slope_bounds = NULL,
                          IFR_multiplier_bounds = NULL,
                          combined_data = NULL,
                          drj_mcmc = NULL,
                          combined_data_week = NULL,
                          frac_reg = 1,
                          log_likelihood = NULL,
                          lld =NULL,
                          k_death = 2,
                          Prior_Rt_rw_unif_lim = NULL,

                          ...) {
# browser()
  ## -----------------------------------------------------------------------------
  ## Step 1 DATA CLEANING AND ORDERING
  ## -----------------------------------------------------------------------------

  # order data
  data <- data[order(data$date),]
  data$date <- as.Date(data$date)

  # and remove the rows with no data up to the first date that a death was reported
  first_report <- which(data$deaths>0)[1]
  missing <- which(data$deaths == 0 | is.na(data$deaths))
  to_remove <- missing[missing<first_report]
  if(length(to_remove) > 0) {
    if(length(to_remove) == (nrow(data)-1)) {
      data <- data[-head(to_remove,-1),]
    } else {
      data <- data[-to_remove,]
    }
  }

  combined_data <- combined_data %>% dplyr::arrange(Age_gr, Week_gr)
  if(length(dim(drj_mcmc)) != 3){stop("Baseline mortality must be an array")}
  if(!(lld %in% c("","no week 5 PM","no scaling no weeks 4-5","no under 5s"))){stop("likelihood not found")}
  if(lld == "no scaling no weeks 4-5"){filter_vec <- !combined_data$Week_gr %in% c(4,5)} else
    if(lld == "no week 5 PM"){filter_vec <- combined_data$Week_gr != 5} else
      if(lld == "no under 5s") {filter_vec <- combined_data$Age_gr != 1} else {filter_vec <- NULL}


  ## -----------------------------------------------------------------------------
  ## Step 2a: PMCMC SETUP
  ## -----------------------------------------------------------------------------

  # dat_0 is just the current date now
  date_0 <- max(data$date)

  # what is the date of first death
  null_na <- function(x) {if(is.null(x)) {NA} else {x}}
  min_death_date <- data$date[which(data$deaths>0)][1]

  # We set the R0_change here to be 1 everywhere to effectively turn off mobility
  R0_change <- rep(1, nrow(data))
  date_R0_change <- data$date
  R0_change <- R0_change[as.Date(date_R0_change) <= date_0]
  date_R0_change <- date_R0_change[as.Date(date_R0_change) <= date_0]

  # pmcmc args
  n_particles <- 2 # we use the deterministic model now so this does nothing (makes your life quicker and easier too)
  n_chains <- 8 # number of chains
  start_adaptation <- as.integer(0.1*n_mcmc) # how long before adapting

  # parallel call
  suppressWarnings(future::plan(future::multisession()))

  # Defualt parameter edges for pmcmc
  R0_min <- 1.6
  R0_max <- 5.6
  last_start_date <- as.Date(null_na(min_death_date))-10
  first_start_date <- as.Date(null_na(min_death_date))-55
  start_date <- as.Date(null_na(min_death_date))-50

  # These 4 parameters do nothign as setting R0_change to 1
  Meff_min <- -2
  Meff_max <- 2
  Meff_pl_min <- 0
  Meff_pl_max <- 1
  Rt_shift_min <- 0
  Rt_shift_max <- 0.001
  Rt_shift_scale_min <- 0.1
  Rt_shift_scale_max <- 10


  ## -----------------------------------------------------------------------------
  ## Step 2b: Sourcing suitable starting conditions
  ## -----------------------------------------------------------------------------

  date_start <- data$date[which(cumsum(data$deaths)>10)[1]] - 50
  R0_start <- 3

  # These are the the initial conditions now loaded from our previous run.
  R0_start <- min(max(R0_start, R0_min), R0_max)
  date_start <- min(max(as.Date(start_date), as.Date(first_start_date)), as.Date(last_start_date))

  # again these all do nothing
  Meff_start <- min(max(0, Meff_min), Meff_max)
  Meff_pl_start <- min(max(0.5, Meff_pl_min), Meff_pl_max)
  Rt_shift_start <- min(max(0.0005, Rt_shift_min), Rt_shift_max)
  Rt_shift_scale_start <- min(max(5, Rt_shift_scale_min), Rt_shift_scale_max)

  # Our random walk parameters start after the Meff change
  # Basically just set this suitably far back in the past
  date_Meff_change <- date_start - 1

  ## -----------------------------------------------------------------------------
  ## Step 2c: Spline set up
  ## -----------------------------------------------------------------------------

  last_shift_date <- as.Date(date_Meff_change) + 1
  remaining_days <- as.Date(date_0) - last_shift_date - 14 # reporting delay in place

  # how many spline pars do we need
  Rt_rw_duration <- rw_duration # i.e. we fit with a 2 week duration for our random walks.
  rw_needed <- as.numeric(ceiling(remaining_days/Rt_rw_duration))

  # set up rw pars
  pars_init_rw <- as.list(rep(0, rw_needed))
  pars_min_rw <- as.list(rep(-5, rw_needed))
  pars_max_rw <- as.list(rep(5, rw_needed))
  pars_discrete_rw <- as.list(rep(FALSE, rw_needed))
  names(pars_init_rw) <- names(pars_min_rw) <- names(pars_max_rw) <- names(pars_discrete_rw) <- paste0("Rt_rw_", seq_len(rw_needed))

  ## -----------------------------------------------------------------------------
  ## Step 2d: PMCMC parameter set up
  ## -----------------------------------------------------------------------------

  # seroconversion data from brazeau report 34
  # sero_sens = 0.9
  prob_conversion <-  cumsum(dgamma(0:300,shape = 5, rate = 1/2))/max(cumsum(dgamma(0:300,shape = 5, rate = 1/2)))
  sero_det <- cumsum(dweibull(0:300, 3.669807, scale = 143.7046))
  sero_det <- prob_conversion-sero_det
  sero_det[sero_det < 0] <- 0
  sero_det <- sero_det/max(sero_det)*sero_sens  # assumed maximum test sensitivitys

  # from Hay et al 2021 Science (actually from preprint)
  # pcr_sens = 0.95
  if(is.null(pcr_det)){
    pcr_det <- c(9.206156e-13, 9.206156e-13, 3.678794e-01, 9.645600e-01,
                         9.575796e-01, 9.492607e-01, 9.393628e-01, 9.276090e-01,
                         9.136834e-01, 8.972309e-01, 8.778578e-01, 8.551374e-01,
                         8.286197e-01, 7.978491e-01, 7.623916e-01, 7.218741e-01,
                         6.760375e-01, 6.248060e-01, 5.683688e-01, 5.072699e-01,
                         4.525317e-01, 4.036538e-01, 3.600134e-01, 3.210533e-01,
                         2.862752e-01, 2.552337e-01, 2.275302e-01, 2.028085e-01,
                         1.807502e-01, 1.610705e-01, 1.435151e-01, 1.278563e-01,
                         1.138910e-01, 1.014375e-01, 9.033344e-02)
    pcr_det <- (pcr_det/max(pcr_det))*pcr_sens
  }

  # PMCMC Parameters
  pars_init = list('start_date' = date_start,
                   'R0' = R0_start,
                   'Meff' = Meff_start,
                   'Meff_pl' = Meff_pl_start,
                   "Rt_shift" = 0,
                   "Rt_shift_scale" = Rt_shift_scale_start)
  pars_min = list('start_date' = first_start_date,
                  'R0' = R0_min,
                  'Meff' = Meff_min,
                  'Meff_pl' = Meff_pl_min,
                  "Rt_shift" = Rt_shift_min,
                  "Rt_shift_scale" = Rt_shift_scale_min)
  pars_max = list('start_date' = last_start_date,
                  'R0' = R0_max,
                  'Meff' = Meff_max,
                  'Meff_pl' = Meff_pl_max,
                  "Rt_shift" = Rt_shift_max,
                  "Rt_shift_scale" = Rt_shift_scale_max)
  pars_discrete = list('start_date' = TRUE, 'R0' = FALSE, 'Meff' = FALSE,
                       'Meff_pl' = FALSE, "Rt_shift" = FALSE, "Rt_shift_scale" = FALSE)
  pars_obs = list(phi_cases = 1, k_cases = 2, phi_death = 1, k_death = k_death, exp_noise = 1e6,
                  sero_det = sero_det, pcr_det = pcr_det,
                  combined_data = combined_data[,c("Samples","PosTests","BurRegs")],
                  combined_data_week = combined_data_week, drj_mcmc = list(drj_mcmc_data_baseline = drj_mcmc[,"Mort_ncd_mcmc",],
                                                                           drj_mcmc_data_baseline_agstd = drj_mcmc[,"ag1std",]),
                  lld = lld, filter_vec = filter_vec,
                  frac_reg = frac_reg)

  # add in the spline list
  pars_init <- append(pars_init, pars_init_rw)
  pars_min <- append(pars_min, pars_min_rw)
  pars_max <- append(pars_max, pars_max_rw)
  pars_discrete <- append(pars_discrete, pars_discrete_rw)

  # add reporting bounds if given
  if(!is.null(reporting_fraction_bounds)){
    pars_init <- append(pars_init, c("rf"=reporting_fraction_bounds[1]))
    pars_min <- append(pars_min, c("rf"=reporting_fraction_bounds[2]))
    pars_max <- append(pars_max, c("rf"=reporting_fraction_bounds[3]))
    pars_discrete <- append(pars_discrete, c("rf"=FALSE))
  }


  if(!is.null(comb_df)){
    pars_obs <- append(pars_obs, c("comb_df" = list(comb_df)))

    pcr_det_sero_len <- c(pcr_det, rep(0, length(sero_det)-length(pcr_det)))
    comb_det <-   1-((1-c(0, 0, 0, 0, head(sero_det, -4))) * (1-pcr_det_sero_len))
    pars_obs <- append(pars_obs, c("comb_det" = list(comb_det)))

  }


  if(!is.null(pcr_df)){
    pars_obs <- append(pars_obs, c("pcr_df" = list(pcr_df)))
  }

  if(!is.null(sero_df)){
    pars_obs <- append(pars_obs, c("sero_df" = list(sero_df)))
  }

######################################################################
  if(!is.null(IFR_slope_bounds)){
    pars_init <- append(pars_init, c("IFR_slope"=IFR_slope_bounds[1]))
    pars_min <- append(pars_min, c("IFR_slope"=IFR_slope_bounds[2]))
    pars_max <- append(pars_max, c("IFR_slope"=IFR_slope_bounds[3]))
    pars_discrete <- append(pars_discrete, c("IFR_slope"=FALSE))
  }

  if(!is.null(IFR_multiplier_bounds)){
    pars_init <- append(pars_init, c("IFR_mult"=IFR_multiplier_bounds[1]))
    pars_min <- append(pars_min, c("IFR_mult"=IFR_multiplier_bounds[2]))
    pars_max <- append(pars_max, c("IFR_mult"=IFR_multiplier_bounds[3]))
    pars_discrete <- append(pars_discrete, c("IFR_mult"=FALSE))
  }
#######################################################################
  # Covariance Matrix
  proposal_kernel <- diag(length(names(pars_init))) * 0.3
  rownames(proposal_kernel) <- colnames(proposal_kernel) <- names(pars_init)
  proposal_kernel["start_date", "start_date"] <- 1.5

  # MCMC Functions - Prior and Likelihood Calculation
  logprior <- function(pars){
    ret <- dunif(x = pars[["start_date"]], min = -55, max = -10, log = TRUE) +
      dnorm(x = pars[["R0"]], mean = 3, sd = 1, log = TRUE) +
      dnorm(x = pars[["Meff"]], mean = 0, sd = 1, log = TRUE) +
      dunif(x = pars[["Meff_pl"]], min = 0, max = 1, log = TRUE) +
      dnorm(x = pars[["Rt_shift"]], mean = 0, sd = 1, log = TRUE) +
      dunif(x = pars[["Rt_shift_scale"]], min = 0.1, max = 10, log = TRUE)

    # get rw spline parameters
    if(any(grepl("Rt_rw", names(pars)))) {
      Rt_rws <- pars[grepl("Rt_rw", names(pars))]
      for (i in seq_along(Rt_rws)) {
        # ret <- ret + dnorm(x = Rt_rws[[i]], mean = 0, sd = 0.2, log = TRUE)
        ret <- ret + dunif(x = Rt_rws[[i]], min = -Prior_Rt_rw_unif_lim, max = Prior_Rt_rw_unif_lim, log = TRUE) # default = 1
      }
    }

    if(!is.null(reporting_fraction_bounds)){
      ret <- ret + dunif(pars[["rf"]], min = pars_min[["rf"]], max = pars_max[["rf"]]) #???+ dnorm(x = Rt_rws[[i]], mean = 0, sd = 0.2, log = TRUE) uniform dist?
    }
##################################################################
    if(!is.null(IFR_slope_bounds)){
      ret <- ret + dunif(pars[["IFR_slope"]], min = pars_min[["IFR_slope"]], max = pars_max[["IFR_slope"]]) #???+ dnorm(x = Rt_rws[[i]], mean = 0, sd = 0.2, log = TRUE) uniform dist?
    }

    if(!is.null(IFR_multiplier_bounds)){
      ret <- ret + dunif(pars[["IFR_mult"]], min = pars_min[["IFR_mult"]], max = pars_max[["IFR_mult"]]) #???+ dnorm(x = Rt_rws[[i]], mean = 0, sd = 0.2, log = TRUE) uniform dist?
    }
###################################################################

    return(ret)
  }

  ## -----------------------------------------------------------------------------
  ## Step 3: Run PMCMC
  ## -----------------------------------------------------------------------------

  # mixing matrix - assume is same as country as whole
  # mix_mat <- squire::get_mixing_matrix(country)

  # run the pmcmc
  res <- squire::pmcmc(data = data,
                       n_mcmc = n_mcmc,
                       log_prior = logprior,
                       n_particles = 1,
                       steps_per_day = 1,
                       log_likelihood = log_likelihood,
                       reporting_fraction = reporting_fraction,
                       # squire_model = squire:::explicit_model(),
                       squire_model = squire:::deterministic_model(),
                       output_proposals = FALSE,
                       n_chains = n_chains,
                       pars_init = pars_init,
                       pars_min = pars_min,
                       pars_max = pars_max,
                       pars_discrete = pars_discrete,
                       pars_obs = pars_obs,
                       proposal_kernel = proposal_kernel,
                       population = population,
                       baseline_contact_matrix = baseline_contact_matrix,
                       R0_change = R0_change,
                       date_R0_change = date_R0_change,
                       Rt_args = squire:::Rt_args_list(
                         date_Meff_change = date_Meff_change,
                         scale_Meff_pl = TRUE,
                         Rt_shift_duration = 1,
                         Rt_rw_duration = Rt_rw_duration),
                       burnin = ceiling(n_mcmc/10),
                       seeding_cases = 5,
                       replicates = replicates,
                       required_acceptance_ratio = 0.20,
                       start_adaptation = start_adaptation,
                       baseline_hosp_bed_capacity = hosp_beds,
                       baseline_ICU_bed_capacity = icu_beds,
                       ...)

  ## remove things so they don't atke up so much memory when you save them :)

  # Add the prior
  res$pmcmc_results$inputs$prior <- as.function(c(formals(logprior),
                                                  body(logprior)),
                                                envir = new.env(parent = environment(stats::acf)))

  # remove states to keep object memory save down
  for(i in seq_along(res$pmcmc_results$chains)) {
    res$pmcmc_results$chains[[i]]$states <- NULL
    res$pmcmc_results$chains[[i]]$covariance_matrix <- tail(res$pmcmc_results$chains$chain1$covariance_matrix,1)
  }

  # Add sero_det and pcr_det to output
  res$pmcmc_results$inputs$pars_obs$sero_sens <- sero_sens
  res$pmcmc_results$inputs$pars_obs$sero_det <- sero_det
  res$pmcmc_results$inputs$pars_obs$pcr_sens <- pcr_sens
  res$pmcmc_results$inputs$pars_obs$pcr_det <- pcr_det

  return(res)

}

#' Extract PCR prevalence and seroprevalence from squire model fit
#'
#' @param res Output of [[squire::pmcmc]]
seroprev_df <- function(res){

  if(is.null(res)){return(NULL)}

  # seroconversion data from brazeay report 34
  sero_sens <- res$pmcmc_results$inputs$pars_obs$sero_sens
  sero_det <- res$pmcmc_results$inputs$pars_obs$sero_det

  # from Hay et al 2021 Science (actually from preprint)
  pcr_sens <- res$pmcmc_results$inputs$pars_obs$pcr_sens
  pcr_det <- res$pmcmc_results$inputs$pars_obs$pcr_det

  # additional_functions for rolling
  roll_func <- function(x, det) {
    l <- length(det)
    ret <- rep(0, length(x))
    for(i in seq_along(ret)) {
      to_sum <- tail(x[seq_len(i)], length(det))
      ret[i] <- sum(rev(to_sum)*head(det, length(to_sum)))
    }
    return(ret)
  }

  # get symptom onset data
  date_0 <- max(res$pmcmc_results$inputs$data$date)
  inf <- squire::format_output(res, c("S"), date_0 = max(res$pmcmc_results$inputs$data$date)) %>%
    na.omit() %>%
    mutate(S = .data$y) %>%
    group_by(replicate) %>%
    mutate(infections = c(0, as.integer(diff(max(.data$S)-.data$S)))) %>%
    select(replicate, t, date, .data$S, .data$infections)

  # correctly format
  inf <- left_join(inf,
                   squire::format_output(
                     res, c("infections"),
                     date_0 = max(res$pmcmc_results$inputs$data$date)
                   ) %>%
                     mutate(symptoms = as.integer(.data$y)) %>%
                     select(replicate, t, .data$date, .data$symptoms),
                   by = c("replicate", "t", "date"))

  pcr_det_sero_len <- c(pcr_det, rep(0, length(sero_det)-length(pcr_det)))
  combined_det <-   1-((1-c(0, 0, 0, 0, head(sero_det, -4))) * (1-pcr_det_sero_len))

  inf <- inf %>%
    group_by(replicate) %>%
    na.omit() %>%
    mutate(pcr_positive = roll_func(.data$infections, pcr_det),
           sero_positive = roll_func(.data$symptoms, sero_det),
           combined_positive = roll_func(.data$infections, combined_det),
           ps_ratio = .data$pcr_positive/.data$sero_positive,
           sero_perc = .data$sero_positive/max(.data$S,na.rm = TRUE),
           pcr_perc = .data$pcr_positive/max(.data$S,na.rm = TRUE),
           combined_perc = .data$combined_positive/max(.data$S,na.rm = TRUE)) %>%
    ungroup

  inf$reporting_fraction <- res$pmcmc_results$inputs$pars_obs$phi_death

  return(inf)

}

Summ_sero_pcr_data <- function(x){
  if(is.null(x)){return(NULL)}
  x %>% group_by(date) %>%
    summarise(median_pcr = median(pcr_perc)*100, min_pcr = min(pcr_perc)*100, max_pcr = max(pcr_perc)*100, ci_low_pcr = 100*bayestestR::ci(pcr_perc)$CI_low, ci_high_pcr = 100*bayestestR::ci(pcr_perc)$CI_high,
              median_sero = median(sero_perc)*100, min_sero = min(sero_perc)*100, max_sero = max(sero_perc)*100, ci_low_sero = 100*bayestestR::ci(sero_perc)$CI_low, ci_high_sero = 100*bayestestR::ci(sero_perc)$CI_high,
              median_combined = median(combined_perc)*100, min_combined = min(combined_perc)*100, max_combined = max(combined_perc)*100)}
