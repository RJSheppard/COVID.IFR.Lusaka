rm(list = ls())
gc()

devtools::install()
devtools::load_all()
library(squire)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(Brobdingnag)

#######################################
##'[Data to run sensitivity analyses]##
#######################################

## Official data for generating initial start date
data <- Off_data <- readRDS(file = "~/Documents/Zambia/covid-mortality-ascertainment/analysis/data/Code-generated-data/00_01_Lusaka_Dist_Deaths_Official.rds")

# Lusaka population
population <- readRDS("~/Documents/Zambia/covid-mortality-ascertainment/analysis/data/Code-generated-data/00_02_Lusaka_Dist_Pop_Str_2020_imp_ests.rds")
# Nyanga contact matrix
baseline_contact_matrix <- as.matrix(readRDS("~/Documents/Zambia/covid-mortality-ascertainment/analysis/data/Code-generated-data/00_11_Nyanga_Mixing_Matrix.rds"))

# Burial registrations and Post-Mortem data
Comb_data <- readRDS("~/Documents/Zambia/covid-mortality-ascertainment/analysis/data/Code-generated-data/00_13_Combined_bur_regs_postmortem_data_complete.rds")
Comb_data <- Comb_data %>% mutate(PosTests = PosTests_Strict)

# Baseline registration estimates
dfj_mcmc_data <- readRDS(file = "~/Documents/Zambia/covid-mortality-ascertainment/analysis/data/Code-generated-data/00_16_03_drj_mcmc_data_new_pop_str_2.rds")

# Population PCR and seroprevalence data
pcr_df <- readRDS("~/Documents/Zambia/covid-mortality-ascertainment/analysis/data/Code-generated-data/00_10_Lancet_Data.rds")$pcr_df
sero_df <- readRDS("~/Documents/Zambia/covid-mortality-ascertainment/analysis/data/Code-generated-data/00_10_Lancet_Data.rds")$sero_df

# probability of hospitalisation and death
probs_hosp_death <- readRDS("~/Documents/Zambia/covid-mortality-ascertainment/analysis/data/Code-generated-data/00_03_IFR_probs_death_hosp.rds")
names(probs_hosp_death) <- paste0("X",1:81)

# Durations until death or survival following infection
Weighted_Durs_Hosp <- readRDS("~/Documents/Zambia/covid-mortality-ascertainment/analysis/data/Code-generated-data/00_04_Weighted_durations_death_survive.rds")
dur_get_ox_survive <- Weighted_Durs_Hosp$Surv_Dur_Weighted
dur_get_ox_die <- Weighted_Durs_Hosp$Death_Dur_Weighted
# Assume that deaths without hospital treatment have half duration of treatment and that 70% die at home
dur_death_default <- dur_get_ox_die*0.3 + (dur_get_ox_die*0.5)*0.7

# PCR prevalence with 100% maximum sensitivity
pcr_det_100 <- readRDS("~/Documents/Zambia/covid-mortality-ascertainment/analysis/data/Code-generated-data/00_15_pcr_det_hall_100.rds")

#########################
###????### Should this be a new package??? I suppose I can create a data folder, just to make sure everything works...
#########################
Llike <- cma:::calc_loglikelihood_10_pois_bin_bin_ag1std_agRR


########################################
##'[Send to cluster, varying severity]##
########################################

didehpc::didehpc_config()
didehpc::web_login()


path <- pkgbuild::build("~/Documents/Zambia/covid-mortality-ascertainment", ".")
src <- conan::conan_sources(path)
packages <- c("reshape2", "Brobdingnag")
ctx_cma_01 <- context::context_save("context_cma_01", packages = packages, package_sources = src, sources = "Covid_Zambia_Scripts/2_drjacoby_MCMC_Full_Model.R")
obj_cma_01 <- didehpc::queue_didehpc(ctx_cma_01)

### Standard sensitivity subset:
IFR_coefficients <- readRDS("~/Documents/Zambia/covid-mortality-ascertainment/analysis/data/Code-generated-data/00_03_IFR_matrix_coefficients_log_scale_new_pop_str_ests.rds") %>%
  mutate(Index = 1:nrow(.))

Select_Runs <- IFR_coefficients %>%
  filter(Slope_x %in% c(0.8,1,1.25) & round(IFR_x,2) %in% c(0.4,0.6,0.8,1,1.25,1.67, 2.5)) %>%  pull(Index)


######################
## Sensitivity 1: remove age group <5
######################
Comb_data_no_age_1 <- Comb_data %>% filter(!Age_gr %in% 1)
dfj_mcmc_data_no_age_1 <- lapply(dfj_mcmc_data, function(x){return(x %>% dplyr::filter(
  Age_gr != 1))})


t_Sensitivity_1_remove_age_1 <- obj_cma_01$lapply(probs_hosp_death[Select_Runs], cma_fit,

                                                        data = data,
                                                        population = population,
                                                        baseline_contact_matrix = baseline_contact_matrix,

                                                        dur_get_ox_survive = dur_get_ox_survive,
                                                        dur_get_ox_die = dur_death_default,

                                                        n_mcmc = 30000, replicates = 100,

                                                        log_likelihood = Llike,
                                                        lld = "",
                                                        Prior_Rt_rw_unif_lim = 1,

                                                        pcr_df = pcr_df,
                                                        sero_df = sero_df,

                                                        combined_data = list(Comb_data = Comb_data_no_age_1,
                                                                             dfj_mcmc_data = dfj_mcmc_data_no_age_1),

                                                        pcr_det = pcr_det_100,
                                                        pcr_det_PM = pcr_det_100,

                                                        frac_reg = 0.9
)

S1_Res <- t_Sensitivity_1_remove_age_1$results()
names(S1_Res) <- paste0("X",1:71)[Select_Runs]
saveRDS(S1_Res, file = "Results/Sensitivity_Analysis_1_Remove_Age_gr_1.rds")


######################
## Sensitivity 2: Unscaled mortality data, and no weeks 4-5
######################
t_Sensitivity_2_no_scaling_no_weeks_4_5 <- obj_cma_01$lapply(probs_hosp_death[Select_Runs], cma_fit,

                                                    data = data,
                                                    population = population,
                                                    baseline_contact_matrix = baseline_contact_matrix,

                                                    dur_get_ox_survive = dur_get_ox_survive,
                                                    dur_get_ox_die = dur_death_default,

                                                    n_mcmc = 30000, replicates = 100,

                                                    log_likelihood = Llike,
                                                    lld = "no scaling no weeks 4-5",
                                                    Prior_Rt_rw_unif_lim = 1,

                                                    pcr_df = pcr_df,
                                                    sero_df = sero_df,

                                                    combined_data = list(Comb_data = Comb_data,
                                                                         dfj_mcmc_data = dfj_mcmc_data),

                                                    pcr_det = pcr_det_100,
                                                    pcr_det_PM = pcr_det_100,

                                                    frac_reg = 0.9
)
S2_Res <- t_Sensitivity_2_no_scaling_no_weeks_4_5$results()
names(S2_Res) <- paste0("X",1:71)[Select_Runs]
saveRDS(S2_Res, file = "Results/Sensitivity_Analysis_2_No_scaling_no_weeks_4_or_5.rds")


######################
## Sensitivity 3: Fit only to burials and prevalence
######################
t_Sensitivity_3_fit_only_to_burials_and_prevalence <- obj_cma_01$lapply(probs_hosp_death[Select_Runs], cma_fit,

                                                               data = data,
                                                               population = population,
                                                               baseline_contact_matrix = baseline_contact_matrix,

                                                               dur_get_ox_survive = dur_get_ox_survive,
                                                               dur_get_ox_die = dur_death_default,

                                                               n_mcmc = 30000, replicates = 100,

                                                               log_likelihood = Llike,
                                                               lld = "no post mortem",
                                                               Prior_Rt_rw_unif_lim = 1,

                                                               pcr_df = pcr_df,
                                                               sero_df = sero_df,

                                                               combined_data = list(Comb_data = Comb_data,
                                                                                    dfj_mcmc_data = dfj_mcmc_data),

                                                               pcr_det = pcr_det_100,
                                                               pcr_det_PM = pcr_det_100,

                                                               frac_reg = 0.9
)
S3_Res <- t_Sensitivity_3_fit_only_to_burials_and_prevalence$results()
names(S3_Res) <- paste0("X",1:71)[Select_Runs]
saveRDS(S3_Res, file = "Results/Sensitivity_Analysis_3_Only_burials_and_pop_prevalence.rds")



######################
## Sensitivity 4a: Change relative rates by 20%: Lower
######################
dfj_mcmc_data_lower_10 <- readRDS(file = "~/Documents/Zambia/covid-mortality-ascertainment/analysis/data/Code-generated-data/00_16_03_drj_mcmc_data_new_pop_str_2_lower_rr_10.rds")
t_Sensitivity_4a_change_relative_rates_by_10pc_lower_failed_run <- obj_cma_01$lapply(
  probs_hosp_death[Select_Runs], cma_fit,

                                                                        data = data,
                                                                        population = population,
                                                                        baseline_contact_matrix = baseline_contact_matrix,

                                                                        dur_get_ox_survive = dur_get_ox_survive,
                                                                        dur_get_ox_die = dur_death_default,

                                                                        n_mcmc = 30000, replicates = 100,

                                                                        log_likelihood = Llike,
                                                                        lld = "",
                                                                        Prior_Rt_rw_unif_lim = 1,

                                                                        pcr_df = pcr_df,
                                                                        sero_df = sero_df,

                                                                        combined_data = list(Comb_data = Comb_data,
                                                                                             dfj_mcmc_data = dfj_mcmc_data_lower_10),

                                                                        pcr_det = pcr_det_100,
                                                                        pcr_det_PM = pcr_det_100,

                                                                        frac_reg = 0.9
)

S4a_Res <- t_Sensitivity_4a_change_relative_rates_by_10pc_lower_failed_run$results()
names(S4a_Res) <- paste0("X",1:71)[Select_Runs]
saveRDS(S4a_Res, file = "Results/Sensitivity_Analysis_4a_Change_age_RR_lower_10.rds")

######################
## Sensitivity 4b: Change relative rates by 20%: Higher
######################
dfj_mcmc_data_higher_10 <- readRDS(file = "~/Documents/Zambia/covid-mortality-ascertainment/analysis/data/Code-generated-data/00_16_03_drj_mcmc_data_new_pop_str_2_higher_rr_10.rds")

t_Sensitivity_4b_change_relative_rates_by_10pc_higher_failed_run <- obj_cma_01$lapply(
  probs_hosp_death[Select_Runs], cma_fit,

  data = data,
  population = population,
  baseline_contact_matrix = baseline_contact_matrix,

  dur_get_ox_survive = dur_get_ox_survive,
  dur_get_ox_die = dur_death_default,

  n_mcmc = 30000, replicates = 100,

  log_likelihood = Llike,
  lld = "",
  Prior_Rt_rw_unif_lim = 1,

  pcr_df = pcr_df,
  sero_df = sero_df,

  combined_data = list(Comb_data = Comb_data,
                       dfj_mcmc_data = dfj_mcmc_data_higher_10),

  pcr_det = pcr_det_100,
  pcr_det_PM = pcr_det_100,

  frac_reg = 0.9
)

S4b_Res <- t_Sensitivity_4b_change_relative_rates_by_10pc_higher_failed_run$results()
names(S4b_Res) <- paste0("X",1:71)[Select_Runs]
saveRDS(S4b_Res, file = "Results/Sensitivity_Analysis_4b_Change_age_RR_higher_10.rds")



######################
## Sensitivity 4c: Change relative rates by 10%: Lower
######################
dfj_mcmc_data_lower_20 <- readRDS(file = "~/Documents/Zambia/covid-mortality-ascertainment/analysis/data/Code-generated-data/00_16_03_drj_mcmc_data_new_pop_str_2_lower_rr_20.rds")

t_Sensitivity_4c_change_relative_rates_by_20pc_lower_failed_run <- obj_cma_01$lapply(
  probs_hosp_death[Select_Runs], cma_fit,

  data = data,
  population = population,
  baseline_contact_matrix = baseline_contact_matrix,

  dur_get_ox_survive = dur_get_ox_survive,
  dur_get_ox_die = dur_death_default,

  n_mcmc = 30000, replicates = 100,

  log_likelihood = Llike,
  lld = "",
  Prior_Rt_rw_unif_lim = 1,

  pcr_df = pcr_df,
  sero_df = sero_df,

  combined_data = list(Comb_data = Comb_data,
                       dfj_mcmc_data = dfj_mcmc_data_lower_20),

  pcr_det = pcr_det_100,
  pcr_det_PM = pcr_det_100,

  frac_reg = 0.9
)

S4c_Res <- t_Sensitivity_4c_change_relative_rates_by_20pc_lower_failed_run$results()
names(S4c_Res) <- paste0("X",1:71)[Select_Runs]
saveRDS(S4c_Res, file = "Results/Sensitivity_Analysis_4c_Change_age_RR_lower_20.rds")




######################
## Sensitivity 4d: Change relative rates by 10%: Higher
######################
dfj_mcmc_data_higher_20 <- readRDS(file = "~/Documents/Zambia/covid-mortality-ascertainment/analysis/data/Code-generated-data/00_16_03_drj_mcmc_data_new_pop_str_2_higher_rr_20.rds")

t_Sensitivity_4d_change_relative_rates_by_20pc_higher_failed_run <- obj_cma_01$lapply(
  probs_hosp_death[Select_Runs], cma_fit,

  data = data,
  population = population,
  baseline_contact_matrix = baseline_contact_matrix,

  dur_get_ox_survive = dur_get_ox_survive,
  dur_get_ox_die = dur_death_default,

  n_mcmc = 30000, replicates = 100,

  log_likelihood = Llike,
  lld = "",
  Prior_Rt_rw_unif_lim = 1,

  pcr_df = pcr_df,
  sero_df = sero_df,

  combined_data = list(Comb_data = Comb_data,
                       dfj_mcmc_data = dfj_mcmc_data_higher_20),

  pcr_det = pcr_det_100,
  pcr_det_PM = pcr_det_100,

  frac_reg = 0.9
)

S4d_Res <- t_Sensitivity_4d_change_relative_rates_by_20pc_higher_failed_run$results()
names(S4d_Res) <- paste0("X",1:71)[Select_Runs]
saveRDS(S4d_Res, file = "Results/Sensitivity_Analysis_4d_Change_age_RR_higher_20.rds")



######################
## Sensitivity 5a: Change durations: duration until death without treatement is 20% of that with treatmentt
######################

dur_death_default_low <- dur_get_ox_die*0.3 + (dur_get_ox_die*0.2)*0.7
t_Sensitivity_5a_change_durations_low <- obj_cma_01$lapply(probs_hosp_death[Select_Runs], cma_fit,

                                                                    data = data,
                                                                    population = population,
                                                                    baseline_contact_matrix = baseline_contact_matrix,

                                                                    dur_get_ox_survive = dur_get_ox_survive,
                                                                    dur_get_ox_die = dur_death_default_low,

                                                                    n_mcmc = 30000, replicates = 100,

                                                                    log_likelihood = Llike,
                                                                    lld = "",
                                                                    Prior_Rt_rw_unif_lim = 1,

                                                                    pcr_df = pcr_df,
                                                                    sero_df = sero_df,

                                                                    combined_data = list(Comb_data = Comb_data,
                                                                                         dfj_mcmc_data = dfj_mcmc_data),

                                                                    pcr_det = pcr_det_100,
                                                                    pcr_det_PM = pcr_det_100,

                                                                    frac_reg = 0.9
)
S5aRes <- t_Sensitivity_5a_change_durations_low$results()
names(S5aRes) <- paste0("X",1:71)[Select_Runs]
saveRDS(S5aRes, file = "Results/Sensitivity_Analysis_5a_Change_duration_to_death_lower.rds")

######################
## Sensitivity 5b: Change durations: duration until death without treatement is equal to that with treatment
######################

dur_death_default_high <- dur_get_ox_die*0.3 + (dur_get_ox_die*1)*0.7
t_Sensitivity_5b_change_durations_high <- obj_cma_01$lapply(probs_hosp_death[Select_Runs], cma_fit,

                                                                   data = data,
                                                                   population = population,
                                                                   baseline_contact_matrix = baseline_contact_matrix,

                                                                   dur_get_ox_survive = dur_get_ox_survive,
                                                                   dur_get_ox_die = dur_death_default_high,

                                                                   n_mcmc = 30000, replicates = 100,

                                                                   log_likelihood = Llike,
                                                                   lld = "",
                                                                   Prior_Rt_rw_unif_lim = 1,

                                                                   pcr_df = pcr_df,
                                                                   sero_df = sero_df,

                                                                   combined_data = list(Comb_data = Comb_data,
                                                                                        dfj_mcmc_data = dfj_mcmc_data),

                                                                   pcr_det = pcr_det_100,
                                                                   pcr_det_PM = pcr_det_100,

                                                                   frac_reg = 0.9
)

S5bRes <- t_Sensitivity_5b_change_durations_low$results()
names(S5bRes) <- paste0("X",1:71)[Select_Runs]
saveRDS(S5bRes, file = "Results/Sensitivity_Analysis_5b_Change_duration_to_death_higher.rds")
