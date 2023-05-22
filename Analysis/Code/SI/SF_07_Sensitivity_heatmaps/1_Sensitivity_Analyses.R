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
data <- Off_data <- readRDS(file = "Analysis/Data/derived_data/04_Lusaka_Dist_Deaths_Official.rds")

# Lusaka population
population <- readRDS("Analysis/Data/derived_data/08_Lusaka_Dist_Pop_Str_2020_imp_ests.rds")
# Nyanga contact matrix
baseline_contact_matrix <- as.matrix(readRDS("Analysis/Data/derived_data/11_Nyanga_Mixing_Matrix.rds"))

# Burial registrations and Post-Mortem data
Comb_data <- readRDS("Analysis/Data/derived_data/12_Combined_bur_regs_postmortem_data_complete.rds")
Comb_data <- Comb_data %>% mutate(PosTests = PosTests_Strict)

# Baseline registration estimates
dfj_mcmc_data <- readRDS(file = "Analysis/Data/derived_data/13_drj_mcmc_data_new_pop_str_2.rds")

# Population PCR and seroprevalence data
pcr_df <- readRDS("Analysis/Data/derived_data/14_Lancet_Data.rds")$pcr_df
sero_df <- readRDS("Analysis/Data/derived_data/14_Lancet_Data.rds")$sero_df

# probability of hospitalisation and death
probs_hosp_death <- readRDS("Analysis/Data/derived_data/15_IFR_probs_death_hosp.rds")
names(probs_hosp_death) <- paste0("X",1:81)

# Durations until death or survival following infection
Weighted_Durs_Hosp <- readRDS("Analysis/Data/derived_data/16_Weighted_durations_death_survive.rds")
dur_get_ox_survive <- Weighted_Durs_Hosp$Surv_Dur_Weighted
dur_get_ox_die <- Weighted_Durs_Hosp$Death_Dur_Weighted
# Assume that deaths without hospital treatment have half duration of treatment and that 70% die at home
dur_death_default <- dur_get_ox_die*0.3 + (dur_get_ox_die*0.5)*0.7

# PCR prevalence with 100% maximum sensitivity
pcr_det_100 <- readRDS("Analysis/Data/derived_data/17_pcr_det_hall_100.rds")

Llike <- cma:::calc_loglikelihood_10_pois_bin_bin_ag1std_agRR


########################################
##'[Send to cluster, varying severity]##
########################################

didehpc::didehpc_config()
didehpc::web_login()


path <- pkgbuild::build("COVID_IFR_Lusaka", ".")
src <- conan::conan_sources(path)
packages <- c("reshape2", "Brobdingnag")
ctx_cma_01 <- context::context_save("context_cma_01", packages = packages, package_sources = src, sources = "Analysis/Code/Code_Functions/drjacoby_MCMC_Full_Model.R")
obj_cma_01 <- didehpc::queue_didehpc(ctx_cma_01)

### Standard sensitivity subset:
IFR_coefficients <- readRDS("Analysis/Data/derived_data/18_IFR_matrix_coefficients_log_scale_new_pop_str_ests.rds") %>%
  mutate(Index = 1:nrow(.))

Select_Runs <- IFR_coefficients %>%
  filter(Slope_x %in% c(0.8,1,1.25) & round(IFR_x,2) %in% c(0.4,0.6,0.8,1,1.25,1.67, 2.5)) %>%  pull(Index)


######################
## Sensitivity A: remove age group <5
######################
Comb_data_no_age_1 <- Comb_data %>% filter(!Age_gr %in% 1)
dfj_mcmc_data_no_age_1 <- lapply(dfj_mcmc_data, function(x){return(x %>% dplyr::filter(
  Age_gr != 1))})


t_Sensitivity_A_remove_age_1 <- obj_cma_01$lapply(probs_hosp_death[Select_Runs], cma_fit,

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

SA_Res <- t_Sensitivity_A_remove_age_1$results()
names(SA_Res) <- paste0("X",1:71)[Select_Runs]
saveRDS(SA_Res, file = "Analysis/Results/SI_Sensitivity_Analysis_A_Remove_Age_gr_1.rds")


######################
## Sensitivity B: Unscaled mortality data, and no weeks 4-5
######################
t_Sensitivity_B_no_scaling_no_weeks_4_5 <- obj_cma_01$lapply(probs_hosp_death[Select_Runs], cma_fit,

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
SB_Res <- t_Sensitivity_B_no_scaling_no_weeks_4_5$results()
names(SB_Res) <- paste0("X",1:71)[Select_Runs]
saveRDS(SB_Res, file = "Analysis/Results/SI_Sensitivity_Analysis_B_No_scaling_no_weeks_4_or_5.rds")


######################
## Sensitivity C: Change durations: duration until death without treatement is 20% of that with treatmentt
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
SCRes <- t_Sensitivity_C_change_durations_low$results()
names(SCRes) <- paste0("X",1:71)[Select_Runs]
saveRDS(SCRes, file = "Analysis/Results/SI_Sensitivity_Analysis_C_Change_duration_to_death_lower.rds")

######################
## Sensitivity D: Change durations: duration until death without treatement is equal to that with treatment
######################

dur_death_default_high <- dur_get_ox_die*0.3 + (dur_get_ox_die*1)*0.7
t_Sensitivity_D_change_durations_high <- obj_cma_01$lapply(probs_hosp_death[Select_Runs], cma_fit,

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

SDRes <- t_Sensitivity_D_change_durations_low$results()
names(SDRes) <- paste0("X",1:71)[Select_Runs]
saveRDS(SDRes, file = "Analysis/Results/SI_Sensitivity_Analysis_D_Change_duration_to_death_higher.rds")



######################
## Sensitivity E: Change relative rates by 10%: Lower
######################
dfj_mcmc_data_lower_10 <- readRDS(file = "Analysis/Data/derived_data/20_drj_mcmc_data_new_pop_str_2_lower_rr_10.rds")
t_Sensitivity_E_change_relative_rates_by_10pc_lower_failed_run <- obj_cma_01$lapply(
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

SE_Res <- t_Sensitivity_E_change_relative_rates_by_10pc_lower_failed_run$results()
names(SE_Res) <- paste0("X",1:71)[Select_Runs]
saveRDS(SE_Res, file = "Analysis/Results/SI_Sensitivity_Analysis_E_Change_age_RR_lower_10.rds")

######################
## Sensitivity F: Change relative rates by 10%: Higher
######################
dfj_mcmc_data_higher_10 <- readRDS(file = "Analysis/Data/derived_data/21_drj_mcmc_data_new_pop_str_2_higher_rr_10.rds")

t_Sensitivity_F_change_relative_rates_by_10pc_higher_failed_run <- obj_cma_01$lapply(
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

SF_Res <- t_Sensitivity_F_change_relative_rates_by_10pc_higher_failed_run$results()
names(SF_Res) <- paste0("X",1:71)[Select_Runs]
saveRDS(SF_Res, file = "Analysis/Results/SI_Sensitivity_Analysis_F_Change_age_RR_higher_10.rds")



######################
## Sensitivity G: Change relative rates by 20%: Lower
######################
dfj_mcmc_data_lower_20 <- readRDS(file = "Analysis/Data/derived_data/22_drj_mcmc_data_new_pop_str_2_lower_rr_20.rds")

t_Sensitivity_G_change_relative_rates_by_20pc_lower_failed_run <- obj_cma_01$lapply(
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

SG_Res <- t_Sensitivity_G_change_relative_rates_by_20pc_lower_failed_run$results()
names(SG_Res) <- paste0("X",1:71)[Select_Runs]
saveRDS(SG_Res, file = "Analysis/Results/SI_Sensitivity_Analysis_G_Change_age_RR_lower_20.rds")




######################
## Sensitivity H: Change relative rates by 20%: Higher
######################
dfj_mcmc_data_higher_20 <- readRDS(file = "Analysis/Data/derived_data/23_drj_mcmc_data_new_pop_str_2_higher_rr_20.rds")

t_Sensitivity_H_change_relative_rates_by_20pc_higher_failed_run <- obj_cma_01$lapply(
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

SH_Res <- t_Sensitivity_H_change_relative_rates_by_20pc_higher_failed_run$results()
names(SH_Res) <- paste0("X",1:71)[Select_Runs]
saveRDS(SH_Res, file = "Analysis/Results/SI_Sensitivity_Analysis_H_Change_age_RR_higher_20.rds")
