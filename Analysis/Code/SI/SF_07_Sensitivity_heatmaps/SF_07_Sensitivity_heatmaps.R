## Get data inputs ready
library(tidyr)
library(dplyr)

# Lusaka population
population <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Population_structure_Lusaka_2020_CDC.rds")

# Nyanga contact matrix
baseline_contact_matrix <- as.matrix(readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Mixing_Matrix_Nyanga.rds"))

# Burial registrations and Post-Mortem data
Comb_data <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/combined_squire_data.rds") %>%
  mutate(PosTests = PosTests_Strict) |>
  arrange(Age_gr, Week_gr)

# Baseline registration estimates
dfj_mcmc_data <- readRDS(file = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/Baseline_Mortality_MCMC_formatted_for_squire.rds")
dfj_mcmc_data_lower_10 <- readRDS(file = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/Baseline_Mortality_MCMC_formatted_for_squire_l10.rds")
dfj_mcmc_data_lower_20 <- readRDS(file = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/Baseline_Mortality_MCMC_formatted_for_squire_l20.rds")
dfj_mcmc_data_higher_10 <- readRDS(file = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/Baseline_Mortality_MCMC_formatted_for_squire_h10.rds")
dfj_mcmc_data_higher_20 <- readRDS(file = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/Baseline_Mortality_MCMC_formatted_for_squire_h20.rds")

# Population PCR and seroprevalence data
pcr_df <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Population_prevalence.rds")$pcr_df
sero_df <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Population_prevalence.rds")$sero_df

# probability of hospitalisation and death
probs_hosp_death <- readRDS("Analysis/Data/derived_data/IFR_phosp_pdeath.rds")
names(probs_hosp_death) <- paste0("X",1:81)

# Durations until death or survival following infection
Weighted_Durs_Hosp <- readRDS("Analysis/Data/derived_data/Weighted_durations_death_survive.rds")
dur_get_ox_survive <- Weighted_Durs_Hosp$Surv_Dur_Weighted
dur_get_ox_die <- Weighted_Durs_Hosp$Death_Dur_Weighted
# Assume that deaths without hospital treatment have half duration of treatment and that 70% die at home
dur_death_default <- dur_get_ox_die*0.3 + (dur_get_ox_die*0.5)*0.7
dur_death_default_low <- dur_get_ox_die*0.3 + (dur_get_ox_die*0.2)*0.7
dur_death_default_high <- dur_get_ox_die*0.3 + (dur_get_ox_die*1)*0.7

# PCR prevalence with 100% maximum sensitivity
pcr_det_100 <- readRDS("Analysis/Data/derived_data/pcr_det_hall_100.rds")

Llike <- calc_loglikelihood_11_fully_vectorised

# Use Official deaths to generate window of start dates for pandemic
data <- readRDS(file = "Analysis/Data/derived_data/Lusaka_Province_Dashboard.rds") %>%
  rename(date = "Dates")

########################
########################
## Sensitivity Analyses
########################
########################

### Sensitivity Analyses run using Imperial College DIDE HPC cluster ###
getwd()
didehpc::didehpc_config()
didehpc::web_login()
batch <- "2023_March_Reviewer_Batch"

# options(didehpc.cluster = "fi--dideclusthn")
path <- pkgbuild::build("COVID.IFR.Lusaka", ".")
src <- conan::conan_sources(path)
packages <- c("reshape2", "Brobdingnag")
ctx_cma_01 <- context::context_save("context_cma_01", packages = packages, package_sources = src, sources = "R/squire_MCMC_Full_Model.R")
obj_cma_01 <- didehpc::queue_didehpc(ctx_cma_01)

IFR_coefficients <- readRDS("Analysis/Data/derived_data/IFR_matrix.rds") %>%
  mutate(Index = 1:nrow(.))

Select_Runs <- IFR_coefficients %>%
  filter(Slope_x %in% c(0.8,1,1.25) & round(IFR_x,2) %in% c(0.4,0.6,0.8,1,1.25,1.67, 2.5)) %>%  pull(Index)

######################
## Sensitivity 1: remove week 5 from PM fit
######################
t_S1 <- obj_cma_shared$lapply(probs_hosp_death[Select_Runs], cma_fit,
                              data = data,
                              population = population,
                              baseline_contact_matrix = baseline_contact_matrix,

                              dur_get_ox_survive = dur_get_ox_survive,
                              dur_get_ox_die = dur_death_default,

                              n_mcmc = 30000, replicates = 100,

                              log_likelihood = Llike,
                              lld = "no week 5 PM",
                              Prior_Rt_rw_unif_lim = 1,

                              pcr_df = pcr_df,
                              sero_df = sero_df,

                              combined_data = Comb_data,
                              drj_mcmc = dfj_mcmc_data,

                              pcr_det = pcr_det_100,

                              frac_reg = 0.9
)

# saveRDS(t_S1$results(), file = paste0("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/S1_remove_week_5_finished.rds"))

######################
## Sensitivity 2: remove age group <5
######################
Select_Runs2 <- IFR_coefficients %>%
  filter(round(Slope_x,2) %in% c(0.8,1,1.25) & round(IFR_x,2) %in% c(0.4,0.6,0.8,1,1.25,1.67,2.5) |
           round(Slope_x,2) %in% c(1.67) & round(IFR_x,2) %in% c(0.4,0.6,0.8,1)) %>%  pull(Index)

t_S2 <- obj_cma_01$lapply(probs_hosp_death[Select_Runs2], cma_fit,
                          data = data,
                          population = population,
                          baseline_contact_matrix = baseline_contact_matrix,

                          dur_get_ox_survive = dur_get_ox_survive,
                          dur_get_ox_die = dur_death_default,

                          n_mcmc = 30000, replicates = 100,

                          log_likelihood = Llike,
                          lld = "no under 5s",
                          Prior_Rt_rw_unif_lim = 1,

                          pcr_df = pcr_df,
                          sero_df = sero_df,

                          combined_data = Comb_data,
                          drj_mcmc = dfj_mcmc_data,

                          pcr_det = pcr_det_100,

                          frac_reg = 0.9
)

saveRDS(t_S2$results(), file = paste0("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/S2_Remove_Age_gr_1.rds"))


######################
## Sensitivity 3: Unscaled mortality data, and no weeks 4-5
######################

Select_Runs3 <- IFR_coefficients %>%
  filter(round(Slope_x,2) %in% c(0.8,1,1.25) & round(IFR_x,2) %in% c(0.4,0.6,0.8,1,1.25,1.67,2.5) |
           round(Slope_x,2) %in% c(1.67) & round(IFR_x,2) %in% c(0.4,0.6,0.8,1)) %>%  pull(Index)


t_S3 <- obj_cma_shared$lapply(probs_hosp_death[Select_Runs3], cma_fit,
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

                              combined_data = Comb_data,
                              drj_mcmc = dfj_mcmc_data,

                              pcr_det = pcr_det_100,

                              frac_reg = 0.9
)

saveRDS(t_S3$results(), file = paste0("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/S3_No_scaling_no_weeks_4_or_5_finished.rds"))

######################
## Sensitivity 4: Duration increased
######################

Select_Runs4 <- IFR_coefficients %>%
  filter(Slope_x %in% c(0.8,1,1.25) & round(IFR_x,2) %in% c(0.4,0.6,0.8,1,1.25,1.67, 2.5) |
           Slope_x %in% c(1,1.25) & round(IFR_x,2) %in% c(0.2) |
           round(Slope_x,2) %in% c(1.67) & round(IFR_x,2) %in% c(0.2,0.4,0.6,0.8,1)) %>%  pull(Index)


t_S4 <- obj_cma_shared$lapply(probs_hosp_death[Select_Runs4], cma_fit,
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

                              combined_data = Comb_data,
                              drj_mcmc = dfj_mcmc_data,

                              pcr_det = pcr_det_100,

                              frac_reg = 0.9
)

saveRDS(t_S4$results(), file = paste0("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/S4_Increase_duration_to_death.rds"))


######################
## Sensitivity 5: Durations decreased
######################
t_S5 <- obj_cma_01$lapply(probs_hosp_death[c(31,48,53)], cma_fit,
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

                          combined_data = Comb_data,
                          drj_mcmc = dfj_mcmc_data,

                          pcr_det = pcr_det_100,

                          frac_reg = 0.9
)

saveRDS(t_S5$results(), file = paste0("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/S5_Decrease_duration_to_death_batch_3.rds"))


######################
## Sensitivity 6: Change relative rates by 10%: Lower
######################
t_S6 <- obj_cma_01$lapply(probs_hosp_death[c(54)], cma_fit,
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

                          combined_data = Comb_data,
                          drj_mcmc = dfj_mcmc_data_lower_10,

                          pcr_det = pcr_det_100,

                          frac_reg = 0.9
)

saveRDS(t_S6$results(partial = T), file = paste0("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/S6_Change_age_RR_10_lower_batch_3.rds"))

######################
## Sensitivity 7: Change relative rates by 10%: Higher
######################
t_S7 <- obj_cma_01$lapply(probs_hosp_death[Select_Runs], cma_fit,
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

                          combined_data = Comb_data,
                          drj_mcmc = dfj_mcmc_data_higher_10,

                          pcr_det = pcr_det_100,

                          frac_reg = 0.9
)

saveRDS(t_S7$results(), file = paste0("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/S7_Change_age_RR_10_higher_better_batch_2.rds"))


######################
## Sensitivity 8: Change relative rates by 20%: Lower
######################

Select_Runs8 <- IFR_coefficients %>%
  filter(Slope_x %in% c(0.8,1,1.25) & round(IFR_x,2) %in% c(0.4,0.6,0.8,1,1.25,1.67, 2.5,5)) %>%  pull(Index)

t_S8 <- obj_cma_01$lapply(probs_hosp_death[Select_Runs8], cma_fit,
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

                          combined_data = Comb_data,
                          drj_mcmc = dfj_mcmc_data_lower_20,

                          pcr_det = pcr_det_100,

                          frac_reg = 0.9
)

saveRDS(t_S8$results(), file = paste0("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/S8_Change_age_RR_20_lower_batch_2.rds"))


######################
## Sensitivity 9: Change relative rates by 20%: Higher
######################
Select_Runs9 <- IFR_coefficients %>%
  filter(Slope_x %in% c(0.8,1,1.25) & round(IFR_x,2) %in% c(0.4,0.6,0.8,1,1.25,1.67, 2.5) |
           Slope_x %in% c(1,1.25) & round(IFR_x,2) %in% c(0.2) |
           round(Slope_x,2) %in% c(1.67) & round(IFR_x,2) %in% c(0.2,0.4,0.6,0.8,1)) %>%  pull(Index)

t_S9 <- obj_cma_shared$lapply(probs_hosp_death[Select_Runs9], cma_fit,
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

                              combined_data = Comb_data,
                              drj_mcmc = dfj_mcmc_data_higher_20,

                              pcr_det = pcr_det_100,

                              frac_reg = 0.9
)

saveRDS(t_S9$results(), file = paste0("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/S9_Change_age_RR_20_higher_better_batch_2.rds"))



######################
## Sensitivity 10: 100% registration
######################
t_S10 <- obj_cma_01$lapply(probs_hosp_death[c(38,48)], cma_fit,
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

                           combined_data = Comb_data,
                           drj_mcmc = dfj_mcmc_data,

                           pcr_det = pcr_det_100,

                           frac_reg = 1
)

saveRDS(t_S10$results(), file = paste0("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/S10_100pc_Registration_batch_2.rds"))


######################
## Sensitivity 11: 80% registration
######################
t_S11 <- obj_cma_shared$lapply(probs_hosp_death[Select_Runs], cma_fit,
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

                               combined_data = Comb_data,
                               drj_mcmc = dfj_mcmc_data,

                               pcr_det = pcr_det_100,

                               frac_reg = 0.8
)

saveRDS(t_S11$results(), file = paste0("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/S11_80pc_registration.rds"))


################################
################################
## Plot Sensitivity Results ####
################################
################################

## Load the data
devtools::load_all()
library(ggplot2)
library(Brobdingnag)

IFR_mat <- readRDS("Analysis/Data/derived_data/IFR_matrix.rds")

# Default
batch <- "2023_March_Reviewer_Batch"

# Res_default_a <- readRDS(paste0("../Bonus Files/",batch,"/01_Full_fit_round_1.rds"))
# Res_default_b <- readRDS(paste0("../Bonus Files/",batch,"/01_Full_fit_round_2.rds"))
# Res_default_a <- c(Res_default_a, Res_default_b)
# Res_default_a <- Res_default_a[order(as.numeric(gsub("X","",names(Res_default_a))))]
# Res_default_a <- Res_default_a[unlist(lapply(Res_default_a, function(x){any(class(x)=="squire_simulation")}))]
# Res_default_a <- Res_default_a[c(29:35,38:44,47:53)]
# HDefault <- Heatmap_Post(Plot_Post(Res_default_a, IFR_mat))
# rm(list = c("Res_default_a","Res_default_b"))
Res_default <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/Full_Set_Default_Setting.rds/Full_Set_Default_Setting.rds")
HDefault <- Heatmap_Post(Plot_Post(Res_default, IFR_mat))
rm(list = c("Res_default"))
gc()

### Heatmap A
Res_Sens_A <- readRDS(paste0("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/S1_remove_week_5_finished.rds"))
HA <- Heatmap_Post(Plot_Post(Res_Sens_A, IFR_mat))
rm(list = c("Res_Sens_A"))
gc()

### Heatmap B
Res_Sens_B <- readRDS(paste0("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/S2_Remove_Age_gr_1.rds"))
HB <- Heatmap_Post(Plot_Post(Res_Sens_B[-c(22:25)], IFR_mat))
rm(list = c("Res_Sens_B"))
gc()

### Heatmap C
Res_Sens_C <- readRDS(paste0("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/S3_No_scaling_no_weeks_4_or_5_finished.rds"))
HC <- Heatmap_Post(Plot_Post(Res_Sens_C[-c(22:25)], IFR_mat))
rm(list = c("Res_Sens_C"))
gc()

### Heatmap D
Res_Sens_D <- readRDS(paste0("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/S4_Increase_duration_to_death.rds"))
HD <- Heatmap_Post(Plot_Post(Res_Sens_D[-c(8,16,24:28)], IFR_mat))
rm(list = c("Res_Sens_D"))
gc()

### Heatmap E (5c)
Res_Sens_E <- readRDS(paste0("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/S5_Decrease_duration_to_death_batch_1.rds"))
Res_Sens_E <- Res_Sens_E[unlist(lapply(Res_Sens_E, function(x){any(class(x)=="squire_simulation")}))]
Res_Sens_E_2 <- readRDS(paste0("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/S5_Decrease_duration_to_death_batch_2.rds"))
Res_Sens_E_3 <- readRDS(paste0("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/S5_Decrease_duration_to_death_batch_3.rds"))
Res_Sens_E <- c(Res_Sens_E,Res_Sens_E_2,Res_Sens_E_3)
Res_Sens_E <- Res_Sens_E[order(names(Res_Sens_E))]
HE <- Heatmap_Post(Plot_Post(Res_Sens_E, IFR_mat))
rm(list = c("Res_Sens_E","Res_Sens_E_2","Res_Sens_E_3"))
gc()

### Heatmap F (5d)
Res_Sens_F <- readRDS(paste0("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/S6_Change_age_RR_10_lower.rds"))
Res_Sens_F_2 <- readRDS(paste0("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/S6_Change_age_RR_10_lower_batch_2.rds"))
Res_Sens_F_3 <- readRDS(paste0("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/S6_Change_age_RR_10_lower_batch_3.rds"))
Res_Sens_F_2 <- Res_Sens_F_2[unlist(lapply(Res_Sens_F_2, function(x){any(class(x)=="squire_simulation")}))]
Res_Sens_F <- c(Res_Sens_F,Res_Sens_F_2,Res_Sens_F_3)
Res_Sens_F <- Res_Sens_F[order(names(Res_Sens_F))]
HF <- Heatmap_Post(Plot_Post(Res_Sens_F, IFR_mat))
rm(list = c("Res_Sens_F","Res_Sens_F_2","Res_Sens_F_3"))
gc()

### Heatmap G
Res_Sens_G <- readRDS(paste0("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/S7_Change_age_RR_10_batch_1.rds"))
Res_Sens_G_2 <- readRDS(paste0("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/S7_Change_age_RR_10_batch_2.rds"))
Res_Sens_G <- c(Res_Sens_G,Res_Sens_G_2)
Res_Sens_G <- Res_Sens_G[order(names(Res_Sens_G))]
HG <- Heatmap_Post(Plot_Post(Res_Sens_G, IFR_mat))
rm(list = c("Res_Sens_G","Res_Sens_G_2"))
gc()

### Heatmap H
Res_Sens_H <- readRDS(paste0("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/S8_Change_age_RR_20_lower.rds"))
Res_Sens_H_2 <- readRDS(paste0("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/S8_Change_age_RR_20_lower_batch_2.rds"))
Res_Sens_H <- c(Res_Sens_H,Res_Sens_H_2)
Res_Sens_H <- Res_Sens_H[order(names(Res_Sens_H))]
HH <- Heatmap_Post(Plot_Post(Res_Sens_H, IFR_mat))
rm(list = c("Res_Sens_H","Res_Sens_H_2"))
gc()

### Heatmap I
Res_Sens_I <- readRDS(paste0("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/S9_Change_age_RR_20_batch_1.rds"))
Res_Sens_I_2 <- readRDS(paste0("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/S9_Change_age_RR_20_batch_2.rds"))
Res_Sens_I <- c(Res_Sens_I,Res_Sens_I_2)
Res_Sens_I <- Res_Sens_I[order(names(Res_Sens_I))]
HI <- Heatmap_Post(Plot_Post(Res_Sens_I[-c(15,23,27)], IFR_mat))
rm(list = c("Res_Sens_I","Res_Sens_I_2"))
gc()

### Heatmap J
Res_Sens_J <- readRDS(paste0("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/S10_100pc_Registration.rds"))
HJ <- Heatmap_Post(Plot_Post(Res_Sens_J, IFR_mat))
rm(list = c("Res_Sens_J"))
gc()

### Heatmap K
Res_Sens_K <- readRDS(paste0("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/S11_80pc_registration.rds"))
HK <- Heatmap_Post(Plot_Post(Res_Sens_K, IFR_mat))
rm(list = c("Res_Sens_K"))
gc()

update_geom_defaults("text", list(size = 6/.pt))
Full_Heatmap_Sensitivities <- cowplot::plot_grid(cowplot::plot_grid(
  HDefault + ggtitle("a: Default") + theme(axis.title.x = element_blank(), axis.title.y = element_blank()),
  HA + ggtitle("b: Rm PM week 5 ") +theme(axis.title.x = element_blank(), axis.title.y = element_blank()),
  HB + ggtitle("c: Rm PM and BR 0-4 age group")+theme(axis.title.x = element_blank()),
  HC + ggtitle("d: Rm scaling, Rm BR weeks 4-5")+theme(axis.title.x = element_blank(), axis.title.y = element_blank()),
  HE + ggtitle("e: Dec. duration until death to 4.41 days")+theme(axis.title.x = element_blank(),axis.title.y = element_blank()),
  HD + ggtitle("f: Inc. duration until death to 10.01 days")+theme(axis.title.x = element_blank(), axis.title.y = element_blank()),
  HF + ggtitle("g: Dec. baseline BR by 10%")+theme(axis.title.x = element_blank(), axis.title.y = element_blank()),
  HG + ggtitle("h: Inc. baseline BR by 10%")+theme(axis.title.x = element_blank(), axis.title.y = element_blank()),
  HH + ggtitle("i: Dec. baseline BR by 20%")+theme(axis.title.x = element_blank()),
  HI + ggtitle("j: Inc. baseline BR  by 20")+theme(axis.title.x = element_blank(), axis.title.y = element_blank()),
  HK + ggtitle("k: Dec. BR mortality capture to 80%")+ theme(axis.title.y = element_blank()),
  HJ + ggtitle("l: Inc. BR mortality capture to 100%") + theme(axis.title.y = element_blank()),
  ncol = 2, align = "hv"),
  ggplot() +
    theme_void(base_size = 7, base_family = "Helvetica") +
    annotate("text",x=6, y = 1, label = "Rm: Remove, PM: Post-mortem, BR: Burial registration",size = 7/.pt, family = "Helvetica"),
  nrow=2, rel_heights = c(1, 0.05))

ggsave("Analysis/Figures/Supplementary_Figures/SFigure_07_Sensitivity_Heatmaps.pdf",
       Full_Heatmap_Sensitivities,
       width = 145, height = 175, units = "mm")

tiff("Analysis/Figures/Supplementary_Figures/SFigure_07_Sensitivity_Heatmaps.tiff", width = 145, height = 175, units = "mm", res = 300)
Full_Heatmap_Sensitivities
dev.off()
