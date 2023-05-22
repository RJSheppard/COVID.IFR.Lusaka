rm(list = ls())
gc()

devtools::install()
devtools::load_all()
# devtools::install_github("RJSheppard/squire")
library(squire)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(Brobdingnag)

##############################################
##'[Model Fitting Under Default Assumptions]##
##############################################
## Load data
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
# dfj_mcmc_data <- turn_list_into_array(readRDS(file = "~/Documents/Zambia/covid-mortality-ascertainment/analysis/data/Code-generated-data/00_16_03_drj_mcmc_data_new_pop_str_2_2023_March.rds"))


# Population PCR and seroprevalence data
pcr_df <- readRDS("Analysis/Data/derived_data/Population_prevalence.rds")$pcr_df
sero_df <- readRDS("Analysis/Data/derived_data/Population_prevalence.rds")$sero_df

# probability of hospitalisation and death
probs_hosp_death <- readRDS("Analysis/Data/derived_data/IFR_phosp_pdeath.rds")
names(probs_hosp_death) <- paste0("X",1:81)

# Durations until death or survival following infection
Weighted_Durs_Hosp <- readRDS("Analysis/Data/derived_data/Weighted_durations_death_survive.rds")
dur_get_ox_survive <- Weighted_Durs_Hosp$Surv_Dur_Weighted
dur_get_ox_die <- Weighted_Durs_Hosp$Death_Dur_Weighted
# Assume that deaths without hospital treatment have half duration of treatment and that 70% die at home
dur_death_default <- dur_get_ox_die*0.3 + (dur_get_ox_die*0.5)*0.7

# PCR prevalence with 100% maximum sensitivity
pcr_det_100 <- readRDS("Analysis/Data/derived_data/pcr_det_hall_100.rds")

Llike <- calc_loglikelihood_11_fully_vectorised

# Use Official deaths to generate window of start dates for pandemic
data <- readRDS(file = "Analysis/Data/derived_data/Lusaka_Province_Dashboard.rds") %>%
  rename(date = "Dates")


##################################################
Model_Fit <- fit_spline_rt(
  country = "Zambia",
  data = data,

  population = population,
  baseline_contact_matrix = baseline_contact_matrix,

  dur_get_ox_survive = dur_get_ox_survive,
  dur_get_ox_die = dur_death_default,


  pcr_df = pcr_df,
  sero_df = sero_df,

  combined_data = Comb_data,
  drj_mcmc = dfj_mcmc_data,
  frac_reg = 0.9,

  pcr_det = pcr_det_100,

  log_likelihood = Llike,
  lld = "",

  prob_hosp = probs_hosp_death$X41$phi_1,
  prob_non_severe_death_treatment = probs_hosp_death$X41$phi_2,
  prob_severe = rep(0,17),
  hosp_beds = 1e10,
  icu_beds = 1e10,
  dur_R = Inf,

  n_mcmc = 30000,
  replicates = 100,
  rw_duration = 14,
  Prior_Rt_rw_unif_lim = 1
)

# saveRDS(object = Model_Fit, file = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/Model_Fit_Default_Assumptions.rds")

###########################
##'[Plot the for 1x1 fit]##
###########################
Model_Fit <- readRDS(file = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/Model_Fit_Default_Assumptions.rds")

IFR_coefficients <- readRDS("Analysis/Data/derived_data/IFR_matrix.rds") %>%
  mutate(Index = 1:nrow(.))

# library(data.table)
library(reshape2)
devtools::load_all()
Model_Fit_1x1_plot <- Diagnostic_Plot(1, fit_Model_list = list(Model_Fit), IFRvals = IFR_coefficients[41,], Return_likelihoods_only = F)
Model_Fit_1x1_plot$Poisson_Figure_weeks$layers[[1]]$aes_params$size <- 1
Model_Fit_1x1_plot$Poisson_Figure_age$layers[[1]]$aes_params$size <- 1
Model_Fit_1x1_plot$Week_prev_plot$layers[[1]]$aes_params$size <- 1
Model_Fit_1x1_plot$Age_prev_plot$layers[[1]]$aes_params$size <- 1
Model_Fit_1x1_plot$PCR_sero_prev_plot$layers[[1]]$aes_params$size <- 1
Model_Fit_1x1_plot$PCR_sero_prev_plot$layers[[4]]$aes_params$size <- 1

Model_Fit_1x1_plot$Poisson_Figure_weeks$layers[[2]]$aes_params$linewidth <- 0.4
Model_Fit_1x1_plot$Poisson_Figure_age$layers[[2]]$aes_params$linewidth <- 0.4
Model_Fit_1x1_plot$Week_prev_plot$layers[[3]]$aes_params$linewidth <- 0.4
Model_Fit_1x1_plot$Week_prev_plot$layers[[5]]$aes_params$linewidth <- 0.4
Model_Fit_1x1_plot$Week_prev_plot$layers[[7]]$aes_params$linewidth <- 0.4

Model_Fit_1x1_plot$Age_prev_plot$layers[[3]]$aes_params$linewidth <- 0.4
Model_Fit_1x1_plot$Age_prev_plot$layers[[5]]$aes_params$linewidth <- 0.4
Model_Fit_1x1_plot$Age_prev_plot$layers[[7]]$aes_params$linewidth <- 0.4



Model_Fit_1x1_plot$PCR_sero_prev_plot$layers[[1]]$aes_params$linewidth <- 0.4
Model_Fit_1x1_plot$PCR_sero_prev_plot$layers[[2]]$aes_params$linewidth <- 0.4
Model_Fit_1x1_plot$PCR_sero_prev_plot$layers[[3]]$aes_params$linewidth <- 0.4
Model_Fit_1x1_plot$PCR_sero_prev_plot$layers[[5]]$aes_params$linewidth <- 0.4
Model_Fit_1x1_plot$PCR_sero_prev_plot$layers[[6]]$aes_params$linewidth <- 0.4
Model_Fit_1x1_plot$PCR_sero_prev_plot$layers[[7]]$aes_params$linewidth <- 0.4
Model_Fit_1x1_plot$PCR_sero_prev_plot$layers[[9]]$aes_params$linewidth <- 0.4
Model_Fit_1x1_plot$PCR_sero_prev_plot$layers[[11]]$aes_params$linewidth <- 0.4

###########################
Figure_5_plot <- cowplot::plot_grid(Model_Fit_1x1_plot$Poisson_Figure_weeks + theme(legend.position = c(0.8,0.9), legend.key.height = unit(4, "mm")) +
                                      guides(color = guide_legend(nrow = 2, override.aes = list(shape = c(20, NA), linetype = c(0,1)))) +
                                      labs(title = "a", y = "Burial registrations"),

                                    Model_Fit_1x1_plot$Poisson_Figure_age + theme(legend.position = c(0.6,0.9),
                                                                                  axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
                                                                                  legend.key.height = unit(4, "mm")) +
                                      guides(color = guide_legend(nrow = 2, override.aes = list(shape = c(20, NA), linetype = c(0,1)))) +
                                      labs(title = "b", y = "Burial registrations") + ylim(0,700),

                                    Model_Fit_1x1_plot$PCR_sero_prev_plot + theme(legend.position = c(0.3,0.7), legend.key.height = unit(4, "mm")) +
                                      coord_cartesian(xlim = as.Date(c("2020-05-01","2020-10-05")), ylim = c(0,0.28)) +
                                      labs(title = "c") +
                                      guides(colour = guide_legend(override.aes = list(linewidth = 0.33,0.4,0.4,NA,NA,
                                                                                       shape = c(NA,NA,NA,19,19),
                                                                                       linetype = c(1,2,2,NA,NA)))),

                                    Model_Fit_1x1_plot$Week_prev_plot + theme(legend.position = c(0.7,0.9), legend.key.height = unit(2.5, "mm")) +
                                      ggtitle("d") + scale_y_continuous(labels = scales::percent),

                                    Model_Fit_1x1_plot$Age_prev_plot +
                                      theme(legend.position = c(0.55,0.9), legend.key.height = unit(2.5, "mm")) +
                                      ggtitle("e"),

                                    Model_Fit_1x1_plot$p_Rt_Reff_a +
                                      coord_cartesian(xlim = as.Date(c("2020-04-24","2020-10-13")), ylim = c(0,5)) + theme_minimal(base_size = 7, base_family = "Helvetica") + xlab("Date") +
                                      theme(legend.position = c(0.8,0.9),
                                            plot.title = element_text(size = 7, face = "bold"),
                                            legend.text = element_text(hjust = 0),
                                            legend.key.size = unit(3, "mm"), axis.title.y = element_text(vjust = -0.1)) +
                                      guides(fill = guide_legend(nrow = 2, reverse = T), alpha = F) +
                                      scale_fill_manual(name = "", labels =c(expression(italic(R)[eff]),expression(italic(R)[italic(0)](italic(t)))), values = c("#48996b", "#3f8da7")) +
                                      ggtitle("f") + ylab("Time-varying reproduction number"),
                                    ncol = 3, nrow = 2, align = "hv")

ggsave("Analysis/Figures/Figure_5_Transmission.pdf",
       Figure_5_plot,
       width = 180, height = 180*6/10, units = "mm")

# saveRDS(Res_10_1x1_plot, "../Bonus Files/2023_March_Reviewer_Batch/Figure_3_plot_data.rds")
# Res_10_1x1_plot <- readRDS("../Bonus Files/2023_March_Reviewer_Batch/Figure_3_plot_data.rds")

