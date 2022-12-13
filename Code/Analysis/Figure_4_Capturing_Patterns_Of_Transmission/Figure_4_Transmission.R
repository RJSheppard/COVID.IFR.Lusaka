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

##############################################
##'[Model Fitting Under Default Assumptions]##
##############################################
## Load data
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

Model_Fit <- fit_spline_rt(
  country = "Zambia",
  data = data,

  population = population,
  baseline_contact_matrix = baseline_contact_matrix,

  combined_data = list(Comb_data = Comb_data,
                       dfj_mcmc_data = dfj_mcmc_data),
  pcr_df = pcr_df,
  sero_df = sero_df,

  prob_hosp = probs_hosp_death$Prob_Hosp$X41,
  prob_non_severe_death_treatment = probs_hosp_death$Prob_death$X41,
  prob_severe = rep(0,17),
  hosp_beds = 1e10,
  icu_beds = 1e10,

  dur_get_ox_survive = dur_get_ox_survive,
  dur_get_ox_die = dur_death_default,

  pcr_det = pcr_det_100,
  pcr_det_PM = pcr_det_100,

  dur_R = Inf,

  reporting_fraction = 1,

  frac_reg = 0.9,

  log_likelihood = Llike,
  lld = "",

  n_mcmc = 30000,
  replicates = 100,
  rw_duration = 14,
  Prior_Rt_rw_unif_lim = 1
)

# saveRDS(object = Model_Fit, file = "../Bonus Files/Model_Fit_Default_Assumptions.rds")

###########################
##'[Plot the for 1x1 fit]##
###########################
IFR_coefficients <- readRDS("~/Documents/Zambia/covid-mortality-ascertainment/analysis/data/Code-generated-data/00_03_IFR_matrix_coefficients_log_scale_new_pop_str_ests.rds") %>%
  mutate(Index = 1:nrow(.))

Model_Fit_1x1_plot <- Diagnostic_Plot(1, fit_Model_list = Model_Fit, IFRvals = IFR_coefficients[41,], Return_likelihoods_only = F)
###########################
Figure_4_plot <- cowplot::plot_grid(Model_Fit_1x1_plot[[1]]$Poisson_Figure_weeks + theme(legend.position = c(0.8,0.9),
                                                                                      axis.text.x = element_text(size = 11),
                                                                                      legend.text = element_text(size = 11)) +
                                      guides(color = guide_legend(nrow = 2, override.aes = list(shape = c(20, NA), linetype = c(0,1)))) +
                                      ggtitle("A") + ylab("Burial registrations"),
                                    Model_Fit_1x1_plot[[1]]$Poisson_Figure_age + theme(legend.position = c(0.6,0.9),
                                                                                    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
                                                                                    legend.text = element_text(size = 11)) +
                                      guides(color = guide_legend(nrow = 2, override.aes = list(shape = c(20, NA), linetype = c(0,1)))) +

                                      ggtitle("B") + ylab("Burial registrations") + ylim(0,700),
                                    Model_Fit_1x1_plot[[1]]$PCR_sero_prev_plot + theme(legend.position = c(0.3,0.8),
                                                                                    axis.text.x = element_text(size = 11),
                                                                                    legend.text = element_text(size = 11)) +
                                      guides(color = guide_legend(nrow = 3, override.aes = list(shape = c(NA,20,20)))) +
                                      coord_cartesian(xlim = as.Date(c("2020-05-01","2020-10-05")), ylim = c(0,0.28)) +
                                      ggtitle("C"),
                                    Model_Fit_1x1_plot[[1]]$Week_prev_plot + theme(legend.position = c(0.7,0.9),
                                                                                legend.text = element_text(size = 11),
                                                                                axis.text.x = element_text(size = 11)) +
                                      guides(color = guide_legend(nrow = 3, override.aes = list(shape = c(20,NA,NA)))) +
                                      ggtitle("D") + scale_y_continuous(labels = scales::percent),
                                    Model_Fit_1x1_plot[[1]]$Age_prev_plot + theme(legend.position = c(0.55,0.9), legend.text = element_text(size = 11)) + guides(colour = guide_legend(nrow = 3, override.aes = list(shape = c(20,NA,NA)))) + ggtitle("E"),
                                    Model_Fit_1x1_plot[[1]]$p_Rt_Reff_a +
                                      coord_cartesian(xlim = as.Date(c("2020-04-24","2020-10-13")), ylim = c(0,5)) + theme_minimal() + xlab("Date") +
                                      theme(legend.position = c(0.9,0.9),
                                            axis.text.x = element_text(size = 11),
                                            legend.text = element_text(size = 11, hjust = 0)) + guides(fill = guide_legend(nrow = 2, reverse = T)) + ggtitle("F"),
                                    ncol = 3, nrow = 2, align = "hv")

pdf(file = "Figure_4.pdf", width = 10, height = 6)
Figure_4_plot
dev.off()

tiff(file = "Figure_4.tiff", units = "in", height = 6, width = 10, res = 250)
Figure_4_plot
dev.off()
