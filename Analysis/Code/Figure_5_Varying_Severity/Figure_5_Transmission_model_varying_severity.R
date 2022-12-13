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
## Official data for generating initial start date
data <- Off_data <- readRDS(file = "Analysis/Data/derived_data/00_01_Lusaka_Dist_Deaths_Official.rds")

# Lusaka population
population <- readRDS("Analysis/Data/derived_data/00_02_Lusaka_Dist_Pop_Str_2020_imp_ests.rds")
# Nyanga contact matrix
baseline_contact_matrix <- as.matrix(readRDS("Analysis/Data/derived_data/00_11_Nyanga_Mixing_Matrix.rds"))

# Burial registrations and Post-Mortem data
Comb_data <- readRDS("Analysis/Data/derived_data/00_13_Combined_bur_regs_postmortem_data_complete.rds")
Comb_data <- Comb_data %>% mutate(PosTests = PosTests_Strict)

# Baseline registration estimates
dfj_mcmc_data <- readRDS(file = "Analysis/Data/derived_data/00_16_03_drj_mcmc_data_new_pop_str_2.rds")

# Population PCR and seroprevalence data
pcr_df <- readRDS("Analysis/Data/derived_data/00_10_Lancet_Data.rds")$pcr_df
sero_df <- readRDS("Analysis/Data/derived_data/00_10_Lancet_Data.rds")$sero_df

# probability of hospitalisation and death
probs_hosp_death <- readRDS("Analysis/Data/derived_data/00_03_IFR_probs_death_hosp.rds")
names(probs_hosp_death) <- paste0("X",1:81)

# Durations until death or survival following infection
Weighted_Durs_Hosp <- readRDS("Analysis/Data/derived_data/00_04_Weighted_durations_death_survive.rds")
dur_get_ox_survive <- Weighted_Durs_Hosp$Surv_Dur_Weighted
dur_get_ox_die <- Weighted_Durs_Hosp$Death_Dur_Weighted
# Assume that deaths without hospital treatment have half duration of treatment and that 70% die at home
dur_death_default <- dur_get_ox_die*0.3 + (dur_get_ox_die*0.5)*0.7

# PCR prevalence with 100% maximum sensitivity
pcr_det_100 <- readRDS("Analysis/Data/derived_data/00_15_pcr_det_hall_100.rds")

#########################
###????### Should this be a new package??? I suppose I can create a data folder, just to make sure everything works...
#########################
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


IFR_coefficients <- readRDS("Analysis/Data/derived_data/00_03_IFR_matrix_coefficients_log_scale_new_pop_str_ests.rds") %>%
  mutate(Index = 1:nrow(.))

t_full_fit <- obj_cma_01$lapply(probs_hosp_death[1:71], cma_fit,

                                data = data,

                                population = population,
                                baseline_contact_matrix = baseline_contact_matrix,

                                combined_data = list(Comb_data = Comb_data,
                                                     dfj_mcmc_data = dfj_mcmc_data),

                                pcr_df = pcr_df,
                                sero_df = sero_df,

                                dur_get_ox_survive = dur_get_ox_survive,
                                dur_get_ox_die = dur_death_default,

                                pcr_det = pcr_det_100,
                                pcr_det_PM = pcr_det_100,

                                frac_reg = 0.9,

                                log_likelihood = Llike,
                                lld = "",
                                Prior_Rt_rw_unif_lim = 1,

                                n_mcmc = 30000, replicates = 100

)

t_full_fit_finished <- t_full_fit$results()
names(t_full_fit_finished) <- paste0("X",1:71)
saveRDS(t_full_fit_finished, file = "Analysis/Results/Full_Set_Default_Setting.rds")

#######################################
##'[Plot heatmap results for fit]##
#######################################
library(ggplot2)
library(dplyr)
library(tidyr)
devtools::load_all()

t_full_fit_finished_LL_only <- lapply(X =1:length(t_full_fit_finished), FUN = Diagnostic_Plot, fit_Model_list = t_full_fit_finished, IFRvals = IFR_coefficients,
                         Return_likelihoods_only = T)

t_full_fit_finished_Heatmaps <- Plot_Heatmaps(Mod_Res = t_full_fit_finished, Res_Figs = t_full_fit_finished_LL_only, Select_Runs = 1:81, Title = "")

############################
##'[Generate infographics]##
############################
## Use IFR slope/intercept to get actual estimates:
Age_groups <- readRDS("Analysis/Data/derived_data/00_03_IFR_values_Brazeau.rds")$IFR_Age_gr
IFR_Ests_age <- apply(t_full_fit_finished_Heatmaps$IFR_mat,1, function(x){exp(x["IFR_abs"] + x["Slope_abs"]*Age_groups)}) %>% t() %>% as.data.frame()
colnames(IFR_Ests_age) <- Age_groups

IFR_df <- readRDS("Analysis/Data/derived_data/00_03_IFR_matrix_coefficients_log_scale.rds") %>%
  filter(IFR_x %in% c(0.2,1,5) & Slope_x == 1 | Slope_x %in% c(0.2,1,2.5) & IFR_x == 1)

IFR_df_by_age <- cbind(IFR_df %>% mutate(IFR_no = c("Age gradient x0.2",
                                                    "Overall severity x0.2",
                                                    "Default",
                                                    "Overall severity x5",
                                                    "Age gradient x2.5")),
                       t(apply(IFR_df, 1, function(x){
                         df_tmp <- exp(x["Int_abs"] + x["Slope_abs"]*Age_groups)
                       })))
colnames(IFR_df_by_age)[-c(1:6)] <- Age_groups

IFR_df_by_age_plot <- IFR_df_by_age %>%
  mutate(Top = `82.5`) %>%
  pivot_longer(cols = 7:24, names_to = "Age", values_to = "IFR") %>%
  mutate(Age = replace(Age, Age == "Top", Inf))

IFR_df_by_age_plot$IFR_no <- factor(IFR_df_by_age_plot$IFR_no, levels = c("Default", "Overall severity x5", "Overall severity x0.2", "Age gradient x2.5", "Age gradient x0.2"))

IFRs_Braz <- readxl::read_xlsx("Analysis/Data/raw_data/Brazeau_et_al.xlsx")

p_info_1 <- ggplot(data = IFR_df_by_age_plot %>% filter(Slope_x != 5) %>% select(IFR_no, Age, IFR),
                   aes(x = as.numeric(Age)-2.5, y = IFR/100, group = IFR_no, linetype = as.factor(IFR_no), color = as.factor(IFR_no))) +
  geom_step() +
  coord_cartesian(ylim = c(0,0.1)) +
  theme_minimal() +
  xlab("Age") + ylab("IFR") +
  scale_y_continuous(labels = scales::percent) +
  scale_linetype_manual(name = "", values = c(1,2,2,4,4)) +
  scale_color_manual(name = "", values = c("black", "darkred","red2","darkgoldenrod4","goldenrod2")) +
  ggtitle("A")

library(scales)
p_info_2 <- ggplot(data = IFR_df_by_age_plot %>% filter(Slope_x != 5) %>% select(IFR_no, Age, IFR), aes(x = as.numeric(Age)-2.5, y = IFR/10, group = as.factor(IFR_no), linetype = as.factor(IFR_no),color = as.factor(IFR_no))) +
  geom_step() +
  coord_cartesian(ylim = c(1E-9,1E3)) +
  theme_minimal() +
  xlab("Age") + ylab("IFR") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = function(x) paste0(x, '%')) +
  scale_linetype_manual(name = "", values = c(1,2,2,4,4)) +
  scale_color_manual(name = "", values = c("black", "darkred","red2","darkgoldenrod4","goldenrod2")) +
  ggtitle("B")

P_info_legend <- cowplot::plot_grid(ggpubr::get_legend(p = p_info_2))

## Colour IFR curves by goodness of fit:
IFR_df_full <- readRDS("Analysis/Data/derived_data/00_03_IFR_matrix_coefficients_log_scale.rds") %>%
  mutate(linetype_col = t_full_fit_finished_Heatmaps$IFR_mat$Post_col_group) %>%
  arrange(linetype_col) %>%
  mutate(IFR_no = 1:81)

IFR_df_by_age_full <- cbind(IFR_df_full,
                            t(apply(IFR_df_full, 1, function(x){df_tmp <- exp(x["Int_abs"] + x["Slope_abs"]*Age_groups)})))
colnames(IFR_df_by_age_full)[-c(1:7)] <- Age_groups

IFR_df_by_age_plot_full <- IFR_df_by_age_full %>%
  mutate(Top = `82.5`) %>%
  pivot_longer(cols = 8:25, names_to = "Age", values_to = "IFR") %>%
  mutate(Age = replace(Age, Age == "Top", Inf))


p_colored_lines <- ggplot(data = IFR_df_by_age_plot_full %>% filter(Slope_x != 5) %>% select(IFR_no, Age, IFR, linetype_col), aes(x = as.numeric(Age)-2.5, y = IFR/100, group = IFR_no, color = as.factor(linetype_col))) +
  geom_step(linewidth = 1.2) +
  geom_step(data = IFR_df_by_age_plot_full %>% filter(Slope_x == 1 & IFR_x == 1) %>% select(IFR_no, Age, IFR, linetype_col), aes(x = as.numeric(Age)-2.5, y = IFR/100), color= "black", linewidth = 0.8) +
  coord_cartesian(ylim = c(0,0.1)) +
  theme_minimal() +
  xlab("Age") + ylab("IFR") +
  scale_y_continuous(labels = scales::percent) +
  viridis::scale_color_viridis(name = "", discrete = T) +
  theme(legend.position = "none") +
  ggtitle("D")

p_colored_lines_log <- ggplot(data = IFR_df_by_age_plot_full %>% filter(Slope_x != 5) %>% select(IFR_no, Age, IFR, linetype_col), aes(x = as.numeric(Age)-2.5, y = IFR, group = IFR_no, color = as.factor(linetype_col))) +
  geom_step(linewidth = 1.2) +
  geom_step(data = IFR_df_by_age_plot_full %>% filter(Slope_x == 1 & IFR_x == 1) %>% select(IFR_no, Age, IFR, linetype_col), aes(x = as.numeric(Age)-2.5, y = IFR), color= "black", linewidth = 0.8) +
  coord_cartesian(ylim = c(1E-9,1E3)) +
  theme_minimal() +
  xlab("Age") + ylab("IFR") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = function(x) paste0(x, '%')) +
  viridis::scale_color_viridis(name = "", discrete = T) +
  theme(legend.position = "none") +
  ggtitle("E")

h <- t_full_fit_finished_Heatmaps$p1+ ggtitle("C") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5))

p_i <- cowplot::plot_grid(p_info_1 + theme(legend.position = c(0.4, 0.75), legend.text = element_text(size = 8),
                                           legend.spacing.y = unit(-0.2, 'cm')) +
                            guides(color = guide_legend(byrow = TRUE),
                                   linetype = guide_legend(byrow = TRUE)),
                          p_info_2 + theme(legend.position = "none"),
                          ncol = 1)

c1 <- p_colored_lines
c2 <- p_colored_lines_log
c_i <- cowplot::plot_grid(c1,c2, ncol = 1)

##########################
tiff("Analysis/Figures/Figure_5_severity_heatmaps.tiff", height = 5, width = 10, units = "in", res = 150)
cowplot::plot_grid(p_i, h, c_i, nrow = 1, rel_widths = c(1,2,1))
dev.off()
