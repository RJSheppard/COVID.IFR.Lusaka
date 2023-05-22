## Get data inputs ready
library(tidyr)
library(dplyr)
library(ggplot2)

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

##############
## Full fit ##
##############
t_full_fit <- obj_cma_01$lapply(probs_hosp_death[Rerun[-c(73:81)]], cma_fit,

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

                                frac_reg = 0.9
)

# saveRDS(t_full_fit$results(), file = paste0("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/01_Full_fit_round_1.rds"))


devtools::load_all()

IFR_mat <- readRDS("analysis/data/Code-generated-data/00_03_IFR_matrix_coefficients_log_scale_new_pop_str_ests.rds")

# Default
# t_full_fit <- readRDS("../Bonus Files/2023_March_Reviewer_Batch/01_Full_fit_round_1.rds")
# t_full_fit_2 <- readRDS("../Bonus Files/2023_March_Reviewer_Batch/01_Full_fit_round_2.rds")
# order(as.numeric(gsub("X", "", names(c(Res_10,Res_10_2)))))
# t_full_fit <- c(t_full_fit,t_full_fit_2)[order(as.numeric(gsub("X", "", names(c(t_full_fit,t_full_fit_2)))))]
# t_full_fit <- t_full_fit[unlist(lapply(t_full_fit, function(x){any(class(x)=="squire_simulation")}))]
# rm(t_full_fit_2)
# gc()

H_default <- ggplot(Plot_Post(t_full_fit, IFR_mat) %>% filter(Slope_x !=5), aes(x = as.factor(round(IFR_x,2)), y = as.factor(round(Slope_x,2)), fill = as.factor(Post_col_group))) +
  geom_tile() +
  geom_text(aes(label = round(AvPost), colour = (Post_col_group >= max(Post_col_group, na.rm=T))), size = 7/.pt) +
  scale_colour_manual(values = c("white", "black")) +
  labs(title = "c", x = "Relative overall severity", y = "Relative IFR age gradient") +
  labs(fill = "Mean Post") +
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme(legend.position = "none",
        axis.text.y = element_text(angle = 90, hjust = 0.5), plot.title = element_text(size = 7, face = "bold")) +
  scale_fill_discrete(type = tail(viridis::viridis(n = 9), length(table(Plot_Post(Res_10, IFR_mat)$Post_col_group)))) +
  coord_cartesian(ylim = c(0.5,8.5), expand = F) +
  scale_x_discrete(expand = c(0,0), labels = paste0(100*c(0.2,0.4,0.6,0.8,1,1.25,1.67,2.5,5),"%")) +
  scale_y_discrete(expand = c(0,0), labels = paste0(100*c(0.2,0.4,0.6,0.8,1,1.25,1.67,2.5,5),"%"))

###########################
## Generate infographics ##
###########################
## Use IFR slope/intercept to get actual estimates:
Age_groups <- seq(2.5, 82.5, by = 5)
Age_groups_IFRs <- as.data.frame(do.call(rbind,lapply(readRDS("Analysis/Data/derived_data/IFR_phosp_pdeath.rds")[c(5,37,41,45,68)], function(x){unlist(x["IFR"])})))
colnames(Age_groups_IFRs) <- Age_groups

IFR_df_by_age <- cbind(Age_groups_IFRs %>% mutate(IFR_no = c("Age gradient 20%",
                                                             "Overall severity 20%",
                                                             "Default",
                                                             "Overall severity 500%",
                                                             "Age gradient 250%")))

IFR_df_by_age_plot <- IFR_df_by_age %>%
  mutate(Top = `82.5`) %>%
  pivot_longer(cols = c(1:17,19), names_to = "Age", values_to = "IFR") %>%
  mutate(Age = replace(Age, Age == "Top", Inf))

IFR_df_by_age_plot$IFR_no <- factor(IFR_df_by_age_plot$IFR_no, levels = c("Default", "Overall severity 500%", "Overall severity 20%", "Age gradient 250%", "Age gradient 20%"))

IFRs_Braz <- readxl::read_xlsx("Analysis/Data/raw_data/Brazeau_et_al.xlsx")

p_info_1 <- ggplot(data = IFR_df_by_age_plot,
                   aes(x = as.numeric(Age)-2.5, y = IFR, group = IFR_no, linetype = as.factor(IFR_no), color = as.factor(IFR_no))) +
  geom_step() +
  coord_cartesian(ylim = c(0,0.1)) +
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  labs(title = "a", x = "Age", y = "IFR") +
  theme(plot.title = element_text(size = 7, face = "bold")) +
  scale_y_continuous(labels = scales::percent) +
  scale_linetype_manual(name = "", values = c(1,2,2,4,4)) +
  scale_color_manual(name = "", values = c("black", "darkred","red2","darkgoldenrod4","goldenrod2"))

library(scales)
p_info_2 <- ggplot(data = IFR_df_by_age_plot,
                   aes(x = as.numeric(Age)-2.5, y = IFR/10, group = as.factor(IFR_no), linetype = as.factor(IFR_no),color = as.factor(IFR_no))) +
  geom_step() +
  coord_cartesian(ylim = c(1E-9,1E3)) +
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme(plot.title = element_text(size = 7, face = "bold")) +
  labs(title = "b", x = "Age", y = "IFR") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = function(x) paste0(x, '%')) +
  scale_linetype_manual(name = "", values = c(1,2,2,4,4)) +
  scale_color_manual(name = "", values = c("black", "darkred","red2","darkgoldenrod4","goldenrod2"))


P_info_legend <- cowplot::plot_grid(ggpubr::get_legend(p = p_info_2))

## Colour IFR curves by goodness of fit:
Age_groups_IFRs_Full <- as.data.frame(do.call(rbind,lapply(readRDS("Analysis/Data/derived_data/IFR_phosp_pdeath.rds"), function(x){unlist(x["IFR"])})))
colnames(Age_groups_IFRs_Full) <- Age_groups
Age_groups_IFRs_Full <- Age_groups_IFRs_Full %>% mutate(linetype_col = Plot_Post(Res_10, IFR_mat)$Post_col_group) %>%
  arrange(linetype_col) %>%
  mutate(IFR_no = 1:81) %>%
  mutate(Top = `82.5`) %>%
  pivot_longer(cols = c(1:17,20), names_to = "Age", values_to = "IFR") %>%
  mutate(Age = replace(Age, Age == "Top", Inf))


p_colored_lines <- ggplot(data = Age_groups_IFRs_Full, aes(x = as.numeric(Age)-2.5, y = IFR, group = IFR_no, color = as.factor(linetype_col))) +
  geom_step(linewidth = 1.2) +
  geom_step(data = IFR_df_by_age_plot %>% filter(IFR_no =="Default"), aes(x = as.numeric(Age)-2.5, y = IFR), color= "black", linewidth = 0.8) +
  coord_cartesian(ylim = c(0,0.1)) +
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme(plot.title = element_text(size = 7, face = "bold"), legend.position = "none") +
  labs(title = "d", x = "Age", y = "IFR") +
  scale_y_continuous(labels = scales::percent) +
  viridis::scale_color_viridis(name = "", discrete = T)

p_colored_lines_log <- ggplot(data = Age_groups_IFRs_Full, aes(x = as.numeric(Age)-2.5, y = IFR, group = IFR_no, color = as.factor(linetype_col))) +
  geom_step(linewidth = 1.2) +
  geom_step(data = IFR_df_by_age_plot %>% filter(IFR_no =="Default"), aes(x = as.numeric(Age)-2.5, y = IFR), color= "black", linewidth = 0.8) +
  coord_cartesian(ylim = c(1E-9,1E3)) +
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme(plot.title = element_text(size = 7, face = "bold"), legend.position = "none") +
  labs(title = "e", x = "Age", y = "IFR") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = function(x) paste0(x, '%')) +
  viridis::scale_color_viridis(name = "", discrete = T)

h <- H_default+ ggtitle("C")

p_i <- cowplot::plot_grid(p_info_1 + theme(legend.position = c(0.2, 0.75), legend.text = element_text(size = 6)) +
                            guides(color = guide_legend(byrow = TRUE),
                                   linetype = guide_legend(byrow = TRUE)),
                          p_info_2 + theme(legend.position = "none"),
                          ncol = 1)

c1 <- p_colored_lines
c2 <- p_colored_lines_log
c_i <- cowplot::plot_grid(c1,c2, ncol = 1)

##########################
ggsave(filename = "Analysis/Figures/Figure_6_Severity.pdf",
       cowplot::plot_grid(p_i, h, c_i, nrow = 1, rel_widths = c(1,2,1)),
       height = 180/2, width = 180, units = "mm")
dev.off()
