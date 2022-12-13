#########################################
#########################################
##'[Plot Supplementary Results figures]##
#########################################
#########################################

### Supplementary plot:
Select_Runs <- IFR_coefficients %>%
  filter(Slope_x == 1 & round(IFR_x,2) %in% c(0.4,0.6,0.8,1,1.25,1.67) |
           Slope_x == 0.8 & round(IFR_x,2) %in% c(0.4,0.6,0.8,1,1.25,1.67) |
           Slope_x == 1.25 & round(IFR_x,2) %in% c(0.4,0.6,0.8,1) |
           round(Slope_x,2) == 1.67 & round(IFR_x,2) %in% c(0.4,0.6,0.8,1)) %>%
  pull(Index)

t_full_fit_finished_supp <- t_full_fit_finished[Select_Runs]
Title <- "Supplememtary: Default Model Fit Diagnostics"
devtools::load_all()
t_full_fit_finished_supp_plots <- lapply(X =1:length(t_full_fit_finished_supp), FUN = Diagnostic_Plot_10, fit_Model_list = t_full_fit_finished_supp, IFRvals = IFR_coefficients[Select_Runs,],
                                         Return_likelihoods_only = F)

saveRDS(t_full_fit_finished_supp_plots, "../Bonus Files/2022-11-10_Supp_plots.rds")
t_full_fit_finished_supp_plots <- readRDS("../Bonus Files/2022-11-10_Supp_plots.rds")

Supp_plots <- lapply(1:length(t_full_fit_finished_supp_plots), function(x){
  title_gg <- ggplot() +
    labs(title = paste0("Default fit (overall severity: ", round(IFR_coefficients[Select_Runs,"IFR_x"][x],2),"x, age gradient: ",round(IFR_coefficients[Select_Runs,"Slope_x"][x],2),"x)" )) +
    theme_minimal()

  cowplot::plot_grid(title_gg,
                     cowplot::plot_grid(t_full_fit_finished_supp_plots[[x]]$Poisson_Figure_weeks + theme(legend.position = "none") + ggtitle("A"),
                                        t_full_fit_finished_supp_plots[[x]]$Poisson_Figure_age + theme(legend.position = "none") + ggtitle("B"),
                                        t_full_fit_finished_supp_plots[[x]]$PCR_sero_prev_plot + theme(legend.position = "none") + ggtitle("C"), ncol = 3),
                     cowplot::plot_grid(ggpubr::get_legend(t_full_fit_finished_supp_plots[[x]]$Poisson_Figure_weeks),
                                        ggpubr::get_legend(t_full_fit_finished_supp_plots[[x]]$PCR_sero_prev_plot + guides(colour = guide_legend(nrow = 1))), rel_widths = c(2,1), ncol = 2),
                     cowplot::plot_grid(t_full_fit_finished_supp_plots[[x]]$Week_prev_plot + theme(legend.position = "none") + ggtitle("D"),
                                        t_full_fit_finished_supp_plots[[x]]$Age_prev_plot + theme(legend.position = "none") + ggtitle("E"),
                                        t_full_fit_finished_supp_plots[[x]]$p_Rt_Reff_a + coord_cartesian(as.Date(c("2020-06-01","2020-10-05"))) + theme_minimal() + theme(legend.position = "none") + ggtitle("F"), ncol = 3),
                     cowplot::plot_grid(ggpubr::get_legend(t_full_fit_finished_supp_plots[[x]]$Week_prev_plot),
                                        ggpubr::get_legend(t_full_fit_finished_supp_plots[[x]]$p_Rt_Reff_a), rel_widths = c(2,1), ncol = 2),
                     ncol = 1, rel_heights = c(0.1,1,0.2,1,0.2))
})


pdf(file = "analysis/figures/39_figure_3_supp.pdf", width = 11)
Supp_plots
dev.off()

tiff(file = "analysis/figures/39_figure_3_supp.tiff", units = "in", height = 7, width = 11, res = 250)
Supp_plots
dev.off()

###### I think that's it then.






devtools::load_all()
library(squire)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(Brobdingnag)

IFR_coefficients <- readRDS("analysis/data/Code-generated-data/00_03_IFR_matrix_coefficients_log_scale_new_pop_str_ests.rds") %>%
  mutate(Index = 1:nrow(.))

#######################################
#######################################
##'[Plot Supplementary Sens 1]##
#######################################
#######################################
Select_Runs_set <- IFR_coefficients %>%
  filter(Slope_x %in% c(0.8,1,1.25) & round(IFR_x,2) %in% c(0.4,0.6,0.8,1,1.25,1.67, 2.5)) %>%  pull(Index)


### Supplementary plot:
Select_Runs <- IFR_coefficients %>%
  filter(Slope_x == 0.8 & round(IFR_x,2) %in% c(0.8,1,1.25,1.67) |
           Slope_x == 1 & round(IFR_x,2) %in% c(0.6,0.8,1,1.25,1.67,2.5) |
           Slope_x == 1.25 & round(IFR_x,2) %in% c(0.4,0.6,0.8,1,1.25,1.67)) %>%
  pull(Index)

# Res_10_supp <- readRDS("../Bonus Files/2022-10-20_L10_default_full_set.rds")[Select_Runs]
Res_10_supp <- Res_10[names(Res_10) %in% paste0("X",Select_Runs)]
Title <- "Supplememtary: Sens 1 Diagnostics"
# Res_10_supp <- lapply(Res_10_supp, function(x){if(class(x)[1] != "squire_simulation"){NULL}else{x}})
# dfj_mcmc_data <- readRDS("analysis/data/Code-generated-data/00_16_03_drj_mcmc_data_new_pop_str_2.rds")
devtools::load_all()
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)
Res_10_supp_plots <- lapply(X =1:length(Res_10_supp), FUN = Diagnostic_Plot_10, fit_Model_list = Res_10_supp, IFRvals = IFR_coefficients[Select_Runs,],
                            Return_likelihoods_only = F)

saveRDS(Res_10_supp_plots, "../Bonus Files/2022-10-30_Supp_plots_S1.rds")

Supp_plots <- lapply(1:length(Res_10_supp_plots), function(x){

  title_gg <- ggplot() +
    labs(title = paste0("Default fit (overall severity: ", round(IFR_coefficients[Select_Runs,"IFR_x"][x],2),"x, age gradient: ",round(IFR_coefficients[Select_Runs,"Slope_x"][x],2),"x)" )) +
    theme_minimal()

  cowplot::plot_grid(title_gg,#cowplot::plot_grid(ggplot() + theme_void() + ggtitle("A") + theme(plot.title = element_text(hjust = 0.5, vjust = - 4)),
                     #                  ggplot() + theme_void() + ggtitle("C") + theme(plot.title = element_text(hjust = 0.5, vjust = - 4)), rel_widths = c(2,1)),
                     cowplot::plot_grid(Res_10_supp_plots[[x]]$Poisson_Figure_weeks + theme(legend.position = "none") + ggtitle("A"),
                                        Res_10_supp_plots[[x]]$Poisson_Figure_age + theme(legend.position = "none") + ggtitle("B"),
                                        Res_10_supp_plots[[x]]$PCR_sero_prev_plot + theme(legend.position = "none") + ggtitle("C"), ncol = 3),
                     cowplot::plot_grid(ggpubr::get_legend(Res_10_supp_plots[[x]]$Poisson_Figure_weeks),
                                        ggpubr::get_legend(Res_10_supp_plots[[x]]$PCR_sero_prev_plot + guides(colour = guide_legend(nrow = 1))), rel_widths = c(2,1), ncol = 2),
                     # cowplot::plot_grid(ggplot() + theme_void() + ggtitle("UTH mortuary post-mortem PCR prevalence") + theme(plot.title = element_text(hjust = 0.5, vjust = - 4)),
                     # ggplot() + theme_void() + ggtitle("Modelled transmissibility") + theme(plot.title = element_text(hjust = 0.5, vjust = - 4)), rel_widths = c(2,1)),
                     cowplot::plot_grid(Res_10_supp_plots[[x]]$Week_prev_plot + theme(legend.position = "none") + ggtitle("D"),
                                        Res_10_supp_plots[[x]]$Age_prev_plot + theme(legend.position = "none") + ggtitle("E"),
                                        #Res_10_1x1_plot[[1]]$p_Rt_Reff_a + theme_minimal() + theme(legend.position = "none") + ggtitle("F"), ncol = 3),
                                        Res_10_supp_plots[[x]]$p_Rt_Reff_a + theme(legend.position = "none") + ggtitle("F"), ncol = 3),
                     cowplot::plot_grid(ggpubr::get_legend(Res_10_supp_plots[[x]]$Week_prev_plot),
                                        ggpubr::get_legend(Res_10_supp_plots[[x]]$p_Rt_Reff_a), rel_widths = c(2,1), ncol = 2),
                     ncol = 1, rel_heights = c(0.1,1,0.2,1,0.2))
})


pdf(file = "analysis/figures/52_supp_plot_S1.pdf", width = 11)
Supp_plots
dev.off()

tiff(file = "analysis/figures/52_supp_plot_S1.tiff", units = "in", height = 7, width = 11, res = 250)
Supp_plots
dev.off()


### Supplementary plot:
dfj_mcmc_data <- readRDS("analysis/data/Code-generated-data/00_16_03_drj_mcmc_data_new_pop_str_2.rds")
devtools::load_all()
library(dplyr)
library(reshape2)
library(ggplot2)
library(tidyr)

Supp_Plot_func <- function(Res_10, Select_Runs, sens_name){
  Res_10_supp_plots <- lapply(X =1:length(Select_Runs), FUN = Diagnostic_Plot_10, fit_Model_list = Res_10[names(Res_10) %in% paste0("X",Select_Runs)], IFRvals = IFR_coefficients[Select_Runs,],
                              Return_likelihoods_only = F)

  saveRDS(Res_10_supp_plots, paste0("../Bonus Files/Supp_plots_",sens_name,".rds"))

  Supp_plots <- lapply(1:length(Res_10_supp_plots), function(x){

    title_gg <- ggplot() +
      labs(title = paste0(sens_name," (overall severity: ", round(IFR_coefficients[Select_Runs,"IFR_x"][x],2),"x, age gradient: ",round(IFR_coefficients[Select_Runs,"Slope_x"][x],2),"x)" )) +
      theme_minimal()

    cowplot::plot_grid(title_gg,#cowplot::plot_grid(ggplot() + theme_void() + ggtitle("A") + theme(plot.title = element_text(hjust = 0.5, vjust = - 4)),
                       #                  ggplot() + theme_void() + ggtitle("C") + theme(plot.title = element_text(hjust = 0.5, vjust = - 4)), rel_widths = c(2,1)),
                       cowplot::plot_grid(Res_10_supp_plots[[x]]$Poisson_Figure_weeks + theme(legend.position = "none") + ggtitle("A"),
                                          Res_10_supp_plots[[x]]$Poisson_Figure_age + theme(legend.position = "none") + ggtitle("B"),
                                          Res_10_supp_plots[[x]]$PCR_sero_prev_plot + theme(legend.position = "none") + ggtitle("C"), ncol = 3),
                       cowplot::plot_grid(ggpubr::get_legend(Res_10_supp_plots[[x]]$Poisson_Figure_weeks),
                                          ggpubr::get_legend(Res_10_supp_plots[[x]]$PCR_sero_prev_plot + guides(colour = guide_legend(nrow = 1))), rel_widths = c(2,1), ncol = 2),
                       # cowplot::plot_grid(ggplot() + theme_void() + ggtitle("UTH mortuary post-mortem PCR prevalence") + theme(plot.title = element_text(hjust = 0.5, vjust = - 4)),
                       # ggplot() + theme_void() + ggtitle("Modelled transmissibility") + theme(plot.title = element_text(hjust = 0.5, vjust = - 4)), rel_widths = c(2,1)),
                       cowplot::plot_grid(Res_10_supp_plots[[x]]$Week_prev_plot + theme(legend.position = "none") + ggtitle("D"),
                                          Res_10_supp_plots[[x]]$Age_prev_plot + theme(legend.position = "none") + ggtitle("E"),
                                          #Res_10_1x1_plot[[1]]$p_Rt_Reff_a + theme_minimal() + theme(legend.position = "none") + ggtitle("F"), ncol = 3),
                                          Res_10_supp_plots[[x]]$p_Rt_Reff_a + theme(legend.position = "none") + ggtitle("F"), ncol = 3),
                       cowplot::plot_grid(ggpubr::get_legend(Res_10_supp_plots[[x]]$Week_prev_plot),
                                          ggpubr::get_legend(Res_10_supp_plots[[x]]$p_Rt_Reff_a), rel_widths = c(2,1), ncol = 2),
                       ncol = 1, rel_heights = c(0.1,1,0.2,1,0.2))
  })
# browser()
  return(Supp_plots)
}


#######################################
#######################################
##'[Plot Supplementary Sens 2]##
#######################################
#######################################

devtools::load_all()
Res_10_Sens_2 <- readRDS("../Bonus Files/2022-11-10_Sensitivity_Analysis_2_Remove_Age_gr_1.rds")
Select_Runs <- IFR_coefficients %>%
  filter(Slope_x %in% c(0.8,1,1.25) & round(IFR_x,2) %in% c(0.6,0.8,1,1.25,1.67)) %>%
  pull(Index)

names(Res_10_Sens_2) <- paste0("X",Select_Runs_set)



Supp_plots_2 <- Supp_Plot_func(Res_10_Sens_2, Select_Runs, "02_Remove_age_group_1")
pdf(file = paste0("analysis/figures/55_supp_plot_02_Remove_age_group_1.pdf"), width = 11)
Supp_plots_2
dev.off()

tiff(file = paste0("analysis/figures/55_supp_plot_02_Remove_age_group_1.tiff"), units = "in", height = 7, width = 11, res = 250)
Supp_plots_2
dev.off()


devtools::load_all()
Res_10_Sens_3 <- readRDS("../Bonus Files/2022-10-30_Sensitivity_Analysis_2_Remove_Age_gr_1.rds")
Select_Runs <- IFR_coefficients %>%
  filter(Slope_x == 0.8 & round(IFR_x,2) %in% c(0.8,1,1.25,1.67) |
           Slope_x == 1 & round(IFR_x,2) %in% c(0.6,0.8,1,1.25,1.67,2.5) |
           Slope_x == 1.25 & round(IFR_x,2) %in% c(0.8,1,1.25,1.67)) %>%
  pull(Index)

Supp_plots_2 <- Supp_Plot_func(Res_10_Sens_2, Select_Runs, "02_Remove_age_group_1")
pdf(file = paste0("analysis/figures/55_supp_plot_02_Remove_age_group_1.pdf"), width = 11)
Supp_plots_2
dev.off()

tiff(file = paste0("analysis/figures/55_supp_plot_02_Remove_age_group_1.tiff"), units = "in", height = 7, width = 11, res = 250)
Supp_plots_2
dev.off()


# Res_10_supp <- readRDS("../Bonus Files/2022-10-20_L10_default_full_set.rds")[Select_Runs]
# Res_10_supp_2 <- Res_10[names(Res_10) %in% paste0("X",Select_Runs)]
# Title <- "Supplememtary: Sens 1 Diagnostics"
# Res_10_supp <- lapply(Res_10_supp, function(x){if(class(x)[1] != "squire_simulation"){NULL}else{x}})







#######################################
#######################################
##'[Plot Supplementary Sens 5a]##
#######################################
#######################################

devtools::load_all()
Res_10_Sens_5a <- readRDS("../Bonus Files/2022-11-19_Sensitivity_Analysis_5a_Change_age_RR_lower_20.rds")
Select_Runs <- IFR_coefficients %>%
  filter(Slope_x %in% c(0.8,1,1.25) & round(IFR_x,2) %in% c(0.8,1,1.25,1.67,2.5)) %>%
  pull(Index)

Supp_plots_5a <- Supp_Plot_func(Res_10_Sens_5a, Select_Runs, "Decrease baseline deaths by 20%")
pdf(file = paste0("analysis/figures/55_supp_plot_5a_Decrease_baseline_deaths_20pc.pdf"), width = 11)
Supp_plots_5a
dev.off()

tiff(file = paste0("analysis/figures/55_supp_plot_5a_Decrease_baseline_deaths_20pc.tiff"), units = "in", height = 7, width = 11, res = 250)
Supp_plots_5a
dev.off()

### 5b
devtools::load_all()
Res_10_Sens_5b <- readRDS("../Bonus Files/2022-11-19_Sensitivity_Analysis_5b_Change_age_RR_higher_20.rds")
Select_Runs <- IFR_coefficients %>%
  filter(Slope_x %in% c(0.8) & round(IFR_x,2) %in% c(0.6,0.8,1) |
           Slope_x %in% c(1) & round(IFR_x,2) %in% c(0.4,0.6,0.8,1) |
           Slope_x %in% c(1.25) & round(IFR_x,2) %in% c(0.4,0.6,0.8,1)) %>%
  pull(Index)

Supp_plots_5b <- Supp_Plot_func(Res_10_Sens_5b, Select_Runs, "Increase baseline deaths by 20%")
pdf(file = paste0("analysis/figures/55_supp_plot_5b_Increase_baseline_deaths_20pc.pdf"), width = 11)
Supp_plots_5b
dev.off()

tiff(file = paste0("analysis/figures/55_supp_plot_5b_Increase_baseline_deaths_20pc.tiff"), units = "in", height = 7, width = 11, res = 250)
Supp_plots_5b
dev.off()

### 5c
devtools::load_all()
Res_10_Sens_5c <- readRDS("../Bonus Files/2022-11-19_Sensitivity_Analysis_5c_Change_age_RR_lower_10.rds")
Select_Runs <- IFR_coefficients %>%
  filter(Slope_x %in% c(0.8,1,1.25) & round(IFR_x,2) %in% c(0.8,1,1.25,1.67,2.5)) %>%
  pull(Index)

Supp_plots_5c <- Supp_Plot_func(Res_10_Sens_5c, Select_Runs, "Decrease baseline deaths by 10%")
pdf(file = paste0("analysis/figures/55_supp_plot_5c_Decrease_baseline_deaths_10pc.pdf"), width = 11)
Supp_plots_5c
dev.off()

tiff(file = paste0("analysis/figures/55_supp_plot_5c_Decrease_baseline_deaths_10pc.tiff"), units = "in", height = 7, width = 11, res = 250)
Supp_plots_5c
dev.off()
