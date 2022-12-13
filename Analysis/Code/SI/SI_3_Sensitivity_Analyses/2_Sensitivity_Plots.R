getwd()
rm(list = ls())
gc()

devtools::load_all()
library(squire)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(Brobdingnag)
library(dplyr)

####
IFR_coefficients <- readRDS("Analysis/Data/derived_data/00_03_IFR_matrix_coefficients_log_scale_new_pop_str_ests.rds") %>%
  mutate(Index = 1:nrow(.))

Select_Runs <- IFR_coefficients %>%
  filter(Slope_x %in% c(0.8,1,1.25) & round(IFR_x,2) %in% c(0.4,0.6,0.8,1,1.25,1.67, 2.5)) %>%  pull(Index)

#######################################
#######################################
##'[Plot Sensitivity heatmap 2: remove age group 1]##
#######################################
#######################################

##########################
devtools::load_all()
Res_Sens_1 <- readRDS("Analysis/Results/SI_Sensitivity_Analysis_1_Remove_Age_gr_1.rds")
Res_LL_only_Sens_1 <- lapply(X =1:length(Res_Sens_1), FUN = Diagnostic_Plot_10, fit_Model_list = Res_Sens_1, IFRvals = IFR_coefficients[Select_Runs,],
                         Return_likelihoods_only = T)
Res_Heatmaps_Sens_1 <- Plot_Heatmaps(Mod_Res = c(Res_Sens_1), Res_Figs = c(Res_LL_only_Sens_1), Select_Runs = c(Select_Runs,as.numeric(gsub("X","",names(Res_Sens_1)))), Title = "")

sp1 <- cowplot::plot_grid(Res_Heatmaps_Sens_1$p1+ ggtitle("Overall") + coord_cartesian(xlim = c(1.5,8.5), ylim = c(3.5,6.5), expand = F) +
                            scale_x_discrete(expand = c(0,0), breaks = c(0.4,0.6,0.8,1,1.25,1.67,2.5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25)) +
                            theme(axis.text.y = element_text(angle = 90, hjust = 0.5)),
                          cowplot::plot_grid(Res_Heatmaps_Sens_1$p2+ ggtitle("Burial registration rate")+
                                               theme(plot.title = element_text(size=8),
                                                     axis.text = element_text(size = 6),
                                                     axis.title = element_text(size = 6),
                                                     plot.title.position = "plot",
                                                     axis.text.y = element_text(angle = 90, hjust = 0.5)) +
                                               scale_x_discrete(expand = c(0,0), breaks = c(0.6,1,1.67)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25)) +
                                               coord_cartesian(xlim = c(1.5,8.5), ylim = c(3.5,6.5), expand = F),
                                             Res_Heatmaps_Sens_1$p3+ ggtitle("C19 post-mortem prev.")+
                                               theme(plot.title = element_text(size=8),
                                                     axis.text = element_text(size = 6),
                                                     axis.title = element_text(size = 6),
                                                     plot.title.position = "plot",
                                                     axis.text.y = element_text(angle = 90, hjust = 0.5)) +
                                               scale_x_discrete(expand = c(0,0), breaks = c(0.6,1,1.67)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25)) +
                                               coord_cartesian(xlim = c(1.5,8.5), ylim = c(3.5,6.5), expand = F),
                                             Res_Heatmaps_Sens_1$p4+ ggtitle("C19 population surveys")+
                                               theme(plot.title = element_text(size=8),
                                                     axis.text = element_text(size = 6),
                                                     axis.title = element_text(size = 6),
                                                     plot.title.position = "plot",
                                                     axis.text.y = element_text(angle = 90, hjust = 0.5)) +
                                               scale_x_discrete(expand = c(0,0), breaks = c(0.6,1,1.67)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25)) +
                                               coord_cartesian(xlim = c(1.5,8.5), ylim = c(3.5,6.5), expand = F), ncol = 1), rel_widths = c(3,1))

tiff("Analysis/Figures/SI_Sens_1_remove_age_1.tiff", height = 3.5, width = 6, units = "in", res = 300)
sp1
dev.off()


#######################################
#######################################
##'[Plot Sensitivity heatmap 2: no scaling, no weeks 4 or 5]##
#######################################
#######################################

##########################
devtools::load_all()
Res_Sens_2 <- readRDS("Analysis/Results/SI_Sensitivity_Analysis_2_No_scaling_no_weeks_4_or_5.rds")
Res_LL_only_Sens_2 <- lapply(X =1:length(Res_Sens_2), FUN = Diagnostic_Plot_10, fit_Model_list = Res_Sens_2, IFRvals = IFR_coefficients[Select_Runs,],
                                Return_likelihoods_only = T)
Res_Heatmaps_Sens_2 <- Plot_Heatmaps(Mod_Res = Res_Sens_2, Res_Figs = Res_LL_only_Sens_2, Select_Runs = Select_Runs, Title = "")

sp2 <- cowplot::plot_grid(Res_Heatmaps_Sens_2$p1+ ggtitle("Overall") + coord_cartesian(xlim = c(1.5,8.5), ylim = c(3.5,6.5), expand = F) +
                            scale_x_discrete(expand = c(0,0), breaks = c(0.4,0.6,0.8,1,1.25,1.67,2.5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25)) +
                            theme(axis.text.y = element_text(angle = 90, hjust = 0.5)),
                          cowplot::plot_grid(Res_Heatmaps_Sens_2$p2+ ggtitle("Burial registration rate")+
                                               theme(plot.title = element_text(size=8),
                                                     axis.text = element_text(size = 6),
                                                     axis.title = element_text(size = 6),
                                                     plot.title.position = "plot",
                                                     axis.text.y = element_text(angle = 90, hjust = 0.5)) +
                                               scale_x_discrete(expand = c(0,0), breaks = c(0.4,0.6,0.8,1,1.25,1.67,2.5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25)) +
                                               coord_cartesian(xlim = c(1.5,8.5), ylim = c(3.5,6.5), expand = F),
                                             Res_Heatmaps_Sens_2$p3+ ggtitle("C19 post-mortem prevalence")+
                                               theme(plot.title = element_text(size=8),
                                                     axis.text = element_text(size = 6),
                                                     axis.title = element_text(size = 6),
                                                     plot.title.position = "plot",
                                                     axis.text.y = element_text(angle = 90, hjust = 0.5)) +
                                               scale_x_discrete(expand = c(0,0), breaks = c(0.4,0.6,0.8,1,1.25,1.67,2.5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25)) +
                                               coord_cartesian(xlim = c(1.5,8.5), ylim = c(3.5,6.5), expand = F),
                                             Res_Heatmaps_Sens_2$p4+ ggtitle("C19 population surveys")+
                                               theme(plot.title = element_text(size=8),
                                                     axis.text = element_text(size = 6),
                                                     axis.title = element_text(size = 6),
                                                     plot.title.position = "plot",
                                                     axis.text.y = element_text(angle = 90, hjust = 0.5)) +
                                               scale_x_discrete(expand = c(0,0), breaks = c(0.4,0.6,0.8,1,1.25,1.67,2.5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25)) +
                                               coord_cartesian(xlim = c(1.5,8.5), ylim = c(3.5,6.5), expand = F), ncol = 1), rel_widths = c(3,1))

tiff("Analysis/Figures/SI_Sens_2_no_scaling_no_weeks_4_5.tiff", height = 3.5, width = 7, units = "in", res = 300)
sp2
dev.off()

#######################################
#######################################
##'[Plot Sensitivity heatmap 3: Fit to burial registrations and population surveys, not post mortem data]##
#######################################
#######################################
devtools::load_all()

Res_Sens_3 <- readRDS("Analysis/Results/SI_Sensitivity_3_Only_burials_and_pop_prevalence.rds")
Res_LL_only_Sens_3 <- lapply(X =1:length(Res_Sens_3), FUN = Diagnostic_Plot_10, fit_Model_list = Res_Sens_3, IFRvals = IFR_coefficients[as.numeric(gsub(pattern = "X","",names(Res_Sens_3))),],
                                 Return_likelihoods_only = T)
Res_Heatmaps_Sens_3 <- Plot_Heatmaps(Mod_Res = Res_Sens_3, Res_Figs = Res_LL_only_Sens_3, Select_Runs = Select_Runs, Title = "")

sp3 <- cowplot::plot_grid(Res_Heatmaps_Sens_3$p1+ ggtitle("Overall") + coord_cartesian(xlim = c(1.5,9.5), ylim = c(3.5,7.5), expand = F) +
                             scale_x_discrete(expand = c(0,0), breaks = c(0.4,0.6,0.8,1,1.25,1.67,2.5,5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25,1.67)) +
                             theme(axis.text.y = element_text(angle = 90, hjust = 0.5)),
                           cowplot::plot_grid(Res_Heatmaps_Sens_3$p2+ ggtitle("Burial registration rate")+
                                                theme(plot.title = element_text(size=8),
                                                      axis.text = element_text(size = 6),
                                                      axis.title = element_text(size = 6),
                                                      plot.title.position = "plot",
                                                      axis.text.y = element_text(angle = 90, hjust = 0.5)) +
                                                scale_x_discrete(expand = c(0,0), breaks = c(0.4,0.6,0.8,1,1.25,1.67,2.5,5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25,1.67)) +
                                                coord_cartesian(xlim = c(1.5,9.5), ylim = c(3.5,7.5), expand = F),
                                              Res_Heatmaps_Sens_3$p4+ ggtitle("C19 population surveys")+
                                                theme(plot.title = element_text(size=8),
                                                      axis.text = element_text(size = 6),
                                                      axis.title = element_text(size = 6),
                                                      plot.title.position = "plot",
                                                      axis.text.y = element_text(angle = 90, hjust = 0.5)) +
                                                scale_x_discrete(expand = c(0,0), breaks = c(0.4,0.6,0.8,1,1.25,1.67,2.5,5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25,1.67)) +
                                                coord_cartesian(xlim = c(1.5,9.5), ylim = c(3.5,7.5), expand = F), ncol = 1), rel_widths = c(2,1))

pdf("Analysis/Figures/SI_Sens_3_Only_burials_and_pop_prevalence.pdf", height = 3, width = 7)
sp3
dev.off()
tiff("Analysis/Figures/SI_Sens_3_Only_burials_and_pop_prevalence.tiff", height = 3, width = 7, units = "in", res = 300)
sp3
dev.off()


#######################################
#######################################
##'[Plot Sensitivity heatmap 4a: Change age relative risk: decrease by 10%]##
#######################################
#######################################
Res_Sens_4a <- readRDS("Analysis/Results/SI_Sensitivity_Analysis_4a_Change_age_RR_lower_10.rds")
Res_LL_only_Sens_4a <- lapply(X =1:length(Res_Sens_4a), FUN = Diagnostic_Plot_10, fit_Model_list = Res_Sens_4a, IFRvals = IFR_coefficients[as.numeric(gsub(pattern = "X","",names(Res_Sens_4a))),],
                                 Return_likelihoods_only = T)
Res_Heatmaps_Sens_4a <- Plot_Heatmaps(Mod_Res = Res_Sens_4a, Res_Figs = Res_LL_only_Sens_4a, Select_Runs = Select_Runs, Title = "")

sp4a <- cowplot::plot_grid(Res_Heatmaps_Sens_4a$p1+ ggtitle("Overall") + coord_cartesian(xlim = c(1.5,9.5), ylim = c(3.5,6.5), expand = F) +
                             scale_x_discrete(expand = c(0,0), breaks = c(0.4,0.6,0.8,1,1.25,1.67,2.5,5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25)) +
                             theme(axis.text.y = element_text(angle = 90, hjust = 0.5)),
                           cowplot::plot_grid(Res_Heatmaps_Sens_4a$p2+ ggtitle("Burial registration rate")+
                                                theme(plot.title = element_text(size=8),
                                                      axis.text = element_text(size = 6),
                                                      axis.title = element_text(size = 6),
                                                      plot.title.position = "plot",
                                                      axis.text.y = element_text(angle = 90, hjust = 0.5)) +
                                                scale_x_discrete(expand = c(0,0), breaks = c(0.4,0.6,0.8,1,1.25,1.67,2.5,5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25)) +
                                                coord_cartesian(xlim = c(1.5,9.5), ylim = c(3.5,6.5), expand = F),
                                              Res_Heatmaps_Sens_4a$p3+ ggtitle("C19 post-mortem prev.")+
                                                theme(plot.title = element_text(size=8),
                                                      axis.text = element_text(size = 6),
                                                      axis.title = element_text(size = 6),
                                                      plot.title.position = "plot",
                                                      axis.text.y = element_text(angle = 90, hjust = 0.5)) +
                                                scale_x_discrete(expand = c(0,0), breaks = c(0.4,0.6,0.8,1,1.25,1.67,2.5,5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25)) +
                                                coord_cartesian(xlim = c(1.5,9.5), ylim = c(3.5,6.5), expand = F),
                                              Res_Heatmaps_Sens_4a$p4+ ggtitle("C19 population surveys")+
                                                theme(plot.title = element_text(size=8),
                                                      axis.text = element_text(size = 6),
                                                      axis.title = element_text(size = 6),
                                                      plot.title.position = "plot",
                                                      axis.text.y = element_text(angle = 90, hjust = 0.5)) +
                                                scale_x_discrete(expand = c(0,0), breaks = c(0.4,0.6,0.8,1,1.25,1.67,2.5,5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25)) +
                                                coord_cartesian(xlim = c(1.5,9.5), ylim = c(3.5,6.5), expand = F), ncol = 1), rel_widths = c(3,1))

tiff("Analysis/Figures/SI_Sens_4a_Change_age_RR_lower_10.tiff", height = 3.5, width = 8, units = "in", res = 300)
sp4a
dev.off()

#######################################
#######################################
##'[Plot Sensitivity heatmap 4b: Change age relative risk: increase by 10%]##
#######################################
#######################################

Res_Sens_4b <- readRDS("Analysis/Results/SI_Sensitivity_Analysis_4b_Change_age_RR_higher_10.rds")
Res_LL_only_Sens_4b <- lapply(X =1:length(Res_Sens_4b), FUN = Diagnostic_Plot_10, fit_Model_list = Res_Sens_4b, IFRvals = IFR_coefficients[as.numeric(gsub(pattern = "X", "", names(Res_Sens_4b))),],
                                Return_likelihoods_only = T)
Res_Heatmaps_Sens_4b <- Plot_Heatmaps(Mod_Res = Res_Sens_4b, Res_Figs = Res_LL_only_Sens_4b, Select_Runs = IFR_coefficients[as.numeric(gsub(pattern = "X", "", names(Res_Sens_4b))),], Title = "")

sp4b <- cowplot::plot_grid(Res_Heatmaps_Sens_4b$p1+ ggtitle("Overall") + coord_cartesian(xlim = c(0.5,8.5), ylim = c(3.5,7.5), expand = F) +
                             scale_x_discrete(expand = c(0,0), breaks = c(0.2,0.4,0.6,0.8,1,1.25,1.67,2.5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25,1.67)) +
                             theme(axis.text.y = element_text(angle = 90, hjust = 0.5)),
                          cowplot::plot_grid(Res_Heatmaps_Sens_4b$p2+ ggtitle("Burial registration rate")+
                                               theme(plot.title = element_text(size=8),
                                                     axis.text = element_text(size = 6),
                                                     axis.title = element_text(size = 6),
                                                     plot.title.position = "plot",
                                                     axis.text.y = element_text(angle = 90, hjust = 0.5)) +
                                               scale_x_discrete(expand = c(0,0), breaks = c(0.2,0.4,0.6,0.8,1,1.25,1.67,2.5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25,1.67)) +
                                               coord_cartesian(xlim = c(0.5,8.5), ylim = c(3.5,7.5), expand = F),
                                             Res_Heatmaps_Sens_4b$p3+ ggtitle("C19 post-mortem prev.")+
                                               theme(plot.title = element_text(size=8),
                                                     axis.text = element_text(size = 6),
                                                     axis.title = element_text(size = 6),
                                                     plot.title.position = "plot",
                                                     axis.text.y = element_text(angle = 90, hjust = 0.5)) +
                                               scale_x_discrete(expand = c(0,0), breaks = c(0.2,0.4,0.6,0.8,1,1.25,1.67,2.5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25,1.67)) +
                                               coord_cartesian(xlim = c(0.5,8.5), ylim = c(3.5,7.5), expand = F),
                                             Res_Heatmaps_Sens_4b$p4+ ggtitle("C19 population surveys")+
                                               theme(plot.title = element_text(size=8),
                                                     axis.text = element_text(size = 6),
                                                     axis.title = element_text(size = 6),
                                                     plot.title.position = "plot",
                                                     axis.text.y = element_text(angle = 90, hjust = 0.5)) +
                                               scale_x_discrete(expand = c(0,0), breaks = c(0.2,0.4,0.6,0.8,1,1.25,1.67,2.5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25,1.67)) +
                                               coord_cartesian(xlim = c(0.5,8.5), ylim = c(3.5,7.5), expand = F), ncol = 1), rel_widths = c(3,1))

tiff("Analysis/Figures/SI_Sens_4b_Change_age_RR_higher_10.tiff", height = 3.5, width = 7.5, units = "in", res = 300)
sp4b
dev.off()

#######################################
#######################################
##'[Plot Sensitivity heatmap 4c: Change age relative risk: decrease by 20%]##
#######################################
#######################################
Res_Sens_4c <- readRDS("Analysis/Results/SI_Sensitivity_Analysis_4c_Change_age_RR_lower_20.rds")
Res_LL_only_Sens_4c <- lapply(X =1:length(Res_Sens_4c), FUN = Diagnostic_Plot_20, fit_Model_list = Res_Sens_4c, IFRvals = IFR_coefficients[Select_Runs,],
                                 Return_likelihoods_only = T)
Res_Heatmaps_Sens_4c <- Plot_Heatmaps(Mod_Res = Res_Sens_4c, Res_Figs = Res_LL_only_Sens_4c, Select_Runs = Select_Runs, Title = "")

sp4c <- cowplot::plot_grid(Res_Heatmaps_Sens_4c$p1+ ggtitle("Overall") + coord_cartesian(xlim = c(1.5,8.5), ylim = c(3.5,6.5), expand = F) +
                             scale_x_discrete(expand = c(0,0), breaks = c(0.4,0.6,0.8,1,1.25,1.67,2.5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25)) +
                             theme(axis.text.y = element_text(angle = 90, hjust = 0.5)),
                           cowplot::plot_grid(Res_Heatmaps_Sens_4c$p2+ ggtitle("Burial registration rate")+
                                                theme(plot.title = element_text(size=8),
                                                      axis.text = element_text(size = 6),
                                                      axis.title = element_text(size = 6),
                                                      plot.title.position = "plot",
                                                      axis.text.y = element_text(angle = 90, hjust = 0.5)) +
                                                scale_x_discrete(expand = c(0,0), breaks = c(0.4,0.6,0.8,1,1.25,1.67,2.5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25)) +
                                                coord_cartesian(xlim = c(1.5,8.5), ylim = c(3.5,6.5), expand = F),
                                              Res_Heatmaps_Sens_4c$p3+ ggtitle("C19 post-mortem prev.")+
                                                theme(plot.title = element_text(size=8),
                                                      axis.text = element_text(size = 6),
                                                      axis.title = element_text(size = 6),
                                                      plot.title.position = "plot",
                                                      axis.text.y = element_text(angle = 90, hjust = 0.5)) +
                                                scale_x_discrete(expand = c(0,0), breaks = c(0.4,0.6,0.8,1,1.25,1.67,2.5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25)) +
                                                coord_cartesian(xlim = c(1.5,8.5), ylim = c(3.5,6.5), expand = F),
                                              Res_Heatmaps_Sens_4c$p4+ ggtitle("C19 population surveys")+
                                                theme(plot.title = element_text(size=8),
                                                      axis.text = element_text(size = 6),
                                                      axis.title = element_text(size = 6),
                                                      plot.title.position = "plot",
                                                      axis.text.y = element_text(angle = 90, hjust = 0.5)) +
                                                scale_x_discrete(expand = c(0,0), breaks = c(0.4,0.6,0.8,1,1.25,1.67,2.5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25)) +
                                                coord_cartesian(xlim = c(1.5,8.5), ylim = c(3.5,6.5), expand = F), ncol = 1), rel_widths = c(3,1))

tiff("Analysis/Figures/SI_Sens_4c_Change_age_RR_lower_20.tiff", height = 3.5, width = 6, units = "in", res = 300)
sp4c
dev.off()

#######################################
#######################################
##'[Plot Sensitivity heatmap 4d: Change age relative risk: increase by 20%]##
#######################################
#######################################
Res_Sens_4d <- readRDS("Analysis/Results/SI_Sensitivity_Analysis_4d_Change_age_RR_higher_20.rds")
Res_LL_only_Sens_4d <- lapply(X =1:length(Res_Sens_4d), FUN = Diagnostic_Plot_20, fit_Model_list = Res_Sens_4d, IFRvals = IFR_coefficients[Select_Runs,],
                                 Return_likelihoods_only = T)
Res_Heatmaps_Sens_4d <- Plot_Heatmaps(Mod_Res = Res_Sens_4d, Res_Figs = Res_LL_only_Sens_4d, Select_Runs = Select_Runs, Title = "")

sp4d <- cowplot::plot_grid(Res_Heatmaps_Sens_4d$p1+ ggtitle("Overall") + coord_cartesian(xlim = c(1.5,8.5), ylim = c(3.5,6.5), expand = F) +
                             scale_x_discrete(expand = c(0,0), breaks = c(0.4,0.6,0.8,1,1.25,1.67,2.5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25)) +
                             theme(axis.text.y = element_text(angle = 90, hjust = 0.5)),
                           cowplot::plot_grid(Res_Heatmaps_Sens_4d$p2+ ggtitle("Burial registration rate")+
                                                theme(plot.title = element_text(size=8),
                                                      axis.text = element_text(size = 6),
                                                      axis.title = element_text(size = 6),
                                                      plot.title.position = "plot",
                                                      axis.text.y = element_text(angle = 90, hjust = 0.5)) +
                                                scale_x_discrete(expand = c(0,0), breaks = c(0.4,0.6,0.8,1,1.25,1.67,2.5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25)) +
                                                coord_cartesian(xlim = c(1.5,8.5), ylim = c(3.5,6.5), expand = F),
                                              Res_Heatmaps_Sens_4d$p3+ ggtitle("C19 post-mortem prevalence")+
                                                theme(plot.title = element_text(size=8),
                                                      axis.text = element_text(size = 6),
                                                      axis.title = element_text(size = 6),
                                                      plot.title.position = "plot",
                                                      axis.text.y = element_text(angle = 90, hjust = 0.5)) +
                                                scale_x_discrete(expand = c(0,0), breaks = c(0.4,0.6,0.8,1,1.25,1.67,2.5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25)) +
                                                coord_cartesian(xlim = c(1.5,8.5), ylim = c(3.5,6.5), expand = F),
                                              Res_Heatmaps_Sens_4d$p4+ ggtitle("C19 population surveys")+
                                                theme(plot.title = element_text(size=8),
                                                      axis.text = element_text(size = 6),
                                                      axis.title = element_text(size = 6),
                                                      plot.title.position = "plot",
                                                      axis.text.y = element_text(angle = 90, hjust = 0.5)) +
                                                scale_x_discrete(expand = c(0,0), breaks = c(0.4,0.6,0.8,1,1.25,1.67,2.5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25)) +
                                                coord_cartesian(xlim = c(1.5,8.5), ylim = c(3.5,6.5), expand = F), ncol = 1), rel_widths = c(3,1))

tiff("Analysis/Figures/SI_Sens_4d_Change_age_RR_higher_20.tiff", height = 3.5, width = 6, units = "in", res = 300)
sp4d
dev.off()


#######################################
#######################################
##'[Plot Sensitivity heatmap 5a: longer duration until death (20%)]##
#######################################
#######################################

Res_Sens_5a <- readRDS("Analysis/Results/SI_Sensitivity_Analysis_5a_Change_duration_to_death_lower.rds")
Res_LL_only_Sens_5a <- lapply(X =1:length(Res_Sens_5a), FUN = Diagnostic_Plot_10, fit_Model_list = Res_Sens_5a, IFRvals = IFR_coefficients[Select_Runs,],
                                 Return_likelihoods_only = T)
Res_Heatmaps_Sens_5a <- Plot_Heatmaps(Mod_Res = Res_Sens_5a, Res_Figs = Res_LL_only_Sens_5a, Select_Runs = Select_Runs, Title = "")

sp5a <- cowplot::plot_grid(Res_Heatmaps_Sens_5a$p1+ ggtitle("Overall") + coord_cartesian(xlim = c(1.5,8.5), ylim = c(3.5,6.5), expand = F) +
                             scale_x_discrete(expand = c(0,0), breaks = c(0.4,0.6,0.8,1,1.25,1.67,2.5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25)) +
                             theme(axis.text.y = element_text(angle = 90, hjust = 0.5)),
                           cowplot::plot_grid(Res_Heatmaps_Sens_5a$p2+ ggtitle("Burial registration rate")+
                                                theme(plot.title = element_text(size=8),
                                                      axis.text = element_text(size = 6),
                                                      axis.title = element_text(size = 6),
                                                      plot.title.position = "plot",
                                                      axis.text.y = element_text(angle = 90, hjust = 0.5)) +
                                                scale_x_discrete(expand = c(0,0), breaks = c(0.4,0.6,0.8,1,1.25,1.67,2.5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25)) +
                                                coord_cartesian(xlim = c(1.5,8.5), ylim = c(3.5,6.5), expand = F),
                                              Res_Heatmaps_Sens_5a$p3+ ggtitle("C19 post-mortem prevalence")+
                                                theme(plot.title = element_text(size=8),
                                                      axis.text = element_text(size = 6),
                                                      axis.title = element_text(size = 6),
                                                      plot.title.position = "plot",
                                                      axis.text.y = element_text(angle = 90, hjust = 0.5)) +
                                                scale_x_discrete(expand = c(0,0), breaks = c(0.4,0.6,0.8,1,1.25,1.67,2.5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25)) +
                                                coord_cartesian(xlim = c(1.5,8.5), ylim = c(3.5,6.5), expand = F),
                                              Res_Heatmaps_Sens_5a$p4+ ggtitle("C19 population surveys")+
                                                theme(plot.title = element_text(size=8),
                                                      axis.text = element_text(size = 6),
                                                      axis.title = element_text(size = 6),
                                                      plot.title.position = "plot",
                                                      axis.text.y = element_text(angle = 90, hjust = 0.5)) +
                                                scale_x_discrete(expand = c(0,0), breaks = c(0.4,0.6,0.8,1,1.25,1.67,2.5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25)) +
                                                coord_cartesian(xlim = c(1.5,8.5), ylim = c(3.5,6.5), expand = F), ncol = 1), rel_widths = c(3,1))

pdf("Analysis/Figures/SI_Sens_5a_Change_duration_to_death_lower.pdf", height = 3, width = 6)
sp5a
dev.off()
tiff("Analysis/Figures/SI_Sens_5a_Change_duration_to_death_lower.tiff", height = 3, width = 6, units = "in", res = 300)
sp5a
dev.off()



#######################################
#######################################
##'[Plot Sensitivity heatmap 5b: longer duration until death (100%)]##
#######################################
#######################################

Res_Sens_5b <- readRDS("Analysis/Results/SI_Sensitivity_Analysis_5b_Change_duration_to_death_higher.rds")
Res_LL_only_Sens_5b <- lapply(X =1:length(Res_Sens_5b), FUN = Diagnostic_Plot_10, fit_Model_list = Res_Sens_5b, IFRvals = IFR_coefficients[Select_Runs,],
                                 Return_likelihoods_only = T)
Res_Heatmaps_Sens_5b <- Plot_Heatmaps(Mod_Res = Res_Sens_5b, Res_Figs = Res_LL_only_Sens_5b, Select_Runs = Select_Runs, Title = "")

sp5b <- cowplot::plot_grid(Res_Heatmaps_Sens_5b$p1+ ggtitle("Overall") + coord_cartesian(xlim = c(1.5,8.5), ylim = c(3.5,6.5), expand = F) +
                             scale_x_discrete(expand = c(0,0), breaks = c(0.4,0.6,0.8,1,1.25,1.67,2.5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25)) +
                             theme(axis.text.y = element_text(angle = 90, hjust = 0.5)),
                           cowplot::plot_grid(Res_Heatmaps_Sens_5b$p2+ ggtitle("Burial registration rate")+
                                                theme(plot.title = element_text(size=8),
                                                      axis.text = element_text(size = 6),
                                                      axis.title = element_text(size = 6),
                                                      plot.title.position = "plot",
                                                      axis.text.y = element_text(angle = 90, hjust = 0.5)) +
                                                scale_x_discrete(expand = c(0,0), breaks = c(0.4,0.6,0.8,1,1.25,1.67,2.5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25)) +
                                                coord_cartesian(xlim = c(1.5,8.5), ylim = c(3.5,6.5), expand = F),
                                              Res_Heatmaps_Sens_5b$p3+ ggtitle("C19 post-mortem prevalence")+
                                                theme(plot.title = element_text(size=8),
                                                      axis.text = element_text(size = 6),
                                                      axis.title = element_text(size = 6),
                                                      plot.title.position = "plot",
                                                      axis.text.y = element_text(angle = 90, hjust = 0.5)) +
                                                scale_x_discrete(expand = c(0,0), breaks = c(0.4,0.6,0.8,1,1.25,1.67,2.5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25)) +
                                                coord_cartesian(xlim = c(1.5,8.5), ylim = c(3.5,6.5), expand = F),
                                              Res_Heatmaps_Sens_5b$p4+ ggtitle("C19 population surveys")+
                                                theme(plot.title = element_text(size=8),
                                                      axis.text = element_text(size = 6),
                                                      axis.title = element_text(size = 6),
                                                      plot.title.position = "plot",
                                                      axis.text.y = element_text(angle = 90, hjust = 0.5)) +
                                                scale_x_discrete(expand = c(0,0), breaks = c(0.4,0.6,0.8,1,1.25,1.67,2.5)) + scale_y_discrete(expand = c(0,0), breaks = c(0.8,1,1.25)) +
                                                coord_cartesian(xlim = c(1.5,8.5), ylim = c(3.5,6.5), expand = F), ncol = 1), rel_widths = c(3,1))

pdf("Analysis/Figures/SI_Sens_5b_Change_duration_to_death_higher.pdf", height = 3, width = 6)
sp5b
dev.off()
tiff("Analysis/Figures/SI_Sens_5b_Change_duration_to_death_higher.tiff", height = 3, width = 6, units = "in", res = 300)
sp5b
dev.off()
