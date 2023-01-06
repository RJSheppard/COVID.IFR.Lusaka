###################################################
###################################################
##'[Plot Rt trends and attack rate for best fits]##
###################################################
###################################################
library(dplyr);library(ggplot2);library(ggh4x)
devtools::load_all()

### Sample plots
IFR_coefficients <- readRDS("Analysis/Data/derived_data/18_IFR_matrix_coefficients_log_scale_new_pop_str_ests.rds") %>%
  mutate(Index = 1:nrow(.))

Select_Runs <- IFR_coefficients %>%
  filter(Slope_x == 1 & round(IFR_x,2) %in% c(0.6,0.8,1,1.25) |
           Slope_x == 1.25 & round(IFR_x,2) %in% c(0.6)) %>%
  pull(Index)

Res_10_supp <- readRDS("Analysis/Results/Model_Fit_Default_Assumptions.rds")[Select_Runs]
Res_10_Rt_supp_plot <- lapply(X =1:length(Res_10_supp), FUN = Diagnostic_Plot_10, fit_Model_list = Res_10_supp, IFRvals = IFR_coefficients[Select_Runs,],
                              Return_likelihoods_only = F, get_data = T)

Sample_Rt_Plots <- Plot_Rt_Samples(Res_10_Rt_supp_plot, pcr_sero_data, Select_Runs)

tiff("analysis/figures/39_Sample_Best_Fits.tiff", height = 3, width = 10, units = "in", res = 300)
cowplot::plot_grid(cowplot::plot_grid(Sample_Rt_Plots$ps_1 + theme(legend.position = "none"),
                                      Sample_Rt_Plots$Rt_plot + theme(legend.position = "none") + ggtitle("B") + xlab("Date"),
                                      Sample_Rt_Plots$Reff_plot + theme(legend.position = "none") + ggtitle("C") + xlab("Date"),
                                      ncol = 3, align = "hv"),
                   ggpubr::get_legend(Sample_Rt_Plots$ps_1), ncol = 2, rel_widths = c(0.85,0.15))
dev.off()


#########################################################
#########################################################
##'[Plot Samples of model fits varying IFR assumptions]##
#########################################################
#########################################################
library(dplyr);library(ggplot2);library(ggh4x)
devtools::load_all()

### Sample plots
IFR_coefficients <- readRDS("analysis/data/Code-generated-data/00_03_IFR_matrix_coefficients_log_scale.rds") %>%
  mutate(Index = 1:nrow(.))

Select_Runs <- IFR_coefficients %>%
  filter(Slope_x == 1 & round(IFR_x,2) %in% c(0.4,1,2.5) |
           Slope_x == 0.4 & round(IFR_x,2) %in% c(0.4,1,2.5) |
           Slope_x == 2.5 & round(IFR_x,2) %in% c(0.4,1,2.5)) %>%
  pull(Index)

Res_10_supp <- readRDS("Analysis/Results/Model_Fit_Default_Assumptions.rds")[Select_Runs]
Res_10_Rt_supp_plot <- lapply(X =1:length(Res_10_supp), FUN = Diagnostic_Plot_10, fit_Model_list = Res_10_supp, IFRvals = IFR_coefficients[Select_Runs,],
                              Return_likelihoods_only = F, get_data = T)

pcr_sero_data <- readRDS("Analysis/Data/derived_data/14_Lancet_Data.rds")
Sample_Plots <- Plot_Samples(Res_10_supp_plots, pcr_sero_data, Select_Runs)


tiff("Analysis/Figures/SI_Sample_Fits.tiff", height = 10, width = 9, units = "in", res = 300)
cowplot::plot_grid(cowplot::plot_grid(Sample_Plots$ps_1 + theme(legend.position = "none") + theme(plot.margin = unit(c(0.1, 0.5, 0.1, 0.2), "cm")),
                                      Sample_Plots$ps_2 + theme(legend.position = "none") + theme(plot.margin = unit(c(0.5, 0.5, 0.1, 0.2), "cm")),
                                      Sample_Plots$ps_3 + theme(legend.position = "none") + theme(plot.margin = unit(c(0.5, 0.5, 0.1, 0.2), "cm")),
                                      Sample_Plots$ps_4 + theme(legend.position = "none") + theme(plot.margin = unit(c(0.5, 0.5, 0.1, 0.2), "cm")),
                                      Sample_Plots$ps_5 + theme(legend.position = "none") + theme(plot.margin = unit(c(0.5, 0.5, 0.1, 0.2), "cm")), ncol = 1, align = "v", rel_heights = c(1.05,1,1,1,1)),
                   cowplot::plot_grid(ggpubr::get_legend(Sample_Plots$ps_1),
                                      ggpubr::get_legend(Sample_Plots$ps_5), ncol = 1, rel_heights = c(4,1)), ncol = 2, rel_widths = c(0.85,0.15))
dev.off()


