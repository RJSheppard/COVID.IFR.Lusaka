library(dplyr);library(ggplot2);library(ggh4x);library(reshape2)
devtools::load_all()

### Sample plots
IFR_coefficients <- readRDS("Analysis/Data/derived_data/IFR_matrix.rds") %>%
  mutate(Index = 1:nrow(.))

Select_Runs <- IFR_coefficients %>%
  filter(Slope_x == 1 & round(IFR_x,2) %in% c(0.8,0.8,1,1.25,1.67)) %>%
  pull(Index)

Res_default <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/Full_Set_Default_Setting.rds")
# Res_default <- readRDS("../Bonus Files/2023_March_Reviewer_Batch/01_Full_fit_round_1.rds")
# Res_default_2 <- readRDS("../Bonus Files/2023_March_Reviewer_Batch/01_Full_fit_round_2.rds")
# Res_default <- c(Res_default,Res_default_2)[order(as.numeric(gsub("X", "", names(c(Res_default,Res_default_2)))))]
# Res_default <- Res_default[Select_Runs]
# rm(Res_default_2)
# gc()

# Res_default <- lapply(Res_default, FUN = function(x){
#   x$pmcmc_results$inputs$pars_obs$combined_data$BurRegs <- x$pmcmc_results$inputs$pars_obs$combined_data$Bur_regs
#   x$pmcmc_results$inputs$pars_obs$combined_data$Bur_regs <- NULL
#   return(x)})

Res_default_Rt_supp_plot <- lapply(X =1:length(Res_default), FUN = Diagnostic_Plot, fit_Model_list = Res_default, IFRvals = IFR_coefficients[Select_Runs,],
                              Return_likelihoods_only = F, get_data = T)


devtools::load_all()
Sample_Rt_Plots <- Plot_Rt_Samples(Res_default_Rt_supp_plot, pcr_sero_data, Select_Runs)

ggsave("Analysis/Figures/Supplementary_Figures/SFigure_08_Spread_and_Transmissibility.pdf",
       cowplot::plot_grid(cowplot::plot_grid(Sample_Rt_Plots$ps_1 + theme(legend.position = "none"),
                                             Sample_Rt_Plots$Rt_plot + theme(legend.position = "none"),
                                             Sample_Rt_Plots$Reff_plot + theme(legend.position = "none"),
                                             ncol = 3, align = "hv"),
                          ggpubr::get_legend(Sample_Rt_Plots$ps_1), ncol = 2, rel_widths = c(0.85,0.15)),
       height = 180*3/10, width = 180, units = "mm")

tiff("Analysis/Figures/Supplementary_Figures/SFigure_08_Spread_and_Transmissibility.tiff", height = 180*3/10, width = 180, units = "mm", res = 300)
cowplot::plot_grid(cowplot::plot_grid(Sample_Rt_Plots$ps_1 + theme(legend.position = "none"),
                                      Sample_Rt_Plots$Rt_plot + theme(legend.position = "none"),
                                      Sample_Rt_Plots$Reff_plot + theme(legend.position = "none"),
                                      ncol = 3, align = "hv"),
                   ggpubr::get_legend(Sample_Rt_Plots$ps_1), ncol = 2, rel_widths = c(0.85,0.15))
dev.off()
