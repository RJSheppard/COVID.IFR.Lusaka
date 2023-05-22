rm(list = ls()); gc()
library(dplyr);library(ggplot2);library(ggh4x);library(reshape2)
devtools::load_all()

### Sample plots
IFR_coefficients <- readRDS("Analysis/Data/derived_data/IFR_matrix.rds") %>%
  mutate(Index = 1:nrow(.))

Select_Runs <- IFR_coefficients %>%
  filter(Slope_x == 1 & round(IFR_x,2) %in% c(0.4,1,2.5) |
           Slope_x == 0.4 & round(IFR_x,2) %in% c(0.4,1,2.5) |
           Slope_x == 2.5 & round(IFR_x,2) %in% c(0.4,1,2.5)) %>%
  pull(Index)

Res_default <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/Full_Set_Default_Setting.rds")
# Res_default <- readRDS("../Bonus Files/2023_March_Reviewer_Batch/01_Full_fit_round_1.rds")
# Res_default_2 <- readRDS("../Bonus Files/2023_March_Reviewer_Batch/01_Full_fit_round_2.rds")
# Res_default <- c(Res_default,Res_default_2)[order(as.numeric(gsub("X", "", names(c(Res_default,Res_default_2)))))]
# Res_supp <- Res_default[Select_Runs]
rm(list =c("Res_default_2","Res_default"))
gc()

# Res_supp <- lapply(Res_supp, FUN = function(x){
#   x$pmcmc_results$inputs$pars_obs$combined_data$BurRegs <- x$pmcmc_results$inputs$pars_obs$combined_data$Bur_regs
#   x$pmcmc_results$inputs$pars_obs$combined_data$Bur_regs <- NULL
#   return(x)})

Res_supp_plots <- lapply(X =1:length(Res_supp), FUN = Diagnostic_Plot, fit_Model_list = Res_supp, IFRvals = IFR_coefficients[Select_Runs,],
                            Return_likelihoods_only = F, get_data = T)


pcr_sero_data <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Population_prevalence.rds")
devtools::load_all()
Sample_Plots <- Plot_Samples(Res_supp_plots, pcr_sero_data, Select_Runs)

Sample_plots_Figure <- cowplot::plot_grid(cowplot::plot_grid(Sample_Plots$ps_1 + theme(legend.position = "none"),
                                                             Sample_Plots$ps_2 + theme(legend.position = "none"),
                                                             Sample_Plots$ps_3 + theme(legend.position = "none"),
                                                             Sample_Plots$ps_4 + theme(legend.position = "none"),
                                                             Sample_Plots$ps_5 + theme(legend.position = "none"), ncol = 1, align = "v", rel_heights = c(1.05,1,1,1,1)),
                                          cowplot::plot_grid(ggpubr::get_legend(Sample_Plots$ps_1),
                                                             ggpubr::get_legend(Sample_Plots$ps_5 + theme(legend.key.size = unit(4, "mm"))), ncol = 1, rel_heights = c(4,1)), ncol = 2, rel_widths = c(0.8,0.2))

ggsave("Analysis/Figures/Supplementary_Figures/SFigure_06_Fit_Varying_Severity.pdf",
       Sample_plots_Figure,
       width = 180, height = 180, units = "mm")

tiff("Analysis/Figures/Supplementary_Figures/SFigure_06_Fit_Varying_Severity.tiff",
     width = 180, height = 180, units = "mm", res = 300)
Sample_plots_Figure
dev.off()

