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
IFR_mat <- readRDS("Analysis/Data/derived_data/18_IFR_matrix_coefficients_log_scale_new_pop_str_ests.rds")

### Heatmap A
Res_A <- readRDS("Analysis/Results/SI_Sensitivity_Analysis_A_Remove_Age_gr_1.rds")
H_A <- Heatmap_Post(Plot_Post(Res_A, IFR_mat))
# saveRDS(H_A, file = "../Bonus Files/Heatmap_SA_Remove_Age_gr_1.rds")


### Heatmap B
Res_B <- readRDS("Analysis/Results/SI_Sensitivity_Analysis_B_No_scaling_no_weeks_4_or_5.rds")
H_B <- Heatmap_Post(Plot_Post(Res_B, IFR_mat))
# saveRDS(H_B, file = "../Bonus Files/Heatmap_SB_No_scaling_no_weeks_4_5.rds")

### Heatmap C
Res_C <- readRDS("Analysis/Results/SI_Sensitivity_Analysis_C_Change_duration_to_death_lower.rds")
H_C <- Heatmap_Post(Plot_Post(Res_C, IFR_mat))
# saveRDS(H_C, file = "../Bonus Files/Heatmap_SC_lower_duration_until_death")


### Heatmap D
Res_D <- readRDS("Analysis/Results/SI_Sensitivity_Analysis_D_Change_age_RR_higher_20.rds")
H_D <- Heatmap_Post(Plot_Post(Res_D, IFR_mat))
# saveRDS(H_D, file = "../Bonus Files/Heatmap_SD_higher_duration_until_death")


### Heatmap E
Res_E <- readRDS("Analysis/Results/SI_Sensitivity_Analysis_E_Change_age_RR_lower_10.rds")
H_E <- Heatmap_Post(Plot_Post(Res_E, IFR_mat))
saveRDS(H5, file = "Analysis/Results/Heatmap_SE_Change_age_RR_lower_10.rds")


### Heatmap F
Res_F <- readRDS("Analysis/Results/SI_Sensitivity_Analysis_F_Change_age_RR_higher_10.rds")
H_F <- Heatmap_Post(Plot_Post(Res_F, IFR_mat))
saveRDS(H_F, file = "../Bonus Files/Heatmap_S_F_Change_age_RR_higher_10.rds")


### Heatmap G
Res_G <- readRDS("Analysis/Results/SI_Sensitivity_Analysis_G_Change_age_RR_lower_20.rds")
H_G <- Heatmap_Post(Plot_Post(Res_G, IFR_mat))
saveRDS(H_G, file = "../Bonus Files/Heatmap_S_G_Change_age_RR_lower_20.rds")


### Heatmap F
Res_H <- readRDS("Analysis/Results/SI_Sensitivity_Analysis_H_Change_age_RR_higher_20.rds")
H_H <- Heatmap_Post(Plot_Post(Res_H, IFR_mat))
saveRDS(H_H, file = "../Bonus Files/Heatmap_S_H_Change_age_RR_higher_20.rds")


Full_Heatmap_Sensitivities <- cowplot::plot_grid(
  H_A + ggtitle("A"),
  H_B + ggtitle("B"),
  H_C + ggtitle("C"),
  H_D + ggtitle("D"),
  H_E + ggtitle("E"),
  H_F + ggtitle("F"),
  H_G + ggtitle("G"),
  H_H + ggtitle("H"),
  ncol = 2)

tiff("analysis/figures/SI_Sensitivity_analyses.tiff", height = 10, width = 8, units = "in", res = 300)
Full_Heatmap_Sensitivities
dev.off()
