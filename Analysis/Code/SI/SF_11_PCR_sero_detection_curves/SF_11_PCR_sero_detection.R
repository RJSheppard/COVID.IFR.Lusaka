#####################################
###### Supplementary Figure 11 ######
#####################################
# PCR and seroprevalence detection probabilities
library(ggplot2);library(dplyr)
# panel a
pcr_det <- readRDS("Analysis/Data/derived_data/pcr_det_hall_100.rds")
p1 <- ggplot(data = data.frame("Days_since_inf" = 1:length(pcr_det),
                               pcr_det), aes(y = pcr_det, x = Days_since_inf, linetype = "100%")) + geom_line() +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "a", x = "Days since infection", y = "Detection probability") +
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme(legend.position = c("none"), plot.title = element_text(size = 7, face = "bold"))

# panel b
sero_sens <- 0.9
prob_conversion <-  cumsum(dgamma(0:300,shape = 5, rate = 1/2))/max(cumsum(dgamma(0:300,shape = 5, rate = 1/2)))
sero_det <- cumsum(dweibull(0:300, 3.669807, scale = 143.7046))
sero_det <- prob_conversion-sero_det
sero_det[sero_det < 0] <- 0
sero_det <- sero_det/max(sero_det)*sero_sens  # assumed maximum test sensitivitys

p2 <- ggplot(data = data.frame("Days_since_symptom_onset" = 1:length(sero_det),
                               sero_det), aes(y = sero_det, x = Days_since_symptom_onset)) + geom_line() +
  labs(title = "b", x = "Days since symptom onset", y = "Detection probability") +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme(legend.position = c("none"), plot.title = element_text(size = 7, face = "bold"))

ggsave("Analysis/Figures/Supplementary_Figures/SFigure_11_PCR_seroprevalence_detection_curves.pdf",
       cowplot::plot_grid(p1,p2, nrow = 1),
       height = 180*3/7, width = 180, units = "mm")

tiff("Analysis/Figures/Supplementary_Figures/SFigure_11_PCR_seroprevalence_detection_curves.tiff",
     height = 180*3/7, width = 180, units = "mm", res = 300)
cowplot::plot_grid(p1,p2, nrow = 1)
dev.off()

# tiff("Analysis/Figures/SI_squire_utils_detection_probabilities.tiff", units = "in", res = 300, height = 3, width = 7)
# cowplot::plot_grid(p1,p2, nrow = 1)
# dev.off()

