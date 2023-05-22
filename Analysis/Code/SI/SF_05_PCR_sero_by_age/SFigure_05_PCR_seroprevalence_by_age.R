####################################
###### Supplementary Figure 5 ######
####################################
library(dplyr);library(tidyr);library(reshape2);library(ggplot2)
# PCR and seroprevalence by age
fit_Model <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/Full_Set_Default_Setting.rds")$X41
# fit_Model <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/01_Full_fit_round_1.rds")$X41
index <- squire:::odin_index(fit_Model$model)
pars_obs <- fit_Model$pmcmc_results$inputs$pars_obs
pcr_det <- pars_obs$pcr_det
Age_gr_conversion <- data.frame(Age_gr = 1:17, Age_gr_2 = c(1,2,2,3,3,4,4,4,5,5,6,6,7,7,8,8,8))


pcr_sero_data <- fit_Model$output[,c(index$S, index$E2),] %>%
  melt(varnames = c("date","var","Replicate"), value.name = "value") %>%
  mutate(Age_gr = as.numeric(gsub(".*?([0-9]+).", '\\1', var)),
         var = substr(var, 1, 1),
         date = as.Date(date)) %>%
  dcast(... ~ var, value.var="value") %>%
  replace_na(list(E = 0)) %>%
  group_by(Age_gr) %>% mutate(S = as.integer(ifelse(is.na(S), max(S, na.rm = T),S))) %>% ungroup() %>%
  dplyr::rename(Sus = "S", Exp = "E") %>%
  merge(Age_gr_conversion) %>%
  group_by(date, Replicate, Age_gr_2) %>% summarise(Exp = sum(Exp), Sus = sum(Sus)) %>% ungroup() %>%
  group_by(Replicate, Age_gr_2) %>%
  mutate(cum_infs = max(Sus, na.rm = T)-Sus,
         attack_rate = cum_infs/sum(max(Sus)),

         infs = c(0, diff(max(Sus, na.rm = T)-Sus)),
         pcr_pos = cma:::roll_func_10(infs, pcr_det),
         pcr_perc = pcr_pos/max(Sus,na.rm = T),

         Symps = as.integer(Exp) * fit_Model$parameters$gamma_E,
         sero_pos = cma:::roll_func_10(Symps, pars_obs$sero_det),
         sero_perc = sero_pos/max(Sus,na.rm = T)) %>%
  group_by(date, Age_gr_2) %>%
  summarise_at(c("attack_rate", "pcr_perc","sero_perc"), list(median = median, ci = bayestestR::ci), na.rm = TRUE)


Prev_data <- readRDS("Analysis/Data/derived_data/Population_prevalence_by_age.rds")
Age_groups.labs <- c(paste0(c(0,5,15,25,40,50,60),"-",c(4,14,24,39,49,59,69)), "70+")
names(Age_groups.labs) <- 1:8

PCR_by_age_plot <- ggplot(Prev_data$pcr_by_age, aes(group= as.factor(Age_gr_2))) +
  facet_wrap(~Age_gr_2, nrow = 1, labeller = labeller(Age_gr_2 = Age_groups.labs)) +
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme(legend.position = "bottom", axis.text.x = element_blank(),axis.title.x = element_blank(), plot.title = element_text(size = 7, face = "bold")) +
  geom_point(aes(x= as.Date(pars_obs$pcr_df$date_start) + 0.5*(as.Date(pars_obs$pcr_df$date_end)-as.Date(pars_obs$pcr_df$date_start)), y=Perc, color = "PCR %")) +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper, x=as.Date(pars_obs$pcr_df$date_start) + 0.5*(as.Date(pars_obs$pcr_df$date_end)-as.Date(pars_obs$pcr_df$date_start)),
                    width=10,color = "PCR %")) +
  geom_errorbarh(aes(xmin=as.Date(pars_obs$pcr_df$date_start),xmax=as.Date(pars_obs$pcr_df$date_end),y=Perc, height=0,color = "PCR %")) +
  geom_ribbon(data = pcr_sero_data, aes(x=date,ymin=pcr_perc_ci$CI_low, ymax=pcr_perc_ci$CI_high), alpha=0.3, fill = "darkgoldenrod1")+
  geom_line(data = pcr_sero_data, aes(x=date, y=pcr_perc_median, color = "Modelled PCR %"), linetype=2) +
  labs(y = "PCR prevalence", x = "Date", title = "a") +
  coord_cartesian(xlim = c(as.Date("2020-04-15"), as.Date("2020-10-01")), ylim = c(0,0.4)) +
  scale_color_manual(name=NULL,
                     breaks = c("Modelled attack rate", "Modelled PCR %", "Modelled Sero %","PCR %","Sero %"),
                     values = c("black","darkgoldenrod2","chartreuse4","darkgoldenrod2","chartreuse4")) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_date(date_breaks = "2 months", date_labels =  "%b") +
  guides(color = guide_legend(override.aes = list(linetype = c(2,0),
                                                  shape = c(NA,19))))



Sero_by_age_plot <- ggplot(Prev_data$sero_by_age, aes(group= as.factor(Age_gr_2))) +
  facet_wrap(~Age_gr_2, nrow = 1, labeller = labeller(Age_gr_2 = Age_groups.labs)) +
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme(legend.position = "bottom", strip.text.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1),axis.title.x = element_blank(), plot.title = element_text(size = 7, face = "bold")) +
  geom_point(aes(x= as.Date(pars_obs$sero_df$date_start) + 0.5*(as.Date(pars_obs$sero_df$date_end)-as.Date(pars_obs$sero_df$date_start)), y=Perc, color = "Sero %")) +
  geom_errorbar(aes(ymin=CI_lower, ymax=CI_upper, x=as.Date(pars_obs$sero_df$date_start) + 0.5*(as.Date(pars_obs$sero_df$date_end)-as.Date(pars_obs$sero_df$date_start)),
                    width=10,color = "Sero %")) +
  geom_errorbarh(aes(xmin=as.Date(pars_obs$sero_df$date_start),xmax=as.Date(pars_obs$sero_df$date_end),y=Perc, height=0,color = "Sero %")) +
  geom_ribbon(data = pcr_sero_data, aes(x=date,ymin=sero_perc_ci$CI_low, ymax=sero_perc_ci$CI_high), alpha=0.3, fill = "chartreuse4")+
  geom_line(data = pcr_sero_data, aes(x=date, y=sero_perc_median, color = "Modelled Sero %"), linetype=2) +
  labs(y = "Seroprevalence", x = "Date", title = "b") +
  coord_cartesian(xlim = c(as.Date("2020-04-15"), as.Date("2020-10-01")), ylim = c(0, 0.4)) +
  scale_color_manual(name=NULL,
                     breaks = c("Modelled attack rate", "Modelled PCR %", "Modelled Sero %","PCR %","Sero %"),
                     values = c("black","darkgoldenrod2","chartreuse4","darkgoldenrod2","chartreuse4")) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_date(date_breaks = "2 months", date_labels =  "%b") +
  guides(color = guide_legend(override.aes = list(linetype = c(2,0),
                                                  shape = c(NA,19))))

ggsave(filename = "Analysis/Figures/Supplementary_Figures/SFigure_5_Population_prevalence_by_age.pdf",
       cowplot::plot_grid(
         cowplot::plot_grid(PCR_by_age_plot + theme(legend.position = "none"), Sero_by_age_plot  + theme(legend.position = "none"), nrow = 2),
         cowplot::plot_grid(ggpubr::get_legend(PCR_by_age_plot), ggpubr::get_legend(Sero_by_age_plot), nrow = 1), nrow = 2, rel_heights = c(1, 0.1)),
       width = 180, height = 180/2, units = "mm")

tiff(filename = "Analysis/Figures/Supplementary_Figures/SFigure_5_Population_prevalence_by_age.tiff",
     width = 180, height = 180/2, units = "mm", res = 300)
cowplot::plot_grid(
  cowplot::plot_grid(PCR_by_age_plot + theme(legend.position = "none"), Sero_by_age_plot  + theme(legend.position = "none"), nrow = 2),
  cowplot::plot_grid(ggpubr::get_legend(PCR_by_age_plot), ggpubr::get_legend(Sero_by_age_plot), nrow = 1), nrow = 2, rel_heights = c(1, 0.1))
dev.off()
