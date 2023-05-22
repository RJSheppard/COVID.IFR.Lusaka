## get a basic estimate of what the figure 2 might look like:
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)

####################################
##'[Panel A: plot the burial data]##
####################################

Burial_df <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Burial_registrations_by_week_2018_2021_Fig4_age_groups.rds")

Baseline_Excess_Burs <- Burial_df |>
  filter(Week >=as.Date("2018-01-01"), Week <as.Date("2019-12-25")) |>
  group_by(Age_gr) %>%
  summarise(Median_Baseline = median(BurRegs),
            Baseline_ci_low = bayestestR::ci(BurRegs)$CI_low,
            Baseline_ci_high = bayestestR::ci(BurRegs)$CI_high) |>
  merge(Burial_df) |>
  mutate(Excess_burs = BurRegs - Median_Baseline,
         Excess_burs_ci_low = BurRegs - Baseline_ci_low,
         Excess_burs_ci_high = BurRegs - Baseline_ci_high)


# Burs_2018_2019 <- Burial_df %>%
#   filter(Week_st >=as.Date("2018-01-01"),
#          Week_st <as.Date("2019-12-25"))
#
# Burs_2018_2021 <- Burial_df %>%
#   filter(Week_st >=as.Date("2018-01-01"),
#          Week_st <as.Date("2021-06-12"))
#
# Burs_2018_2021_Total <- Burs_2018_2021 %>% group_by(Week_st) %>%
#   summarise(Total_deaths = sum(Total_deaths))
#
# Burs_2018_2021_5plus_Total <- Burs_2018_2021 %>% filter(Age_gr_fig2 != 1) %>%group_by(Week_st) %>%
#   summarise(Total_deaths = sum(Total_deaths))
#
# Burs_2020_2021 <- Burial_df %>%
#   filter(Week_st >=as.Date("2019-12-30"),
#          Week_st <as.Date("2021-06-12"))


## Get median all cause burial registrations for 2018-2019:
# Median_Baseline <- Burs_2018_2019 %>% group_by(Age_gr_fig2) %>%
#   summarise(Median_Baseline = median(Total_deaths),
#             Baseline_ci_low = bayestestR::ci(Total_deaths)$CI_low,
#             Baseline_ci_high = bayestestR::ci(Total_deaths)$CI_high)
#
# Burs_Excess <- Burs_2018_2021 %>%
#   merge(Median_Baseline) %>%
#   mutate(Excess_burs = Total_deaths - Median_Baseline,
#          Excess_burs_ci_low = Total_deaths - Baseline_ci_low,
#          Excess_burs_ci_high = Total_deaths - Baseline_ci_high,
#          Age_gr_fig2 = as.factor(Age_gr_fig2))
#
#
# Median_Baseline_total <- Burs_2018_2019 %>% group_by(Week_st) %>% summarise(Total_deaths = sum(Total_deaths)) %>%
#   summarise(Median_Baseline = median(Total_deaths),
#             Baseline_ci_low = bayestestR::ci(Total_deaths)$CI_low,
#             Baseline_ci_high = bayestestR::ci(Total_deaths)$CI_high)
#
#
# Burs_Excess_total <- Burs_2018_2021 %>%
#   group_by(Week_st) %>% summarise(Total_deaths = sum(Total_deaths)) %>%
#   merge(Median_Baseline_total) %>%
#   mutate(Excess_burs = Total_deaths - Median_Baseline,
#          Excess_burs_ci_low = Total_deaths - Baseline_ci_low,
#          Excess_burs_ci_high = Total_deaths - Baseline_ci_high)

Age_gr_labs <- unique(Baseline_Excess_Burs$Age_gr_labs)
names(Age_gr_labs) <- 1:8

p1 <- ggplot(Baseline_Excess_Burs, aes(x = Week, y = BurRegs, color = Age_gr_labs)) +
  geom_line() +
  geom_hline(aes(yintercept = Median_Baseline), linetype = 2) +
  facet_wrap(~Age_gr, strip.position = "left", ncol = 1, labeller = labeller(Age_gr = Age_gr_labs), scales = "free_y") +
  scale_x_date(expand = c(0,0)) +
  labs(title = "a",
       y = "Burial registrations compared with Jan 2018-Dec 2019 median") +
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme(axis.title.x = element_blank(),
        strip.placement = "outside",
        plot.title = element_text(size = 7, face = "bold")) +
  viridis::scale_color_viridis("Age group", discrete = T, option = "H",
                               labels = c(paste0(c(0,5, 15,25, 40,50,60),"-",c(5, 15,25,40,50,60,70)),"70+"))



#######################################
##'[Calculate the model plot results]##
#######################################

### mcmc results
mcmc <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/Baseline_Mortality_MCMC.rds")
# mcmc <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/2023-03-07_Baseline_Mortality_mcmc.rds")
mcmc_samples <- mcmc$output %>% filter(phase =="sampling")

AG1_mcmc <- mcmc_samples[, c(paste0("Week_rate_0to5_",c(1:180)))]

# Get mortality in <5 for use in scaling factor
Ag1std <- Burial_df %>%
  filter(Age_gr ==1, Week>="2018-01-01") %>%
  rename("Ag1std" = BurRegs) %>%
  mutate(Week_gr = 1:180)

# Get week to date conversion
Dates_df <- Ag1std %>% select(Week, Week_gr)

# Calculate mean 2018-2019 baseline for use in scaling factor
AG1_pre2020_mcmc <- mcmc_samples[, c(paste0("Week_rate_0to5_",1:104))]
U5_Baseline_2018_2019 <- rowMeans(AG1_pre2020_mcmc)
# Apply <5 mortality patterns to baseline and 90% registration rate
WeeklyStandardise <- apply(AG1_mcmc, 2, function(x){x/U5_Baseline_2018_2019})

colnames(WeeklyStandardise) <- 1:180
WeeklyStandardise <- WeeklyStandardise %>%
  reshape2::melt(value.name = "Standard", varnames = c("list_names","Week_gr"))

# Get population structure
Pop_Str <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Population_structure_Lusaka_2020_CDC.rds")

BurRegs_2018_2021_5yr_age_groups <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Burial_registrations_by_week_2018_2021_5yr_age_groups.rds") %>%
  merge(Dates_df)

# Generate baseline registration estimates, compare with actual data and apply scaling factor
Mort_excess_deaths <- lapply(1:17, function(y){
  tmp_Mort_excess <- lapply(1:180, function(x){
    if(y == 1){tmp_df <- AG1_mcmc[,paste0("Week_rate_0to5_",x)]} else {tmp_df <- AG1_mcmc[,paste0("Week_rate_0to5_",x)] * mcmc_samples[,paste0("RR",y)]}
    Excess <- BurRegs_2018_2021_5yr_age_groups %>% filter(Week_gr == x, Age_gr == y) %>% pull(BurRegs) - tmp_df
    Excess_std <- Excess/(WeeklyStandardise %>% filter(Week_gr == x) %>% pull(Standard))
    return(data.frame(Sample = as.integer(rownames(AG1_mcmc)),
                      Week_gr = x,
                      Age_gr = y,
                      Total = tmp_df,
                      Excess = Excess,
                      Excess_std = Excess_std,
                      Pop = Pop_Str[y]))
  })
  tmp_Mort_excess <- do.call(rbind, tmp_Mort_excess)
})

Mort_excess_deaths <- do.call(rbind, Mort_excess_deaths)

# Group age groups
Mort_excess_deaths_grouped <- Mort_excess_deaths %>% mutate(Age_gr_fig4 = case_when(Age_gr == 1 ~ 1,
                                                                                    Age_gr %in% 2:3 ~ 2,
                                                                                    Age_gr %in% 4:5 ~ 3,
                                                                                    Age_gr %in% 6:8 ~ 4,
                                                                                    Age_gr %in% 9:10 ~ 5,
                                                                                    Age_gr %in% 11:12 ~ 6,
                                                                                    Age_gr %in% 13:14 ~ 7,
                                                                                    Age_gr %in% 15:17 ~ 8)) %>%
  group_by(Sample, Week_gr, Age_gr_fig4) %>% summarise_at(c("Total","Excess","Excess_std", "Pop"), sum) %>%
  merge(Dates_df) |>
  rename(Age_gr = Age_gr_fig4)


# saveRDS(Mort_excess_deaths_grouped, "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/Excess_Deaths_for_figure_4.rds")
Mort_excess_deaths_grouped <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/Excess_Deaths_for_figure_4.rds")

# Get median and confidence intervals
Panel_B_C_plot_data <- Mort_excess_deaths_grouped %>% group_by(Week, Age_gr) %>%
  summarise_at(c("Total","Excess","Excess_std","Pop"),funs(median = median(.), CI_low = bayestestR::ci(.)$CI_low, CI_high = bayestestR::ci(.)$CI_high))  %>% merge(Burial_df |> filter(Week >="2018-01-01")) %>%
  mutate(Age_gr = as.factor(Age_gr)) %>%
  select(-Pop_CI_low,-Pop_CI_high) %>% rename(Pop = Pop_median)


#########################
##'[Panel B: Model Fit]##
#########################

p2 <- ggplot(Panel_B_C_plot_data, aes(x = Week, color = Age_gr, fill = Age_gr)) +
  geom_line(aes(y = Total_median)) +
  geom_point(aes(y = BurRegs), col = "black", size = 0.3) +
  geom_ribbon(aes(ymin = Total_CI_low, ymax = Total_CI_high), alpha = 0.4, color = NA) +
  facet_wrap(~Age_gr, ncol = 1, labeller = labeller(Age_gr = Age_gr_labs), scales = "free_y") +
  scale_x_date(expand = c(0,0)) +
  labs(title = "b",
       y = "Burial registrations with model fit and predictions") +
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme(plot.title = element_text(size=7, face = "bold"),
        axis.title.x = element_blank(),
        strip.text.x = element_blank()) +
  viridis::scale_color_viridis("Age group", discrete = T, option = "H",
                               labels = c(paste0(c(0, 5, 15,25, 40,50,60),"-",c(5,15,25,40,50,60,70)),"70+")) +
  viridis::scale_fill_viridis("Age group", discrete = T, option = "H",
                              labels = c(paste0(c(0, 5, 15,25, 40,50,60),"-",c(5,15,25,40,50,60,70)),"70+"))


####################################
##'[Panel C: Excess Registrations per capita]##
####################################

p3 <- ggplot(Panel_B_C_plot_data %>% filter(Week >=as.Date("2019-12-30")), aes(x = Week)) +
  geom_line(aes(y = Excess_median*1000/Pop, color = Age_gr), linetype = 1) +
  geom_ribbon(aes(ymin = Excess_CI_low*1000/Pop, ymax = Excess_CI_high*1000/Pop, fill = Age_gr), alpha = 0.4, color = NA) +
  facet_wrap(~Age_gr, strip.position = "right", ncol = 1, labeller = labeller(Age_gr = Age_gr_labs)) +
  geom_hline(aes(yintercept = 0)) +
  scale_x_date(position = "bottom", breaks = seq(as.Date("2020-01-01"), as.Date("2021-07-01"), by="6 months"), date_labels = "%b-%y", expand = c(0,0)) +
  labs(title = "c",
       y = "Burial registrations - pre-pandemic model predictions, per thousand population") +
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size = 7, face = "bold")) +
  viridis::scale_color_viridis("Age group", discrete = T, option = "H",
                               labels = c(paste0(c(0,5,15,25,40,50,60),"-",c(5,15,25,40,50,60,70)),"70+")) +
  viridis::scale_fill_viridis("Age group", discrete = T, option = "H",
                              labels = c(paste0(c(0,5,15,25,40,50,60),"-",c(5,15,25,40,50,60,70)),"70+"))

######################################
##'[Panel D-E: Scaling total excess]##
######################################

Panel_D_plot_data <- Mort_excess_deaths_grouped %>% filter(Week >="2020-01-01", Age_gr != 1) %>% group_by(Sample, Week) %>%
  summarise_at(c("Excess","Excess_std","Pop"),sum) %>%
  merge(WeeklyStandardise %>% rename(Sample = list_names) %>% merge(Dates_df) %>% filter(Week >="2020-01-01") %>% select(Week, Sample, Standard)) %>%
  group_by(Week) %>%
  summarise_at(c("Excess","Excess_std","Pop","Standard"),funs(median = median(.), CI_low = bayestestR::ci(.)$CI_low, CI_high = bayestestR::ci(.)$CI_high))  %>%
  merge(Burial_df |> filter(Week>="2020-01-01") |> group_by(Week) |> summarise(BurRegs = sum(BurRegs))) %>%
  select(-Pop_CI_low,-Pop_CI_high) %>% rename(Pop = Pop_median) %>%
  filter(Week >=as.Date("2020-01-01")) |>
  mutate(Excess_std_median_90 = Excess_std_median/0.9,
         Excess_std_median_80 = Excess_std_median/0.8) |>
  select(-Excess_std_CI_low, -Excess_std_CI_high)


p4 <- ggplot(Panel_D_plot_data, aes(x = Week)) +
  geom_ribbon(aes(ymin = 1/Standard_CI_low, ymax = 1/Standard_CI_high), alpha = 0.2, linetype = 1, linewidth = 0.01, color = NA) +
  geom_line(aes(y = 1/Standard_median), linetype = 1) +
  geom_hline(yintercept = 1, linetype = 2) +
  scale_x_date(expand = c(0,0), limits = as.Date(c("2019-12-30","2021-07-01")), date_labels = "%m-%Y") +
  coord_cartesian(xlim = as.Date(c("2020-01-01","2021-06-07")), ylim = c(0.5,40)) +
  scale_y_continuous(limits = c(0.5,100), trans = "log10") +
  labs(title = "d",
       y = "Scaling\nfactor") +
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme(axis.text.x =element_blank(), axis.title.x = element_blank(),
        plot.title = element_text(size = 7, face = "bold"))


p5 <- ggplot(Panel_D_plot_data, aes(x = Week)) +
  geom_ribbon(aes(ymin = Excess_CI_low, ymax = Excess_CI_high), alpha = 0.2, linetype = 1, linewidth = 0.01, color = NA) +
  geom_line(aes(y = Excess_median), linetype = 1) +
  geom_line(aes(y = Excess_std_median), linetype = 1, color = "blue", linewidth = 0.3) +
  geom_line(aes(y = Excess_std_median_90), linetype = 2, color = "blue", linewidth = 0.3) +
  geom_line(aes(y = Excess_std_median_80), linetype = 3, color = "blue", linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_x_date(expand = c(0,0), limits = as.Date(c("2019-12-30","2021-07-01")), date_labels = "%b-%Y") +
  coord_cartesian(xlim = as.Date(c("2020-01-01","2021-06-07"))) +
  labs(title = "e",
       y = "Excess burial\nregistrations and mortality") +
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size = 7, face = "bold"))


####################################
##'[Panel F: Cumulative Excess]##
####################################

Panel_F_plot_data <- Mort_excess_deaths_grouped %>% filter(Week >="2020-01-01", Age_gr != 1) %>% group_by(Sample, Week) %>%
  summarise_at(c("Excess","Excess_std","Pop"),sum) %>%
  group_by(Sample) %>%
  mutate_at(c("Excess","Excess_std"),cumsum) %>%
  group_by(Week) %>%
  mutate(Excess_std_90 = Excess_std/0.9,
         Excess_std_80 = Excess_std/0.8) %>%
  summarise_at(c("Excess","Excess_std","Excess_std_90","Excess_std_80"),funs(median = median(.), CI_low = bayestestR::ci(.)$CI_low, CI_high = bayestestR::ci(.)$CI_high))  %>%
  merge(Burial_df |> filter(Week>="2020-01-01") |> group_by(Week) |> summarise(BurRegs = sum(BurRegs)))

Panel_F_plot_data_long <- merge(merge(Panel_F_plot_data %>% select(Week, ends_with("median")) %>% setNames(gsub("_median","",names(.))) %>% pivot_longer(cols = 2:5, names_to = "Line", values_to = "median"),
                                      Panel_F_plot_data %>% select(Week, ends_with("CI_low")) %>% setNames(gsub("_CI_low","",names(.))) %>% pivot_longer(cols = 2:5, names_to = "Line", values_to = "CI_low")),
                                Panel_F_plot_data %>% select(Week, ends_with("CI_high")) %>% setNames(gsub("_CI_high","",names(.))) %>% pivot_longer(cols = 2:5, names_to = "Line", values_to = "CI_high"))


p6 <- ggplot(Panel_F_plot_data_long, aes(x = Week, fill = Line, linetype = Line, group = Line)) +
  geom_ribbon(data = Panel_F_plot_data_long %>% filter(Line =="Excess"), aes(ymin = CI_low, ymax = CI_high), inherit.aes = T, alpha = 0.2, linewidth = 0) +
  geom_line(aes(y = median, color = Line)) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_x_date(expand = c(0,0), limits = as.Date(c("2019-12-30","2021-07-01")), date_labels = "%b %y", date_breaks = "6 months", date_minor_breaks = "3 month") +  # scale_x_date(expand = c(0,0)) +
  coord_cartesian(xlim = as.Date(c("2020-01-01","2021-06-07"))) +
  labs(title = "f",
       y = "Cumulative excess burial\nregistrations and mortality") +
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme(axis.title.x = element_blank(),
        plot.title = element_text(size = 7, face = "bold")) +
  scale_fill_manual(name = "", labels = c("Burial registrations","Mortality (scaled)","Mortality (scaled, 90% capture)","Mortality (scaled, 80% capture)"), values = c("black","white","white","white")) +
  scale_color_manual(name = "", labels = c("Burial registrations","Mortality (scaled)","Mortality (scaled, 90% capture)","Mortality (scaled, 80% capture)"), values = c("black","blue","blue","blue")) +
  scale_linetype_manual(name = "", labels = c("Burial registrations","Mortality (scaled)","Mortality (scaled, 90% capture)","Mortality (scaled, 80% capture)"), values = c(1,1,3,2))


####################################
##'[Panel H: DMVI]##
####################################
## Get overall population IFR for Zambia and Lusaka
IFR_weighted_Zamb <- readRDS("Analysis/Data/derived_data/IFR_Zambia.rds") %>% pull(IFR_Braz)
IFR_weighted_Lus <- readRDS("Analysis/Data/derived_data/IFR_Lusaka.rds")
un_demog_Z <- readRDS("Analysis/Data/derived_data/UN_Demography_Zambia.rds")

# Get WHO excess mortality for Zambia by month and calculate DVWI based on Zambia IFR
WHO_Covid_deaths <- readxl::read_xlsx("Analysis/Data/raw_data/WHO_COVID_Excess_Deaths_EstimatesByCountry.xlsx", sheet = "Country by year and month", skip =12) %>%
  filter(country =="Zambia",year == 2020 | year == 2021 & month %in% 1:6) %>%
  select(year, month, cumul.excess.mean, cumul.excess.low, cumul.excess.high) %>%
  add_row(.before = 1, year = 2020, month = 0, cumul.excess.mean = 0, cumul.excess.low = 0, cumul.excess.high = 0) %>%
  mutate(Date = seq.Date(from = as.Date("2020-01-01"), to = as.Date("2021-07-01"), by = "1 month")-1) %>%
  mutate(Excess_deaths_pc = cumul.excess.mean*100000/(sum(un_demog_Z$pop)*1E3)) %>%
  mutate(DMVI = (cumul.excess.mean/(IFR_weighted_Zamb/100))/(sum(un_demog_Z$pop)*1E3),
         DMVI_low = (cumul.excess.low/(IFR_weighted_Zamb/100))/(sum(un_demog_Z$pop)*1E3),
         DMVI_high = (cumul.excess.high/(IFR_weighted_Zamb/100))/(sum(un_demog_Z$pop)*1E3))

# Calculate DVWI for Lusaka
Panel_G_plot_data_long <-
  Panel_F_plot_data_long %>% group_by(Line, Week) %>%
  mutate_at(c("median","CI_low","CI_high"), function(x)(10000*x/(IFR_weighted_Lus))/sum(Pop_Str))

# plot
p7 <- ggplot(Panel_G_plot_data_long, aes(x = Week, color = Line, fill = Line, linetype = Line)) +
  geom_ribbon(data = WHO_Covid_deaths,
              aes(x = as.Date(Date), ymin = DMVI_low, ymax = DMVI_high, fill = "WHO Zambia", linetype = "WHO Zambia"), alpha = 0.2, inherit.aes = F, linewidth = 0) +
  geom_ribbon(data = Panel_G_plot_data_long %>% filter(Line =="Excess"),
              aes(ymin = CI_low/100, ymax = CI_high/100), alpha = 0.1, linewidth = 0, color = NA) +
  geom_line(data = WHO_Covid_deaths, aes(x = as.Date(Date), y = DMVI, color = "WHO Zambia", linetype = "WHO Zambia"), inherit.aes = F) +
  geom_line(aes(y = median/100), alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_x_date(expand = c(0,0), limits = as.Date(c("2019-12-30","2021-07-01")), date_labels = "%b %y", date_breaks = "6 months", date_minor_breaks = "3 month") +  # scale_x_date(expand = c(0,0)) +
  coord_cartesian(xlim = as.Date(c("2020-01-01","2021-06-07"))) +
  labs(title = "g",
       y = "DVWI")+
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme(legend.box.just = "top",
        axis.title.x = element_blank(),
        plot.title = element_text(size = 7, face = "bold")) +
  scale_fill_manual(name = "", limits = c("Excess","Excess_std","Excess_std_90","Excess_std_80","WHO Zambia"), labels = c("Burial registrations","Mortality (scaled)","Mortality (scaled, 90% capture)","Mortality (scaled, 80% capture)","WHO Zambia"), values = c("black","white","white","white","darkgreen")) +
  scale_color_manual(name = "", limits = c("Excess","Excess_std","Excess_std_90","Excess_std_80","WHO Zambia"), labels = c("Burial registrations","Mortality (scaled)","Mortality (scaled, 90% capture)","Mortality (scaled, 80% capture)","WHO Zambia"), values = c("black","blue","blue","blue","darkgreen")) +
  scale_linetype_manual(name = "", limits = c("Excess","Excess_std","Excess_std_90","Excess_std_80","WHO Zambia"), labels = c("Burial registrations","Mortality (scaled)","Mortality (scaled, 90% capture)","Mortality (scaled, 80% capture)","WHO Zambia"), values = c(1,1,2,3,4)) +
  guides(colour = guide_legend(nrow = 2),
         fill = guide_legend(nrow = 2),
         linetype = guide_legend(nrow = 2))

ggsave("Analysis/Figures/Figure_4_Excess_Mortality.pdf",
       cowplot::plot_grid(cowplot::plot_grid(p1 + theme(legend.position = "none",
                                                        plot.title = element_text(vjust = 1)),
                                             p2 + theme(legend.position = "none",
                                                        plot.title = element_text(vjust = 1)),
                                             p3 + theme(legend.position = "none",
                                                        plot.title = element_text(vjust = 1)), nrow = 1, rel_widths = c(1.05,1,1,1.05)),
                          cowplot::plot_grid(
                            cowplot::plot_grid(p4 + theme(legend.position = "none"),
                                               p5 + theme(legend.position = "none"),
                                               nrow = 2, align = "v", rel_heights = c(1,1.5)),
                            cowplot::plot_grid(
                              cowplot::plot_grid(
                                p6 + theme(legend.position = "none"),
                                p7 + theme(legend.position = "none"), nrow = 1, rel_widths = c(1,1)),
                              ggpubr::get_legend(p7 + theme(legend.position = "top",
                                                            # legend.spacing.x = unit(0.1, "cm"),
                                                            legend.text.align = 0,
                                                            legend.key.size = unit(4, units = "mm"),
                                                            legend.text = element_text(margin = margin(l = 0, r = 5)))), nrow = 2, rel_heights = c(1,0.3))), nrow = 2, rel_heights = c(7,3)),

       width = 180, height = 180, units = "mm")


#######################################################
#######################################################
## Figures for text
# Results for use in text
Excess_Mortality_for_paper <- cbind(Panel_F_plot_data_long %>% filter(Week =="2020-12-28") %>% ungroup() %>%  select(-Week) |>
                                      rename(median_2020 = median, CI_low_2020 = CI_low, CI_high_2020 = CI_high),
                                    Panel_F_plot_data_long %>% filter(Week =="2021-06-07") %>% ungroup() %>% select(-Week, -Line) |>
                                      rename(median_2021 = median, CI_low_2021 = CI_low, CI_high_2021 = CI_high)) %>%
  column_to_rownames(var = "Line") %>% round()

Excess_per_million_for_paper <- round(1000000*Excess_Mortality_for_paper/sum(Pop_Str), 1)

DVIW_for_paper <- cbind(Panel_G_plot_data_long %>% filter(Week =="2020-12-28") %>% ungroup() %>%  select(-Week) %>% mutate(median = round(median/100,3), CI_low = round(CI_low/100,3), CI_high = round(CI_high/100,3)),
      Panel_G_plot_data_long %>% filter(Week =="2021-06-07") %>% ungroup() %>% select(-Week, -Line) %>% mutate(median = round(median/100,3), CI_low = round(CI_low/100,3), CI_high = round(CI_high/100,3)))


##############
Perc_mortality

Baseline_Excess_Burs_2018_2019 <- Baseline_Excess_Burs |>
  filter(Week>=as.Date("2018-01-01"),
         Week<as.Date("2019-12-25")) |>
  ungroup() |>
  group_by(Week) |>
  summarise(BurRegs = sum(BurRegs)) |>
  pull(BurRegs) |>
  median()

Total_Burial_Registrations_2020 <- sum(Baseline_Excess_Burs %>% filter(Week>= "2019-12-30", Week <"2021-01-01") %>% pull(BurRegs))
Total_Burial_Registrations_2020_2021 <- sum(Baseline_Excess_Burs %>% filter(Week>= "2019-12-30") %>% pull(BurRegs))


Percentage_of_registrations_2020_median_1819 <- 100*(Mort_excess_deaths_grouped %>%
                                                  filter(Age_gr != 1, Week_gr >105, Week <"2021-01-01") %>%
                                                  group_by(Sample, Week) %>%
                                                  summarise(Excess = sum(Excess),
                                                            Excess_std = sum(Excess_std)) |>
                                                  ungroup() |>
                                                  group_by(Sample) |>
                                                  summarise(Excess = mean(Excess)/Baseline_Excess_Burs_2018_2019,
                                                            Excess_std = mean(Excess_std)/Baseline_Excess_Burs_2018_2019)) |>
  ungroup() |>
  summarise(median_2020_Excess = median(Excess),
            CI_low_2020_Excess = bayestestR::ci(Excess)$CI_low,
            CI_high_median_2020_Excess = bayestestR::ci(Excess)$CI_high,
            median_2020_Excess_std = median(Excess_std),
            CI_low_2020_Excess_std = bayestestR::ci(Excess_std)$CI_low,
            CI_high_median_2020_Excess_std = bayestestR::ci(Excess_std)$CI_high,
  )


Percentage_of_registrations_2021_median_1819 <- 100*(Mort_excess_deaths_grouped %>%
                                                  filter(Age_gr != 1, Week_gr >105) %>%
                                                  group_by(Sample, Week) %>%
                                                  summarise(Excess = sum(Excess),
                                                            Excess_std = sum(Excess_std)) |>
                                                  ungroup() |>
                                                  group_by(Sample) |>
                                                  summarise(Excess = mean(Excess)/Baseline_Excess_Burs_2018_2019,
                                                            Excess_std = mean(Excess_std)/Baseline_Excess_Burs_2018_2019)) |>
  ungroup() |>
  summarise(median_2021_Excess = median(Excess),
            CI_low_2021_Excess = bayestestR::ci(Excess)$CI_low,
            CI_high_median_2021_Excess = bayestestR::ci(Excess)$CI_high,
            median_2021_Excess_std = median(Excess_std),
            CI_low_2021_Excess_std = bayestestR::ci(Excess_std)$CI_low,
            CI_high_median_2021_Excess_std = bayestestR::ci(Excess_std)$CI_high,
  )
