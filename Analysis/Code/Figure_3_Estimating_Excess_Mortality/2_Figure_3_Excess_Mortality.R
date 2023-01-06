## get a basic estimate of what the figure 2 might look like:
library(dplyr)
library(ggplot2)
library(tidyr)

####################################
##'[Panel A: plot the burial data]##
####################################

Burial_df <- readRDS("Analysis/Data/derived_data/07_Burial_registrations_2017_2021.rds") %>%
  rename(Age_gr_fig2 = Age_gr_fig_2b) %>%
  group_by(Age_gr_fig2, Week_st) %>%
  summarise(Total_deaths = n()) %>%
  ungroup() %>% complete(Age_gr_fig2, Week_st, fill = list(Total_deaths = 0))

Burs_2018_2019 <- Burial_df %>%
  filter(Week_st >=as.Date("2018-01-01"),
         Week_st <as.Date("2019-12-25"))

Burs_2018_2021 <- Burial_df %>%
  filter(Week_st >=as.Date("2018-01-01"),
         Week_st <as.Date("2021-06-12"))

Burs_2018_2021_Total <- Burs_2018_2021 %>% group_by(Week_st) %>%
  summarise(Total_deaths = sum(Total_deaths))

Burs_2018_2021_5plus_Total <- Burs_2018_2021 %>% filter(Age_gr_fig2 != 1) %>%group_by(Week_st) %>%
  summarise(Total_deaths = sum(Total_deaths))

Burs_2020_2021 <- Burial_df %>%
  filter(Week_st >=as.Date("2019-12-30"),
         Week_st <as.Date("2021-06-12"))


## Get median all cause burial registrations for 2018-2019:
Median_Baseline <- Burs_2018_2019 %>% group_by(Age_gr_fig2) %>%
  summarise(Median_Baseline = median(Total_deaths),
            Baseline_ci_low = bayestestR::ci(Total_deaths)$CI_low,
            Baseline_ci_high = bayestestR::ci(Total_deaths)$CI_high)

Burs_Excess <- Burs_2018_2021 %>%
  merge(Median_Baseline) %>%
  mutate(Excess_burs = Total_deaths - Median_Baseline,
         Excess_burs_ci_low = Total_deaths - Baseline_ci_low,
         Excess_burs_ci_high = Total_deaths - Baseline_ci_high,
         Age_gr_fig2 = as.factor(Age_gr_fig2))


Median_Baseline_total <- Burs_2018_2019 %>% group_by(Week_st) %>% summarise(Total_deaths = sum(Total_deaths)) %>%
  summarise(Median_Baseline = median(Total_deaths),
            Baseline_ci_low = bayestestR::ci(Total_deaths)$CI_low,
            Baseline_ci_high = bayestestR::ci(Total_deaths)$CI_high)


Burs_Excess_total <- Burs_2018_2021 %>%
  group_by(Week_st) %>% summarise(Total_deaths = sum(Total_deaths)) %>%
  merge(Median_Baseline_total) %>%
  mutate(Excess_burs = Total_deaths - Median_Baseline,
         Excess_burs_ci_low = Total_deaths - Baseline_ci_low,
         Excess_burs_ci_high = Total_deaths - Baseline_ci_high)

Age_gr_fig2.labs <- c(paste0(c(0,5, 15,25, 40,50,60),"-",c(4, 14,24,39,49,59,69)),"70+")
names(Age_gr_fig2.labs) <- 1:8

p1 <- ggplot(Burs_Excess, aes(x = Week_st, y = Total_deaths, color = Age_gr_fig2)) +
  geom_line() +
  facet_wrap(~Age_gr_fig2, strip.position = "left", ncol = 1, labeller = labeller(Age_gr_fig2 = Age_gr_fig2.labs), scales = "free_y") +
  theme_minimal() +
  xlab("") + ylab("Burial registrations compared with Jan 2018-Dec 2019 median") +
  viridis::scale_color_viridis("Age group", discrete = T, option = "H",
                               labels = c(paste0(c(0,5, 15,25, 40,50,60),"-",c(5, 15,25,40,50,60,70)),"70+")) +
  geom_hline(aes(yintercept = Median_Baseline), linetype = 2) +
  scale_x_date(expand = c(0,0)) +
  ggtitle("A") +
  theme(strip.placement = "outside", strip.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))



#######################################
##'[Calculate the model plot results]##
#######################################

### mcmc results
mcmc <- readRDS("Analysis/Data/derived_data/06_Baseline_Mortality_MCMC_Gamma_Prior_inc_Feb_2021.rds")
mcmc_samples <- mcmc$output %>% filter(phase =="sampling")

AG1_mcmc <- mcmc_samples[, c(paste0("Week_rate_0to5_",c(1:180)))]

# Get mortality in <5 for use in scaling factor
Ag1std <- Burs_2018_2021 %>%
  rename("Ag1std" = Total_deaths) %>%
  filter(Age_gr_fig2 ==1) %>%
  mutate(Week_gr = 1:180)

# Get week to date conversion
Dates_df <- Ag1std %>% select(Week_st, Week_gr)

# Calculate mean 2018-2019 baseline for use in scaling factor
AG1_pre2020_mcmc <- mcmc_samples[, c(paste0("Week_rate_0to5_",1:104))]
U5_Baseline_2018_2019 <- rowMeans(AG1_pre2020_mcmc)
# Apply <5 mortality patterns to baseline and 90% registration rate
WeeklyStandardise <- apply(AG1_mcmc, 2, function(x){0.9*x/U5_Baseline_2018_2019})

colnames(WeeklyStandardise) <- 1:180
WeeklyStandardise <- WeeklyStandardise %>%
  reshape2::melt(value.name = "Standard", varnames = c("list_names","Week_gr"))

# Get population structure
Pop_Str <- readRDS("Analysis/Data/derived_data/08_Lusaka_Dist_Pop_Str_2020_imp_ests.rds")

Burs_2018_2021_all_age_groups <- readRDS("Analysis/Data/derived_data/07_Burial_registrations_2017_2021.rds") %>%
  group_by(Age_gr, Week_st) %>%
  summarise(Total_deaths = n()) %>%
  ungroup() %>% complete(Age_gr, Week_st, fill = list(Total_deaths = 0)) %>%
  filter(Week_st >=as.Date("2018-01-01"),
         Week_st <as.Date("2021-06-12")) %>%
  merge(Dates_df)

# Generate baseline mortality estimates, compare with actual data and apply scaling factor
Mort_excess_deaths <- lapply(1:17, function(y){
  lapply(1:180, function(x){
    if(y == 1){tmp_df <- AG1_mcmc[,paste0("Week_rate_0to5_",x)]} else {tmp_df <- AG1_mcmc[,paste0("Week_rate_0to5_",x)] * mcmc_samples[,paste0("RR",y)]}
    Excess <- Burs_2018_2021_all_age_groups %>% filter(Week_gr == x, Age_gr == y) %>% pull(Total_deaths) - tmp_df
    Excess_std <- Excess/(WeeklyStandardise %>% filter(Week_gr == x) %>% pull(Standard))
    return(data.frame(Sample = as.integer(rownames(AG1_mcmc)), Week_gr = x, Age_gr = y, Total = tmp_df, Excess = Excess, Excess_std = Excess_std, Pop = Pop_Str[y]))
  }) %>% str2str::ld2d() %>% select(-row_names)
}) %>% str2str::ld2d() %>% select(-row_names, -list_names)

# Group age groups
Mort_excess_deaths_fig2_age_gr <- Mort_excess_deaths %>% mutate(Age_gr_fig2 = case_when(Age_gr == 1 ~ 1,
                                                                                        Age_gr %in% 2:3 ~ 2,
                                                                                        Age_gr %in% 4:5 ~ 3,
                                                                                        Age_gr %in% 6:8 ~ 4,
                                                                                        Age_gr %in% 9:10 ~ 5,
                                                                                        Age_gr %in% 11:12 ~ 6,
                                                                                        Age_gr %in% 13:14 ~ 7,
                                                                                        Age_gr %in% 15:17 ~ 8)) %>%
  group_by(Sample, Week_gr, Age_gr_fig2) %>% summarise_at(c("Total","Excess","Excess_std", "Pop"), sum) %>%
  merge(Dates_df)

# Get median and confidence intervals
Panel_B_C_D_plot_data <- Mort_excess_deaths_fig2_age_gr %>% group_by(Week_st, Age_gr_fig2) %>%
  summarise_at(c("Total","Excess","Excess_std","Pop"),funs(median = median(.), CI_low = bayestestR::ci(.)$CI_low, CI_high = bayestestR::ci(.)$CI_high))  %>% merge(Burs_2018_2021) %>%
  mutate(Age_gr_fig2 = as.factor(Age_gr_fig2)) %>%
  select(-Pop_CI_low,-Pop_CI_high) %>% rename(Pop = Pop_median)


#########################
##'[Panel B: Model Fit]##
#########################

p2 <- ggplot(Panel_B_C_D_plot_data, aes(x = Week_st, y = Total_median, color = Age_gr_fig2, fill = Age_gr_fig2)) +
  geom_line() +
  geom_point(aes(y = Total_deaths), col = "black", size = 0.3) +
  facet_wrap(~Age_gr_fig2, ncol = 1, labeller = labeller(Age_gr_fig2 = Age_gr_fig2.labs), scales = "free_y") +
  geom_ribbon(aes(ymin = Total_CI_low, ymax = Total_CI_high), alpha = 0.4, color = NA) +
  theme_minimal() +
  xlab("") + ylab("Burial registrations with model fit and predictions") +
  viridis::scale_color_viridis("Age group", discrete = T, option = "H",
                               labels = c(paste0(c(0, 5, 15,25, 40,50,60),"-",c(5,15,25,40,50,60,70)),"70+")) +
  viridis::scale_fill_viridis("Age group", discrete = T, option = "H",
                              labels = c(paste0(c(0, 5, 15,25, 40,50,60),"-",c(5,15,25,40,50,60,70)),"70+")) +
  scale_x_date(expand = c(0,0)) +
  ggtitle("B") +
  theme(strip.text.x = element_blank(),
        axis.text.x = element_text(size = 10))


####################################
##'[Panel C: Excess Registrations]##
####################################
p3 <- ggplot(Panel_B_C_D_plot_data %>% filter(Week_st >=as.Date("2019-12-30")) %>% mutate(Excess_median = ifelse(Age_gr_fig2==1, 0, Excess_median),
                                                                                          Excess_CI_low = ifelse(Age_gr_fig2==1, 0, Excess_CI_low),
                                                                                          Excess_CI_high = ifelse(Age_gr_fig2==1, 0, Excess_CI_high)), aes(x = Week_st)) +
  geom_line(aes(y = Excess_median, color = Age_gr_fig2), linetype = 1) +
  geom_ribbon(aes(ymin = Excess_CI_low, ymax = Excess_CI_high, fill = Age_gr_fig2), alpha = 0.4, color = NA) +
  facet_wrap(~Age_gr_fig2, ncol = 1, labeller = labeller(Age_gr_fig2 = Age_gr_fig2.labs), scales = "free_y") +
  theme_minimal() +
  xlab("") + ylab("Burial registrations - pre-pandemic model predictions") +
  viridis::scale_color_viridis("Age group", discrete = T, option = "H",
                               labels = c(paste0(c(0,5,15,25,40,50,60),"-",c(5,15,25,40,50,60,70)),"70+")) +
  viridis::scale_fill_viridis("Age group", discrete = T, option = "H",
                              labels = c(paste0(c(0,5,15,25,40,50,60),"-",c(5,15,25,40,50,60,70)),"70+")) +
  scale_y_continuous(limits = c(-30,100)) +
  ggh4x::facetted_pos_scales(y = list(Age_gr_fig2 %in% as.character(1:8) ~ scale_y_continuous(limits = c(-30,100)))) +
  scale_x_date(position = "bottom", breaks = seq(as.Date("2020-01-01"), as.Date("2021-07-01"), by="6 months"), date_labels = "%b-%y", expand = c(0,0)) +
  geom_hline(yintercept = 0, linetype = 2) +
  ggtitle("C") +
  theme(strip.text.x = element_blank(),
        axis.text.x = element_text(size = 10))

####################################
##'[Panel D: Excess Registrations per capita]##
####################################

p4 <- ggplot(Panel_B_C_D_plot_data %>% filter(Week_st >=as.Date("2019-12-30")), aes(x = Week_st)) +
  geom_line(aes(y = Excess_median*1000/Pop, color = Age_gr_fig2), linetype = 1) +
  geom_ribbon(aes(ymin = Excess_CI_low*1000/Pop, ymax = Excess_CI_high*1000/Pop, fill = Age_gr_fig2), alpha = 0.4, color = NA) +
  facet_wrap(~Age_gr_fig2, strip.position = "right", ncol = 1, labeller = labeller(Age_gr_fig2 = Age_gr_fig2.labs)) +
  theme_minimal() +
  xlab("") + ylab("Burial registrations - pre-pandemic model predictions, per capita (/1000)") +
  viridis::scale_color_viridis("Age group", discrete = T, option = "H",
                               labels = c(paste0(c(0,5,15,25,40,50,60),"-",c(5,15,25,40,50,60,70)),"70+")) +
  viridis::scale_fill_viridis("Age group", discrete = T, option = "H",
                              labels = c(paste0(c(0,5,15,25,40,50,60),"-",c(5,15,25,40,50,60,70)),"70+")) +
  scale_x_date(position = "bottom", breaks = seq(as.Date("2020-01-01"), as.Date("2021-07-01"), by="6 months"), date_labels = "%b-%y", expand = c(0,0)) +
  ggtitle("D") +
  geom_hline(aes(yintercept = 0)) +
  theme(strip.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10))

######################################
##'[Panel E-F: Scaling total excess]##
######################################

Panel_E_plot_data <- Mort_excess_deaths_fig2_age_gr %>% filter(Week_st >="2020-01-01", Week_st <"2022-07-01") %>% group_by(Sample, Week_st) %>%
  summarise_at(c("Excess","Excess_std","Pop"),sum) %>%
  merge(WeeklyStandardise %>% rename(Sample = list_names) %>% merge(Dates_df) %>% filter(Week_st >="2020-01-01") %>% select(Week_st, Sample, Standard)) %>%
  group_by(Week_st) %>%
  summarise_at(c("Excess","Excess_std","Pop","Standard"),funs(median = median(.), CI_low = bayestestR::ci(.)$CI_low, CI_high = bayestestR::ci(.)$CI_high))  %>% merge(Burs_2018_2021_Total) %>%
  select(-Pop_CI_low,-Pop_CI_high) %>% rename(Pop = Pop_median) %>%
  filter(Week_st >=as.Date("2020-01-01"))

Panel_E_plot_data_long <- merge(merge(Panel_E_plot_data %>% select(Week_st, ends_with("median")) %>% setNames(gsub("_median","",names(.))) %>% pivot_longer(cols = 2:4, names_to = "Line", values_to = "median"),
                                      Panel_E_plot_data %>% select(Week_st, ends_with("CI_low")) %>% setNames(gsub("_CI_low","",names(.))) %>% pivot_longer(cols = 2:4, names_to = "Line", values_to = "CI_low")),
                                Panel_E_plot_data %>% select(Week_st, ends_with("CI_high")) %>% setNames(gsub("_CI_high","",names(.))) %>% pivot_longer(cols = 2:4, names_to = "Line", values_to = "CI_high")) %>%
  mutate(Plot_order = case_when(Line == "Excess" ~ 1,
                                Line == "Excess_std" ~ 2,
                                Line == "Standard" ~ 3)) %>%
  mutate(Plot_order = as.factor(Plot_order))


p5 <- ggplot(Panel_E_plot_data_long %>% filter(Line == "Standard"), aes(x = Week_st, color = Plot_order, fill = Plot_order, alpha = Line)) +
  geom_ribbon(aes(ymin = 1/CI_low, ymax = 1/CI_high), alpha = 0.2, linetype = 1, linewidth = 0.01, color = NA) +
  geom_line(aes(y = 1/median), linewidth = 0.8, linetype = 1) +
  scale_y_continuous(name = "Scaling\nfactor",limits = c(0.5,100), trans = "log10") +
  theme_minimal() +
  xlab("") +
  geom_hline(yintercept = 1, linetype = 2) +
  scale_x_date(expand = c(0,0), limits = as.Date(c("2019-12-30","2021-07-01")), date_labels = "%m-%Y") +
  ggtitle("E") +
  scale_fill_manual(name = "", labels = c("Scaling\nfactor"), values = c("black")) +
  scale_color_manual(name = "", labels = c("Scaling\nfactor"), values = c("black")) +
  scale_alpha_manual(name = "", labels = c("Scaling\nfactor"), values = c(1)) +
  theme(axis.text.x =element_blank()) +
  coord_cartesian(xlim = as.Date(c("2020-01-01","2021-06-07")), ylim = c(0.5,40))


p6 <- ggplot(Panel_E_plot_data_long %>% filter(Line != "Standard"), aes(x = Week_st, color = Plot_order, fill = Plot_order, alpha = Line)) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.2, linetype = 1, linewidth = 0.01, color = NA) +
  geom_line(aes(y = median), linetype = 1, linewidth = 0.8) +
  scale_y_continuous(name = "Excess burial\nregistrations/mortality") +
  theme_minimal() +
  xlab("") +
  coord_cartesian(ylim = c(-200,500), xlim = as.Date(c("2020-01-01","2021-06-07"))) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_x_date(expand = c(0,0), limits = as.Date(c("2019-12-30","2021-07-01")), date_labels = "%b-%Y") +
  ggtitle("F") +
  scale_fill_manual(name = "", labels = c("Burial registrations","Mortality","Scaling factor"), values = c("black","darkblue","darkred")) +
  scale_color_manual(name = "", labels = c("Burial registrations","Mortality","Scaling factor"), values = c("black","darkblue","darkred")) +
  scale_alpha_manual(name = "", labels = c("Burial registrations","Mortality","Scaling factor"), values = c(0.8,0.8,1)) +
  theme(axis.text.x = element_text(size = 10))


####################################
##'[Panel G: Cumulative Excess]##
####################################

Panel_G_plot_data <- Mort_excess_deaths_fig2_age_gr %>% filter(Week_st >="2020-01-01") %>% group_by(Sample, Week_st) %>%
  summarise_at(c("Excess","Excess_std","Pop"),sum) %>%
  group_by(Sample) %>%
  mutate_at(c("Excess","Excess_std"),cumsum) %>%
  group_by(Week_st) %>%
  summarise_at(c("Excess","Excess_std"),funs(median = median(.), CI_low = bayestestR::ci(.)$CI_low, CI_high = bayestestR::ci(.)$CI_high))  %>% merge(Burs_2018_2021_Total)

Panel_G_plot_data_long <- merge(merge(Panel_G_plot_data %>% select(Week_st, ends_with("median")) %>% setNames(gsub("_median","",names(.))) %>% pivot_longer(cols = 2:3, names_to = "Line", values_to = "median"),
                                      Panel_G_plot_data %>% select(Week_st, ends_with("CI_low")) %>% setNames(gsub("_CI_low","",names(.))) %>% pivot_longer(cols = 2:3, names_to = "Line", values_to = "CI_low")),
                                Panel_G_plot_data %>% select(Week_st, ends_with("CI_high")) %>% setNames(gsub("_CI_high","",names(.))) %>% pivot_longer(cols = 2:3, names_to = "Line", values_to = "CI_high"))

Pandemic_Excess_Deaths <- Panel_G_plot_data_long %>% filter(Week_st =="2020-05-25" | Week_st == "2020-09-28", Line == "Excess_std") %>% select(median, CI_low, CI_high)

############################
# Results for use in text
Pandemic_Excess_Deaths <- Panel_G_plot_data_long %>% filter(Week_st =="2020-12-21" | Week_st =="2021-06-07")
Pandemic_Excess_Deaths %>% mutate(pc_Excess_median = 1E6*median/sum(Pop_Str),
                                  pc_Excess_CI_low = 1E6*CI_low/sum(Pop_Str),
                                  pc_Excess_CI_high = 1E6*CI_high/sum(Pop_Str))

## As a percentage of total deaths
Total_Deaths_in_2020 <- Burs_Excess_total %>% filter(Week_st <="2020-12-21" & Week_st >="2020-01-01") %>% pull(Total_deaths) %>% sum()
Pandemic_Excess_Deaths %>% filter(Line == "Excess", Week_st == "2020-12-21") %>% select(-Line, -Week_st) %>%
  mutate(100*median/(Total_Deaths_in_2020/0.9),
         100*CI_low/(Total_Deaths_in_2020/0.9),
         100*CI_high/(Total_Deaths_in_2020/0.9))


p7 <- ggplot(Panel_G_plot_data_long %>% merge(Dates_df), aes(x = Week_st, color = Line, fill = Line)) +
  geom_line(aes(y = median), linewidth = 0.8) +
  geom_ribbon(aes(ymin = CI_low, ymax = CI_high), alpha = 0.4, linewidth = 0.2, color = NA) +
  ggtitle("G") +
  scale_y_continuous(name = "Cumulative excess burial\nregistrations/mortality") +
  theme_minimal() +
  xlab("") +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_x_date(expand = c(0,0), limits = as.Date(c("2019-12-30","2021-07-01")), date_labels = "%b %y", date_breaks = "6 months", date_minor_breaks = "3 month") +  # scale_x_date(expand = c(0,0)) +
  scale_fill_manual(name = "", labels = c("Burial registrations","Mortality","Scaling factor"), values = c("black","darkblue","red")) +
  scale_color_manual(name = "", labels = c("Burial registrations","Mortality","Scaling factor"), values = c("black","darkblue","red")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10)) +
  coord_cartesian(xlim = as.Date(c("2020-01-01","2021-06-07")))

####################################
##'[Panel H: DMVI]##
####################################
## Get overall population IFR for Zambia and Lusaka
IFR_weighted_Zamb <- readRDS("Analysis/Data/derived_data/09_IFR_Zambia.rds") %>% pull(IFR_Braz)
IFR_weighted_Lus <- readRDS("Analysis/Data/derived_data/10_Lusaka_Overall_IFR.rds")

# Get WHO excess mortality for Zambia by month and calculate DVWI based on Zambia IFR
WHO_Covid_deaths <- readxl::read_xlsx("Analysis/Data/raw_data/02_WHO_COVID_Excess_Deaths_EstimatesByCountry.xlsx", sheet = "Country by year and month", skip =12) %>%
  filter(country =="Zambia",year == 2020 | year == 2021 & month %in% 1:6) %>%
  select(year, month, cumul.excess.mean, cumul.excess.low, cumul.excess.high) %>%
  add_row(.before = 1, year = 2020, month = 0, cumul.excess.mean = 0, cumul.excess.low = 0, cumul.excess.high = 0) %>%
  mutate(Date = seq.Date(from = as.Date("2020-01-01"), to = as.Date("2021-07-01"), by = "1 month")-1) %>%
  mutate(Excess_deaths_pc = cumul.excess.mean*100000/(sum(un_demog_Z$pop)*1E3)) %>%
  mutate(DMVI = (cumul.excess.mean/(IFR_weighted_Zamb/100))/(sum(un_demog_Z$pop)*1E3),
         DMVI_low = (cumul.excess.low/(IFR_weighted_Zamb/100))/(sum(un_demog_Z$pop)*1E3),
         DMVI_high = (cumul.excess.high/(IFR_weighted_Zamb/100))/(sum(un_demog_Z$pop)*1E3))


# Calculate DVWI for Lusaka
Panel_H_plot_data_long <-
  Panel_G_plot_data_long %>% group_by(Line, Week_st) %>%
  mutate_at(c("median","CI_low","CI_high"), function(x)(10000*x/(IFR_weighted_Lus))/sum(Pop_Str))

######################
# Get results for text
Pandemic_DVWI <- Panel_G_plot_data_long %>% filter(Week_st =="2020-12-21" | Week_st =="2021-06-07")
Pandemic_DVWI

######################
## As a percentage of total deaths
Total_Deaths_in_2020 <- Burs_Excess_total %>% filter(Week_st <="2020-12-21" & Week_st >="2020-01-01") %>% pull(Total_deaths) %>% sum()
Pandemic_Excess_Deaths %>% filter(Line == "Excess", Week_st == "2020-12-21") %>% select(-Line, -Week_st) %>%
  mutate(100*median/Total_Deaths_in_2020,
         100*CI_low/Total_Deaths_in_2020,
         100*CI_high/Total_Deaths_in_2020)


p7 <- ggplot(Panel_G_plot_data_long %>% merge(Dates_df), aes(x = Week_st, color = Line, fill = Line, linetype = Line)) +
  geom_ribbon(data = WHO_Covid_deaths, aes(x = as.Date(Date), ymin = DMVI_low, ymax = DMVI_high, color = "WHO Zambia", fill = "WHO Zambia", linetype = "WHO Zambia"), alpha = 0.2, inherit.aes = F) +
  geom_ribbon(aes(ymin = CI_low/100, ymax = CI_high/100), alpha = 0.1, linewidth = 0.1, color = NA) +
  geom_line(data = WHO_Covid_deaths, aes(x = as.Date(Date), y = DMVI, color = "WHO Zambia", fill = "WHO Zambia", linetype = "WHO Zambia"), inherit.aes = F, linewidth = 0.8) +
  geom_line(aes(y = median/100), linewidth = 0.8, alpha = 0.3) +
  ggtitle("H") +
  scale_y_continuous(labels = scales::percent, name = stringr::str_wrap("DVWI", width = 40)) +
  theme_minimal() +
  xlab("") +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_x_date(expand = c(0,0), limits = as.Date(c("2019-12-30","2021-07-01")), date_labels = "%b %y", date_breaks = "6 months", date_minor_breaks = "3 month") +  # scale_x_date(expand = c(0,0)) +
  scale_fill_manual(name = "", labels = c("Burial registrations","Mortality","WHO Zambia"), values = c("black","darkblue","darkgreen")) +
  scale_color_manual(name = "", labels = c("Burial registrations","Mortality","WHO Zambia"), values = c("black","darkblue","darkgreen")) +
  scale_linetype_manual(name = "", labels = c("Burial registrations","Mortality","WHO Zambia"), values = c(1,1,4)) +
  coord_cartesian(xlim = as.Date(c("2020-01-01","2021-06-07"))) +
  guides(colour = guide_legend(nrow = 1),
         fill = guide_legend(nrow = 1),
         linetype = guide_legend(nrow = 1)) +
  theme(legend.box.just = "top",
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10))


tiff("Analysis/Figures/Figure_3_Excess_Mortality.tiff", width = 12, height = 12, units = "in", res = 150)
cowplot::plot_grid(cowplot::plot_grid(p1 + theme(legend.position = "none",
                                                 plot.title = element_text(vjust = 1, hjust = -0.09)),
                                      p2 + theme(legend.position = "none",
                                                 plot.title = element_text(vjust = 1, hjust = -0.09)),
                                      p3 + theme(legend.position = "none",
                                                 plot.title = element_text(vjust = 1, hjust = -0.09)),
                                      p4 + theme(legend.position = "none",
                                                 plot.title = element_text(vjust = 1, hjust = -0.09)), nrow = 1, rel_widths = c(1.05,1,1,1.05)),
                   cowplot::plot_grid(
                     cowplot::plot_grid(p5a + theme(legend.position = "none"),
                                        p5b + theme(legend.position = "none"),
                                        nrow = 2, align = "v", rel_heights = c(1,1.5)),
                     cowplot::plot_grid(
                       cowplot::plot_grid(
                         p6 + theme(legend.position = "none"),
                         p7 + theme(legend.position = "none"), nrow = 1, rel_widths = c(1,1)),
                       ggpubr::get_legend(p7 + theme(legend.position = "top",
                                                     legend.spacing.x = unit(0.3, "cm"),
                                                     legend.text.align = 0,
                                                     legend.text = element_text(margin = margin(l = 0, r = 40)))), nrow = 2, rel_heights = c(1,0.2))), nrow = 2, rel_heights = c(7,3))

dev.off()
