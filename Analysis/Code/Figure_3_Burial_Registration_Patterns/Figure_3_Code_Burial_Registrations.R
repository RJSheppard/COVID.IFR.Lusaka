### Figure 2 plot:
library(shadowtext)
library(ggplot2)
library(dplyr)
library(tidyr)

#########################################
### Panel A: Official COVID-19 Deaths ###
#########################################

# Read and format official COVID_19 deaths in Lusaka Province
# Data from Lusaka Province Dashboard
Off_data_prov_daily <- readRDS("Analysis/Data/derived_data/03_Lusaka_Province_Dashboard.rds") |>
  mutate(Dates = as.Date(Dates, format = "%y/%m/%d"),
         deaths = as.integer(deaths)) |>
  complete(Dates = seq.Date(min(Dates, na.rm = T), max(Dates, na.rm = T), by = "day"), fill = list(deaths = 0))

Off_data_prov <- Off_data_prov_daily |>
  mutate(Week = lubridate::floor_date(Dates, unit = "week", week_start = 1)) |>
  group_by(Week) %>% summarise(Off_deaths = sum(deaths)) %>%
  mutate(Roll_Av =zoo::rollapply(Off_deaths,3,mean,fill=NA))

p1 <- ggplot() +
  geom_line(data = Off_data_prov, aes(y = Off_deaths, x = Week, linetype = stringr::str_wrap("Lusaka Province", width = 25)),
            linewidth = 0.6, inherit.aes = F, color = "black") +
  ylab("Confirmed\nCOVID-19\ndeaths") +
  scale_x_date(position = "bottom", date_breaks = "1 month", date_labels = "%b %Y",
               sec.axis = dup_axis(breaks = as.Date(c("2020-03-17",
                                                      "2020-04-24",
                                                      "2020-06-01",
                                                      "2020-10-05")),
                                   labels = c("NPIs implemented",
                                              "Some restrictions lifted",
                                              "Further restrictions lifted",
                                              "Restrictions fully lifted")), expand = c(0,0)) +
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.text.x.top = element_text(angle = 45, hjust = 0),
        plot.title = element_text(vjust = -35, hjust = 0, face = "bold", size = 7),
        legend.text.align = 0,
        legend.justification = c("left","center")) +
  coord_cartesian(xlim = as.Date(c("2017-12-25","2021-06-15"))) +
  geom_vline(xintercept = as.Date(c("2020-03-17",
                                    "2020-04-24",
                                    "2020-06-01",
                                    "2020-10-05")), color = "black", linetype = "dashed", linewidth = 0.6) +
  ggtitle("a") +
  scale_linetype_manual("", values = c(1,3))


#####################################
### Panel B: Burial Registrations ###
#####################################
Total_burial_registrations <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Burial_registrations_by_week_2018_2021_total.rds")

p2 <- ggplot(Total_burial_registrations, aes(x = Week)) +
  geom_line(aes(y = BurRegs, alpha = "Burial\nregistrations", linewidth = "Burial\nregistrations")) +
  geom_line(aes(y = Roll_Av, alpha = "5-week\nRolling average", linewidth = "5-week\nRolling average")) +
  geom_vline(xintercept = as.Date(c("2020-03-17",
                                    "2020-04-24",
                                    "2020-06-01",
                                    "2020-10-05")), color = "black", linetype = "dashed", linewidth = 0.6) +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y", expand = c(0,0)) +
  coord_cartesian(xlim = as.Date(c("2017-12-25","2021-06-15")), ylim = c(0,500)) +
  labs(title = "b",
       y = "Burial\nregistrations") +
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title=element_text(face = "bold", hjust = 0, vjust = 0, size = 7),
        legend.text.align = 0,
        legend.justification = c("left","center"),
        legend.key.height = unit(6, 'mm')) +
  scale_alpha_manual("", values = c(0.5,1), breaks = c("Burial\nregistrations","5-week\nRolling average")) +
  scale_linewidth_manual("", values = c(0.4,0.6), breaks = c("Burial\nregistrations","5-week\nRolling average")) +
  guides(linewidth = guide_legend(override.aes = list(alpha = c(0.3,1))))


############################################################################
### Panel C: burial registrations by age relative to pre-pandemic median ###
############################################################################
Relative_BurRegs_Week_Age_Fig3c <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Burial_registrations_by_week_2018_2021_relative_to_pre_pandemic_mean.rds")

p3 <- ggplot(Relative_BurRegs_Week_Age_Fig3c) +
  geom_line(aes(x = Week, y = Rel_BurRegs, color = Age_gr_labs), alpha = 1, linewidth = 0.4) +
  geom_vline(xintercept = as.Date(c("2020-03-17",
                                    "2020-04-24",
                                    "2020-06-01",
                                    "2020-10-05")), color = "black", linetype = "dashed", linewidth = 0.6) +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y", expand = c(0,0)) +
  coord_cartesian(xlim = as.Date(c("2017-12-25","2021-06-15"))) +
  labs(title = "c",
       y = "Relative\nregistrations") +
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0, vjust = 0, size = 7),
        legend.text.align = 0,
        legend.justification = c("left","center"),
        legend.key.height = unit(4, 'mm')) +
  scale_color_manual("Age group", values = viridis::turbo(n = 5), breaks = c("0-4","5-14","15-24","25-54","55+"))


###########################################################
### Panel D: Average Age at death, burial registrations ###
###########################################################

# Read and format burial registrations to calculate average age
Av_Age <- readRDS(file = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Burial_registrations_by_week_2018_2021_average_age.rds")

p4 <- ggplot() +
  geom_line(data = Av_Age, aes(y = av_age, x = Week, color = "Average age"), linewidth = 0.4, inherit.aes = F) +
  geom_line(data = Av_Age, aes(y = Roll_av, x = Week, color = "5-week\nRolling average"), linewidth = 0.6, inherit.aes = F) +
  ylab("Average age\nat death (years)") +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y", expand = c(0,0)) +
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0, vjust = 0, size = 7),
        legend.text = element_text(),
        legend.text.align = 0,
        legend.justification = c("left","center"),
        legend.spacing.y = unit(-0.5, 'cm'),
        legend.key.height = unit(5, 'mm')) +
  coord_cartesian(xlim = as.Date(c("2017-12-25","2021-06-15"))) +
  geom_vline(xintercept = as.Date(c("2020-03-17","2020-04-24","2020-06-01","2020-10-05")), color = "black", linetype = "dashed", linewidth = 0.6) +
  ggtitle("d") +
  scale_color_manual("", values = c("darkgrey", "black"), breaks = c("Average age","5-week\nRolling average"))



####################################################
### Panel E: Age-structured burial registrations ###
####################################################
Age_Proportion_BurRegs <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Burial_registrations_by_week_2018_2021_proportion_of_total.rds")

p5 <- ggplot(Age_Proportion_BurRegs, aes(x = Week, y = Weekly_Prop_deaths, group = Age_gr_labs, fill = Age_gr_labs)) +
  geom_area() +
  geom_vline(xintercept = as.Date(c("2020-03-17","2020-04-24","2020-06-01","2020-10-05")), color = "white", linetype = "dashed", linewidth = 0.6) +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y", expand = c(0,0)) +
  coord_cartesian(as.Date(c("2017-12-25","2021-06-15"))) +
  scale_y_continuous(expand = c(0,0)) +
  labs(title = "e",
       x = "Date",
       y = "Proportion of\nregistrations") +
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme(axis.text.x = element_text(angle = 90),
        plot.title = element_text(hjust = 0, vjust = 0, face = "bold", size = 7),
        legend.text.align = 0,
        legend.justification = c("left","center"),
        legend.key.size = unit(4,"mm")) +
  viridis::scale_fill_viridis(discrete = T, option = "H") +
  guides(fill=guide_legend(title="Age group", label.hjust = 0),
         color = guide_legend(title = "Age group",
                              override.aes = list(color = c("black","white"),
                                                  bordercolor = c("white", viridis::turbo(n = 9, begin = 0, end =1)[9]),
                                                  bordersize = c(NULL, 2.5))))


ggsave(filename = "Analysis/Figures/Figure_3_Prop_Burial_Registrations.pdf",
       plot = cowplot::plot_grid(p1,p2,p3,p4,p5, ncol = 1, rel_heights = c(2.1,1.5,1.5,1.5,3), align = "v"),
       width = 180, height = 180, units = "mm")



#######################################################
## Results for text
Wave_1_official_deaths <- Off_data_prov_daily %>%
  filter(Dates>=as.Date("2020-06-01") & Dates < as.Date("2020-11-01")) %>%
  pull(deaths) %>% sum()

Wave_2_official_deaths <- Off_data_prov_daily %>%
  filter(Dates>=as.Date("2021-01-01") & Dates < as.Date("2021-04-01")) |>
  pull(deaths) %>% sum()

Wave_3_official_deaths <- Off_data_prov_daily %>%
  filter(Dates>=as.Date("2021-04-01")) |>
  pull(deaths) %>% sum()

## Results for text
## Total deaths in 2018

Total_2018_deaths <- read.csv(file = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/raw_data/mortuary_records_v2.csv") |>
  filter(dod !=".") %>%
  mutate(date = as.Date(dod, "%m/%d/%y")) |>
  filter(date >= "2018-01-01",date < "2019-01-01") |>
  nrow()

Total_2019_deaths <- read.csv(file = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/raw_data/mortuary_records_v2.csv") |>
  filter(dod !=".") %>%
  mutate(date = as.Date(dod, "%m/%d/%y")) |>
  filter(date >= "2019-01-01",date < "2020-01-01") |>
  nrow()

1000*Total_2018_deaths/2526102
1000*Total_2019_deaths/2627716
5.82/6.7
5.72/6.6

Total_burial_registrations %>% filter(Week < "2020-01-01") %>% pull(BurRegs) %>% median()

Total_burial_registrations %>% filter(BurRegs >500)
Total_burial_registrations %>% filter(Week < "2020-01-01", BurRegs<100)

Av_Age %>% filter(Week >= "2018-01-01", Week < "2019-12-25") %>% pull(av_age) %>% median()
Av_Age %>% filter(Week %in% as.Date(c("2020-07-20", "2021-01-04","2021-06-07")))
Av_Age %>% filter(Week >= "2018-01-01", Week < "2019-12-25") %>% pull(av_age) %>% max()

Total_burial_registrations %>% filter(Week %in% as.Date(c("2021-01-18","2021-06-07")))


Total_burial_registrations %>% filter(Week == "2020-03-16")
BurRegs %>% filter(Week == "2020-07-20")
BurRegs %>% filter(Week > "2020-08-17" & Week < "2020-10-20")
Total_burial_registrations %>% filter(Week > "2021-01-01" & Week < "2021-02-01")
Total_burial_registrations %>% filter(Week > "2021-06-01")

## Results for text
