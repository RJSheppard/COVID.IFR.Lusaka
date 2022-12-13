### Figure 1 plot:
library(shadowtext)
library(ggplot2)
library(dplyr)

#########################################
### Panel A: Official COVID-19 Deaths ###
#########################################

# Read and format official COVID_19 deaths in Lusaka
Off_data_prov <- readRDS("Data/derived_data/00_01_Lusaka_Prov_Deaths_Official.rds") %>%
  mutate(Week = lubridate::floor_date(date, unit = "week", week_start = 1)) %>%
  group_by(Week) %>% summarise(Off_deaths = sum(deaths)) %>%
  mutate(Roll_Av =zoo::rollapply(Off_deaths,3,mean,fill=NA))

Off_data_dist <- readRDS("Data/derived_data/00_01_Lusaka_Dist_Deaths_Official.rds") %>%
  mutate(Week = lubridate::floor_date(date, unit = "week", week_start = 1)) %>%
  group_by(Week) %>% summarise(Off_deaths = sum(deaths)) %>%
  mutate(Roll_Av =zoo::rollapply(Off_deaths,3,mean,fill=NA))

p1 <- ggplot() +
  geom_line(data = Off_data_dist, aes(y = Off_deaths, x = Week, linetype = stringr::str_wrap("Lusaka district", width = 25)), linewidth = 1.2, inherit.aes = F, color = "black") +
  geom_line(data = Off_data_prov, aes(y = Off_deaths, x = Week, linetype = stringr::str_wrap("Lusaka province", width = 25)), linewidth = 1, inherit.aes = F, color = "black") +
  theme_minimal() +
  ylab("Confirmed\nCOVID-19 deaths") +
  scale_x_date(position = "bottom", date_breaks = "1 month", date_labels = "%b %Y",
               sec.axis = dup_axis(breaks = as.Date(c("2020-03-17","2020-04-24","2020-06-01","2020-10-05")),
                                   labels = c("NPIs implemented","Some restrictions lifted","Further restrictions lifted","Restrictions fully lifted")), expand = c(0,0)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.text.x.top = element_text(angle = 45, hjust = 0),
        plot.title = element_text(vjust = -25),
        legend.spacing.y = unit(0.3, 'cm'),
        axis.title.y = element_text(size = 12), legend.text = element_text(size = 11), legend.text.align = 0, legend.justification = c("left","center")) +
  coord_cartesian(xlim = as.Date(c("2017-12-25","2021-06-15"))) +
  geom_vline(xintercept = as.Date(c("2020-03-17","2020-04-24","2020-06-01","2020-10-05")), color = "black", linetype = "dashed", linewidth = 1) +
  ggtitle("A") +
  scale_linetype_manual("", values = c(1,3))

## Results for text
readRDS("Data/derived_data/00_01_Lusaka_Prov_Deaths_Official.rds") %>% filter(date>=as.Date("2020-06-01") & date < as.Date("2020-09-01")) %>%
  pull(deaths) %>% sum()


#####################################
### Panel B: Burial Registrations ###
#####################################

# Get and format total burial registrations
BurRegs <- readRDS("Data/derived_data/00_07_Burial_registrations_by_week_2017_to_2021.rds") %>%
  rename(Week = Week_st,
         BurRegs = Total_deaths) %>%
  filter(Week >= as.Date("2017-12-25"),
         Week < as.Date("2021-06-12")) %>%
  mutate(Roll_Av =zoo::rollapply(BurRegs,5,mean,fill=NA))

p2 <- ggplot() +
  geom_line(data = BurRegs, aes(y = BurRegs, x = Week, color = "Burial registrations"), linewidth = 0.8, inherit.aes = F) +
  geom_line(data = BurRegs, aes(y = Roll_Av, x = Week, color = "Rolling average"), linewidth = 1, inherit.aes = F) +
  theme_minimal() +
  ylab("Burial\nregistrations") +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y", expand = c(0,0)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 12), legend.text = element_text(size = 11), legend.text.align = 0, legend.justification = c("left","center")) +
  guides(linetype = guide_legend(byrow = TRUE))+
  coord_cartesian(xlim = as.Date(c("2017-12-25","2021-06-15")), ylim = c(0,500)) +
  geom_vline(xintercept = as.Date(c("2020-03-17","2020-04-24","2020-06-01","2020-10-05")), color = "black", linetype = "dashed", linewidth = 1) +
  ggtitle(label = "B") +
  scale_color_manual("", values = c("darkgrey", "black"), breaks = c("Burial registrations","Rolling average"))


## Results for text
BurRegs %>% filter(Week < "2020-01-01") %>% pull(BurRegs) %>% median()
BurRegs %>% filter(BurRegs >500)
BurRegs %>% filter(BurRegs<100)
BurRegs %>% filter(Week == "2020-03-16")
BurRegs %>% filter(Week == "2020-07-20")
BurRegs %>% filter(Week > "2020-08-17" & Week < "2020-10-20")
BurRegs %>% filter(Week > "2021-01-01" & Week < "2021-02-01")
BurRegs %>% filter(Week > "2021-06-01")


###########################################################
### Panel C: Average Age at death, burial registrations ###
###########################################################

# Read and format burial registrations to calculate average age
Av_Age <- read.csv(file = "Data/raw_data/mortuary_records_v2.csv") %>%
  filter(age_years !=".",dod !=".") %>%
  mutate(date = as.Date(dod, "%m/%d/%y")) %>%
  mutate(Week = lubridate::floor_date(date, unit = "week", week_start = 1)) %>%
  group_by(Week) %>% summarise(av_age = mean(as.numeric(age_years))) %>%
  filter(Week >= as.Date("2017-12-25"),
         Week < as.Date("2021-06-15")) %>%
  mutate(Roll_av =zoo::rollapply(av_age,5,mean,fill=NA))

p3 <- ggplot() +
  geom_line(data = Av_Age, aes(y = av_age, x = Week, color = "Average age"), linewidth = 0.8, inherit.aes = F) +
  geom_line(data = Av_Age, aes(y = Roll_av, x = Week, color = "Rolling average"), linewidth = 1, inherit.aes = F) +
  theme_minimal() +
  ylab("Average age\nat death") +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y", expand = c(0,0)) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_text(size = 12), legend.text = element_text(size = 11), legend.text.align = 0, legend.justification = c("left","center")) +
  coord_cartesian(xlim = as.Date(c("2017-12-25","2021-06-15"))) +
  geom_vline(xintercept = as.Date(c("2020-03-17","2020-04-24","2020-06-01","2020-10-05")), color = "black", linetype = "dashed", linewidth = 1) +
  ggtitle("C") +
  scale_color_manual("", values = c("darkgrey", "black"), breaks = c("Average age","Rolling average"))

## Results for text
Av_Age %>% filter(Week >= "2018-01-01", Week < "2019-12-25") %>% pull(av_age) %>% median()
Av_Age %>% filter(Week >= "2018-01-01", Week < "2019-12-25") %>% pull(av_age) %>% quantile(c(0.025,0.975))
Av_Age %>% filter(Week %in% as.Date(c("2020-07-20", "2021-01-04","2021-06-14")))
Av_Age %>% filter(Week >= "2018-01-01", Week < "2019-12-25") %>% pull(av_age) %>% max()


####################################################
### Panel D: Age-structured burial registrations ###
####################################################

# Read and format age-structured burial registrations
UTH_Mortality_Total <- read.csv(file = "Data/raw_data/mortuary_records_v2.csv") %>%
  filter(age_years !=".",
         dod !=".") %>%
  mutate(date = as.Date(dod, "%m/%d/%y"),
         Age_gr = cut(as.numeric(age_years), c(seq(0,80,by = 5),Inf), right = F, labels = F)) %>%
  select(-sex, -dod, -age_years)

PlotData <- UTH_Mortality_Total %>%
  mutate(Week = lubridate::floor_date(date, unit = "week", week_start = 1),
         Age_gr = cut(Age_gr, c(seq(0,18,by = 2)), right = T, labels = F)) %>%
  group_by(Week, Age_gr) %>%
  summarise(Weekly_Age_deaths = n()) %>%
  ungroup() %>% tidyr::complete(Week, Age_gr, fill = list(Weekly_Age_deaths = 0)) %>%
  group_by(Week) %>%
  mutate(Weekly_Total_deaths = sum(Weekly_Age_deaths),
         Weekly_Prop_deaths = Weekly_Age_deaths/Weekly_Total_deaths) %>%
  filter(Week >= as.Date("2017-12-25"),
         Week < as.Date("2021-06-15"))

p4 <- ggplot(data = PlotData, aes(x = Week, y = Weekly_Prop_deaths, group = as.factor(Age_gr), fill = as.factor(Age_gr))) +
  geom_area() +
  viridis::scale_fill_viridis(discrete = T, option = "H", labels = c("0-4",paste0(c(1:7)*10-5,"-",c(1:7)*10+4),"75+")) +
  guides(fill=guide_legend(title="", label.hjust = 0)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90), axis.title.y = element_text(size = 12), legend.text = element_text(size = 11), legend.text.align = 0, legend.justification = c("left","center")) +
  ylab("Proportion") + xlab("Date") +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y", expand = c(0,0)) +
  geom_vline(xintercept = as.Date(c("2020-03-17","2020-04-24","2020-06-01","2020-10-05")), color = "white", linetype = "dashed", linewidth = 1) +
  scale_y_continuous(name = "Proportion", expand = c(0,0)) +
  guides(color = guide_legend(title = "Age group",override.aes = list(color = c("black","white"),
                                                                      bordercolor = c("white", viridis::turbo(n = 9, begin = 0, end =1)[9]),
                                                                      bordersize = c(NULL, 2.5)))) +
  coord_cartesian(as.Date(c("2017-12-25","2021-06-15"))) +
  ggtitle(label = "D")


#####################
### Format figure ###
#####################

g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)
g3 <- ggplotGrob(p3)
g4 <- ggplotGrob(p4)

# set the same widths for both blots
library(grid)
g1$widths <- unit.pmax(g1$widths, g2$widths, g3$widths, g4$widths)
g2$widths <- unit.pmax(g1$widths, g2$widths, g3$widths, g4$widths)
g3$widths <- unit.pmax(g1$widths, g2$widths, g3$widths, g4$widths)
g4$widths <- unit.pmax(g1$widths, g2$widths, g3$widths, g4$widths)

# stack them afterwards
g <- rbind(g1, g2, g3, g4, size="first")
g$heights[7] <- unit(1,"null")
g$heights[19] <- unit(2,"null")
g$heights[31] <- unit(1,"null")
g$heights[43] <- unit(4,"null")
grid.newpage()
grid.draw(g)

tiff("Figures/Figure_2_Prop_Burial_Registrations.tiff", width = 10, height = 10, units = "in", res = 150)
grid.draw(g)
dev.off()
