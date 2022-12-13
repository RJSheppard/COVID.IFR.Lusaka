library(dplyr)
library(ggpubr)
library(tidyverse)
library(viridisLite)
library(apyramid)
library(forcats)
library(gridExtra)
library(readxl)

## Data is from:
# Excess death estimates by region (e.g. global, AFR, EMR etc) by age and gender for each year (By age is important for IFR estimates), also gives population sizes
# Excess death estimates by region (e.g. global, AFR, EMR etc) by year and month
# Country death estimates by country: by age and gender, also gives population sizes
# Country death estimates for each year, with population sizes

## WHO estimates of excess mortality
WHO_estimate_region_age_year <- read_xlsx(path = "Data/raw_data/WHO_COVID_Excess_Deaths_EstimatesByRegion.xlsx",sheet="Region by year, sex and age", skip = 10)
WHO_estimate_region_year <- read_xlsx("Data/raw_data/WHO_COVID_Excess_Deaths_EstimatesByRegion.xlsx",sheet="Region by year and month", skip = 11)
WHO_estimate_ctry_age_year <- read_xlsx("Data/raw_data/WHO_COVID_Excess_Deaths_EstimatesByCountry.xlsx",sheet="Country by year, sex and age", skip = 10)
WHO_estimate_ctry_year <- read_xlsx("Data/raw_data/WHO_COVID_Excess_Deaths_EstimatesByCountry.xlsx",sheet="Country rate by year", skip = 8)

# UN population demographics
un_demog <- read_xlsx(path = "Data/raw_data/unpopulation_dataportal_20220810111500.xlsx",sheet="Data",skip = 4)
names(un_demog) <- c("country","year",paste0("age_cat",1:21))
un_demog <- un_demog %>% mutate(country=replace(country,country=="Dem. Rep. of the Congo","Democratic Republic of the Congo"))

## Un data for ISO codes (not found in un demog)
un_for_iso <- read_xlsx("Data/raw_data/WPP2022_GEN_F01_DEMOGRAPHIC_INDICATORS_COMPACT_REV1.xlsx",sheet="Estimates",skip = 15)
names(un_for_iso) <- un_for_iso[1,]
un_for_iso <- un_for_iso[-1,]

# Create ISO3 conversion
un_for_iso <- un_for_iso %>%
  filter(Year==2020)%>%
  rename(country='Region, subregion, country or area *',
         iso3='ISO3 Alpha-code') %>%
  select(country,iso3) %>%
  mutate(country=replace(country,country=="Dem. Republic of the Congo","Democratic Republic of the Congo"))

# Get WHO region to country level ISO codes conversion
WHO_regions<-read.csv("Data/raw_data/who-regions.csv")%>%
  rename(iso3=Code)

# Get Brazeau et al. estimates of COVID-19 IFR (including CIs)
brazeau <- read_xlsx("Data/raw_data/Brazeau_et_al.xlsx")


#########################
#### Excess mortality ###
#########################

# Formatting for consistency
WHO_estimate_ctry_year <- WHO_estimate_ctry_year %>% rename(country=Country)

# Match with ISO code
WHO_estimate_ctry_year_region <- merge(WHO_estimate_ctry_year,WHO_regions,by="iso3")

# Get WHO estimates for African countries
WHO_estimate_ctry_year_africa <- WHO_estimate_ctry_year_region %>%
  filter(WHO.region=="Africa",year==2020) %>%
  rename(type=WHO.region,area=country) %>%
  select(type,area,excess.mean,excess.low,excess.high)

# Get population sizes for regions
WHO_pop <- WHO_estimate_ctry_year_region %>%
  filter(year==2020) %>%
  group_by(WHO.region) %>%
  summarise(population=sum(pop.e5))

# Get WHO region estimates
WHO_estimate_region_2020 <- WHO_estimate_region_year %>% filter(month==12,year==2020)
WHO_estimate_region_2020$area <- c("Global",WHO_pop$WHO.region)
WHO_estimate_region_2020$population <- c(sum(WHO_pop$population),WHO_pop$population)

# Divide by population to get comparative metric, where high and low are 95% confidence intervals
WHO_estimate_region_2020 <- WHO_estimate_region_2020 %>%
  mutate(excess.mean=cumul.excess.mean/population,
         excess.high=cumul.excess.high/population,
         excess.low=cumul.excess.low/population,
         type=c("Global",rep("WHO Region",6))) %>%
  select(type,area,excess.mean,excess.low,excess.high) %>%
  mutate(type=factor(type,levels=c("Global","WHO Region","Africa")))

# Join region excess with country excess.
pc_excess_plot_df <- rbind(WHO_estimate_region_2020,WHO_estimate_ctry_year_africa)


################################
#### Population weighted IFR ###
################################

# Combine older age categories into 19, 5-year age groups usable for Brazeau data
un_demog_long <- un_demog %>%
  filter(year==2020) %>%
  mutate(age_cat19=age_cat19+age_cat20+age_cat21) %>%
  select(!c(age_cat20,age_cat21)) %>%
  pivot_longer(!c(country,year),names_to = "age_cat",names_prefix="age_cat",values_to = "pop") %>%
  mutate(country=replace(country,country=="Dem. Republic of the Congo","Democratic Republic of the Congo"),
         country=replace(country,country=="Dem. People's Rep. of Korea","Dem. People's Republic of Korea"),
         country=replace(country,country=="Lao People's Dem. Republic","Lao People's Democratic Republic"))


# Rename country for consistency and merge with iso codes
un_demog_long<-merge(un_demog_long, un_for_iso,by="country", all = F)

# Merge with UN age structured deographies
IFR_merge<-merge(un_demog_long,brazeau,by="age_cat")

# Merge with WHO region info
IFR_merge<-merge(IFR_merge,WHO_regions,by="iso3")

# Get global, regional and country age-structured populations
Global_demog_i<-un_demog_long%>%
  group_by(age_cat)%>%
  summarise(pop=sum(pop)) %>%
  mutate(area = "global")

Region_demog_ii <-IFR_merge%>%
  group_by(WHO.region, age_cat)%>%
  summarise(n=sum(pop))%>%
  rename(area=WHO.region)%>%
  select(area,age_cat, n)#,IFR_low,IFR_up)

Country_demog_iii <- IFR_merge%>%
  filter(WHO.region=="Africa")

Complete_Demog <- Country_demog_iii %>% rename(area = country, n = pop) %>%
  select(area, age_cat, n) %>%
  merge(Region_demog_ii, all = T) %>%
  merge(Global_demog_i %>% rename(n = pop), all = T)

# Save complete age-specific demographies for use in "2_Overall_COVID_IFR_CIs"
saveRDS(Complete_Demog, "Data/derived_data/Complete_Population_demog.rds")

# Calculate IFR values for each area
IFR_global<-merge(brazeau,Global_demog_i,by="age_cat")%>%
  summarise(IFR_med=sum(median*pop)/sum(pop))%>%
  mutate(area="Global")%>%
  select(area,IFR_med)

IFR_country<-IFR_merge%>%
  filter(WHO.region=="Africa")%>%
  group_by(country,iso3)%>%
  summarise(IFR_med=sum(median*pop)/sum(pop))%>%
  rename(area=country)%>%
  select(area,IFR_med)

IFR_Region<-IFR_merge%>%
  group_by(WHO.region)%>%
  summarise(IFR_med=sum(median*pop)/sum(pop))%>%
  rename(area=WHO.region)%>%
  select(area,IFR_med)

IFR_plot_df<-rbind(IFR_global,IFR_Region,IFR_country)

# Graphics designation
africa_pal<-viridis(length(WHO_estimate_ctry_year_africa$area), alpha = 1, begin = 0, end = 1, direction = -1, option = "D")
WHO_pal<-viridis(6, alpha = 1, begin = 0.25, end = 0.75, direction = -1, option = "A")
region_ord<-WHO_estimate_region_2020$area[2:7][order(-WHO_estimate_region_2020$excess.mean[2:7])]
Africa_ord<-WHO_estimate_ctry_year_africa$area[order(-WHO_estimate_ctry_year_africa$excess.mean)]
Africa_ord[which(Africa_ord=="Democratic Republic of the Congo")]="Dem. Rep. Congo"
Africa_ord[which(Africa_ord=="United Republic of Tanzania")]="Tanzania"
col_scheme<-c("black",WHO_pal,africa_pal)
names(col_scheme)<-c("Global",region_ord,Africa_ord)

# For each country: calculate confirmed deaths per million in 2020
Confirmed_deaths_WHO <- read.csv(file = "Data/raw_data/WHO-COVID-19-global-data.csv", header = T)
Confirmed_deaths_pc_countries <- Confirmed_deaths_WHO %>%
  select(Country, WHO_region, Cumulative_deaths, Date_reported) %>%
  filter(Date_reported =="2020-12-31",
         Country %in% pc_excess_plot_df$area) %>%
  rename(area = Country) %>% select(-WHO_region) %>%
  merge(WHO_estimate_ctry_year_region %>% filter(year ==2020) %>% select(country, pop.e5) %>%
          rename(area = country)) %>%
  mutate(Confirmed_deaths_pc = Cumulative_deaths*10/pop.e5) %>%
  select(area, Confirmed_deaths_pc)

# For each region: calculate confirmed deaths per million in 2020
Confirmed_deaths_pc_regions <- Confirmed_deaths_WHO %>% filter(Date_reported =="2020-12-31") %>% select(Country, WHO_region, Cumulative_deaths) %>%
  filter(WHO_region != "Other") %>%
  group_by(WHO_region) %>%
  summarise(Confirmed_deaths = sum(Cumulative_deaths)) %>%
  cbind(WHO_pop) %>%
  add_row(WHO_region = "Global", Confirmed_deaths = sum(.$Confirmed_deaths),WHO.region = "Global", population = sum(.$population)) %>%
  mutate(Confirmed_deaths_pc = Confirmed_deaths*10/population) %>%
  rename(area = WHO.region) %>%
  select(area, Confirmed_deaths_pc)

Confirmed_deaths_pc_all <- rbind(Confirmed_deaths_pc_regions, Confirmed_deaths_pc_countries)

############################################
# Excess mortality results described in text
pc_excess_plot_df %>% filter(area %in% c("Algeria","South Africa","Zambia", "Seychelles","Mauritius", "Kenya","Togo", "Americas","Europe"))

###############################################
### Figure 1: Panel 1: WHO Excess Mortality ###
###############################################

excess_plot<-ggplot(
  pc_excess_plot_df %>%
    merge(Confirmed_deaths_pc_all, all.x = T) %>%
    mutate(area=replace(area,area=="Democratic Republic of the Congo","Dem. Rep. Congo"))%>%
    mutate(area=replace(area,area=="United Republic of Tanzania","Tanzania"))%>%
    mutate(area_order=fct_reorder(area,-excess.mean)),
  aes(x=area_order,excess.mean*10))+ ### *10 to make it per million: standard is per 100,000
  geom_vline(xintercept='Africa',color="deepskyblue",lwd=4,alpha=0.1)+
  geom_vline(xintercept='Zambia',color="deepskyblue",lwd=4,alpha=0.1)+
  geom_hline(yintercept = 0)+
  geom_point(position = position_dodge(0.3))+
  geom_errorbar(
    aes(ymin = excess.low*10, ymax = excess.high*10, color = area_order),
    position = position_dodge(0.3), width = 0.2
  )+
  geom_point(aes(color=area_order), position = position_dodge(0.3))+
  geom_point(aes(y = Confirmed_deaths_pc), shape = 4) + # Confirmed deaths
  theme_minimal()+
  theme(legend.position = "none", axis.text.x = element_text(angle = 65, hjust = 1),
        strip.text.x = element_blank(),
        plot.margin=margin(t = 0, r = 20, b = -10, l = 0, unit = "pt"))+
  xlab("")+ylab("Excess deaths per million")+ggtitle("A") +
  facet_grid(~type,scale="free",space="free")+
  scale_color_manual(values=col_scheme)+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 1))


# Get and manage IFR CrIs and calculate DVWI values
IFR_vals_Braz <- readRDS("Data/derived_data/Brazeau_IFR_Post_predictive_CrIs.rds") %>%
  rename(area = georegion,
         IFR_mean = mean,
         IFR_med = Q50) %>%
  filter(sero == "reg") %>%
  mutate(area = ifelse(area == "global", "Global", area)) %>%
  merge(IFR_plot_df %>% rename(IFR_Braz = IFR_med)) %>%
  # mutate(type = ifelse(area == "Global", 1, ifelse(area %in% c("Europe","South-East Asia","Africa","Eastern Mediterranean", "Western Pacific", "Americas"), 2, 3))) %>%
  merge(pc_excess_plot_df,by="area")%>%
  mutate(excess_weight_IFR_mid=excess.mean/IFR_Braz,
         excess_weight_IFR_high=excess.high/IFR_Braz,
         excess_weight_IFR_low=excess.low/IFR_Braz)%>%
  mutate(excess_weight_IFR_low=replace(excess_weight_IFR_low,excess_weight_IFR_low<0,0),
         excess_weight_IFR_mid=replace(excess_weight_IFR_mid,excess_weight_IFR_mid<0,-10),
         excess_weight_IFR_high=replace(excess_weight_IFR_high,excess_weight_IFR_high<0,0))

# Graphics designation (order)
region_ord_IFR<-IFR_vals_Braz %>% filter(type == "WHO Region") %>% arrange(-IFR_med) %>% pull(area)
Africa_ord_IFR<-IFR_vals_Braz %>% filter(type == "Africa") %>% arrange(-IFR_med) %>% pull(area)
Africa_ord_IFR[which(Africa_ord_IFR=="Democratic Republic of the Congo")]="Dem. Rep. Congo"
Africa_ord_IFR[which(Africa_ord_IFR=="United Republic of Tanzania")]="Tanzania"
col_scheme_IFR<-c("black",WHO_pal,africa_pal)
names(col_scheme_IFR)<-c("Global",region_ord_IFR,Africa_ord_IFR)

###############################################
#### Figure 1: Panel 2: Overall IFR by area ###
###############################################

IFR_plot<-ggplot(
  IFR_vals_Braz %>%
    mutate(area=replace(area,area=="Democratic Republic of the Congo","Dem. Rep. Congo"))%>%
    mutate(area=replace(area,area=="United Republic of Tanzania","Tanzania"))%>%
    mutate(area_order=fct_reorder(area,-IFR_med)),aes(x=area_order,IFR_med/100))+
  geom_vline(xintercept='Africa',color="deepskyblue",lwd=4,alpha=0.1)+
  geom_vline(xintercept='Zambia',color="deepskyblue",lwd=4,alpha=0.1)+
  geom_hline(yintercept = 0)+
  geom_point(aes(color=area_order), position = position_dodge(0.3))+
  geom_errorbar(aes(ymin = Q025/100, ymax = Q975/100, color = area_order), width = 0.2) +
  theme_minimal()+
  theme(legend.position = "none", axis.text.x = element_text(angle = 65, hjust = 1),
        strip.text.x = element_blank(),
        plot.margin=margin(t = 0, r = 20, b = -10, l = 0, unit = "pt")
  )+
  xlab("")+ylab("Infection-fatality ratio")+ggtitle("B") +
  facet_grid(~type,scale="free",space="free")+
  scale_color_manual(values=col_scheme_IFR)+
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(labels = scales::percent,limits=c(0,0.016))


# Annotations for Panel 3
ann_text <- data.frame(area_order = c("Western Pacific","Kenya","Rwanda","Mauritius","Seychelles","Togo"),excess_weight_IFR_mid=0,
                       type = factor(c("WHO Region",rep("Africa",5)),levels = levels(IFR_vals_Braz$type)))

ann_text2 <- data.frame(area_order = c("Liberia"),excess_weight_IFR_mid=1020,
                      type = factor("Africa",levels = levels(IFR_vals_Braz$type)))

# Graphics designation
region_ord_excess_IFR <- IFR_vals_Braz %>% filter(type == "WHO Region") %>% arrange(-excess_weight_IFR_mid) %>% pull(area)
Africa_ord_excess_IFR <- IFR_vals_Braz %>% filter(type == "Africa") %>% arrange(-excess_weight_IFR_mid) %>% pull(area)
Africa_ord_excess_IFR[which(Africa_ord_excess_IFR=="Democratic Republic of the Congo")]="Dem. Rep. Congo"
Africa_ord_excess_IFR[which(Africa_ord_excess_IFR=="United Republic of Tanzania")]="Tanzania"
col_scheme_excess_IFR<-c("black",WHO_pal,africa_pal)
names(col_scheme_excess_IFR)<-c("Global",region_ord_excess_IFR,Africa_ord_excess_IFR)


###################################
# IFR and DVWI results used in text
IFR_vals_Braz %>% filter(area %in% c("Zambia", "South Africa", "Algeria", "Europe", "Americas"))

#####################################
#### Figure 1: Panel 3: DVWI Plot ###
#####################################

DVWI_plot<-ggplot(
  IFR_vals_Braz%>%
    mutate(area=replace(area,area=="Democratic Republic of the Congo","Dem. Rep. Congo"))%>%
    mutate(area=replace(area,area=="United Republic of Tanzania","Tanzania"))%>%
    mutate(area_order=fct_reorder(area,-excess_weight_IFR_mid)),aes(x=area_order,excess_weight_IFR_mid/1000))+
  geom_vline(xintercept='Africa',color="deepskyblue",lwd=4,alpha=0.1)+
  geom_vline(xintercept='Zambia',color="deepskyblue",lwd=4,alpha=0.1)+
  geom_point(aes(color=area_order), position = position_dodge(0.3))+
  geom_text(data = ann_text,label="*",vjust = 0.95,color="darkgrey",size=6)+
  geom_text(data = ann_text2,label="*mean estimated excess deaths<0",vjust = 1,color="darkgrey",size=5,hjust=0)+
  geom_errorbar(aes(ymin = excess_weight_IFR_low/1000, ymax = excess_weight_IFR_high/1000, color = area_order), width = 0.2) +
  geom_hline(yintercept = 0)+
  theme_minimal()+
  theme(legend.position = "none", axis.text.x = element_text(angle = 65, hjust = 1),
        strip.text.x = element_blank(),
        plot.margin=margin(t = 0, r = 20, b = -10, l = 0, unit = "pt"))+
  xlab("")+ylab("DVWI")+ ggtitle("C") +
  facet_grid(~type,scale="free",space="free")+
  scale_color_manual(values=col_scheme_excess_IFR)+
  scale_y_continuous(labels = scales::percent, breaks = seq(0,1,by =0.25)) +
  coord_cartesian(ylim = c(0,1.1))


############################################
# Save plot
tiff("Figures/Figure_1_WHO_Excess_Mortality_DVWI.tiff",height=10,width=10,unit="in",res=200)
grid.arrange(excess_plot,IFR_plot,DVWI_plot,ncol=1)
dev.off()
############################################
############################################
############################################
############################################
############################################
