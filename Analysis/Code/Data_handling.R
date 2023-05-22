### Data

##############################
## Burial Registration Data ##
##############################

# Total Burial registrations by week
Burial_registrations_cleaned <- read.csv(file = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/raw_data/mortuary_records_v2.csv") |>
  filter(age_years !=".",
         dod !=".") %>%
  mutate(date = as.Date(dod, "%m/%d/%y"),
         Age_gr = cut(as.numeric(age_years), c(seq(0,80,by = 5),Inf), right = F, labels = F),
         Age_gr_fig_3c = cut(as.numeric(age_years), c(0,seq(5,25,by = 10),55,Inf), right = F, labels = F),
         Age_gr_fig_3e = cut(as.numeric(age_years), c(0,seq(5,75,by = 10),Inf), right = F, labels = F),
         Age_gr_fig_4 = cut(as.numeric(age_years), c(0,5,15,25,40,50,60,70,Inf), right = F, labels = F),
         Week = lubridate::floor_date(date, unit = "week", week_start = 1)) |>
  filter(Week >= as.Date("2017-12-25"),Week < as.Date("2021-06-12")) |>
  select(-dod)

saveRDS(Burial_registrations_cleaned, file = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Burial_registrations_cleaned.rds")

Total_burial_registration_by_week_2018_2021 <- Burial_registrations_cleaned |>
  group_by(Week) %>%
  summarise(BurRegs = n()) |>
  mutate(Roll_Av =zoo::rollapply(BurRegs,5,mean,fill=NA))

saveRDS(Total_burial_registration_by_week_2018_2021, file = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Burial_registrations_by_week_2018_2021_total")

# Relative age rates of burial registrations
BurRegs_Week_Age_Fig3c <- Burial_registrations_cleaned |>
  group_by(Week, Age_gr_fig_3c) |>
  summarise(BurRegs = n()) |>
  rename(Age_gr = Age_gr_fig_3c) |>
  ungroup() %>% tidyr::complete(Week, Age_gr, fill = list(BurRegs = 0))

Rel_BurRegs_Week_Age_Fig3c <- BurRegs_Week_Age_Fig3c |>
  filter(Week>="2018-01-01", Week<"2020-01-01") |>
  group_by(Age_gr) |>
  summarise(PrePan_av_BurRegs = mean(BurRegs)) |>
  mutate(Age_gr_labs = c(paste0(c(0,5,15,25),"-",c(4,14,24,54)),"55+")) |>
  merge(BurRegs_Week_Age_Fig3c) |>
  mutate(Rel_BurRegs = BurRegs/PrePan_av_BurRegs)

saveRDS(Rel_BurRegs_Week_Age_Fig3c, file = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Burial_registrations_by_week_2018_2021_relative_to_pre_pandemic_mean.rds")

## Average Age
Av_Age <- Burial_registrations_cleaned |>
  group_by(Week) %>% summarise(av_age = mean(as.numeric(age_years))) %>%
  mutate(Roll_av =zoo::rollapply(av_age,5,mean,fill=NA))

saveRDS(Av_Age, file = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Burial_registrations_by_week_2018_2021_average_age.rds")

## Weekly age proportion of burial registrations
Age_Proportion_BurRegs <- Burial_registrations_cleaned %>%
  group_by(Week, Age_gr_fig_3e) %>%
  summarise(BurRegs = n()) %>%
  rename(Age_gr = Age_gr_fig_3e) |>
  ungroup() %>% tidyr::complete(Week, Age_gr, fill = list(BurRegs = 0)) %>%
  group_by(Week) %>%
  mutate(Weekly_Total_deaths = sum(BurRegs),
         Weekly_Prop_deaths = BurRegs/Weekly_Total_deaths) |>
  merge(data.frame(Age_gr = 1:9, Age_gr_labs = c(paste0(c(0,seq(5,65,by=10)),"-",c(seq(4,74,by=10))),"75+")))

saveRDS(Age_Proportion_BurRegs, file = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Burial_registrations_by_week_2018_2021_proportion_of_total.rds")

## Weekly registrations  by 5-year age group.
Age_structured_Burial_Registrations <- Burial_registrations_cleaned %>%
  filter(date >= "2018-01-01") %>%
  group_by(Age_gr, Week) %>%
  summarise(BurRegs = n()) %>%
  ungroup() %>%
  complete(Age_gr, Week, fill = list(BurRegs = 0)) %>%
  arrange(Age_gr, Week)

saveRDS(Age_structured_Burial_Registrations, file = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Burial_registrations_by_week_2018_2021_5yr_age_groups.rds")

## Weekly registrations  by Fig4 age groups.
Age_structured_Burial_Registrations_Fig4 <- Burial_registrations_cleaned |>
  group_by(Age_gr_fig_4, Week) |>
  summarise(BurRegs = n()) |>
  rename(Age_gr = Age_gr_fig_4) |>
  ungroup() %>%
  complete(Age_gr, Week, fill = list(BurRegs = 0)) %>%
  merge(data.frame(Age_gr = 1:8, Age_gr_labs = c(paste0(c(0,5,15,25,40,50,60),"-",c(4,14,24,39,49,59,69)),"70+")))

saveRDS(Age_structured_Burial_Registrations_Fig4, file = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Burial_registrations_by_week_2018_2021_Fig4_age_groups.rds")

########################
## UN demography data ##
########################
un_demog<-readxl::read_xlsx(path = "Analysis/Data/raw_data/unpopulation_dataportal_20220810111500.xlsx",sheet="Data",skip = 4)
names(un_demog)<-c("country","year",paste0("age_cat",1:21))
un_demog_Z <- un_demog %>% filter(country =="Zambia", year == 2020) %>%
  mutate(age_cat19=age_cat19+age_cat20+age_cat21)%>%
  select(!c(age_cat20,age_cat21))%>%
  pivot_longer(!c(country,year),names_to = "age_cat",names_prefix="age_cat",values_to = "pop")

saveRDS(un_demog_Z, file = "Analysis/Data/derived_data/UN_Demography_Zambia.rds")



##########################
## Post_Mortem_PCR_data ##
##########################
post_mortem <- readxl::read_xlsx(path = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/raw_data/post_mortem_results.xlsx") %>%
  mutate(covid19_combined = ifelse(sampleid %in% c(6429,6148,6726),1,covid19_combined)) %>%
  select(deceased_date, covid19_combined, age_death) %>%
  rename(date = deceased_date) %>%
  mutate(date = as.Date(date),
         Age_gr = cut(as.numeric(age_death), c(seq(0,80,by = 5),Inf), right = F, labels = F)) %>%
  select(-age_death) %>%
  filter(covid19_combined != ".") %>%
  group_by(date, Age_gr) %>%
  summarise(Samples = length(date),
            PosTests = sum(ifelse(covid19_combined %in% c(1,2),T,F)),
            PosTests_Strict = sum(ifelse(covid19_combined ==2,T,F))) %>%
  arrange(date) %>% ungroup() %>%
  tidyr::complete(Age_gr, date = seq.Date(min(date), max(date), by="day"),
                  fill = list(Samples = 0, PosTests = 0, PosTests_Strict = 0))

saveRDS(post_mortem, file = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Post_Mortem_cleaned.rds")

Age_structured_Burial_Registrations <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Burial_registrations_by_week_2018_2021_5yr_age_groups.rds")

Combined_data <- post_mortem %>%
  filter(date>="2020-06-15") %>%
  mutate(Week_gr = cut.Date(x = date, breaks = "weeks", labels = F),
         Week = lubridate::floor_date(date, unit = "week", week_start = 1)) %>%
  ungroup() %>%
  select(-date) |>
  group_by(Week_gr, Age_gr) %>%
  summarise(Week = Week[1],
            Samples = sum(Samples),
            PosTests = sum(PosTests),
            PosTests_Strict = sum(PosTests_Strict)) %>%
  merge(Age_structured_Burial_Registrations |> filter(Week>="2020-06-15", Week<"2020-10-05"), all = T)

saveRDS(Combined_data, file = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/combined_squire_data.rds")


#######################################
## Population PCR and seroprevalence ##
#######################################
df_pcr <- read.csv("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/raw_data/S2PS dataset for IFR model.csv") |>
  filter(!is.na(pcr)) %>%
  mutate(Pos = ifelse(pcr ==1, 1, 0)) %>%
  summarise(Samples = length(Pos),
            Pos_tests = sum(Pos),
            Pos_prev = 100*sum(Pos)/length(Pos))

df_sero <- read.csv("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/raw_data/S2PS dataset for IFR model.csv") |>
  filter(!is.na(elisa)) %>%
  mutate(Pos = ifelse(elisa ==1, 1, 0)) %>%
  summarise(Samples = length(Pos),
            Pos_tests = sum(Pos),
            Pos_prev = 100*sum(Pos)/length(Pos))

df_dates <- read.csv("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/raw_data/S2PS dataset for IFR model.csv") |>
  filter(date != "") %>%
  pull(date) %>%
  list(st_date = min(as.Date(., format = "%m/%d/%y")), end_date = max(as.Date(., format = "%m/%d/%y")))

pcr_df <- data.frame(date_start = df_dates$st_date, date_end = df_dates$end_date,
                     pos_tests = df_pcr$Pos_tests, samples = df_pcr$Samples)
sero_df <- data.frame(date_start = df_dates$st_date, date_end = df_dates$end_date,
                      pos_tests = df_sero$Pos_tests, samples = df_sero$Samples)

Lancet_Data <- list(pcr_df=pcr_df, sero_df=sero_df)
saveRDS(Lancet_Data, "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Population_prevalence.rds")


# PCR prevalence by age
PCR_by_age <- read.csv("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/raw_data/S2PS dataset for IFR model.csv") |>
  filter(!is.na(pcr), !is.na(age)) %>%
  mutate(Age_gr_2 = cut(as.numeric(age), c(0,5,15,25,40,50,60,70,Inf), right = F, labels = F)) %>%
  group_by(Age_gr_2) %>% summarise(positive = sum(pcr), total = n(),
                                   Perc = positive/total,
                                   CI_lower = drop(Hmisc::binconf(positive, total))["Lower"],
                                   CI_upper = drop(Hmisc::binconf(positive, total))["Upper"])

## Seroprevalence by age
Sero_by_age <- read.csv("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/raw_data/S2PS dataset for IFR model.csv") |>
  filter(!is.na(elisa), !is.na(age)) %>%
  mutate(Age_gr_2 = cut(as.numeric(age), c(0,5,15,25,40,50,60,70,Inf), right = F, labels = F)) %>%
  group_by(Age_gr_2) %>% summarise(positive = sum(elisa), total = n(),
                                 Perc = positive/total,
                                 CI_lower = drop(Hmisc::binconf(positive, total))["Lower"],
                                 CI_upper = drop(Hmisc::binconf(positive, total))["Upper"])

saveRDS(list(pcr_by_age = PCR_by_age, sero_by_age = Sero_by_age), "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Population_prevalence_by_age.rds")

#########################
## Severity estimataes ##
#########################
Zmb_p <-squire::parameters_explicit_SEEIR("Zambia")

# Calculate the probability of death by weighting probabilities of death given treatment by probability of severe
prob_death_tot <- Zmb_p$prob_severe * Zmb_p$prob_severe_death_treatment +
  (1-Zmb_p$prob_severe) * Zmb_p$prob_non_severe_death_treatment

# Calculate IFR for each age group: multiply probability of death of a case by probability hospitalised
IFR_Age <- Zmb_p$prob_hosp * prob_death_tot

# 9x9 matrix of IFR values:
IFR_vec <- seq(0.2,1, by = 0.2) # IFR vector
IFR_vec_lsc <- sort(unique(c(IFR_vec,1/IFR_vec)))
Slope_vec <- seq(0.2,1, by = 0.2) # IFR vector
Slope_vec_lsc <- sort(unique(c(Slope_vec,1/Slope_vec)))

# Expand IFR/slope
IFR_mat <- expand.grid("IFR_x" = IFR_vec_lsc, "Slope_x" = Slope_vec_lsc)

pop_st_lu <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Population_structure_Lusaka_2020_CDC.rds")
pop_index <- pop_st_lu/sum(pop_st_lu)

## Take IFR_Age, scale and gradient.
IFR_sim_func <- function(IFR_Age, Scale, Slope, pop_index){
  Scale*exp(Slope*log(IFR_Age)) *
    sum(IFR_Age*pop_index) /
    sum(exp(Slope*log(IFR_Age))*pop_index)
}

# Calculate IFR for each slope/scale combination
IFR_sim <- lapply(1:81, function(x){
  IFR_sim_func(IFR_Age,IFR_mat[x,1], IFR_mat[x,2], pop_index)
})

# Calculate Prob of death, given hospitalisation for each slope/scale combination
phi_2 <- lapply(IFR_sim, function(x){
  sapply(1:17, function(y){min(x[y]/Zmb_p$prob_hosp[y], 1)})
})

# Calculate prob of hospitalisation, for each slope/scale combination
phi_1 <- lapply(IFR_sim, function(x){
  sapply(1:17, function(y){max(Zmb_p$prob_hosp[y], x[y])})
})

## State which can be run (that do not exceed 100%):
Run_List <- lapply(phi_1, function(x){ifelse(any(x>1),F,T)})

Model_Fit_Inputs <- lapply(1:81, function(x){
  return(list(IFR = IFR_sim[[x]], phi_1 = phi_1[[x]], phi_2 = phi_2[[x]], OK_to_run = Run_List[[x]]))
})

names(Model_Fit_Inputs) <- paste0("X",1:81)
saveRDS(Model_Fit_Inputs, "Analysis/Data/derived_data/IFR_phosp_pdeath.rds")
saveRDS(IFR_mat, "Analysis/Data/derived_data/IFR_matrix.rds")

##########################
## Duration until death ##
##########################
Death_Dur_WeightedByAge <-
  Zmb_p$dur_get_mv_die *
  Zmb_p$prob_severe * Zmb_p$prob_severe_death_treatment/
  (Zmb_p$prob_severe * Zmb_p$prob_severe_death_treatment + (1-Zmb_p$prob_severe)*Zmb_p$prob_non_severe_death_treatment) +

  Zmb_p$dur_get_ox_die *
  (1-Zmb_p$prob_severe) * Zmb_p$prob_non_severe_death_treatment/
  (Zmb_p$prob_severe * Zmb_p$prob_severe_death_treatment + (1-Zmb_p$prob_severe)*Zmb_p$prob_non_severe_death_treatment)

Surv_Dur_WeightedByAge <-
  Zmb_p$dur_get_mv_survive *
  Zmb_p$prob_severe*(1-Zmb_p$prob_severe_death_treatment)/
  (Zmb_p$prob_severe*(1-Zmb_p$prob_severe_death_treatment) + (1-Zmb_p$prob_severe)*(1-Zmb_p$prob_non_severe_death_treatment)) +

  Zmb_p$dur_get_ox_survive *
  (1-Zmb_p$prob_severe) * (1-Zmb_p$prob_non_severe_death_treatment)/
  (Zmb_p$prob_severe*(1-Zmb_p$prob_severe_death_treatment) + (1-Zmb_p$prob_severe)*(1-Zmb_p$prob_non_severe_death_treatment))

## Weight durations per age group according to age prevalence in Lusaka and probability they need treatment by age
Death_Dur_Weighted <- sum(Death_Dur_WeightedByAge * Zmb_p$prob_hosp/sum(Zmb_p$prob_hosp))
Surv_Dur_Weighted <- sum(Surv_Dur_WeightedByAge * Zmb_p$prob_hosp/sum(Zmb_p$prob_hosp))

saveRDS(list(Death_Dur_Weighted = Death_Dur_Weighted, Surv_Dur_Weighted = Surv_Dur_Weighted),
        "Analysis/Data/derived_data/Weighted_durations_death_survive.rds")

##########################
## PCR detection ##
##########################
pcr_sens = 1
# pcr_sens = 0.8
pcr_det <- c(9.206156e-13, 9.206156e-13, 3.678794e-01, 9.645600e-01,
             9.575796e-01, 9.492607e-01, 9.393628e-01, 9.276090e-01,
             9.136834e-01, 8.972309e-01, 8.778578e-01, 8.551374e-01,
             8.286197e-01, 7.978491e-01, 7.623916e-01, 7.218741e-01,
             6.760375e-01, 6.248060e-01, 5.683688e-01, 5.072699e-01,
             4.525317e-01, 4.036538e-01, 3.600134e-01, 3.210533e-01,
             2.862752e-01, 2.552337e-01, 2.275302e-01, 2.028085e-01,
             1.807502e-01, 1.610705e-01, 1.435151e-01, 1.278563e-01,
             1.138910e-01, 1.014375e-01, 9.033344e-02)
pcr_det_100 <- (pcr_det/max(pcr_det))*1
pcr_det_95 <- (pcr_det/max(pcr_det))*0.95
pcr_det_90 <- (pcr_det/max(pcr_det))*0.9
plot(Time_from_infection,pcr_det, type = "l")
points(Time_from_infection,pcr_det_mort, type = "l")

saveRDS(pcr_det_100, file = "Analysis/Data/derived_data/pcr_det_hall_100.rds")
