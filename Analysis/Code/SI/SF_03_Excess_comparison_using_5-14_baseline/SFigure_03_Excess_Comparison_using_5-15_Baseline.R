library(ggplot2)
library(dplyr)
library(tidyr)
library(drjacoby)



### 5-15 mortality estimates

library(drjacoby)
library(dplyr)
library(tidyr)

## Format inputs
# Data
Burial_registrations <- readRDS(file = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Burial_registrations_by_week_2018_2021_5yr_age_groups.rds") |>
  filter(Week >= "2018-01-01", Week <"2021-06-21", Age_gr != 1) |>
  mutate(Age_gr = ifelse(Age_gr ==2, 3, Age_gr)) |>
  group_by(Week, Age_gr) |>
  summarise(BurRegs = sum(BurRegs)) |>
  mutate(Age_gr-2)

weekly_deaths_list_5to15 <- Burial_registrations |>
  filter(Age_gr==1) |> pull(BurRegs)

weekly_deaths_list_15_plus <- Burial_registrations |>
  filter(Age_gr != 1) |>
  pivot_wider(names_from = c(Week), values_from = BurRegs, values_fill = 0) %>%
  tibble::column_to_rownames(var="Age_gr")

data_list <- list(weekly_deaths_list_5to15 = weekly_deaths_list_5to15, weekly_deaths_list_15_plus = weekly_deaths_list_15_plus)

#### Parameters
age_cats <- 15

# define parameters for each of the age_rates
df_params <- define_params(name = c(paste0("Week_rate_5to15_",1:180), paste0("RR",2:age_cats)),
                           min =c(rep(0,length(weekly_deaths_list_5to15)),
                                  rep(0,age_cats-1)),
                           max = c(rep(200,length(weekly_deaths_list_5to15)),
                                   rep(10,age_cats-1)))

# define log-likelihood function
r_loglike <- function(params, data, misc) {

  age_cats <- misc$age_cats
  # Split data
  data_base <- data$weekly_deaths_list_5to15
  data_age_str <- data$weekly_deaths_list_15_plus

  # Split parameters: under 5 rate and age category variables
  base_w_rates <- as.numeric(params[1:180])
  rel_rates <- as.numeric(params[181:194])

  ret<-0
  ret <- ret + sum(dpois(x = data_base, lambda = base_w_rates, log=TRUE))
  for(j in 1:(age_cats-1)){
    ret <- ret + sum(dpois(x = as.numeric(data_age_str[j,1:104]), lambda = rel_rates[j]*base_w_rates[1:104], log=TRUE))
  }
  return(ret)
}


## define prior
r_logprior <- function(params, misc) {
  age_cats <- misc$age_cats

  # extract parameter values
  base_w_rates <- as.numeric(params[1:180])
  rel_rates <- as.numeric(params[181:194])

  # calculate log-prior
  ret <- 0

  # Add a prior for each of the age group relative risks
  ret <- ret + sum(dlnorm(x = rel_rates, meanlog = 0, sdlog = 50, log = T))
  ret <- ret + sum(dunif(x = base_w_rates, min = 0, max = 200, log = T))

  # return
  return(ret)
}

didehpc::web_login()
setwd("/home/rjs11/net/home/Cluster")
ctx_bl <- context::context_save("context_mcmc_baseline_5_15", package_sources = conan::conan_sources("mrc-ide/drjacoby"))
obj_bl <- didehpc::queue_didehpc(ctx_bl, config = didehpc::didehpc_config(cores = 5))

t_bl <- obj_bl$enqueue(drjacoby::run_mcmc(data = data_list,
                                          df_params = df_params,
                                          loglike = r_loglike,
                                          logprior = r_logprior,
                                          burnin = 5000,
                                          samples = 3000,
                                          pb_markdown = TRUE,
                                          chains = 5,
                                          cluster = parallel::makeCluster(5),
                                          misc = list(age_cats=age_cats)))


t_bl <- obj_bl$task_get(obj_bl$task_list())

# saveRDS(t_bl$result(),"~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/Baseline_Mortality_mcmc_5_15.rds")

###################################
###################################

### Compare 5-15 baseline predictions with 0-5 baseline predictsions
BM_0_5 <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/Baseline_Mortality_mcmc.rds")
BM_5_15 <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/Baseline_Mortality_mcmc_5_15.rds")

mcmc_samples <- BM_0_5$output %>% filter(phase =="sampling")

AG1_2020_mcmc <- mcmc_samples[, c(paste0("Week_rate_0to5_",c(1:180)))]
AG1_mcmc <- mcmc_samples[, c(paste0("Week_rate_0to5_",c(1:180)))]

Burial_df <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Burial_registrations_by_week_2018_2021_Fig4_age_groups.rds")

Burs_2018_2021 <- Burial_df %>%
  filter(Week >=as.Date("2018-01-01"),
         Week <as.Date("2021-06-12"))

Ag1std <- Burs_2018_2021 %>%
  rename("Ag1std" = BurRegs) %>%
  filter(Age_gr == 1) %>%
  mutate(Week_gr = 1:180)

Dates_df <- Ag1std %>% select(Week, Week_gr)

AG1_pre2020_mcmc <- mcmc_samples[, c(paste0("Week_rate_0to5_",1:104))]
U5_Baseline_2018_2019 <- rowMeans(AG1_pre2020_mcmc)
WeeklyStandardise <- apply(AG1_mcmc, 2, function(x){0.9*x/U5_Baseline_2018_2019})

colnames(WeeklyStandardise) <- 1:180
WeeklyStandardise <- WeeklyStandardise %>%
  reshape2::melt(value.name = "Standard", varnames = c("list_names","Week_gr"))

Pop_Str <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Population_structure_Lusaka_2020_CDC.rds")

Burs_2018_2021_all_age_groups <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Burial_registrations_by_week_2018_2021_5yr_age_groups.rds") %>%
  filter(Week >=as.Date("2018-01-01"),
         Week <as.Date("2021-06-12")) %>%
  merge(Dates_df)


Mort_excess_deaths <- lapply(1:17, function(y){
  print(y)
  tmp_Mort_excess <- lapply(1:180, function(x){
    if(y == 1){tmp_df <- AG1_mcmc[,paste0("Week_rate_0to5_",x)]} else {tmp_df <- AG1_mcmc[,paste0("Week_rate_0to5_",x)] * mcmc_samples[,paste0("RR",y)]}
    # browser()
    Excess <- Burs_2018_2021_all_age_groups %>% filter(Week_gr == x, Age_gr == y) %>% pull(Total_deaths) - tmp_df
    Excess_std <- Excess/(WeeklyStandardise %>% filter(Week_gr == x) %>% pull(Standard))

    return(data.frame(Sample = as.integer(rownames(AG1_mcmc)), Week_gr = x, Age_gr = y, Total = tmp_df, Excess = Excess, Excess_std = Excess_std, Pop = Pop_Str[y]))

  })
  tmp_Mort_excess <- do.call(rbind, tmp_Mort_excess)
})

Mort_excess_deaths <- do.call(rbind, Mort_excess_deaths)



#### 5-15
mcmc_samples_2 <- BM_5_15$output %>% filter(phase =="sampling")

AG1_2020_mcmc_2 <- mcmc_samples_2[, c(paste0("Week_rate_5to15_",c(1:180)))]
AG1_mcmc_2 <- mcmc_samples_2[, c(paste0("Week_rate_5to15_",c(1:180)))]

Ag1std_2 <- Burs_2018_2021 %>%
  rename("Ag1std" = BurRegs) %>%
  filter(Age_gr ==2) %>%
  mutate(Week_gr = 1:180)

Dates_df_2 <- Ag1std_2 %>% select(Week, Week_gr)

AG1_pre2020_mcmc_2 <- mcmc_samples_2[, c(paste0("Week_rate_5to15_",1:104))]
U5_Baseline_2018_2019_2 <- rowMeans(AG1_pre2020_mcmc_2)
WeeklyStandardise_2 <- apply(AG1_mcmc_2, 2, function(x){0.9*x/U5_Baseline_2018_2019_2})

colnames(WeeklyStandardise_2) <- 1:180
WeeklyStandardise_2 <- WeeklyStandardise_2 %>%
  reshape2::melt(value.name = "Standard", varnames = c("list_names","Week_gr"))

Pop_Str <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Population_structure_Lusaka_2020_CDC.rds")
Pop_Str_2 <- c(sum(Pop_Str[2:3]), Pop_Str[4:17])

Burs_2018_2021_all_age_groups_2 <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Burial_registrations_by_week_2018_2021_5yr_age_groups.rds") %>%
  filter(Age_gr != 1) %>%
  mutate(Age_gr = ifelse(Age_gr == 2, 3, Age_gr)) %>%
  group_by(Age_gr, Week) %>%
  summarise(BurRegs = sum(BurRegs)) %>%
  ungroup() %>%
  mutate(Age_gr = Age_gr - 2) %>%
  filter(Week >=as.Date("2018-01-01"),
         Week <as.Date("2021-06-12")) %>%
  merge(Dates_df_2)

Mort_excess_deaths_2 <- lapply(1:15, function(y){
  print(y)
  tmp_Mort_excess <- lapply(1:180, function(x){
    if(y == 1){tmp_df <- AG1_mcmc_2[,paste0("Week_rate_5to15_",x)]} else {tmp_df <- AG1_mcmc_2[,paste0("Week_rate_5to15_",x)] * mcmc_samples_2[,paste0("RR",y)]}

    Excess <- Burs_2018_2021_all_age_groups_2 %>% filter(Week_gr == x, Age_gr == y) %>% pull(Total_deaths) - tmp_df
    Excess_std <- Excess/(WeeklyStandardise_2 %>% filter(Week_gr == x) %>% pull(Standard))

    return(data.frame(Sample = as.integer(rownames(AG1_mcmc_2)), Week_gr = x, Age_gr = y, Total = tmp_df, Excess = Excess, Excess_std = Excess_std, Pop = Pop_Str_2[y]))

  })
  tmp_Mort_excess <- do.call(rbind, tmp_Mort_excess)
})

Mort_excess_deaths_2 <- do.call(rbind, Mort_excess_deaths_2)

Combined_excess <- Mort_excess_deaths_2 %>% rename(Excess_2 = Excess,
                                                   Excess_std_2 = Excess_std,
                                                   Pop_2 = Pop,
                                                   Total_2 = Total) %>%
  mutate(Age_gr = ifelse(Age_gr>1, Age_gr+2, Age_gr)) %>%
  filter(Week_gr >104) %>%
  merge(Mort_excess_deaths %>% filter(!Age_gr %in% 2:3, Week_gr >104))

saveRDS(list(Combined_excess=Combined_excess,
             Mort_excess_deaths_0_5 = Mort_excess_deaths,
             Mort_excess_deaths_5_15=Mort_excess_deaths_2), "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/5_15_comp_data.rds")


############################
############################
rm(list = ls())
gc()

Combined_excess <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/5_15_comp_data.rds")$Combined_excess %>%
  select(Sample, Week_gr, Age_gr, Excess, Excess_2, Pop, Pop_2) %>%
  mutate(Excess = 10000*Excess/Pop,
         Excess_2 = 10000*Excess_2/Pop_2) %>%
  select(-Pop, -Pop_2) %>%
  pivot_longer(cols = 4:5, names_to = "Age_base", values_to = "Excess") %>%
  group_by(Age_gr, Week_gr, Age_base) %>%
  summarise(Median_Excess = median(Excess),
            LCI_Excess = quantile(Excess, 0.025),
            HCI_Excess = quantile(Excess, 0.975))

Age_groups_labels <- c(paste0("Age group: ",seq(0,75,by =5),"-",seq(4,80,by =5)),"Age: 80+")
names(Age_groups_labels) <- 1:17

Dates_df <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Burial_registrations_by_week_2018_2021_5yr_age_groups.rds") |>
  select(Week) |>  unique() |>
  mutate(Week_gr =1:length(Week))

Comparison_Figure <- ggplot(Combined_excess %>% merge(Dates_df) %>% filter(Age_gr !=1), aes(x = Week, color = Age_base)) + geom_line(aes(y = Median_Excess)) +
  geom_ribbon(aes(ymin = LCI_Excess, ymax = HCI_Excess, fill = Age_base, linetype = NA), alpha = 0.4) +
  facet_wrap(~Age_gr, labeller = labeller(Age_gr =Age_groups_labels)) +
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  labs(x = "Date", y = "Excess burial registrations") +
  viridis::scale_color_viridis(name = "", discrete = T, option = "B", begin = 0.2, end = 0.6, labels = c("Using 0-4 baseline","Using 5-14 baseline")) +
  viridis::scale_fill_viridis(name = "", discrete = T, option = "B", begin = 0.2, end = 0.6, labels = c("Using 0-4 baseline","Using 5-14 baseline")) +
  scale_x_date(date_labels = "%b-%y") +
  theme(legend.position = c(0.75, 0.1), plot.title = element_text(face = "bold", size = 7)) +
  coord_cartesian(xlim = as.Date(c("2020-01-01", "2021-06-01")))

ggsave("Analysis/Figures/Supplementary_Figures/SFigure_03_Comparison_with_5_15_baseline.pdf",
       Comparison_Figure,
       width = 180, height = 180*6/7, units = "mm")

tiff("Analysis/Figures/Supplementary_Figures/SFigure_03_Comparison_with_5_15_baseline.tiff",
       width = 180, height = 180*6/7, units = "mm", res = 300)
Comparison_Figure
dev.off()
