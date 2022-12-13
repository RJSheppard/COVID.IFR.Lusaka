library(drjacoby)
library(dplyr)
library(tidyr)

## Format inputs
# Burial Registration Data
UTH_Mortality_Total <- read.csv(file = "~/Documents/Zambia/covid-mortality-ascertainment/analysis/data/raw/BMJ_UTH_excess_mortality/mortuary_records_v2.csv")

weekly_deaths_list <- UTH_Mortality_Total %>%
  mutate(date = as.Date(dod, "%m/%d/%y")) %>%
  filter(dod != ".") %>%
  filter(date >= "2018-01-01", date<"2021-06-21", age_years !=".") %>%
  mutate(age_years = as.numeric(age_years)) %>%
  mutate(Age_gr = cut(age_years, c(seq(0,80,by = 5),Inf), right = F, labels = c(1:17)),
         week = cut.Date(date, breaks = "1 week", start.on.monday = T, labels = F)) %>%
  group_by(Age_gr, week) %>%
  summarise(total_deaths = n()) %>% ungroup() %>%
  complete(., Age_gr, week, fill = list(total_deaths = 0)) %>%
  arrange(Age_gr, week) %>%
  # filter(Age_gr !=1) %>%
  ungroup() %>%
  pivot_wider(names_from = c(week), values_from = total_deaths, values_fill = 0) %>%
  tibble::column_to_rownames(var="Age_gr")

weekly_deaths_list_0to5 <- unlist(weekly_deaths_list[1,])
weekly_deaths_list_5_plus <- weekly_deaths_list[-1,1:104]

data_list <- list(weekly_deaths_list_0to5 = weekly_deaths_list_0to5, weekly_deaths_list_5_plus = weekly_deaths_list_5_plus)

#### Parameters
age_cats <- 17
# define parameters for each of the age_rates
df_params <- define_params(name = c(paste0("Week_rate_0to5_",1:181), paste0("RR",2:age_cats)),
                           min =c(rep(0,length(weekly_deaths_list_0to5)),
                                  rep(0,age_cats-1)),
                           max = c(rep(200,length(weekly_deaths_list_0to5)),
                                   rep(10,age_cats-1)))

# define log-likelihood function
r_loglike <- function(params, data, misc) {

  age_cats <- misc$age_cats

  # Split data
  data_base <- data$weekly_deaths_list_0to5
  data_age_str <- data$weekly_deaths_list_5_plus

  # Split parameters: under 5 rate and age category variables
  base_w_rates <- as.numeric(params[1:181])
  rel_rates <- as.numeric(params[182:197])

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
  base_w_rates <- as.numeric(params[1:181])
  rel_rates <- as.numeric(params[182:197])

  # calculate log-prior
  ret <- 0

  # Add a prior for each of the age group relative risks
  ret <- ret + sum(dlnorm(x = rel_rates, meanlog = 1, sdlog = 10, log = T))
  ret <- ret + sum(dunif(x = base_w_rates, min = 0, max = 200, log = T))

  # return
  return(ret)
}



#### HPC Cluster run code
# didehpc::web_login()
# ctx_bl <- context::context_save("context_mcmc_baseline_05", package_sources = conan::conan_sources("mrc-ide/drjacoby"))
# obj_bl <- didehpc::queue_didehpc(ctx_bl)
# t_bl <- obj_bl$enqueue(drjacoby::run_mcmc(data = data_list,
#                                           df_params = df_params,
#                                           loglike = r_loglike,
#                                           logprior = r_logprior,
#                                           burnin = 5e2,
#                                           samples = 2.5e3,
#                                           pb_markdown = TRUE,
#                                           chains = 5,
#                                           misc = list(age_cats=age_cats)))
# obj_bl$task_list()
# saveRDS(t_bl$result(),"Baseline_Mortality_MCMC_Gamma_Prior_inc_Feb_2021.rds")

MCMC_Baseline_Excess_Burials <- drjacoby::run_mcmc(data = data_list,
                                          df_params = df_params,
                                          loglike = r_loglike,
                                          logprior = r_logprior,
                                          burnin = 5e2,
                                          samples = 2.5e3,
                                          pb_markdown = TRUE,
                                          chains = 5,
                                          misc = list(age_cats=age_cats))

saveRDS(MCMC_Baseline_Excess_Burials,"Baseline_Mortality_MCMC_Gamma_Prior_inc_Feb_2021.rds")
