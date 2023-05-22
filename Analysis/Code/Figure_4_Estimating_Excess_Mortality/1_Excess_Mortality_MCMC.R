library(drjacoby)
library(dplyr)
library(tidyr)

## Format inputs
# Burial Registration Data
Burial_registrations <- readRDS(file = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Burial_registrations_by_week_2018_2021_5yr_age_groups.rds")

weekly_BurRegs_list_0to4 <- Burial_registrations |>
  filter(Age_gr==1) |> pull(BurRegs)

weekly_BurRegs_list_5_plus <- Burial_registrations |>
  filter(Age_gr != 1) |>
  pivot_wider(names_from = c(Week), values_from = BurRegs, values_fill = 0) %>%
  tibble::column_to_rownames(var="Age_gr")

data_list <- list(weekly_BurRegs_list_0to4 = weekly_BurRegs_list_0to4, weekly_BurRegs_list_5_plus = weekly_BurRegs_list_5_plus)

#### Parameters
age_cats <- 17

# define parameters for each of the age_rates
df_params <- define_params(name = c(paste0("Week_rate_0to4_",1:180), paste0("RR",2:age_cats)),
                           min =c(rep(0,length(weekly_BurRegs_list_0to4)),
                                  rep(0,age_cats-1)),
                           max = c(rep(200,length(weekly_BurRegs_list_0to4)),
                                   rep(10,age_cats-1)))

# define log-likelihood function
r_loglike <- function(params, data, misc) {

  age_cats <- misc$age_cats

  # Split data
  data_base <- data$weekly_BurRegs_list_0to4
  data_age_str <- data$weekly_BurRegs_list_5_plus

  # Split parameters: under 5 rate and age category variables
  base_w_rates <- as.numeric(params[1:180])
  rel_rates <- as.numeric(params[181:196])

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
  rel_rates <- as.numeric(params[181:196])

  # calculate log-prior
  ret <- 0

  # Add a prior for each of the age group relative risks
  ret <- ret + sum(dlnorm(x = rel_rates, meanlog = 0, sdlog = 50, log = T))
  ret <- ret + sum(dunif(x = base_w_rates, min = 0, max = 200, log = T))

  # return
  return(ret)
}



#### HPC Cluster run code
# didehpc::web_login()
# ctx_bl <- context::context_save("context_mcmc_baseline", package_sources = conan::conan_sources("mrc-ide/drjacoby"))
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
# saveRDS(t_bl$result(),"Data/derived_data/Baseline_Mortality_MCMC_Gamma_Prior_inc_Feb_2021.rds")

## Running just a subsample for testing
MCMC_Baseline_Excess_Burials <- drjacoby::run_mcmc(data = data_list,
                                          df_params = df_params,
                                          loglike = r_loglike,
                                          logprior = r_logprior,
                                          burnin = 50,
                                          samples = 50,
                                          pb_markdown = TRUE,
                                          chains = 5,
                                          misc = list(age_cats=age_cats))

saveRDS(MCMC_Baseline_Excess_Burials,"~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/Baseline_Mortality_MCMC.rds")

## Prepare baseline estimates for squire fitting
devtools::load_all()
MCMC_Baseline_Excess_Burials <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/Baseline_Mortality_mcmc.rds")
# MCMC_Baseline_Excess_Burials <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/2023-03-07_Baseline_Mortality_mcmc.rds")
mcmc_samples <-MCMC_Baseline_Excess_Burials$output %>% filter(phase =="sampling")

set.seed(85)
sample_index <- sample(x = 1:nrow(mcmc_samples), size = 10000, replace = F)
dfj_mcmc_data <- Get_baselines(mcmc_samples[sample_index[1:4000],], "Week_rate_0to5_","RR")
saveRDS(turn_list_into_array(dfj_mcmc_data), file = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/Baseline_Mortality_MCMC_formatted_for_squire.rds")
# saveRDS(dfj_mcmc_data, file = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/Code-generated-data/00_16_03_drj_mcmc_data_new_pop_str_2_2023_March.rds")


## Scale baseline estimates for squire fitting sensitivity analyses
dfj_mcmc_data_lower_10 <- lapply(dfj_mcmc_data, function(x){x %>%
    mutate(Bg_dr = ifelse(Age_gr != 1, 0.9 * Bg_dr, Bg_dr),
           Mort_ncd_mcmc = ifelse(Age_gr != 1, 0.9 * Mort_ncd_mcmc, Mort_ncd_mcmc))})

dfj_mcmc_data_lower_20 <- lapply(dfj_mcmc_data, function(x){x %>%
    mutate(Bg_dr = ifelse(Age_gr != 1, 0.8 * Bg_dr, Bg_dr),
           Mort_ncd_mcmc = ifelse(Age_gr != 1, 0.8 * Mort_ncd_mcmc, Mort_ncd_mcmc))})

dfj_mcmc_data_higher_10 <- lapply(dfj_mcmc_data, function(x){x %>%
    mutate(Bg_dr = ifelse(Age_gr != 1, 1.1*Bg_dr, Bg_dr),
           Mort_ncd_mcmc = ifelse(Age_gr != 1, 1.1*Mort_ncd_mcmc, Mort_ncd_mcmc))})

dfj_mcmc_data_higher_20 <- lapply(dfj_mcmc_data, function(x){x %>%
    mutate(Bg_dr = ifelse(Age_gr != 1, 1.2*Bg_dr, Bg_dr),
           Mort_ncd_mcmc = ifelse(Age_gr != 1, 1.2*Mort_ncd_mcmc, Mort_ncd_mcmc))})

saveRDS(turn_list_into_array(dfj_mcmc_data_lower_10), file = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/Baseline_Mortality_MCMC_formatted_for_squire_l10.rds")
saveRDS(turn_list_into_array(dfj_mcmc_data_higher_10), file = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/Baseline_Mortality_MCMC_formatted_for_squire_h10.rds")
saveRDS(turn_list_into_array(dfj_mcmc_data_lower_20), file = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/Baseline_Mortality_MCMC_formatted_for_squire_l20.rds")
saveRDS(turn_list_into_array(dfj_mcmc_data_higher_20), file = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/Baseline_Mortality_MCMC_formatted_for_squire_h20.rds")
