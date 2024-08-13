orderly2::orderly_description(display = "Perform excess mortality MCMC")

# orderly2::orderly_artefact(description = "MCMC results",
#                            files = c(".rds"))

orderly2::orderly_parameters(derived_data_loc = NULL,
                             baseline_group = NULL)

library(drjacoby)
library(dplyr)
library(tidyr)

## Format inputs
# Burial Registration Data
Burial_registrations_cleaned <- readRDS(file = paste0(derived_data_loc, "Burial_registrations_cleaned.rds"))

Age_structured_Burial_Registrations <- Burial_registrations_cleaned %>%
  group_by(Age_gr, Week) %>%
  summarise(BurRegs = n()) %>%
  ungroup() %>%
  complete(Age_gr, Week, fill = list(BurRegs = 0)) %>%
  arrange(Age_gr, Week) |>
  filter(Week >= "2018-01-01")

if(baseline_group == "0-4"){

  # 0-4
  baseline_group_df <- Age_structured_Burial_Registrations |>
    filter(Age_gr==1) |> pull(BurRegs)

  # 5+
  older_ages <- Age_structured_Burial_Registrations |>
    filter(Age_gr != 1) |>
    pivot_wider(names_from = c(Week), values_from = BurRegs, values_fill = 0) %>%
    tibble::column_to_rownames(var="Age_gr")


} else if (baseline_group == "5-14"){

  # 5-14
  baseline_group_df <- Burial_registrations |>
    filter(Age_gr %in% 2:3) |>
    group_by(Week) |>
    summarise(BurRegs = sum(BurRegs)) |>
    pull(BurRegs)

  # 15+
  older_ages <- Age_structured_Burial_Registrations |>
    filter(!Age_gr %in% 1:3) |>
    pivot_wider(names_from = c(Week), values_from = BurRegs, values_fill = 0) %>%
    tibble::column_to_rownames(var="Age_gr")

}

data_list <- list(baseline_group_df = baseline_group_df, older_ages = older_ages)

#### Parameters
age_cats <- nrow(older_ages) + 1
n_weeks <- length(baseline_group_df)
pp_weeks <- 104

# define parameters for each of the age_rates
df_params <- define_params(name = c(paste0("baseline_rate_week_",1:n_weeks), paste0("RR",2:age_cats)),
                           min =c(rep(0,length(baseline_group_df)),
                                  rep(0,age_cats-1)),
                           max = c(rep(200,length(baseline_group_df)),
                                   rep(10,age_cats-1)))

# define log-likelihood function
r_loglike <- function(params, data, misc) {
browser()
  age_cats <- misc$age_cats
  n_weeks <- misc$n_weeks
  pp_weeks <- misc$pp_weeks

  # Split data
  data_base <- data$baseline_group_df
  data_age_str <- data$older_ages

  # Split parameters: under 5 rate and age category variables
  base_w_rates <- as.numeric(params[1:n_weeks])
  rel_rates <- as.numeric(params[(n_weeks+1):(n_weeks+age_cats-1)])

  ret<-0
  ret <- ret + sum(dpois(x = data_base, lambda = base_w_rates, log=TRUE))
  for(j in 1:(age_cats-1)){
    ret <- ret + sum(dpois(x = as.numeric(data_age_str[j,1:pp_weeks]), lambda = rel_rates[j]*base_w_rates[1:pp_weeks], log=TRUE))
  }
  return(ret)
}


## define prior
r_logprior <- function(params, misc) {
browser()
  age_cats <- misc$age_cats
  n_weeks <- misc$n_weeks

  # extract parameter values
  base_w_rates <- as.numeric(params[1:n_weeks])
  rel_rates <- as.numeric(params[(n_weeks+1):(n_weeks+age_cats-1)])

  # calculate log-prior
  ret <- 0

  # Add a prior for each of the age group relative risks
  ret <- ret + sum(dlnorm(x = rel_rates, meanlog = 0, sdlog = 50, log = T))
  ret <- ret + sum(dunif(x = base_w_rates, min = 0, max = 200, log = T))

  # return
  return(ret)
}

#
MCMC_Baseline_Excess_Burials <- drjacoby::run_mcmc(data = data_list,
                                                   df_params = df_params,
                                                   loglike = r_loglike,
                                                   logprior = r_logprior,
                                                   burnin = 5e3,
                                                   samples = 2.5e3,
                                                   pb_markdown = TRUE,
                                                   chains = 5,
                                                   misc = list(age_cats=age_cats,
                                                               n_weeks=n_weeks,
                                                               pp_weeks = pp_weeks
                                                   ))

saveRDS(MCMC_Baseline_Excess_Burials, paste0(derived_data_loc, "Baseline_Mortality_MCMC_2020_2022.rds"))
