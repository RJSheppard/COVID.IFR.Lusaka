####################################
###### Supplementary Figure 9 ######
####################################
# LL of baseline burial regsitration distribution assumptions
library(drjacoby);library(dplyr);library(tidyr)

## Format inputs
# Data
UTH_Mortality_Total <- read.csv(file = "~/Documents/Zambia/covid-mortality-ascertainment/analysis/data/raw/BMJ_UTH_excess_mortality/mortuary_records_v2.csv")

# UTH_Mortality_Total%>%   mutate(date = as.Date(dod, "%m/%d/%y")) %>%
# filter(date >= "2018-08-13", date <= "2018-08-19")

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
df_params <- define_params(name = c(paste0("Week_rate_0to5_",1:180), paste0("RR",2:age_cats)),
                           min =c(rep(0,length(weekly_deaths_list_0to5)), rep(0,age_cats-1)),
                           max = c(rep(200,length(weekly_deaths_list_0to5)), rep(10,age_cats-1)))

# define log-likelihood function
r_loglike <- function(params, data, misc) {

  age_cats <- misc$age_cats
  # Split data
  data_base <- data$weekly_deaths_list_0to5
  data_age_str <- data$weekly_deaths_list_5_plus

  # Split parameters: under 5 rate and age category variables
  base_w_rates <- as.numeric(params[1:180])
  rel_rates <- as.numeric(params[181:196])

  ret<-0
  ret <- ret + sum(dnbinom(x = data_base, size = misc$size_parm, mu = base_w_rates, log=TRUE))
  for(j in 1:(age_cats-1)){
    ret <- ret + sum(dnbinom(x = as.numeric(data_age_str[j,1:104]), size = misc$size_parm, mu = rel_rates[j]*base_w_rates[1:104], log=TRUE))
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
  ret <- ret + sum(dlnorm(x = rel_rates, meanlog = 0, sdlog = 50, log = T))
  ret <- ret + sum(dunif(x = base_w_rates, min = 0, max = 200, log = T))

  # return
  return(ret)
}


Test <- drjacoby::run_mcmc(data = data_list,
                           df_params = df_params,
                           loglike = r_loglike,
                           logprior = r_logprior,
                           burnin = 5,
                           samples = 5,
                           pb_markdown = TRUE,
                           chains = 5,
                           cluster = parallel::makeCluster(5),
                           misc = list(age_cats = age_cats,
                                       size_parm = 20))


# didehpc::didehpc_config()
didehpc::web_login()

# setwd("/home/rjs11/net/home/Cluster")
ctx_bl <- context::context_save("context_mcmc_baseline_poisson", package_sources = conan::conan_sources("mrc-ide/drjacoby"))
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
                                          misc = list(age_cats=age_cats,
                                                      size_parm = 2)))

t_bl_5 <- obj_bl$enqueue(drjacoby::run_mcmc(data = data_list,
                                            df_params = df_params,
                                            loglike = r_loglike,
                                            logprior = r_logprior,
                                            burnin = 5000,
                                            samples = 3000,
                                            pb_markdown = TRUE,
                                            chains = 5,
                                            cluster = parallel::makeCluster(5),
                                            misc = list(age_cats=age_cats,
                                                        size_parm = 5)))

t_bl_20 <- obj_bl$enqueue(drjacoby::run_mcmc(data = data_list,
                                             df_params = df_params,
                                             loglike = r_loglike,
                                             logprior = r_logprior,
                                             burnin = 5000,
                                             samples = 3000,
                                             pb_markdown = TRUE,
                                             chains = 5,
                                             cluster = parallel::makeCluster(5),
                                             misc = list(age_cats=age_cats,
                                                         size_parm = 20)))

t_bl_50 <- obj_bl$enqueue(drjacoby::run_mcmc(data = data_list,
                                             df_params = df_params,
                                             loglike = r_loglike,
                                             logprior = r_logprior,
                                             burnin = 5000,
                                             samples = 3000,
                                             pb_markdown = TRUE,
                                             chains = 5,
                                             cluster = parallel::makeCluster(5),
                                             misc = list(age_cats=age_cats,
                                                         size_parm = 50)))

t_bl_1 <- obj_bl$task_get(obj_bl$task_list()[1])
t_bl_2 <- obj_bl$task_get(obj_bl$task_list()[2])
t_bl_3 <- obj_bl$task_get(obj_bl$task_list()[3])
t_bl_4 <- obj_bl$task_get(obj_bl$task_list()[4])


saveRDS(t_bl_1$result(),"~/Documents/Zambia/Bonus Files/2023_March_Reviewer_Batch/2023-03-07_Baseline_Mortality_mcmc_negbin_like_size_2.rds")
saveRDS(t_bl_2$result(),"~/Documents/Zambia/Bonus Files/2023_March_Reviewer_Batch/2023-03-07_Baseline_Mortality_mcmc_negbin_like_size_5.rds")
saveRDS(t_bl_3$result(),"~/Documents/Zambia/Bonus Files/2023_March_Reviewer_Batch/2023-03-07_Baseline_Mortality_mcmc_negbin_like_size_20.rds")
saveRDS(t_bl_4$result(),"~/Documents/Zambia/Bonus Files/2023_March_Reviewer_Batch/2023-03-07_Baseline_Mortality_mcmc_negbin_like_size_50.rds")


Res_2 <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/2023-03-07_Baseline_Mortality_mcmc_negbin_like_size_2.rds")
Res_5 <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/2023-03-07_Baseline_Mortality_mcmc_negbin_like_size_5.rds")
Res_20 <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/2023-03-07_Baseline_Mortality_mcmc_negbin_like_size_20.rds")
Res_50 <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/2023-03-07_Baseline_Mortality_mcmc_negbin_like_size_50.rds")
Res_pois <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/2023-03-07_Baseline_Mortality_mcmc.rds")


dispersion_df <- data.frame(
  Median = unlist(lapply(list(Res_2,Res_5,Res_20,Res_50,Res_pois), function(x){median(x$output[Res_2$output$phase!="burnin","loglikelihood"])})),
  disp = c("NB: Size = 2","NB: Size = 5","NB: Size = 20","NB: Size = 50","Poisson"),
  ui=unlist(lapply(list(Res_2,Res_5,Res_20,Res_50,Res_pois), function(x){quantile(x$output[x$output$phase!="burnin","loglikelihood"], 0.975)})),
  li=unlist(lapply(list(Res_2,Res_5,Res_20,Res_50,Res_pois), function(x){quantile(x$output[x$output$phase!="burnin","loglikelihood"], 0.025)}))) |>
  mutate(disp=factor(disp, levels=disp))   # This trick update the factor levels


pois_nbinom_fig <- ggplot(dispersion_df, aes(x = disp, y = Median)) +
  geom_point(size = 0.8) +
  geom_errorbar(aes(x = disp, ymin = li, ymax = ui), width = 0.1) +
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme(plot.title = element_text(size = 7, face = "bold")) +
  labs(x = element_blank(),
       y = "Log likelihood")

ggsave("Analysis/Figures/Supplementary_Figures/SFigure_09_Burial_registration_distribution_assumptions.pdf",
       pois_nbinom_fig,
       width = 180, height = 180*4/10, units = "mm")


tiff("Analysis/Figures/Supplementary_Figures/SFigure_09_Burial_registration_distribution_assumptions.tiff", res = 300, width = 180, height = 180*4/10, units = "mm")
pois_nbinom_fig
dev.off()
