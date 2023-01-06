library(reshape2)
library(ggplot2)
library(dplyr)

############################################
##'[drjacoby MCMC K-fold Cross validation]##
############################################
UTH_Mortality_Total <- read.csv(file = "Analysis/Data/raw_data/11_mortuary_records_v2.csv")

weekly_deaths_list_5_plus_loocv <- UTH_Mortality_Total %>%
  mutate(date = as.Date(dod, "%m/%d/%y")) %>%
  filter(dod != ".") %>%
  filter(date >= "2018-01-01",date < "2019-12-30", age_years !=".") %>%
  mutate(Age_gr = cut(as.numeric(age_years), c(seq(0,80,by = 5),Inf), right = F, labels = c(1:17)),
         week = cut.Date(date, breaks = "1 week", start.on.monday = T, labels = c(1:104))) %>%
  group_by(Age_gr, week) %>%
  summarise(total_deaths = length(date)) %>% ungroup() %>%
  complete(., Age_gr, week, fill = list(total_deaths = 0)) %>%
  arrange(Age_gr, week) %>%
  filter(Age_gr !=1) %>%
  ungroup() %>%
  pivot_wider(names_from = c(week), values_from = total_deaths, values_fill = 0) %>%
  tibble::column_to_rownames(var="Age_gr")

Under_5s_deaths_loocv <- UTH_Mortality_Total %>%
  mutate(date = as.Date(dod, "%m/%d/%y")) %>%
  filter(dod != ".") %>%
  filter(date >= "2018-01-01",date < "2019-12-30", age_years !=".") %>%
  mutate(Age_gr = cut(as.numeric(age_years), c(seq(0,80,by = 5),Inf), right = F, labels = F),
         year = as.numeric(format(date,"%Y")),
         week = as.numeric(cut.Date(date, breaks = "1 week", start.on.monday = T, labels = F))) %>%
  complete(expand(., week, Age_gr)) %>%
  ungroup() %>% group_by(Age_gr, week) %>%
  filter(Age_gr ==1) %>%
  summarise(total_deaths = as.numeric(sum(!is.na(date))), st_date = min(date)) %>%
  ungroup() %>%
  select(-Age_gr, -st_date) %>%
  arrange(week)

data_list_loocv <- list(Under_5s_deaths = weekly_deaths_list_0to5[1:104], weekly_deaths_list_5_plus = weekly_deaths_list_5_plus[,1:104])

df_params_loocv <- define_params(name = c(paste0("U_5_Rate_Week",1:nrow(Under_5s_deaths_loocv)), paste0("RR",2:age_cats)),
                                 min =c(rep(0,nrow(Under_5s_deaths_loocv)), rep(0,age_cats-1)), max = c(rep(200,nrow(Under_5s_deaths_loocv)),rep(10,age_cats-1)))


r_loglike_loocv <- function(params, data, misc) {

  age_cats <- misc$age_cats
  Test_weeks <- misc$Test_weeks

  # Split data
  data_base <- data$Under_5s_deaths
  data_age_str <- data$weekly_deaths_list_5_plus[,-Test_weeks]

  # browser()
  # Split parameters: under 5 rate and age category variables
  base_w_rates <- as.numeric(params[1:104])
  rel_rates <- as.numeric(params[105:120])

  ret<-0
  ret <- ret + sum(dpois(x = data_base, lambda = base_w_rates, log=TRUE))
  for(j in 1:(age_cats-1)){
    ret <- ret + sum(dpois(x = as.numeric(data_age_str[j,]), lambda = rel_rates[j]*base_w_rates[-Test_weeks], log=TRUE))
  }

  return(ret)
}


r_logprior_loocv <- function(params, misc) {

  age_cats <- misc$age_cats

  # extract parameter values
  rel_rates <- as.numeric(params[105:120])
  base_w_rates <- as.numeric(params[1:104])

  # calculate log-prior
  # Add a prior for each of the age group relative risks
  ret <- 0
  ret <- ret + sum(dlnorm(x = rel_rates, meanlog = 1, sdlog = 10, log = T))
  ret <- ret + sum(dunif(x = base_w_rates, min = 0, max = 200, log = T))

  # return
  return(ret)
}

set.seed(130)
Test_weeks_loocv <- as.list(as.data.frame(matrix(sample(1:104, 104, replace = F), nrow = 4)))
names(Test_weeks_loocv) <- NULL

didehpc::web_login()
ctx_bl_loocv <- context::context_save("context_mcmc_baseline_kfcv", sources = "Analysis/Code/Code_Functions/drjacoby_MCMC_K-fold_Cross_validaion.R", package_sources = conan::conan_sources("mrc-ide/drjacoby"))
obj_bl_loocv <- didehpc::queue_didehpc(ctx_bl_loocv)
t_loocv <- obj_bl_loocv$lapply(Test_weeks_loocv, dj_loocv,
                               df_params = df_params_loocv,
                               r_loglike = r_loglike_loocv,
                               r_logprior = r_logprior_loocv,
                               data_list = data_list_loocv,
                               age_cats = age_cats,
                               burnin = 500,
                               samples = 2500,
                               chains = 5)

# saveRDS(t_loocv$results(),"Results~/Documents/Zambia/Bonus Files/2022-12-08_Baseline_Mortality_mcmc_loocv_Gamma_Prior_randomised.rds")

LOOCV <- t_loocv$results()

###################
##'[Plot Results]##
###################

LOOCV_Ests <- lapply(1:nrow(LOOCV[[1]]$output[LOOCV[[1]]$output$phase=="sampling",]), function(y){
  List_rep <- lapply(1:length(Test_weeks), function(x){
    Lamdas <- rbind(t(as.matrix(unlist(LOOCV[[x]]$output[LOOCV[[x]]$output$phase=="sampling", paste0("U_5_Rate_Week",Test_weeks[[x]])][y,]))),
                    as.matrix(unlist(LOOCV[[x]]$output[LOOCV[[x]]$output$phase=="sampling", paste0("RR",2:17)][y,])) %*%
                      t(as.matrix(unlist(LOOCV[[x]]$output[LOOCV[[x]]$output$phase=="sampling", paste0("U_5_Rate_Week",Test_weeks[[x]])][y,]))))
    rownames(Lamdas) <- 1:17
    samples_from_lambdas <- apply(Lamdas, MARGIN = c(1,2), FUN = function(x){rpois(n = 1, lambda = x)})
  })
  do.call(cbind.data.frame, List_rep)
})

LOOCV_Ests_summed <- lapply(LOOCV_Ests, function(x){
  str2str::lv2d(lapply(Test_weeks, function(y){rowSums(x[,paste0("U_5_Rate_Week",y)])}), along = 1)
})

LOOCV_Ests_summed_comb <- str2str::ld2a(LOOCV_Ests_summed)

LOOCV_Ests_Av <- apply(LOOCV_Ests_summed_comb, 2, rowMeans)

LOOCV_Ests_median <- apply(LOOCV_Ests_summed_comb, 2, function(x){
  apply(x, 1, median)
})

LOOCV_Ests_CI_low <- apply(LOOCV_Ests_summed_comb, 2, function(x){
  apply(x, 1, function(y){bayestestR::ci(y)$CI_low})
})

LOOCV_Ests_CI_high <- apply(LOOCV_Ests_summed_comb, 2, function(x){
  apply(x, 1, function(y){bayestestR::ci(y)$CI_high})
})

Test_weeks_match <- str2str::ld2d(lapply(1:26, function(x){
  as.data.frame(cbind(x,Test_weeks[[x]]))
})) %>% select(x, V2) %>% rename(Week_gr = V2, Test_weeks_set = x)

UTH_Mortality <- read.csv(file = "Analysis/Data/raw_data/11_mortuary_records.csv") %>%
  mutate(date = as.Date(dod, "%m/%d/%y")) %>%
  filter(date >= "2018-01-01" & date < "2019-12-30", age_years !=".", dod != ".") %>%
  mutate(Age_gr = cut(as.numeric(age_years), c(seq(0,80,by = 5),Inf), right = F, labels = F),
         Week_gr = cut.Date(date, breaks = "1 week", labels = FALSE, start.on.monday = T),
         Week_gr_date_start = lubridate::floor_date(date, "week", week_start = 1)) %>%
  tidyr::complete(Age_gr, Week_gr, fill = list(Mort_deaths = 0)) %>%
  group_by(Week_gr,Age_gr,Week_gr_date_start) %>%
  summarise(Mort_deaths = length(date)) %>%
  ungroup %>%
  tidyr::complete(Week_gr, Age_gr, fill = list(Mort_deaths = 0)) %>%
  group_by(Week_gr) %>% mutate(Week_gr_date_start = na.omit(unique(Week_gr_date_start)))

UTH_Mortality_Total <- UTH_Mortality %>% merge(Test_weeks_match) %>%
  group_by(Test_weeks_set, Age_gr) %>%
  summarise(Mort_deaths = sum(Mort_deaths)) %>%
  rename(Week_gr = Test_weeks_set)

Plot_Res <- melt(LOOCV_Ests_median, varnames = c("Week_gr", "Age_gr"), value.name = "Deaths_median") %>%
  mutate(Age_gr = as.numeric(gsub(Age_gr, pattern = "RR", replacement = "")),
         Week_gr = as.numeric(gsub(x = Week_gr, "U_5_Rate_Week", ""))) %>%
  merge(melt(LOOCV_Ests_CI_low, varnames = c("Week_gr", "Age_gr"), value.name = "Deaths_low_CI") %>%
          mutate(Age_gr = as.numeric(Age_gr),
                 Week_gr = as.numeric(gsub(x = Week_gr, "U_5_Rate_Week", "")))) %>%
  merge(melt(LOOCV_Ests_CI_high, varnames = c("Week_gr", "Age_gr"), value.name = "Deaths_high_CI") %>%
          mutate(Age_gr = as.numeric(Age_gr),
                 Week_gr = as.numeric(gsub(x = Week_gr, "U_5_Rate_Week", "")))) %>%
  merge(UTH_Mortality_Total)


Age_groups.labs <- c(paste0(c("0-4","5-9","10-14","15-29","20-24","25-29",
                              "30-34","35-39","40-44","45-49","50-54","55-59",
                              "60-64","65-69","70-74","75-79","80+")))
names(Age_groups.labs) <- 1:17

####################################
## Get Medians from data
Median_Baseline_2018_2019 <- UTH_Mortality %>% group_by(Age_gr) %>%
  summarise(median_baseline = median(Mort_deaths)*4)


tiff("Analysis/Figures/SI_KFCV_combined_split.tiff", units = "in", res = 300, height = 6, width = 7)
p3 + #xlim(c(0,165)) + ylim(c(0,160)) +
  geom_hline(data = Median_Baseline_2018_2019 %>% filter(Age_gr !=1), aes(yintercept = median_baseline), color = "darkgrey", size = 0.35) +
  geom_point(color = "black", size = 0.5) +
  geom_point(aes(color = as.factor(Age_gr)), size = 0.3) +
  facet_wrap(~Age_gr, scales = "free", labeller = labeller(Age_gr = Age_groups.labs)) +
  theme(legend.position = "none", title = element_blank(), axis.title = element_text(), plot.margin = margin(5,10,5,5, "points"))
dev.off()
