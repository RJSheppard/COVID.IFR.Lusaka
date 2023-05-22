library(reshape2)
library(ggplot2)
library(dplyr)
library(drjacoby)

############################################
##'[drjacoby MCMC K-fold Cross validation]##
############################################
Burial_registrations <- readRDS(file = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Burial_registrations_by_week_2018_2021_5yr_age_groups.rds") |>
  filter(Week < "2019-12-30")

weekly_BurRegs_list_0to4 <- Burial_registrations |>
  filter(Age_gr==1) |> pull(BurRegs)

weekly_BurRegs_list_5_plus <- Burial_registrations |>
  filter(Age_gr != 1) |>
  pivot_wider(names_from = c(Week), values_from = BurRegs, values_fill = 0) %>%
  tibble::column_to_rownames(var="Age_gr")

data_list <- list(weekly_BurRegs_list_0to4 = weekly_BurRegs_list_0to4, weekly_BurRegs_list_5_plus = weekly_BurRegs_list_5_plus)

age_cats <- 17

df_params_loocv <- define_params(name = c(paste0("U_5_Rate_Week",1:length(weekly_BurRegs_list_0to4)), paste0("RR",2:age_cats)),
                                 min =c(rep(0,length(weekly_BurRegs_list_0to4)), rep(0,age_cats-1)), max = c(rep(200,length(weekly_BurRegs_list_0to4)),rep(10,age_cats-1)))


r_loglike_loocv <- function(params, data, misc) {

  age_cats <- misc$age_cats
  Test_weeks <- misc$Test_weeks

  # Split data
  data_base <- data$Under_5s_deaths
  data_age_str <- data$weekly_deaths_list_5_plus[,-Test_weeks]

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
  ret <- ret + sum(dlnorm(x = rel_rates, meanlog = 0, sdlog = 10, log = T))
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

# saveRDS(t_loocv$results(),"~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/Baseline_Mortality_K_Fold_CV.rds")

LOOCV <- t_loocv$results()

###################
##'[Plot Results]##
###################

library(reshape2)
library(ggplot2)
library(dplyr)

LOOCV <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/Baseline_Mortality_K_fold_CV.rds")

set.seed(130)
Test_weeks <- as.list(as.data.frame(matrix(sample(1:104, 104, replace = F), nrow = 4)))
names(Test_weeks) <- NULL

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

Age_Week_Registrations <- readRDS(file = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Burial_registrations_by_week_2018_2021_5yr_age_groups.rds") |>
  filter(Week < "2019-12-30") |>
  group_by(Age_gr) |>
  mutate(Week_gr = 1:104) |>
  merge(Test_weeks_match) %>%
  group_by(Test_weeks_set, Age_gr) %>%
  summarise(BurRegs = sum(BurRegs)) %>%
  rename(Week_gr = Test_weeks_set)
#read.csv(file = "analysis/data/raw/BMJ_UTH_excess_mortality/mortuary_records_v2.csv") %>%
  # mutate(date = as.Date(dod, "%m/%d/%y")) %>%
  # filter(date >= "2018-01-01" & date < "2019-12-30", age_years !=".", dod != ".") %>%
  # mutate(Age_gr = cut(as.numeric(age_years), c(seq(0,80,by = 5),Inf), right = F, labels = F),
  #        Week_gr = cut.Date(date, breaks = "1 week", labels = FALSE, start.on.monday = T),
  #        Week_gr_date_start = lubridate::floor_date(date, "week", week_start = 1)) %>%
  # tidyr::complete(Age_gr, Week_gr, fill = list(Mort_deaths = 0)) %>%
  # group_by(Week_gr,Age_gr,Week_gr_date_start) %>%
  # summarise(Mort_deaths = length(date)) %>%
  # ungroup %>%
  # tidyr::complete(Week_gr, Age_gr, fill = list(Mort_deaths = 0)) %>%
  # group_by(Week_gr) %>% mutate(Week_gr_date_start = na.omit(unique(Week_gr_date_start)))

# UTH_Mortality_Total <- UTH_Mortality %>% merge(Test_weeks_match) %>%
#   group_by(Test_weeks_set, Age_gr) %>%
#   summarise(Mort_deaths = sum(Mort_deaths)) %>%
#   rename(Week_gr = Test_weeks_set)

Plot_Res <- melt(LOOCV_Ests_median, varnames = c("Week_gr", "Age_gr"), value.name = "Deaths_median") %>%
  mutate(Age_gr = as.numeric(gsub(Age_gr, pattern = "RR", replacement = "")),
         Week_gr = as.numeric(gsub(x = Week_gr, "U_5_Rate_Week", ""))) %>%
  merge(melt(LOOCV_Ests_CI_low, varnames = c("Week_gr", "Age_gr"), value.name = "Deaths_low_CI") %>%
          mutate(Age_gr = as.numeric(Age_gr),
                 Week_gr = as.numeric(gsub(x = Week_gr, "U_5_Rate_Week", "")))) %>%
  merge(melt(LOOCV_Ests_CI_high, varnames = c("Week_gr", "Age_gr"), value.name = "Deaths_high_CI") %>%
          mutate(Age_gr = as.numeric(Age_gr),
                 Week_gr = as.numeric(gsub(x = Week_gr, "U_5_Rate_Week", "")))) %>%
  merge(Age_Week_Registrations)



Age_groups.labs <- c(paste0("Age group: ",c("0-4","5-9","10-14","15-29","20-24","25-29",
                              "30-34","35-39","40-44","45-49","50-54","55-59",
                              "60-64","65-69","70-74","75-79","80+")))
names(Age_groups.labs) <- 1:17


Median_Baseline_2018_2019 <- readRDS(file = "~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Burial_registrations_by_week_2018_2021_5yr_age_groups.rds") |>
  group_by(Age_gr) %>%
  summarise(median_baseline = median(BurRegs)*4)


pK <- ggplot(data = merge(Age_Week_Registrations, Plot_Res %>% filter(Age_gr !=1)),
             aes(x = BurRegs, y = Deaths_median, col = as.factor(Age_gr))) +
  geom_errorbar(aes(ymin = Deaths_low_CI, ymax = Deaths_high_CI), linewidth =0.3) +
  geom_point(size = 0.5) +
  theme_minimal(base_size = 7, base_family = "Helvetica")+
  geom_abline(slope=1, intercept=0, linetype = 2)+
  viridis::scale_color_viridis(discrete = T, name = "Age group",
                               breaks = 1:17,
                               labels= Age_groups.labs) +
  labs(x ="Burial registrations", y = "Predictions") +
  theme(plot.title = element_text(size = 7, face = "bold"),
        legend.position = "none",
        title = element_blank(),
        axis.title = element_text(),
        plot.margin = margin(5,10,5,5, "points")) +
  geom_hline(data = Median_Baseline_2018_2019 %>% filter(Age_gr !=1), aes(yintercept = median_baseline), color = "darkgrey", size = 0.35) +
  geom_point(color = "black", size = 0.5) +
  geom_point(aes(color = as.factor(Age_gr)), size = 0.3) +
  facet_wrap(~Age_gr, scales = "free", labeller = labeller(Age_gr = Age_groups.labs))



#################################### compare median estimates with actual data:
## Get Medians from data
ggsave("Analysis/Figures/Supplementary_Figures/SFigure_02_KFCV.pdf",
       pK, #xlim(c(0,165)) + ylim(c(0,160)) +
       units = "mm", height = 180*6/7, width = 180)

tiff("Analysis/Figures/Supplementary_Figures/SFigure_02_KFCV.tiff", units = "mm", res = 300, height = 180*6/7, width = 180)
pK
dev.off()
