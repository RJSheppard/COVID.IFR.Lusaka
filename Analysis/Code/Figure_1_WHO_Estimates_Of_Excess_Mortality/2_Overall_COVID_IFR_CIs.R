## .................................................................................
## Purpose: Overall IFR from the Brazeau 2020 from the log-normal model that was fit to the
## posterior country-age-level estimates
##
## Author: Nick Brazeau
##
## Date: 20 October, 2022
##
## Notes: These are prediction intervals
## .................................................................................
set.seed(48)
library(tidyverse)

#............................................................
# Import/Tidy
#...........................................................
# ifr dat from lognormal model in scrip `05-analyze_age_IFRs_v2.R`
ifrdat <- readRDS("for_richard/data/ifrdat.RDS")

ifrdat %>% filter(sero == "reg") %>% pull(bestestmod)
ifrdat %>% filter(sero == "reg") %>% pull(linear_predints)

# read in Richard's countries
cntry <- readr::read_tsv("for_richard/data/countries.txt", col_names = F)
cntry <- unlist(cntry)

# tidy up UN World Population Prospects
wpp <- readxl::read_excel("for_richard/data/wpp_un_agepopulations.xlsx")
wpp <- wpp[wpp$`Reference date (as of 1 July)` == 2020,]
wpp <- wpp %>%
  dplyr::select(c("Region, subregion, country or area *",
                  "0-4", "5-9",
                  "10-14", "15-19", "20-24", "25-29",
                  "30-34", "35-39", "40-44", "45-49",
                  "50-54", "55-59", "60-64", "65-69",
                  "70-74", "75-79", "80-84", "85-89",
                  "90-94", "95-99", "100+"), Type)
colnames(wpp)[1] <- "georegion"

wpp <- wpp %>%
  dplyr::filter(georegion %in% cntry) %>%
  tidyr::pivot_longer(., cols = -c("georegion"), names_to = "ageband", values_to = "popN") %>%
  dplyr::mutate(popN = as.numeric(popN)*1e3) # wpp adjustment


# cuts
wpp <- wpp %>%
  dplyr::mutate(
    ageband = ifelse(ageband == "100+", "95-99", ageband),
    age_high = as.numeric(stringr::str_split_fixed(ageband, "-", n=2)[,2]),
    ageband = cut(age_high,
                  breaks = c(0, seq(4, 89, by = 5), 999)),
    ageband = as.character(ageband)) %>%
  dplyr::group_by(georegion, ageband) %>%
  dplyr::summarise(pop_size = sum(popN)) %>%
  dplyr::mutate(age_high = as.numeric(stringr::str_extract(ageband, "[0-9]+?(?=])"))) %>%
  dplyr::arrange(age_high) %>%
  dplyr::select(-c("age_high")) %>%
  dplyr::group_by(georegion) %>%
  tidyr::nest(.) %>%
  dplyr::rename(popN = data)

Complete_Population_demog <- readRDS(file = "Complete_Population_demog.rds") %>% rename(georegion = area,
                                                                                        ageband = age_cat,
                                                                                        popN = n) %>%
  mutate(ageband = as.numeric(ageband),
         ageband_Check = ageband,
         ageband = cut(ageband*5-2,
                       breaks = c(0, seq(4, 89, by = 5), 999)),
         ageband = as.character(ageband)) %>%
  dplyr::group_by(georegion, ageband) %>%
  dplyr::summarise(pop_size = sum(popN)) %>%
  dplyr::group_by(georegion) %>%
  tidyr::nest(.) %>%
  dplyr::rename(popN = data)

#........................
# Monte Carlo Calculations
#........................
calc_monte_carlo_overall_IFRs <- function(new_dat, popN, reps = 1e5) {
  if (!all(c("ageband", "pop_size") %in% colnames(popN))) {
    stop()
  }
  # bring in demog
  new_dat <- dplyr::left_join(new_dat, popN)
  # calculate
  z <- mapply(function(i) {
    sum(new_dat$pop_size * rlnorm(nrow(new_dat), meanlog = new_dat$fit, sdlog = sqrt(new_dat$var_fit)))
  }, seq_len(reps))
  z <- z / sum(new_dat$pop_size)

  # get prediction intervals
  Q025 <- round(quantile(z, 0.025) * 100, 3)
  Q50 <- round(quantile(z, 0.50) * 100, 3)
  Q975 <- round(quantile(z, 0.975) * 100, 3)
  # get mean
  mean <- exp(new_dat$fit + new_dat$var_fit/2)
  mean <- sum( mean * (new_dat$pop_size/sum(new_dat$pop_size)) )
  mean <- round(mean * 100, 3)

  # out
  ret <- tibble::tibble(Q025 = Q025,
                        Q50 = Q50,
                        mean = mean,
                        Q975 = Q975)
  return(ret)
}


# tidy up pieces
overall_IFR_best_est <- tidyr::expand_grid(ifrdat, Complete_Population_demog) %>%
  dplyr::select(c("sero", "linear_predints", "georegion", "popN")) %>%
  dplyr::rename(new_dat = linear_predints)
overall_IFR_best_est$bestest <- purrr::pmap(overall_IFR_best_est[, c("new_dat", "popN")],
                                            calc_monte_carlo_overall_IFRs,
                                            reps = 1e5)

# out
out <- overall_IFR_best_est %>%
  dplyr::select(c("sero", "georegion", "bestest")) %>%
  tidyr::unnest(cols = "bestest")

############################################
############################################
saveRDS(out, "Brazeau_IFR_Post_predictive_CrIs")
############################################
############################################
