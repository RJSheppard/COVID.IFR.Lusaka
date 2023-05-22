Get_baselines <- function(mcmc_samples, Week_names, RR_names){
  # Week_index <- grepl(pattern = Week_names, x = names(mcmc_samples))
  Pre_2020_index <- names(mcmc_samples) %in% paste0(Week_names, 1:104)
  Week_index <- names(mcmc_samples) %in% paste0(Week_names, 129:144)
  RR_index <- grepl(pattern = RR_names, x = names(mcmc_samples))

  ## Pre_2020
  Pre_2020_average_age_u5 <- rowMeans(mcmc_samples[,Pre_2020_index])

  mcmc_baselines <- lapply(1:nrow(mcmc_samples), function(x){
    mcmc_baselines_tmp <- cbind(as.numeric(mcmc_samples[x,Week_index]),
                                as.numeric(mcmc_samples[x,Week_index]) %*% t(as.numeric(mcmc_samples[x,RR_index])))
    rownames(mcmc_baselines_tmp) <- 1:16
    colnames(mcmc_baselines_tmp) <- 1:17
    mcmc_baselines_tmp <- mcmc_baselines_tmp %>% melt(value.name = "Mort_ncd_mcmc", varnames = c("Week_gr", "Age_gr"))

    Pre_2020 <- data.frame(Age_gr = 1:17,
                           Bg_dr = unlist(c(Pre_2020_average_age_u5[x], Pre_2020_average_age_u5[x] * mcmc_samples[x,RR_index])))
    merge(mcmc_baselines_tmp,Pre_2020) %>% mutate(ag1std = Mort_ncd_mcmc/Bg_dr)
  })
}

turn_list_into_array <- function(arr){
  array(data = unlist(arr),
        dim = c(dim(arr[[1]]),
                length(arr)),
        dimnames = list(NULL,colnames(arr[[1]]),NULL))
}
