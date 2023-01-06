####
get_Posterior <- function(model_fit){
  # Select all sampled posteriors
  PosC1 <- model_fit$pmcmc_results$chains$chain1$results$log_posterior[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain1", rownames(model_fit$replicate_parameters), value = T))))]
  PosC2 <- model_fit$pmcmc_results$chains$chain2$results$log_posterior[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain2", rownames(model_fit$replicate_parameters), value = T))))]
  PosC3 <- model_fit$pmcmc_results$chains$chain3$results$log_posterior[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain3", rownames(model_fit$replicate_parameters), value = T))))]
  PosC4 <- model_fit$pmcmc_results$chains$chain4$results$log_posterior[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain4", rownames(model_fit$replicate_parameters), value = T))))]
  PosC5 <- model_fit$pmcmc_results$chains$chain5$results$log_posterior[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain5", rownames(model_fit$replicate_parameters), value = T))))]
  PosC6 <- model_fit$pmcmc_results$chains$chain6$results$log_posterior[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain6", rownames(model_fit$replicate_parameters), value = T))))]
  PosC7 <- model_fit$pmcmc_results$chains$chain7$results$log_posterior[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain7", rownames(model_fit$replicate_parameters), value = T))))]
  PosC8 <- model_fit$pmcmc_results$chains$chain8$results$log_posterior[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain8", rownames(model_fit$replicate_parameters), value = T))))]
  # Average out posteriors
  log(unlist(Brobdingnag::sum(exp(Brobdingnag::cbrob(Brobdingnag::as.brob(PosC1),
                                                     Brobdingnag::as.brob(PosC2),
                                                     Brobdingnag::as.brob(PosC3),
                                                     Brobdingnag::as.brob(PosC4),
                                                     Brobdingnag::as.brob(PosC5),
                                                     Brobdingnag::as.brob(PosC6),
                                                     Brobdingnag::as.brob(PosC7),
                                                     Brobdingnag::as.brob(PosC8)
  )))/nrow(model_fit$replicate_parameters)))

}

get_Likelihood <- function(model_fit){
  # Select all sampled posteriors
  LikC1 <- model_fit$pmcmc_results$chains$chain1$results$log_likelihood[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain1", rownames(model_fit$replicate_parameters), value = T))))]
  LikC2 <- model_fit$pmcmc_results$chains$chain2$results$log_likelihood[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain2", rownames(model_fit$replicate_parameters), value = T))))]
  LikC3 <- model_fit$pmcmc_results$chains$chain3$results$log_likelihood[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain3", rownames(model_fit$replicate_parameters), value = T))))]
  LikC4 <- model_fit$pmcmc_results$chains$chain4$results$log_likelihood[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain4", rownames(model_fit$replicate_parameters), value = T))))]
  LikC5 <- model_fit$pmcmc_results$chains$chain5$results$log_likelihood[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain5", rownames(model_fit$replicate_parameters), value = T))))]
  LikC6 <- model_fit$pmcmc_results$chains$chain6$results$log_likelihood[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain6", rownames(model_fit$replicate_parameters), value = T))))]
  LikC7 <- model_fit$pmcmc_results$chains$chain7$results$log_likelihood[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain7", rownames(model_fit$replicate_parameters), value = T))))]
  LikC8 <- model_fit$pmcmc_results$chains$chain8$results$log_likelihood[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain8", rownames(model_fit$replicate_parameters), value = T))))]
  # Average out posteriors
  log(unlist(Brobdingnag::sum(exp(Brobdingnag::cbrob(Brobdingnag::as.brob(LikC1),
                                                     Brobdingnag::as.brob(LikC2),
                                                     Brobdingnag::as.brob(LikC3),
                                                     Brobdingnag::as.brob(LikC4),
                                                     Brobdingnag::as.brob(LikC5),
                                                     Brobdingnag::as.brob(LikC6),
                                                     Brobdingnag::as.brob(LikC7),
                                                     Brobdingnag::as.brob(LikC8)
  )))/nrow(model_fit$replicate_parameters)))
}

get_Prior <- function(model_fit){
  # Select all sampled posteriors
  PriC1 <- model_fit$pmcmc_results$chains$chain1$results$log_prior[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain1", rownames(model_fit$replicate_parameters), value = T))))]
  PriC2 <- model_fit$pmcmc_results$chains$chain2$results$log_prior[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain2", rownames(model_fit$replicate_parameters), value = T))))]
  PriC3 <- model_fit$pmcmc_results$chains$chain3$results$log_prior[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain3", rownames(model_fit$replicate_parameters), value = T))))]
  PriC4 <- model_fit$pmcmc_results$chains$chain4$results$log_prior[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain4", rownames(model_fit$replicate_parameters), value = T))))]
  PriC5 <- model_fit$pmcmc_results$chains$chain5$results$log_prior[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain5", rownames(model_fit$replicate_parameters), value = T))))]
  PriC6 <- model_fit$pmcmc_results$chains$chain6$results$log_prior[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain6", rownames(model_fit$replicate_parameters), value = T))))]
  PriC7 <- model_fit$pmcmc_results$chains$chain7$results$log_prior[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain7", rownames(model_fit$replicate_parameters), value = T))))]
  PriC8 <- model_fit$pmcmc_results$chains$chain8$results$log_prior[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain8", rownames(model_fit$replicate_parameters), value = T))))]
  # Average out posteriors
  log(mean(exp(as.numeric(unlist(list(PriC1,PriC2,PriC3,PriC4,PriC5,PriC6,PriC7,PriC8))))))
}
