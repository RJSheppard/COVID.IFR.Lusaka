####

# get_Posterior <- function(model_fit){
#   # browser()
#   # Select all sampled posteriors
#   PosC1 <- model_fit$pmcmc_results$chains$chain1$results$log_posterior[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain1", rownames(model_fit$replicate_parameters), value = T))))]
#   PosC2 <- model_fit$pmcmc_results$chains$chain2$results$log_posterior[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain2", rownames(model_fit$replicate_parameters), value = T))))]
#   PosC3 <- model_fit$pmcmc_results$chains$chain3$results$log_posterior[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain3", rownames(model_fit$replicate_parameters), value = T))))]
#   # Average out posteriors
#   log(mean(exp(as.numeric(unlist(list(PosC1,PosC2,PosC3))))))
# }

get_Posterior <- function(model_fit){
  # Select all sampled posteriors
  PosC1 <- model_fit$pmcmc_results$chains$chain1$results$log_posterior[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain1", rownames(model_fit$replicate_parameters), value = T))))]
  PosC2 <- model_fit$pmcmc_results$chains$chain2$results$log_posterior[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain2", rownames(model_fit$replicate_parameters), value = T))))]
  PosC3 <- model_fit$pmcmc_results$chains$chain3$results$log_posterior[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain3", rownames(model_fit$replicate_parameters), value = T))))]
  PosC4 <- model_fit$pmcmc_results$chains$chain4$results$log_posterior[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain4", rownames(model_fit$replicate_parameters), value = T))))]
  # Average out posteriors
  # browser()
  log(unlist(Brobdingnag::sum(exp(Brobdingnag::cbrob(Brobdingnag::as.brob(PosC1),
                                                     Brobdingnag::as.brob(PosC2),
                                                     Brobdingnag::as.brob(PosC3),
                                                     Brobdingnag::as.brob(PosC4))))/nrow(model_fit$replicate_parameters)))
  # log(unlist(Brobdingnag::sum(exp(Brobdingnag::cbrob(Brobdingnag::as.brob(PosC1),Brobdingnag::as.brob(PosC2),Brobdingnag::as.brob(PosC3))))/nrow(model_fit$replicate_parameters)))

}

get_Likelihood <- function(model_fit){
  # Select all sampled posteriors
  LikC1 <- model_fit$pmcmc_results$chains$chain1$results$log_likelihood[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain1", rownames(model_fit$replicate_parameters), value = T))))]
  LikC2 <- model_fit$pmcmc_results$chains$chain2$results$log_likelihood[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain2", rownames(model_fit$replicate_parameters), value = T))))]
  LikC3 <- model_fit$pmcmc_results$chains$chain3$results$log_likelihood[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain3", rownames(model_fit$replicate_parameters), value = T))))]
  LikC4 <- model_fit$pmcmc_results$chains$chain4$results$log_likelihood[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain4", rownames(model_fit$replicate_parameters), value = T))))]
  # Average out posteriors
  # browser()
  log(unlist(Brobdingnag::sum(exp(Brobdingnag::cbrob(Brobdingnag::as.brob(LikC1),
                                                     Brobdingnag::as.brob(LikC2),
                                                     Brobdingnag::as.brob(LikC3),
                                                     Brobdingnag::as.brob(LikC4))))/nrow(model_fit$replicate_parameters)))
  # log(unlist(Brobdingnag::sum(exp(Brobdingnag::cbrob(Brobdingnag::as.brob(LikC3),Brobdingnag::as.brob(LikC1),Brobdingnag::as.brob(LikC2))))/nrow(model_fit$replicate_parameters)))

}

3.602107+138.0567+742.3249
get_Prior <- function(model_fit){
  # browser()
  # Select all sampled posteriors
  PriC1 <- model_fit$pmcmc_results$chains$chain1$results$log_prior[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain1", rownames(model_fit$replicate_parameters), value = T))))]
  PriC2 <- model_fit$pmcmc_results$chains$chain2$results$log_prior[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain2", rownames(model_fit$replicate_parameters), value = T))))]
  PriC3 <- model_fit$pmcmc_results$chains$chain3$results$log_prior[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain3", rownames(model_fit$replicate_parameters), value = T))))]
  PriC4 <- model_fit$pmcmc_results$chains$chain4$results$log_prior[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain4", rownames(model_fit$replicate_parameters), value = T))))]
  # Average out posteriors
  log(mean(exp(as.numeric(unlist(list(PriC1,PriC2,PriC3,PriC4))))))
}

# log(unlist(Brobdingnag::sum(exp(Brobdingnag::cbrob(Brobdingnag::as.brob(LikC1+PriC1),Brobdingnag::as.brob(LikC2+PriC2),Brobdingnag::as.brob(LikC3+PriC3))))/nrow(model_fit$replicate_parameters))) ==
#   log(unlist(Brobdingnag::sum(exp(Brobdingnag::cbrob(Brobdingnag::as.brob(PosC1),Brobdingnag::as.brob(PosC2),Brobdingnag::as.brob(PosC3))))/nrow(model_fit$replicate_parameters)))

#

# unlist(Brobdingnag::sum(exp(Brobdingnag::cbrob(Brobdingnag::as.brob(LikC1),Brobdingnag::as.brob(LikC2),Brobdingnag::as.brob(LikC3))))/nrow(model_fit$replicate_parameters)) *
#   mean(exp(as.numeric(unlist(list(PriC1,PriC2,PriC3))))) ==
#   unlist(Brobdingnag::sum(exp(Brobdingnag::cbrob(Brobdingnag::as.brob(PosC1),Brobdingnag::as.brob(PosC2),Brobdingnag::as.brob(PosC3))))/nrow(model_fit$replicate_parameters))



# log(mean(exp(-c(2,4,6,8))) * mean(exp(-c(1,3,5,7))))
# log(mean(exp(-c(2,4,6,8)))) + log(mean(exp(-c(1,3,5,7))))
#
# log(mean(exp(-c(2,4,6,8) -c(1,3,5,7))))
#
#
# mean(exp(-c(2,4,6,8))) * mean(exp(-c(1,3,5,7)))
#
# mean(exp(-c(2,4,6,8) -c(1,3,5,7)))

# exp(-90)
# exp(-45) * exp(-2)
#
# ## log(2*5) = log(2) + log(5)
# log(0.5) + log(0.2)
# log(0.1)# -90
# exp(-90/2)
#
# exp(-9)
# -9
#
# exp(-4.5) * exp(-2)
#
# PosC1 <- model_fit$pmcmc_results$chains$chain1$results$log_prior[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain1", rownames(model_fit$replicate_parameters), value = T))))]
# PosC2 <- model_fit$pmcmc_results$chains$chain2$results$log_likelihood[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain2", rownames(model_fit$replicate_parameters), value = T))))]
# PosC3 <- model_fit$pmcmc_results$chains$chain3$results$log_posterior[floor(as.numeric(gsub("^\\D*\\d.",replacement = "", grep("chain3", rownames(model_fit$replicate_parameters), value = T))))]
