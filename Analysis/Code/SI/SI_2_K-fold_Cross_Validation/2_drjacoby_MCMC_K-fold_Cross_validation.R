dj_loocv <- function(Test_weeks, df_params, r_loglike, r_logprior, data_list, age_cats,
                     burnin, samples, chains){
  drjacoby::run_mcmc(data = data_list,
                     df_params = df_params,
                     loglike = r_loglike,
                     logprior = r_logprior,
                     burnin = burnin,
                     samples = samples,
                     pb_markdown = TRUE,
                     chains = chains,
                     misc = list(age_cats = age_cats,
                                 Test_weeks = Test_weeks))}
