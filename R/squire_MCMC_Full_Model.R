cma_fit <- function(probs_hosp_death,
                    data,

                    combined_data = NULL,
                    drj_mcmc = NULL,
                    lld,

                    pcr_df,
                    sero_df,

                    population,
                    baseline_contact_matrix,
                    dur_get_ox_survive,
                    dur_get_ox_die,

                    n_mcmc,
                    replicates,

                    log_likelihood=NULL,
                    Prior_Rt_rw_unif_lim = NULL,

                    pcr_det = NULL,
                    frac_reg
){
  Test <- cma::fit_spline_rt(prob_non_severe_death_treatment = probs_hosp_death$phi_2,

                             population = population,
                             baseline_contact_matrix = baseline_contact_matrix,
                             dur_get_ox_survive = dur_get_ox_survive,
                             dur_get_ox_die = dur_get_ox_die,
                             country = "Zambia",
                             reporting_fraction = 1,
                             n_mcmc = n_mcmc,
                             replicates = replicates,
                             rw_duration = 14,
                             hosp_beds = 1e10,
                             icu_beds = 1e10,
                             prob_severe = rep(0,17),
                             prob_hosp = probs_hosp_death$phi_1,
                             dur_R = Inf,

                             data = data,

                             combined_data = combined_data,
                             drj_mcmc = drj_mcmc,

                             log_likelihood = log_likelihood,

                             pcr_df = pcr_df,
                             sero_df = sero_df,


                             lld = lld,
                             Prior_Rt_rw_unif_lim = Prior_Rt_rw_unif_lim,

                             pcr_det = pcr_det,
                             frac_reg = frac_reg)
}
