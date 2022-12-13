#' @noRd
# save figure
save_figs <- function(name,
                      fig,
                      width = 6,
                      height = 6,
                      root = file.path(here::here(), "analysis/plots")) {

  dir.create(root, showWarnings = FALSE)
  fig_path <- function(name) {paste0(root, "/", name)}

  cowplot::save_plot(filename = fig_path(paste0(name,".png")),
                     plot = fig,
                     base_height = height,
                     base_width = width)

  pdf(file = fig_path(paste0(name,".pdf")), width = width, height = height)
  print(fig)
  dev.off()

}

mortality_plot <- function(data, fit){
  plot(fit, particle_fit = TRUE) +
    geom_bar(data = data, stat = "identity", aes(x=date, y=deaths), alpha = 0.6) +
    annotate(geom = "segment", x=as.Date("2020-02-28"), xend = as.Date("2020-03-05"), y = max(fit$pmcmc_results$inputs$data$deaths)/mean(fit$replicate_parameters$rf), yend = max(fit$pmcmc_results$inputs$data$deaths)/mean(fit$replicate_parameters$rf), alpha = 0.5, size = 1.5) +
    annotate(geom = "text", x=as.Date("2020-03-07"), y = max(fit$pmcmc_results$inputs$data$deaths)/mean(fit$replicate_parameters$rf), hjust = 0, label = "Gov. reported deaths") +
    annotate(geom = "point", x= as.Date("2020-03-02"), y = max(fit$pmcmc_results$inputs$data$deaths)*0.95/mean(fit$replicate_parameters$rf)) +
    annotate(geom = "text", x=as.Date("2020-03-07"), y = max(fit$pmcmc_results$inputs$data$deaths)*0.95/mean(fit$replicate_parameters$rf), hjust=0, label = "Gov. reported deaths (scaled)") +
    annotate(geom = "segment", x=as.Date("2020-02-28"), xend = as.Date("2020-03-05"), y = max(fit$pmcmc_results$inputs$data$deaths)*0.9/mean(fit$replicate_parameters$rf), yend = max(fit$pmcmc_results$inputs$data$deaths)*0.9/mean(fit$replicate_parameters$rf), col = "red", alpha = 0.5, size = 1.5) +
    annotate(geom = "rect", xmin=as.Date("2020-02-28"), xmax = as.Date("2020-03-05"), ymin = max(fit$pmcmc_results$inputs$data$deaths)*0.89/mean(fit$replicate_parameters$rf), ymax = max(fit$pmcmc_results$inputs$data$deaths)*0.91/mean(fit$replicate_parameters$rf), fill = NA, color = "black", linetype = 2) +
    annotate(geom = "text", x=as.Date("2020-03-07"), y = max(fit$pmcmc_results$inputs$data$deaths)*0.9/mean(fit$replicate_parameters$rf), hjust = 0, label = "Model fit") +
    theme_bw() + theme(legend.position = "none")
}

pcr_plot <- function(PCR_Sero_PlotData, fit, Lancet_Data){
  ggplot(PCR_Sero_PlotData, aes(x = date, y = mean_pcr)) + geom_line(aes(x=date, y=mean_pcr),linetype="dashed") +
    geom_ribbon(aes(x=date,ymin=min_pcr, ymax=max_pcr), alpha=0.3)+
    annotate("text", x=as.Date("2020-10-01"), y=max(tail(PCR_Sero_PlotData$max_pcr))*2, label="Model fit", alpha=0.8) +

    geom_point(aes(x= as.Date("2020-07-15"),y=Lancet_Data$PCR_prev["val"])) +
    geom_errorbar(aes(ymin=Lancet_Data$PCR_prev["ci95l"],ymax=Lancet_Data$PCR_prev["ci95h"],x=as.Date("2020-07-15"), width=10)) +
    geom_errorbarh(aes(xmin=as.Date("2020-07-04"),xmax=as.Date("2020-07-27"),y=Lancet_Data$PCR_prev["val"], height=0)) +
    annotate("text", x=as.Date("2020-07-15")+5, y=Lancet_Data$PCR_prev["val"]+0.5, label="Mulenga et al.", hjust=0) +

    ylab(paste0("PCR+ % with ",100*round(mean(fit$replicate_parameters$rf),3),"% death reporting")) +
    xlab("Date") +
    theme_bw()
}

sero_plot <- function(PCR_Sero_PlotData, fit, Lancet_Data){
  ggplot(PCR_Sero_PlotData, aes(x = date, y = mean_sero)) + geom_line(aes(x=date, y=mean_sero),linetype="dashed") +
    geom_ribbon(aes(x=date,ymin=min_sero, ymax=max_sero), alpha=0.3)+
    annotate("text", x=as.Date("2020-10-01"), y=max(tail(PCR_Sero_PlotData$max_sero))*1.12, label="Model fit", alpha=0.8) +

    geom_point(aes(x= as.Date("2020-07-15"),y=Lancet_Data$Sero_prev["val"])) +
    geom_errorbar(aes(ymin=Lancet_Data$Sero_prev["ci95l"],ymax=Lancet_Data$Sero_prev["ci95h"],x=as.Date("2020-07-15"), width=10)) +
    geom_errorbarh(aes(xmin=as.Date("2020-07-04"),xmax=as.Date("2020-07-27"),y=Lancet_Data$Sero_prev["val"], height=0)) +
    annotate("text", x=as.Date("2020-07-15")+5, y=Lancet_Data$Sero_prev["val"]+0.5, label="Mulenga et al.", hjust =0) +

    geom_point(aes(x= as.Date("2020-07-29"),y=Lancet_Data$Tot_cov_prev["val"]), color="darkblue") +
    geom_errorbar(aes(ymin=Lancet_Data$Tot_cov_prev["ci95l"],ymax=Lancet_Data$Tot_cov_prev["ci95h"],x=as.Date("2020-07-29"), width=10), color="darkblue") +
    geom_errorbarh(aes(xmin=as.Date("2020-07-18"),xmax=as.Date("2020-08-10"),y=Lancet_Data$Tot_cov_prev["val"], height=0), color="darkblue") +
    annotate("text", x=as.Date("2020-07-29")+15, y=Lancet_Data$Tot_cov_prev["val"], label="Mulenga et al. \n (combined\n  measure)", color = "darkblue", hjust =0) +

    ylab(paste0("Sero+% with ",100*round(mean(fit$replicate_parameters$rf),3),"% death reporting")) +
    xlab("Date") +
    theme_bw()
}





######
part_fit_plot <- function(fit){plot(fit, particle_fit = T) +
    annotate(geom = "point", x= as.Date("2020-04-18"), y = 50) +
    annotate(geom = "text", x=as.Date("2020-04-22"), y = 50, hjust=0, label = "BMJ data estimated deaths") +
    annotate(geom = "segment", x=as.Date("2020-04-15"), xend = as.Date("2020-04-20"), y = 50*0.92, yend = 50*0.92, col = "red", alpha = 0.7) +
    annotate(geom = "rect", xmin=as.Date("2020-04-15"), xmax = as.Date("2020-04-20"), ymin = 50*0.915, ymax = 50*0.925, fill = NA, color = "black", linetype = 2) +
    annotate(geom = "text", x=as.Date("2020-04-22"), y = 50*0.92, hjust = 0, label = "Model fit") +
    xlim(as.Date("2020-04-15"), as.Date("2020-10-01")) + coord_cartesian(ylim=c(0, 50)) +
    theme(legend.position = "none")
}
pcr_fit_plot <- function(Summ_sero_fit){ggplot(Summ_sero_fit, aes(x = date, y = mean_pcr)) +
    geom_line(aes(x=date, y=mean_pcr, color = "Model fit"),linetype="dashed") +
    geom_ribbon(aes(x=date,ymin=min_pcr, ymax=max_pcr), alpha=0.3)+
    geom_point(aes(x= as.Date("2020-07-15"),y=7.6, color="\nData from\nMulenga et al.\n")) +
    # annotate("text", x=as.Date("2020-07-15")+5, y=7.6+1, label="Mulenga et al.", hjust =0) +
    # annotate("text", x=as.Date("2020-09-01"), y=2, label="Model fit", alpha=0.8) +
    geom_errorbar(aes(ymin=4.7,ymax=10.6,x=as.Date("2020-07-15"), width=10)) +
    geom_errorbarh(aes(xmin=as.Date("2020-07-04"),xmax=as.Date("2020-07-27"),y=7.6, height=0)) +
    ylab(paste0("PCR +ve %")) +
    xlim(as.Date("2020-04-15"), as.Date("2020-10-01")) +# coord_cartesian(ylim=c(0, 15)) +
    xlab("Date")+
    theme_bw()} #+
    # ggtitle(paste("bin ll", round(mean(dbinom(x = round(0.076*2990), size = 2990, prob = filter(pcr_sero_data, date %in% as.Date(c("2020-07-15")))$mean_pcr/100,log = T)),1)))}
sero_fit_plot <- function(Summ_sero_fit){ggplot(Summ_sero_fit, aes(x = date, y = mean_sero)) +
    geom_line(aes(x=date, y=mean_sero),linetype="dashed") +
    geom_ribbon(aes(x=date,ymin=min_sero, ymax=max_sero), alpha=0.3)+
    geom_point(aes(x= as.Date("2020-07-15"),y=2.1)) +
    # annotate("text", x=as.Date("2020-07-15")+5, y=2.1+1, label="Mulenga et al.", hjust =0) +
    # annotate("text", x=as.Date("2020-09-01"), y=max(Summ_sero_fit$max_sero), label="Model fit", alpha=0.8, hjust=0, vjust = 2) +
    geom_errorbar(aes(ymin=1.1,ymax=3.1,x=as.Date("2020-07-15"), width=10)) +
    geom_errorbarh(aes(xmin=as.Date("2020-07-04"),xmax=as.Date("2020-07-27"),y=2.1, height=0)) +
    # geom_point(aes(x= as.Date("2020-07-29"),y=10.6), color="darkblue") +
    # annotate("text", x=as.Date("2020-07-29")+15, y=10.6, label="Mulenga et al. \n (extrapolated \n from PCR)", color = "darkblue", hjust =0) +
    # geom_text(aes(label = "Mulenga et al. \n (extrapolated \n from PCR)", x= as.Date("2020-07-29")+15,y=0.106), color = "darkblue", hjust= 0) +
    # geom_errorbar(aes(ymin=7.3,ymax=13.9,x=as.Date("2020-07-29"), width=10), color="darkblue") +
    # geom_errorbarh(aes(xmin=as.Date("2020-07-18"),xmax=as.Date("2020-08-10"),y=10.6, height=0), color="darkblue") +
    ylab(paste0("Sero +ve %")) +
    annotate("text", x=as.Date("2020-05-01"), y=1.8, label=paste(""), color = "darkblue", hjust =0) +
    xlim(as.Date("2020-04-15"), as.Date("2020-10-01")) + #coord_cartesian(ylim=c(0, 20)) +
    xlab("Date")+
    theme_bw()}

combined_fit_plot <- function(Summ_sero_fit){ggplot(Summ_sero_fit, aes(x = date, y = mean_combined)) + geom_line(aes(x=date, y=mean_combined),linetype="dashed") +
    geom_ribbon(aes(x=date,ymin=min_combined, ymax=max_combined), alpha=0.3)+
    geom_point(aes(x= as.Date("2020-07-11"),y=9.1)) +
    # annotate("text", x=as.Date("2020-07-11")+5, y=9.1+1, label="Mulenga et al.", hjust =0) +
    # annotate("text", x=as.Date("2020-09-01"), y=max(Summ_sero_fit$max_combined), label="Model fit", alpha=0.8, hjust=0, vjust = 2) +
    geom_errorbar(aes(ymin=2.6,ymax=15.7,x=as.Date("2020-07-11"), width=10)) +
    geom_errorbarh(aes(xmin=as.Date("2020-07-04"),xmax=as.Date("2020-07-19"),y=9.1, height=0)) +
    # geom_point(aes(x= as.Date("2020-07-29"),y=10.6), color="darkblue") +
    # annotate("text", x=as.Date("2020-07-29")+15, y=10.6, label="Mulenga et al. \n (extrapolated \n from PCR)", color = "darkblue", hjust =0) +
    # geom_text(aes(label = "Mulenga et al. \n (extrapolated \n from PCR)", x= as.Date("2020-07-29")+15,y=0.106), color = "darkblue", hjust= 0) +
    # geom_errorbar(aes(ymin=7.3,ymax=13.9,x=as.Date("2020-07-29"), width=10), color="darkblue") +
    # geom_errorbarh(aes(xmin=as.Date("2020-07-18"),xmax=as.Date("2020-08-10"),y=10.6, height=0), color="darkblue") +
    ylab(paste0("Combined +ve %")) +
    # annotate("text", x=as.Date("2020-05-01"), y=1.8, label=paste(""), color = "darkblue", hjust =0) +
    xlim(as.Date("2020-04-15"), as.Date("2020-10-01")) + #coord_cartesian(ylim=c(0, 20)) +
    xlab("Date")+
    theme_bw()}


Plot_Heatmaps <- function(Mod_Res, Res_Figs, Select_Runs, Title){
# browser()
  IFR_mat <- readRDS("analysis/data/Code-generated-data/00_03_IFR_matrix_coefficients_log_scale_new_pop_str_ests.rds")
  # IFR_vec <- 1:nrow(IFR_mat) %in% readRDS("analysis/data/Code-generated-data/00_03_Prob_Index_Vector_log_sc.rds")
  # IFR_mat_fil <- IFR_mat[IFR_vec,]
  # rownames(IFR_mat_fil) <- 1:nrow(IFR_mat_fil)
  # IFR_mat_fil <- IFR_mat_fil %>% select(IFR_x, Slope_x)

  # IFR_vec <- 1:nrow(IFR_mat) %in% Select_Runs  &  1:nrow(IFR_mat) %in% as.numeric(rownames(readRDS("~/Documents/Imperial/PostDoc/Zambia/covid-mortality-ascertainment/analysis/data/Code-generated-data/00_03_IFR_matrix_coefficients_log_scale.rds")[readRDS("~/Documents/Imperial/PostDoc/Zambia/covid-mortality-ascertainment/analysis/data/Code-generated-data/00_03_Prob_Index_Vector_log_sc.rds"),]))
  # IFR_vec <- 1:81 %in% Select_Runs
  Res_Figs <- Res_Figs[!unlist(lapply(Mod_Res, is.null))]
  Mod_Res <- Mod_Res[!unlist(lapply(Mod_Res, is.null))]
  IFR_vec <- 1:81 %in% as.numeric(gsub(pattern = "X", "",names(Mod_Res)))


  # readRDS("analysis/data/Code-generated-data/00_03_Prob_death_logical_log_sc.rds")


  # Res_Figs[[1]]$ll_tot
# browser()
  # log(Brobdingnag::sum(exp(Brobdingnag::as.brob(Res_Figs[[1]]$ll_all$ll_total)))/length(Res_Figs[[1]]$ll_all$ll_total))
  # get_Likelihood(Mod_Res[[1]])
  # Res_Figs[[1]]$ll_tot

  # IFR_mat$AvPost <- NA

  # sort(unlist(lapply(Res_Figs, function(x){
    # log(Brobdingnag::sum(exp(Brobdingnag::as.brob(x$ll_all$ll_total)))/length(x$ll_all$ll_total))})))

  # log(Brobdingnag::sum(exp(Brobdingnag::as.brob(Res_Figs[[1]]$ll_all$ll_total)))/length(Res_Figs[[1]]$ll_all$ll_total))
  # get_Likelihood(Mod_Res[[1]])
  # sort(unlist(lapply(Mod_Res, get_Likelihood)))

IFR_mat$AvPost[IFR_vec] <- unlist(lapply(Mod_Res, get_Posterior))
IFR_mat$AvLike[IFR_vec] <- unlist(lapply(Mod_Res, get_Likelihood))
IFR_mat$AvPrior[IFR_vec] <- unlist(lapply(Mod_Res, get_Prior))

IFR_mat$ll_binom[IFR_vec] <- unlist(lapply(Res_Figs, function(x){x$ll_bin}))
IFR_mat$ll_pois[IFR_vec] <- unlist(lapply(Res_Figs, function(x){x$ll_pois}))
IFR_mat$ll_pcr[IFR_vec] <- unlist(lapply(Res_Figs, function(x){x$ll_pcr}))
IFR_mat$ll_sero[IFR_vec] <- unlist(lapply(Res_Figs, function(x){x$ll_sero}))
IFR_mat$ll_tot[IFR_vec] <- unlist(lapply(Res_Figs, function(x){x$ll_tot}))
IFR_mat$ll_pcr_sero_sum <- IFR_mat$ll_pcr + IFR_mat$ll_sero

# IFR_mat$AvLike == IFR_mat$ll_tot
# cbind(IFR_mat$AvLike, IFR_mat$ll_tot)
# cbind(IFR_mat$AvPost, IFR_mat$ll_tot)

  # IFR_mat$AvPost[IFR_vec] <- unlist(lapply(Mod_Res[IFR_vec], get_Posterior))
  # IFR_mat$AvLike[IFR_vec] <- unlist(lapply(Mod_Res[IFR_vec], get_Likelihood))
  # IFR_mat$AvPrior[IFR_vec] <- unlist(lapply(Mod_Res[IFR_vec], get_Prior))

  # IFR_mat$ll_bin[IFR_vec] <- unlist(lapply(Res_Figs, function(x){x$lls$ll_bin}))
  # IFR_mat$ll_pois[IFR_vec] <- unlist(lapply(Res_Figs, function(x){x$lls$ll_pois}))
  # IFR_mat$ll_pcr[IFR_vec] <- unlist(lapply(Res_Figs, function(x){x$lls$ll_pcr}))
  # IFR_mat$ll_sero[IFR_vec] <- unlist(lapply(Res_Figs, function(x){x$lls$ll_sero}))
  # IFR_mat$ll_comb[IFR_vec] <- unlist(lapply(Res_Figs, function(x){x$lls$ll_comb}))

  # IFR_mat$ll_binom[IFR_vec] <- unlist(lapply(Res_Figs[IFR_vec], function(x){x$ll_bin}))
  # IFR_mat$ll_pois[IFR_vec] <- unlist(lapply(Res_Figs[IFR_vec], function(x){x$ll_pois}))
  # IFR_mat$ll_pcr[IFR_vec] <- unlist(lapply(Res_Figs[IFR_vec], function(x){x$ll_pcr}))
  # IFR_mat$ll_sero[IFR_vec] <- unlist(lapply(Res_Figs[IFR_vec], function(x){x$ll_sero}))
  # IFR_mat$ll_pcr_sero_sum[IFR_vec] <- IFR_mat$ll_pcr[IFR_vec] + IFR_mat$ll_sero[IFR_vec]

  # IFR_mat$ll_comb[IFR_vec] <- unlist(lapply(Res_Figs, function(x){x$lls$ll_comb}))


  # IFR_mat <- IFR_mat %>% mutate(AvPost = ll_binom + ll_pois + ll_comb + AvPrior)

  IFR_mat <- IFR_mat %>% mutate(
    Prior_col_group = as.numeric(cut(round(AvPrior),
                                    breaks = c(-Inf, round(max(AvPrior,na.rm = T))- c(500,200,100,50,20,12,8,4,0)),
                                    labels = as.character(c(1:9)))),
    Like_col_group = as.numeric(cut(round(AvLike),
                                     breaks = c(-Inf, round(max(AvLike,na.rm = T))- c(500,200,100,50,20,12,8,4,0)),
                                     labels = as.character(c(1:9)))),
    Post_col_group = as.numeric(cut(round(AvPost),
                                                                breaks = c(-Inf, round(max(AvPost,na.rm = T))- c(500,200,100,50,20,12,8,4,0)),
                                                                labels = as.character(c(1:9)))),
                                Bin_col_group = as.numeric(cut(round(ll_binom),
                                                               breaks = c(-Inf, round(max(ll_binom,na.rm = T))- c(500,200,100,50,20,12,8,4,0)),
                                                               labels = as.character(c(1:9)))),
                                Pois_col_group = as.numeric(cut(round(ll_pois),
                                                                breaks = c(-Inf, round(max(ll_pois,na.rm = T))- c(500,200,100,50,20,12,8,4,0)),
                                                                labels = as.character(c(1:9)))),
                                pcr_col_group = as.numeric(cut(round(ll_pcr),
                                                               breaks = c(-Inf, round(max(ll_pcr,na.rm = T))- c(500,200,100,50,20,12,8,4,0)),
                                                               labels = as.character(c(1:9)))),
                                sero_col_group = as.numeric(cut(round(ll_sero),
                                                                breaks = c(-Inf, round(max(ll_sero,na.rm = T))- c(500,200,100,50,20,12,8,4,0)),
                                                                labels = as.character(c(1:9)))),
    pcr_sero_sum_col_group = as.numeric(cut(round(ll_pcr_sero_sum),
                                    breaks = c(-Inf, round(max(ll_pcr_sero_sum,na.rm = T))- c(500,200,100,50,20,12,8,4,0)),
                                    labels = as.character(c(1:9)))))
                                # comb_col_group = as.numeric(cut(round(ll_comb),
                                                                # breaks = c(-Inf, round(max(ll_comb,na.rm = T))- c(500,200,100,50,20,12,8,4,0)),
                                                                # labels = as.character(c(1:9)))))

  # p0_prior <- ggplot(IFR_mat, aes(x = as.factor(round(IFR_x,2)), y = as.factor(round(Slope_x,2)), fill = Prior_col_group)) + geom_tile() +
  #   geom_text(aes(label = round(AvPrior), colour = (Prior_col_group >= max(Prior_col_group, na.rm=T))), size = 4) +
  #   scale_colour_manual(values = c("white", "black")) +
  #   ggtitle("Prior") + xlab("Overall severity") + ylab("IFR age disparity") +
  #   labs(fill = "Mean Prior") + theme(legend.position = "none") +
  #   # xlim(0,2) + ylim(0,2) +
  #   scale_fill_viridis_c()

  p1_like <- ggplot(IFR_mat %>% filter(Slope_x != 5), aes(x = as.factor(round(IFR_x,2)), y = as.factor(round(Slope_x,2)), fill = as.factor(Like_col_group))) + geom_tile() +
    geom_text(aes(label = round(AvLike), colour = (Like_col_group >= max(Like_col_group, na.rm=T))), size = 3) +
    scale_colour_manual(values = c("white", "black")) +
    ggtitle("Likelihood") + xlab("Overall severity") + ylab("IFR age gradient") +
    labs(fill = "Mean Likelihood") + theme(legend.position = "none") +
    # xlim(0,2) + ylim(0,2) +
    scale_fill_discrete(type = tail(viridis::viridis(n = 9), length(table(IFR_mat$Like_col_group))))
    # scale_fill_viridis_c()



   p1_post <- ggplot(IFR_mat %>% filter(Slope_x != 5), aes(x = as.factor(round(IFR_x,2)), y = as.factor(round(Slope_x,2)), fill = as.factor(Post_col_group))) + geom_tile() +
    geom_text(aes(label = round(AvPost), colour = (Post_col_group >= max(Post_col_group, na.rm=T)))) +
    scale_colour_manual(values = c("white", "black")) +
    ggtitle("Overall posterior") + xlab("Overall severity") + ylab("IFR age gradient") +
    labs(fill = "Mean Post") + theme(legend.position = "none") +
    # xlim(0,2) + ylim(0,2) +
     scale_fill_discrete(type = tail(viridis::viridis(n = 9), length(table(IFR_mat$Post_col_group))))
    # scale_fill_viridis_c()


  p2_bin <- ggplot(IFR_mat %>% filter(Slope_x != 5), aes(x = as.factor(round(IFR_x,2)), y = as.factor(round(Slope_x,2)), fill = as.factor(Bin_col_group))) + geom_tile() +
    geom_text(aes(label = round(ll_binom), colour = (Bin_col_group >= max(Bin_col_group, na.rm=T))), size = 1.5) +
    scale_colour_manual(values = c("white", "black")) +
    ggtitle("C19 prevalence in post-mortem samples") + xlab("Overall severity") + ylab("IFR age gradient") +
    labs(fill = "Mean ll bin") + theme(legend.position = "none") +
    # xlim(0,2) + ylim(0,2) +
    scale_fill_discrete(type = tail(viridis::viridis(n = 9), length(table(IFR_mat$Bin_col_group))))
    # scale_fill_viridis_c()

  p3_pois <- ggplot(IFR_mat %>% filter(Slope_x != 5), aes(x = as.factor(round(IFR_x,2)), y = as.factor(round(Slope_x,2)), fill = as.factor(Pois_col_group))) + geom_tile() +
    geom_text(aes(label = round(ll_pois), colour = (Pois_col_group >= max(Pois_col_group, na.rm=T))), size = 1.5) +
    scale_colour_manual(values = c("white", "black")) +
    ggtitle("Burial registration rate") + xlab("Overall severity") + ylab("IFR age gradient") +
    labs(fill = "Mean ll pois") + theme(legend.position = "none") +
    # xlim(0,2) + ylim(0,2) +
    scale_fill_discrete(type = tail(viridis::viridis(n = 9), length(table(IFR_mat$Pois_col_group))))
    # scale_fill_viridis_c()

  p4_pcr <- ggplot(IFR_mat %>% filter(Slope_x != 5), aes(x = as.factor(round(IFR_x,2)), y = as.factor(round(Slope_x,2)), fill = as.factor(pcr_col_group))) + geom_tile() +
    geom_text(aes(label = round(ll_pcr), colour = (pcr_col_group >= max(pcr_col_group, na.rm=T))), size = 1.5) +
    scale_colour_manual(values = c("white", "black")) +
    ggtitle("PCR prevalence") + xlab("Overall severity") + ylab("IFR age gradient") +
    labs(fill = "Mean ll bin pcr") + theme(legend.position = "none") +
    # xlim(0,2) + ylim(0,2) +
    scale_fill_discrete(type = tail(viridis::viridis(n = 9), length(table(IFR_mat$pcr_col_group))))
  # scale_fill_viridis_c()

  p5_sero <- ggplot(IFR_mat %>% filter(Slope_x != 5), aes(x = as.factor(round(IFR_x,2)), y = as.factor(round(Slope_x,2)), fill = as.factor(sero_col_group))) + geom_tile() +
    geom_text(aes(label = round(ll_sero), colour = (sero_col_group >= max(sero_col_group, na.rm=T))), size = 1.5) +
    scale_colour_manual(values = c("white", "black")) +
    ggtitle("Seroprevalence") + xlab("Overall severity") + ylab("IFR age gradient") +
    labs(fill = "Mean ll bin sero") + theme(legend.position = "none") +
    # xlim(0,2) + ylim(0,2) +
    scale_fill_discrete(type = tail(viridis::viridis(n = 9), length(table(IFR_mat$sero_col_group))))
  # scale_fill_viridis_c()

  p6_pcr_sero_sum <- ggplot(IFR_mat %>% filter(Slope_x != 5), aes(x = as.factor(round(IFR_x,2)), y = as.factor(round(Slope_x,2)), fill = as.factor(pcr_sero_sum_col_group))) + geom_tile() +
    geom_text(aes(label = round(ll_pcr_sero_sum), colour = (pcr_sero_sum_col_group >= max(sero_col_group, na.rm=T))), size = 1.5) +
    scale_colour_manual(values = c("white", "black")) +
    ggtitle("C19 population surveys") + xlab("Overall severity") + ylab("IFR age gradient") +
    labs(fill = "Mean ll bin sero") + theme(legend.position = "none") +
    # xlim(0,2) + ylim(0,2) +
  scale_fill_discrete(type = tail(viridis::viridis(n = 9), length(table(IFR_mat$pcr_sero_sum_col_group))))
    # scale_fill_viridis_c()


  # p6_comb <- ggplot(IFR_mat, aes(x = as.factor(round(IFR_x,2)), y = as.factor(round(Slope_x,2)), fill = comb_col_group)) + geom_tile() +
  #   geom_text(aes(label = round(ll_comb), colour = (comb_col_group >= max(comb_col_group, na.rm=T))), size = 4) +
  #   scale_colour_manual(values = c("white", "black")) +
  #   ggtitle("Population prevalence survey") + xlab("Overall severity") + ylab("IFR age disparity") +
  #   labs(fill = "Mean ll bin comb") + theme(legend.position = "none") +
  #   # xlim(0,2) + ylim(0,2) +
  #   scale_fill_viridis_c()

  # title_gg <- ggplot() +
    # labs(title = Title)
  # browser()
  # Figure1 <- cowplot::plot_grid(title_gg, cowplot::plot_grid(p1_post,p1_like, p3_pois, p2_bin, p4_pcr, p5_sero), ncol = 1, rel_heights = c(0.03, 1))+
  # Figure1 <- cowplot::plot_grid(title_gg, cowplot::plot_grid(p1_post, p3_pois, p2_bin, p6_pcr_sero_sum),#p4_pcr, p5_sero
                                                             # ncol = 1, rel_heights = c(0.03, 1))+
  # Figure1 <- cowplot::plot_grid(p1_post, p3_pois, p2_bin, p6_pcr_sero_sum) +
    # theme_minimal()
  return(list(p1 = p1_post,p2 = p3_pois, p3 = p2_bin, p4 = p6_pcr_sero_sum, IFR_mat = IFR_mat))
# browser()
# pdf("analysis/figures/39_heatmap_overall.pdf", width = 4.5, height = 4.5)
# p1_post + ggtitle("")
# dev.off()

# tiff("analysis/figures/39_heatmap_overall.tiff", units = "in", res = 300, width = 4.5, height = 4.5)
# p1_post + ggtitle("")
# dev.off()

  # return(list(p1_post, p2_bin, p3_pois, p4_pcr, p5_sero, p6_comb))
}

