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

Plot_Post <- function(Mod_Res, IFR_mat){
  IFR_vec <- 1:81 %in% as.numeric(gsub(pattern = "X", "",names(Mod_Res)))
  IFR_mat$AvPost[IFR_vec] <- unlist(lapply(Mod_Res, get_Posterior))
  IFR_mat$Post_col_group <- as.numeric(cut(round(IFR_mat$AvPost),
                                           breaks = c(-Inf, round(max(IFR_mat$AvPost,na.rm = T))- c(500,200,100,50,20,12,8,4,0)),
                                           labels = as.character(c(1:9))))

  return(IFR_mat)
}

Heatmap_Post <- function(x){
  ggplot(x, aes(x = as.factor(round(IFR_x,2)), y = as.factor(round(Slope_x,2)), fill = as.factor(Post_col_group))) +
    geom_tile() +
    geom_text(aes(label = round(AvPost), colour = (Post_col_group >= max(Post_col_group, na.rm=T)))) +
    scale_colour_manual(values = c("white", "black")) +
    labs(x = "Relative overall severity", y = "Relative IFR age gradient") +
    labs(fill = "Mean Post") +
    scale_fill_discrete(type = tail(viridis::viridis(n = 9), length(table(x$Post_col_group)))) +
    coord_cartesian(xlim = c(1.5,9.5), ylim = c(2.5,7.5), expand = F) +
    scale_x_discrete(expand = c(0,0), breaks = c(0.4,0.6,0.8,1,1.25,1.67,2.5,5), labels = paste0(100*c(0.4,0.6,0.8,1,1.25,1.67,2.5,5),"%")) +
    scale_y_discrete(expand = c(0,0), breaks = c(0.6,0.8,1,1.25,1.67), labels = paste0(100*c(0.6,0.8,1,1.25,1.67),"%")) +
    theme_minimal(base_size = 7, base_family = "Helvetica") +
    theme(axis.text.y = element_text(angle = 0, hjust = 0.5),
          legend.position = "none", plot.title = element_text(size = 7, face = "bold"))
}



Plot_Samples <- function(Res_supp_data, pcr_sero_data, Select_Runs){

  # For supplementary Figure 6

  library(ggh4x)

  ## Get data for plotting by week
  Res_supp_data_combined_week <- str2str::a2d(str2str::ld2a(lapply(1:length(Res_supp_data), function(x){
    df_tmp <- cbind(IFR_coefficients[Select_Runs[x],c("IFR_x","Slope_x")],merge(Res_supp_data[[x]]$Plotting_data_integrated_summarised_Weeks,Res_supp_data[[x]]$Mod_Age_Deaths_Lus_Av_Week))
    rownames(df_tmp) <- NULL
    df_tmp})), col = 2) %>%
    mutate(IFR_x = paste0(100*round(as.numeric(IFR_x),2),"%"),
           Slope_x = paste0(100*round(as.numeric(Slope_x),2),"%"),
           Facet_Label = "IFR Age Gradient")

  ## Get data for plotting by age group
  Res_supp_data_combined_age <- str2str::a2d(str2str::ld2a(lapply(1:length(Res_supp_data), function(x){
    df_tmp <- cbind(IFR_coefficients[Select_Runs[x],c("IFR_x","Slope_x")],merge(Res_supp_data[[x]]$Plotting_data_integrated_summarised_Ages,Res_supp_data[[x]]$Mod_Age_Deaths_Lus_Av_Age[,1:4]))
    rownames(df_tmp) <- NULL
    df_tmp})), col = 2) %>%
    mutate(IFR_x = paste0(100*round(as.numeric(IFR_x),2),"%"),
           Slope_x = paste0(100*round(as.numeric(Slope_x),2),"%"),
           Facet_Label = "IFR Age Gradient")

  ## Get data for population PCR prevalence and seroprevalence
  Res_supp_data_combined_pcr_sero <- do.call("rbind",lapply(1:length(Res_supp_data), function(x){
    df_tmp <- cbind(x, IFR_coefficients[Select_Runs[x],c("IFR_x","Slope_x")],Res_supp_data[[x]][[2]])
    rownames(df_tmp) <- NULL
    df_tmp}))


  ## Burial registrations by week (Row 1)
  ps_1 <- ggplot(Res_supp_data_combined_week %>% filter(Measure =="Mod_tot_ds_morgue"), aes(x = as.Date(date.x), group = X3)) +
    geom_line(aes(y = as.numeric(Median), color = IFR_x)) +
    geom_point(aes(y = as.numeric(BurRegs)), size = 0.9) +
    facet_nested(~Facet_Label + factor(Slope_x, levels = c("40%","100%","250%"))) +
    theme_minimal(base_size = 7, base_family = "Helvetica") +
    labs(title = "a", x = "Date", y = "Burial\nRegistrations") +
    theme(legend.key = element_rect(fill = "white", linetype = 0)) +
    viridis::scale_color_viridis(option = "A", discrete = T, end = 0.8, breaks =c("40%","60%","80%","100%","125%","167%","250%")) +
    labs(color = "Overall severity") +
    theme(plot.margin = margin(b=0, unit = "cm"),
          plot.title = element_text(size= 7, face = "bold", vjust = -12))


  ## Burial registrations by age (Row 2)
  Res_supp_data_combined_age$Age_gr_label <- factor(Res_supp_data_combined_age$Age_gr_label, levels = unique(Res_supp_data_combined_age$Age_gr_label))
  ps_2 <- ggplot(Res_supp_data_combined_age %>% filter(Measure =="Mod_tot_ds_morgue") , aes(x = Age_gr_label, group = X3)) +
    geom_line(aes(y = as.numeric(Median), color = as.factor(IFR_x))) +
    geom_point(aes(y = as.numeric(BurRegs)), size = 0.9) +
    theme_minimal(base_size = 7, base_family = "Helvetica") +
    labs(title = "b", x = "Age", y = "Burial\nRegistrations") +
    theme(legend.key = element_rect(fill = "white", linetype = 0),
          axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
          strip.text.x = element_blank(),
          plot.margin = margin(t=0, b=0, unit = "cm"),
          plot.title = element_text(size= 7, face = "bold")) +
    facet_nested(~factor(Slope_x, levels = c("40%","100%","250%"))) +
    viridis::scale_color_viridis(option = "A", discrete = T, end = 0.8,breaks =c("40%","100%","250%")) +
    labs(color = "Overall severity")

  ## Post-mortem PCR prevalence by week (Row 3)
  ps_3 <- ggplot(data = Res_supp_data_combined_week %>% filter(Measure =="Pos_prev"), aes(x = as.Date(date.x))) +
    geom_line(aes(y = as.numeric(Median), color = as.factor(IFR_x))) +
    geom_point(aes(y = as.numeric(PosTests)/as.numeric(Samples)), size = 0.9) +
    geom_errorbar(aes(ymin = Hmisc::binconf(as.numeric(PosTests),as.numeric(Samples))[,"Lower"],
                      ymax = Hmisc::binconf(as.numeric(PosTests),as.numeric(Samples))[,"Upper"]), linewidth = 0.3) +
    theme_minimal(base_size = 7, base_family = "Helvetica") +
    labs(title = "c", x = "Date", y = "Post-mortem\nprevalence") +
    coord_cartesian(ylim = c(0,0.90)) +
    theme(plot.title = element_text(size = 7, face = "bold"),
          strip.text.x = element_blank(),
          plot.margin = margin(t=0, b=0, unit = "cm")) +
    scale_y_continuous(labels = scales::percent) +
    facet_nested(~factor(Slope_x, levels = c("40%","100%","250%"))) +
    viridis::scale_color_viridis(option = "A", discrete = T, end = 0.8,breaks =c("40%","100%","250%")) +
    labs(color = "Overall severity")

  ## Post-mortem PCR prevalence by age (Row 4)
  ps_4 <- ggplot(Res_supp_data_combined_age %>% filter(Measure == "Pos_prev"), aes(x = Age_gr_label, group = X3)) +
    geom_line(aes(y = as.numeric(Median), color = as.factor(IFR_x))) +
    geom_point(aes(y = as.numeric(PosTests)/as.numeric(Samples)), size = 0.9) +
    geom_errorbar(aes(ymin = Hmisc::binconf(as.numeric(PosTests),as.numeric(Samples))[,"Lower"],
                      ymax = Hmisc::binconf(as.numeric(PosTests),as.numeric(Samples))[,"Upper"]), linewidth = 0.3) +
    theme_minimal(base_size = 7, base_family = "Helvetica") +
    labs(title = "d", x = "Age", y = "Post-mortem\nprevalence") +
    scale_x_discrete(limits = c(paste0(1:16*5-5, "-",1:16*5-1),"80+")) +
    theme(plot.title = element_text(size = 7, face = "bold"),
          axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
          strip.text.x = element_blank(),
          plot.margin = margin(t=0, b=0, unit = "cm")) +
    scale_y_continuous(labels = scales::percent) +
    facet_nested(~factor(Slope_x, levels = c("40%","100%","250%"))) +
    viridis::scale_color_viridis(option = "A", discrete = T, end = 0.8,breaks =c("40%","100%","250%")) +
    labs(color = "Overall severity")

  ## Population PCR and seroprevalence (Row 5)
  ps_5 <- ggplot(Res_supp_data_combined_pcr_sero %>%
                   mutate(IFR_x = paste0(100*round(as.numeric(IFR_x),2),"%"),
                          Slope_x = paste0(100*round(as.numeric(Slope_x),2),"%")),
                 aes(x = date, y = pcr_perc_median, group = x, color = as.factor(IFR_x))) +
    coord_cartesian(xlim = c(as.Date("2020-06-15"), as.Date("2020-08-05")),
                    ylim = c(0, 0.10)) +
    theme_minimal(base_size = 7, base_family = "Helvetica") +
    labs(title = "e", x = "Date", y = "Population\nprevalence") +
    scale_linetype_manual(name="Model",
                          breaks = c("Attack rate", "PCR %", "Sero %"),
                          values = c(1,1,2)) +
    scale_y_continuous(labels = scales::percent) +
    geom_line(aes(x=date, y=pcr_perc_median, linetype = "PCR %")) +
    geom_line(aes(x=date, y=sero_perc_median, linetype = "Sero %")) +
    geom_point(data = pcr_sero_data$pcr_df,
               aes(x= as.Date(date_start) + 0.5*(as.Date(date_end)-as.Date(date_start)),
                   y=pos_tests/samples, shape = "PCR %"), color = "black", size = 1.4, inherit.aes = F) +
    geom_errorbar(data = pcr_sero_data$pcr_df,
                  aes(ymin=Hmisc::binconf(pos_tests,samples)[,"Lower"],
                      ymax=Hmisc::binconf(pos_tests,samples)[,"Upper"],
                      x=as.Date(date_start) + 0.5*(as.Date(date_end)-as.Date(date_start)),
                      width=5,linetype = "PCR %"), color = "black", inherit.aes = F) +
    geom_errorbarh(data = pcr_sero_data$pcr_df,
                   aes(xmin=as.Date(date_start),xmax=as.Date(date_end),y=pos_tests/samples, height=0,linetype = "PCR %"),
                   color = "black", inherit.aes = F) +
    geom_point(data = pcr_sero_data$sero_df,
               aes(x= as.Date(date_start) + 0.5*(as.Date(date_end)-as.Date(date_start)),
                   y=pos_tests/samples, shape = "Sero %"), size = 1.4, color = "black", inherit.aes = F) +
    geom_errorbar(data = pcr_sero_data$sero_df,
                  aes(ymin=Hmisc::binconf(pos_tests,samples)[,"Lower"],
                      ymax=Hmisc::binconf(pos_tests,samples)[,"Upper"],
                      x=as.Date(date_start) + 0.5*(as.Date(date_end)-as.Date(date_start)),
                      width=5,linetype = "Sero %"), color = "black", inherit.aes = F) +
    geom_errorbarh(data = pcr_sero_data$sero_df,
                   aes(xmin=as.Date(date_start),xmax=as.Date(date_end),y=pos_tests/samples, height=0,linetype = "Sero %"),
                   color = "black", inherit.aes = F) +
    facet_nested(~factor(Slope_x, levels = c("40%","100%","250%"))) +
    viridis::scale_color_viridis(option = "A", discrete = T, end = 0.8,breaks =c("40%","100%","250%"), guide = "none") +
    labs(color = "Overall severity",
         shape = "Data") +
    theme(plot.title = element_text(size = 7, face = "bold"),
          strip.text.x = element_blank(),
          plot.margin = margin(t=0, unit = "cm"))


  return(list(ps_1=ps_1, ps_2=ps_2, ps_3=ps_3, ps_4=ps_4, ps_5=ps_5))
}





Plot_Rt_Samples <- function(Res_supp_data, pcr_sero_data, Select_Runs){

  ## For supplementary Figure 8

  # Get modelled attack rate
  Res_supp_data_combined_pcr_sero <- do.call("rbind",lapply(1:length(Res_supp_data), function(x){
    df_tmp <- cbind(x, IFR_coefficients[Select_Runs[x],c("IFR_x","Slope_x")],Res_supp_data[[x]][[2]])
    rownames(df_tmp) <- NULL
    df_tmp})) %>%
    mutate(IFR_x = paste0(100*round(IFR_x,2),"%"),
           Slope_x = paste0(100*round(Slope_x,2),"%"))


  # Attack rates plot
  ps_1 <- ggplot(Res_supp_data_combined_pcr_sero,
                 aes(x = date, y = attack_rate_median, group = x, fill = as.factor(IFR_x), color = as.factor(IFR_x))) +
    labs(title = "a", y = "Attack rate", x = "Date", color = "Overall severity") +
    coord_cartesian(xlim = c(as.Date("2020-06-01"), as.Date("2020-10-05")),
                    ylim = c(0, 0.40)) +
    theme_minimal(base_size = 7, base_family = "Helvetica") +
    scale_y_continuous(labels = scales::percent) +
    geom_line() +
    geom_ribbon(aes(ymin=attack_rate_ci$CI_low, ymax=attack_rate_ci$CI_high), alpha=0.2, linewidth = 0.2, show.legend = F) +
    viridis::scale_color_viridis(option = "A", discrete = T, end = 0.8, breaks =c("80%","100%","125%","167%")) +
    viridis::scale_fill_viridis(option = "A", discrete = T, end = 0.8, breaks =c("80%","100%","125%","167%")) +
    theme(strip.text.x = element_blank(),
          plot.margin = margin(t=0, unit = "cm"),
          plot.title = element_text(size = 7, face = "bold"))


  # Get Rt trend data
  Res_supp_rt_data <- lapply(Res_supp_data, rt_data_immunity_supp)

  Res_supp_data_rt <- do.call("rbind",lapply(1:length(Res_supp_rt_data), function(x){
    df_tmp <- cbind(x, IFR_coefficients[Select_Runs[x],c("IFR_x","Slope_x")],Res_supp_rt_data[[x]])
    rownames(df_tmp) <- NULL
    df_tmp})) %>%
    mutate(IFR_x = paste0(100*round(IFR_x,2),"%"),
           Slope_x = paste0(100*round(Slope_x,2),"%"))


  # R_eff plot
  Reff_plot <- ggplot(Res_supp_data_rt %>% filter(
    date > as.Date("2020-06-01") & date <= as.Date("2020-10-05")), aes(color = as.factor(IFR_x), fill = as.factor(IFR_x))) +
    geom_ribbon(mapping = aes(x = date, ymin = Reff_min, ymax = Reff_max), alpha = 0.2, linewidth = 0.2, show.legend = F) + #, fill = "#48996b") +
    geom_line(mapping = aes(x = date, y = Reff_median)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    viridis::scale_color_viridis(option = "A", discrete = T, end = 0.8, breaks =c("80%","100%","125%","167%")) +
    viridis::scale_fill_viridis(option = "A", discrete = T, end = 0.8, breaks =c("80%","100%","125%","167%")) +
    theme_minimal(base_size = 7, base_family = "Helvetica") +
    labs(title = "c", y = expression(italic(R)[eff]), x = "Date", color = "Overall severity") +
    scale_x_date(limits = as.Date(c(as.character("2020-06-01"),
                                    as.character("2020-10-05"))),expand = c(0,0)) +
    labs(color = "Overall severity") +
    theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black"),
          plot.title = element_text(size = 7, face = "bold")
    )

  # R_t plot
  Rt_plot <- ggplot(Res_supp_data_rt %>% filter(
    date > as.Date("2020-06-01") & date <= as.Date("2020-10-05")), aes(color = as.factor(IFR_x), fill = as.factor(IFR_x))) +
    geom_ribbon(mapping = aes(x = date, ymin = Rt_min, ymax = Rt_max), alpha = 0.2, linewidth = 0.2, show.legend = F) + #, fill = "#48996b") +
    geom_line(mapping = aes(x = date, y = Rt_median)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    viridis::scale_color_viridis(option = "A", discrete = T, end = 0.8, breaks =c("80%","100%","125%","167%")) +
    viridis::scale_fill_viridis(option = "A", discrete = T, end = 0.8, breaks =c("80%","100%","125%","167%")) +
    theme_minimal(base_size = 7, base_family = "Helvetica") +
    labs(title = "b", y = expression(italic(R)[0](italic(t))), x = "Date", color = "Overall severity") +
    scale_x_date(limits = as.Date(c(as.character("2020-06-01"),
                                    as.character("2020-10-05"))),expand = c(0,0)) +
    labs(color = "Overall severity") +
    theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, colour = "black"),
          plot.title = element_text(size = 7, face = "bold")
    )


  return(list(ps_1=ps_1, Rt_plot=Rt_plot, Reff_plot=Reff_plot))
}

