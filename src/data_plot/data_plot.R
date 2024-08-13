orderly2::orderly_description(display = "Plot burial registration data")

orderly2::orderly_artefact(description = "Burial registration figures",
                           files = c("Burial_Registrations.pdf"))

orderly2::orderly_parameters(derived_data_loc = NULL)

# orderly2::orderly_dependency(
#   name = "data_handle",
#   query = "latest",
#   files = "Burial_registrations_cleaned.rds")

library(ggplot2)

Burial_registrations_cleaned <- readRDS(paste0(derived_data_loc, "Burial_registrations_cleaned.rds"))

# Plot 1: total registrations
Total_burial_registration_by_week_2018_2022 <- Burial_registrations_cleaned |>
  group_by(Week) %>%
  summarise(BurRegs = n()) |>
  filter(Week >= "2018-01-01") |>
  mutate(Roll_Av =zoo::rollapply(BurRegs,5,mean,fill=NA))

plot_basic <- ggplot(mapping = aes(x = Week)) +
  geom_vline(xintercept = as.Date(c("2020-03-17",
                                    "2020-04-24",
                                    "2020-06-01",
                                    "2020-10-05")),
             color = "black", linetype = "dashed", linewidth = 0.6) +
  scale_x_date(date_breaks = "1 month", date_labels =  "%b %Y", expand = c(0,0)) +
  coord_cartesian(xlim = as.Date(c("2017-12-25","2022-12-31"))) +
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        plot.title=element_text(face = "bold", hjust = 0, vjust = 0, size = 7),
        legend.text = element_text(hjust = 0),
        legend.justification = c("left","center"))

p1 <- plot_basic +
  geom_line(data = Total_burial_registration_by_week_2018_2022, aes(y = BurRegs, alpha = "Burial\nregistrations", linewidth = "Burial\nregistrations")) +
  geom_line(data = Total_burial_registration_by_week_2018_2022, aes(y = Roll_Av, alpha = "5-week\nRolling average", linewidth = "5-week\nRolling average")) +
  coord_cartesian(ylim = c(0,1000)) +
  labs(title = "a",
       y = "Burial\nregistrations") +
  scale_alpha_manual("", values = c(0.5,1), breaks = c("Burial\nregistrations","5-week\nRolling average")) +
  scale_linewidth_manual("", values = c(0.4,0.6), breaks = c("Burial\nregistrations","5-week\nRolling average")) +
  guides(linewidth = guide_legend(override.aes = list(alpha = c(0.3,1)))) +
  theme(legend.key.height = unit(6, 'mm'))

# Plot 2: registration relative age rates
BurRegs_Week_Age_Fig3c <- Burial_registrations_cleaned |>
  group_by(Week, Age_gr_fig_3c) |>
  summarise(BurRegs = n()) |>
  rename(Age_gr = Age_gr_fig_3c) |>
  ungroup() %>%
  tidyr::complete(Week, Age_gr, fill = list(BurRegs = 0))

# Calculates and merges the mean baseline for each age group
Rel_BurRegs_Week_Age_Fig3c <- BurRegs_Week_Age_Fig3c |>
  filter(Week>="2018-01-01", Week<"2020-01-01") |>
  group_by(Age_gr) |>
  summarise(PrePan_av_BurRegs = mean(BurRegs)) |>
  mutate(Age_gr_labs = c(paste0(c(0,5,15,25),"-",c(4,14,24,54)),"55+")) |>
  merge(BurRegs_Week_Age_Fig3c) |>
  mutate(Rel_BurRegs = BurRegs/PrePan_av_BurRegs)

p2 <- plot_basic +
  geom_line(data = Rel_BurRegs_Week_Age_Fig3c, aes(y = Rel_BurRegs, color = Age_gr_labs), alpha = 1, linewidth = 0.4) +
  labs(title = "b",
       y = "Relative\nregistrations") +
  theme(legend.key.height = unit(4, 'mm')) +
  scale_color_manual("Age group", values = viridis::turbo(n = 5), breaks = c("0-4","5-14","15-24","25-54","55+"))

# Plot 3: registration average age
Av_Age <- Burial_registrations_cleaned |>
  group_by(Week) %>%
  summarise(av_age = mean(as.numeric(age_years))) %>%
  mutate(Roll_av =zoo::rollapply(av_age,5,mean,fill=NA))

p3 <- plot_basic +
  geom_line(data = Av_Age, aes(y = av_age, x = Week, color = "Average age"), linewidth = 0.4, inherit.aes = F) +
  geom_line(data = Av_Age, aes(y = Roll_av, x = Week, color = "5-week\nRolling average"), linewidth = 0.6, inherit.aes = F) +
  labs(y = "Average age\nat death (years)",
       title = "c") +
  theme(legend.spacing.y = unit(-0.5, 'cm'),
        legend.key.height = unit(5, 'mm')) +
  scale_color_manual("", values = c("darkgrey", "black"), breaks = c("Average age","5-week\nRolling average"))

# Plot 4: registration age proportions
Age_Proportion_BurRegs <- Burial_registrations_cleaned %>%
  group_by(Week, Age_gr_fig_3e) %>%
  summarise(BurRegs = n()) %>%
  rename(Age_gr = Age_gr_fig_3e) |>
  ungroup() %>% tidyr::complete(Week, Age_gr, fill = list(BurRegs = 0)) %>%
  group_by(Week) %>%
  mutate(Weekly_Total_deaths = sum(BurRegs),
         Weekly_Prop_deaths = BurRegs/Weekly_Total_deaths) |>
  merge(data.frame(Age_gr = 1:9, Age_gr_labs = c(paste0(c(0,seq(5,65,by=10)),"-",c(seq(4,74,by=10))),"75+")))

Age_Proportion_BurRegs$Age_gr_labs <- factor(x = Age_Proportion_BurRegs$Age_gr_labs,
                                             levels = unique(Age_Proportion_BurRegs$Age_gr_labs))

p4 <- plot_basic +
  geom_area(data = Age_Proportion_BurRegs, aes(y = Weekly_Prop_deaths, group = Age_gr_labs, fill = Age_gr_labs)) +
  scale_y_continuous(expand = c(0,0)) +
  geom_vline(xintercept = as.Date(c("2020-03-17",
                                    "2020-04-24",
                                    "2020-06-01",
                                    "2020-10-05")),
             color = "white", linetype = "dashed", linewidth = 0.6) +
  labs(title = "d",
       x = "Date",
       y = "Proportion of\nregistrations") +
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.key.size = unit(4,"mm")) +
  viridis::scale_fill_viridis(discrete = T, option = "H") +
  guides(fill=guide_legend(title="Age group", label.hjust = 0),
         color = guide_legend(title = "Age group",
                              override.aes = list(color = c("black","white"),
                                                  bordercolor = c("white", viridis::turbo(n = 9, begin = 0, end =1)[9]),
                                                  bordersize = c(NULL, 2.5))))


ggsave(filename = "Burial_Registrations.pdf",
       plot = cowplot::plot_grid(p1,p2,p3,p4, ncol = 1, rel_heights = c(1.5,1.5,1.5,3), align = "v"),
       width = 180, height = 180, units = "mm")
