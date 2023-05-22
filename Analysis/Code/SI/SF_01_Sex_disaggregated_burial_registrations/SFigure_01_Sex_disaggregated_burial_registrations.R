####################################
###### Supplementary Figure 1 ######
####################################
# Burial registrations by gender and age
Total_burial_registrations_gender_clean <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Burial_registrations_cleaned.rds") %>%
  filter(sex != "U", Week >="2018-01-01",Week < "2021-06-10")

# Panel a
Total_burial_registrations_gender_by_week <- Total_burial_registrations_gender_clean |>
  group_by(Week,sex) %>%
  summarise(Total_deaths = n())

p_tot_bur_reg_all_years_gender <- ggplot(Total_burial_registrations_gender_by_week, aes(x = Week, y = Total_deaths, color = sex)) +
  geom_line() +
  labs(x = "Date",
       y = "Total burial\nregistrations",
       title = "a") +
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme(legend.position = "none", plot.title = element_text(size = 7, face = "bold")) +
  viridis::scale_colour_viridis("", discrete = T, begin = 0.3, end = 0.65, option = "B") +
  scale_x_date(expand = c(0.01,0))

# Panel b
Total_burial_registrations_gender_medians <- Total_burial_registrations_gender_clean |>
  group_by(Age_gr_fig_4, Week, sex) |>
  summarise(Total_deaths = n()) |>
  group_by(Age_gr_fig_4, sex) |>
  summarise(Median_deaths = median(Total_deaths))

Relative_gender_registrations <- Total_burial_registrations_gender_clean %>%
  group_by(Age_gr_fig_4, sex, Week) %>%
  summarise(Total_deaths = n()) |>
  merge(Total_burial_registrations_gender_medians) |>
  mutate(Rel_deaths =Total_deaths/Median_deaths)

Age_labels <- c(paste0("Age: ",c(0,5, 15,25,40,50,60),"-",c(5, 15,25,40,50,60,70)-1),"Age: 70+")
names(Age_labels) <- 1:8

p_tot_bur_reg_all_years_gender_age <- ggplot(Relative_gender_registrations, aes(x = Week, y = Rel_deaths, color = sex)) +
  geom_line() +
  facet_wrap(~Age_gr_fig_4, nrow = 4, labeller = labeller(Age_gr_fig_4 = Age_labels)) +
  theme_minimal() +
  labs(x = "Date",
       y = "Burial registrations relative to 2018-2019 median",
       title = "b") +
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme(legend.position = "bottom", plot.title = element_text(size = 7, face = "bold")) +
  viridis::scale_colour_viridis("", discrete = T, begin = 0.3, end = 0.65, option = "B", labels = c("Female","Male")) +
  scale_x_date(expand = c(0.01,0))

ggsave("Analysis/Figures/Supplementary_Figures/SFigure_1_Burial_registrations_by_age.pdf",
       cowplot::plot_grid(p_tot_bur_reg_all_years_gender,p_tot_bur_reg_all_years_gender_age, nrow = 2, rel_heights = c(1,4), align = "v", axis = "lr"),
       width = 180, height = 180, units = "mm")

tiff(filename = "Analysis/Figures/Supplementary_Figures/SFigure_1_Burial_registrations_by_age.tiff",
     width = 180, height = 180, units = "mm", res = 300)
cowplot::plot_grid(p_tot_bur_reg_all_years_gender,p_tot_bur_reg_all_years_gender_age, nrow = 2, rel_heights = c(1,4), align = "v", axis = "lr")
dev.off()
