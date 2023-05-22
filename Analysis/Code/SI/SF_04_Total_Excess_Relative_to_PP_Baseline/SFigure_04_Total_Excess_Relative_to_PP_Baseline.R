## Load Data
Mort_excess_deaths_grouped <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/results/Excess_Deaths_for_figure_4.rds")

Median_18_19  <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Burial_registrations_by_week_2018_2021_Fig4_age_groups.rds") |>
  group_by(Week) |>
  summarise(BurRegs = sum(BurRegs)) |>
  pull(BurRegs) |>
  median()

Median_18_19_50p <- readRDS("~/Desktop/COVID_IFR_Lusaka_Restricted_Access/derived_data/Burial_registrations_by_week_2018_2021_Fig4_age_groups.rds") |>
  filter(Age_gr %in% 6:8) |>
  group_by(Week) |>
  summarise(BurRegs = sum(BurRegs)) |>
  pull(BurRegs) |>
  median()

### sum excess for age groups of relevance
Proportion_Age_18_19 <- Mort_excess_deaths_grouped |>
  group_by(Week, Sample) |>
  summarise(Excess_std=sum(Excess_std)) |>
  group_by(Week) |>
  summarise_at(c("Excess_std"),funs(median = median(.), CI_low = bayestestR::ci(.)$CI_low, CI_high = bayestestR::ci(.)$CI_high)) |>
  mutate(Proportion = median/Median_18_19,
         Proportion_lower = CI_low/Median_18_19,
         Proportion_higher = CI_high/Median_18_19)

Proportion_Age_18_19_50plus <- Mort_excess_deaths_grouped |>
  filter(Age_gr %in% 6:8) |>
  group_by(Week, Sample) |>
  summarise(Excess_std=sum(Excess_std)) |>
  group_by(Week) |>
  summarise_at(c("Excess_std"),funs(median = median(.), CI_low = bayestestR::ci(.)$CI_low, CI_high = bayestestR::ci(.)$CI_high)) |>
  mutate(Proportion = median/Median_18_19_50p,
         Proportion_lower = CI_low/Median_18_19_50p,
         Proportion_higher = CI_high/Median_18_19_50p)

Proportion_Age_18_19_combined <- rbind(Proportion_Age_18_19 |> mutate(Age = "a"),
                                       Proportion_Age_18_19_50plus |>  mutate(Age = "b"))


Total_relative_plot <- ggplot(data = Proportion_Age_18_19_combined, aes(x = Week)) +
  geom_line(aes(y = Proportion)) +
  geom_ribbon(aes(ymin = Proportion_lower, ymax = Proportion_higher), alpha = 0.2) +
  facet_wrap(~Age, nrow = 2) +
  theme_minimal(base_size = 7, base_family = "Helvetica") +
  theme(plot.title = element_text(size = 7, face = "bold"), strip.text = element_text(size = 7, face = "bold", hjust = 0)) +
  labs(x = "Date",
       y = "Excess deaths relative to median 2018-2019 baseline") +
  coord_cartesian(ylim = c(-1, 2)) +
  scale_x_date(expand = c(0,0)) +
  scale_y_continuous(labels = scales::percent)


ggsave("Analysis/Figures/Supplementary_Figures/SFigure_04_Total_Excess_relative_to_PP_baseline.pdf", Total_relative_plot,
       width = 180, height = 180*5/10, units = "mm")

tiff("Analysis/Figures/Supplementary_Figures/SFigure_04_Total_Excess_relative_to_PP_baseline.tiff",
       width = 180, height = 180*5/10, units = "mm", res = 300)
Total_relative_plot
dev.off()
