orderly2::orderly_description(display = "Prepare data for use in analysis")

# orderly2::orderly_artefact(description = "Burial registration dataframes for use in plotting and analyses",
#                            files = c("Burial_registrations_cleaned.rds"))

orderly2::orderly_parameters(raw_data_loc = NULL,
                             derived_data_loc = NULL)

##############################
## Burial Registration Data ##
##############################

# The original dataset has data until mid-June 2021
# The new dataset begins in June 2021
# We will therefore end the first dataset at the end of May 2021

# Total Burial registrations by week
# Full dataset
Burial_registrations_cleaned <-
  # Original data set
  read.csv(paste0(raw_data_loc,"mortuary_records_v2.csv")) |>
  mutate(date = as.Date(dod, "%m/%d/%y")) |>
  filter(date >= as.Date("2017-12-25"), date < as.Date("2021-06-01")) |>
  # filter(date >= as.Date("2021-06-01")) |>
  select(-dod) |>
  # New data
  rbind(read.csv(file = paste0(raw_data_loc, "burial_register_cleaned_2021_22.csv")) |>
          mutate(date = as.Date(date_death, "%m/%d/%y")) |>
          group_by(date) |>
          select(-PLACE_OF_DEATH, -date_death) |>
          rename(sex = GENDER)) |>
  filter(age_years != ".") %>%
  mutate(Age_gr = cut(as.numeric(age_years), c(seq(0,80,by = 5),Inf), right = F, labels = F),
         Age_gr_fig_3c = cut(as.numeric(age_years), c(0,seq(5,25,by = 10),55,Inf), right = F, labels = F),
         Age_gr_fig_3e = cut(as.numeric(age_years), c(0,seq(5,75,by = 10),Inf), right = F, labels = F),
         Age_gr_fig_4 = cut(as.numeric(age_years), c(0,5,15,25,40,50,60,70,Inf), right = F, labels = F),
         Week = lubridate::floor_date(date, unit = "week", week_start = 1))

saveRDS(Burial_registrations_cleaned, file = paste0(derived_data_loc, "Burial_registrations_cleaned.rds"))
