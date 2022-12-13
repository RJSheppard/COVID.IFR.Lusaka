df <- read.csv("Analysis/Data/raw_data/SeroTracker_ Serosurveys Reporting Prevalence-Grid view_filtered.csv")
names(df)
library(dplyr)
unique(df$Country)

unique(df$Grade.of.Estimate.Scope)
unique(df$Sample.Frame..age.)
unique(df$Sample.Frame..groups.of.interest.)
unique(df$Overall.Risk.of.Bias..JBI.)
unique(df$Data.Quality.Status)
unique(df$Sampling.Method)


df_filtered <- df %>% select(Publication.Date, Country, Test.Type, Grade.of.Estimate.Scope, Sampling.End.Date, Sample.Frame..age., Sample.Frame..groups.of.interest.,Sampling.Method, Overall.Risk.of.Bias..JBI., Source.Type, Data.Quality.Status) %>%
  mutate(Sampling.End.Date = as.Date(Sampling.End.Date, format = "%d-%b-%y")) %>%
  filter(Sampling.End.Date <as.Date("2021-01-01"),
         Overall.Risk.of.Bias..JBI. %in% c("Low","Moderate"),
         Data.Quality.Status == "Verified",
         Grade.of.Estimate.Scope != "Local",
         Sample.Frame..groups.of.interest. %in% c("Household and community samples"),
         !Sampling.Method %in% c("Self-referral", "Convenience","Unclear", "Sequential",""))

df_filtered %>% summarise(length(Country), length(unique(Country)))

df_filtered %>% filter(!Country %in% c("Viet Nam","United States of America",
                         "United Kingdom of Great Britain and Northern Ireland","Switzerland","Sweden",
                         "Spain","Slovenia","Russian Federation",
                         "Peru","Pakistan","Oman",
                         "occupied Palestinian territory - including east Jerusalem","Netherlands",
                         "Nepal","Mongolia","Mexico",
                         "Luxembourg","Lithuania","Lao People's Democratic Republic",
                         "Jordan","Jersey","Japan",
                         "Italy","Israel","Iran (Islamic Republic of)",
                         "Indonesia","India","Iceland",
                         "Hungary","Honduras","Germany",
                         "France","Finland","Faroe Islands",
                         "Denmark","Czechia","Colombia",
                         "China","Chile",
                         "Brazil","Andorra","Afghanistan")) %>%
  summarise(length(Country), length(unique(Country)))
