library(tidyverse)

# Calculating HC50_EC50 for the HESTIA_envirotox Database for backwards compatibility with previous LCA-tox versions.

# function for calculating geometric mean
gmean <- function(x) exp(mean(log(x)))

# function for creating a data frame with HC50 calculations and some quality metrics
HC50_calc_function <- function(input_df) {

  HC50_output <- {{input_df}} %>% 
  # selecting only EC50 data
  filter(Endpoint_conv == "EC50") %>%  
  # Standard USEtox 2.1 ACR conversions (Acute/2 = Chronic) apart from pesticides, where ACR = 2.2 according to Fantke et al., 2017
  mutate(
    Chronic_EC50_mg.l = case_when(
      AcuteChronic == "Acute" & Group == "Pesticide" ~ Value.mg_l/2.2,
      AcuteChronic == "Acute" & Group != "Pesticide" ~ Value.mg_l/2,
      TRUE ~ Value.mg_l),
    Chronic_EC50_mg.l_ac = case_when(
      AcuteChronic == "Acute" ~ Chronic_EC50_mg.l,
      TRUE ~ as.numeric(NA)),
    Chronic_EC50_mg.l_chr = case_when(
      AcuteChronic == "Chronic" ~ Chronic_EC50_mg.l,
      TRUE ~ as.numeric(NA)),
  ) %>% 
  group_by(CAS.Number, AcuteChronic, Taxonomy.Group, Species) %>% 
  # Performing HC50EC50 calculations 
  # 1. species averages 
  summarize(
    n_all = n(),
    sd_all = sd(Chronic_EC50_mg.l, na.rm = T),
    avg_sp_EC50 = as.numeric(gmean(Chronic_EC50_mg.l)),
    avg_sp_EC50_acute = as.numeric(gmean(Chronic_EC50_mg.l_ac)), 
    avg_sp_EC50_chronic = as.numeric(gmean(Chronic_EC50_mg.l_chr)),
    n_chronic = sum(!is.na(Chronic_EC50_mg.l_chr)), # makes a logical operation to evaluate if NA or not NA, sum(all that is true)
  ) %>% 
  ungroup() %>% 
  # 2. Substance averages
  group_by(CAS.Number) %>% 
  summarize(
    n_data_points_HC50 = sum(n_all, na.rm = T),
    n_chronic = sum(n_chronic),
    avlog10_HC50_acute = as.numeric(mean(log10(avg_sp_EC50_acute/1000), na.rm = T)),
    n_species_acute_HC50 = sum(!is.na(avg_sp_EC50_acute), na.rm = T),
    avlog10_HC50_chronic = mean(log10(avg_sp_EC50_chronic/1000), na.rm = T),
    n_species_Chronic_HC50 = sum(!is.na(avg_sp_EC50_chronic), na.rm = T),
    n_Taxonomy.Group_all = n_distinct(Taxonomy.Group[!is.na(avg_sp_EC50)]),
    n_Taxonomy.Group_cr = n_distinct(Taxonomy.Group[!is.na(avg_sp_EC50_chronic)]),
  ) %>% 
  ungroup() %>% 
  # Prioritizing Chronic HC50 values, but using Acute if no Chronic are present
  mutate(
    # transforming NaN into NA
    avlog10_HC50_chronic = case_when(is.nan(avlog10_HC50_chronic) ~ as.numeric(NA), TRUE ~ avlog10_HC50_chronic),
    avlog10_HC50_EC50 = case_when(
      # Placing this operation at the top, else it masks the HC50_chronic values
      is.na(avlog10_HC50_chronic) & !is.na(avlog10_HC50_acute) ~ avlog10_HC50_acute, 
      !is.na(avlog10_HC50_chronic)  ~ avlog10_HC50_chronic,
      TRUE ~ as.numeric(NA)),
    # Annotating data whether Chronic or Acute data is forming the HC50 value.
    HC50_AcuteChronic = case_when(
      is.na(avlog10_HC50_chronic) & !is.na(avlog10_HC50_acute) ~ "Acute",
      !is.na(avlog10_HC50_chronic)  ~ "Chronic",
      TRUE ~ as.character(NA)),
    # Annotating where data comes from, and if it is recommended (>3 USEtox categories) 
    USEtox_HC50_recommend = case_when(
      n_Taxonomy.Group_cr >= 3 ~ "Chronic-Recommended",
      n_Taxonomy.Group_cr < 3 & n_chronic > 0 ~ "Chronic",
      n_Taxonomy.Group_cr == 0 & n_Taxonomy.Group_cr >= 3 ~ "Acute-Recommended",
      TRUE ~ "Acute")
  )
return(HC50_output)
}


