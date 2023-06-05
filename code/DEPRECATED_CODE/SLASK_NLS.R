EC10eq <- Plot_data_one_CAS %>% 
  mutate(Li = log(EC10eq)) %>% 
  group_by(CAS.Number, Taxonomy.Group, Species) %>% 
  # Log10-transforming effect concentrations, then calculating the arithmetic mean per species' log10 effect concentration
  summarise(
    n_samples = n(),
    n_extrapolations = sum(No.Extrapolations),
    # sp_mean = Per-species average
    sp_mean = mean(Li),
    # sd_samples = per-species st.dev, this uses denominator n - 1.
    sd_samples = sd(Li, na.rm = T),
    # sp_error = per-species standard error (CI?)
    sp_error = sd_samples/(sqrt(n_samples)),
    #cum_norm_dist = cum_norm_dist_function(x = Li, M = sp_mean, S = sd_samples) # Måste kanske köra en for-loop???
  ) %>% 
  ungroup()

HC20EC10eq <- EC10eq %>% 
  group_by(CAS.Number) %>% 
  summarise(
    n_data_recs = sum(n_samples),
    n_species = n(),
    sum_extrapolations = sum(n_extrapolations),
    n_Taxonomy.Group = n_distinct(Taxonomy.Group[!is.nan(sp_mean)]),
    mean_EC10eq = mean(sp_mean, na.rm = T),
    Sd_EC10eq = sd(sp_mean, na.rm = T),
    #least_squares = least_square_min_res_function(sd_samples, sp_mean, n_species, cum_norm_dist)  # Måste kanske köra en for-loop??????
  )  %>% 
  mutate(HC20EC10eq = mean_EC10eq +(-0.842 * Sd_EC10eq),
         CRF = 0.2/(10^HC20EC10eq)) 

rm(plot_CAS)
plot_CAS <- HESTIA_BASE_EnviroTox_FILL %>% 
  filter(CAS.Number == "63-25-2") %>% 
  mutate(Li = log(EC10eq),
         mean_Li = mean(log(EC10eq))) #%>% 
plot(plot_CAS$Li, mean(plot_CAS$Li), type = "p")


# testing the ranking function
test <- c(4, 2, 8, 5, 4, 12)
cum_norm_dist_function <- function(x, y, z){
  y <- mean(x)
  z <- sd(x)
  cum_norm_dist <- 0.5 + (0.5*((x - y)/(z*sqrt(2))))
  return(cum_norm_dist)
}
plot(cum_norm_dist_function)

(rank(test, na.last = NA)-0.5)/3

head(HESTIA_BASE_EnviroTox_FILL)


df_out <- HESTIA_BASE_EnviroTox_FILL %>% 
  mutate(Li = log10(as.numeric(EC10eq)), na.rm = T) %>% 
  group_by(CAS.Number, Taxonomy.Group, Species) %>% 
  # Log10-transforming effect concentrations, then calculating the arithmetic mean per species' log10 effect concentration
  summarise(
    n_samples = n(),
    n_extrapolations = sum(No.Extrapolations),
    # sp_mean = Per-species average
    sp_mean = mean(Li),
    # sd_samples = per-species st.dev, this uses denominator n - 1.
    sd_samples = sd(Li, na.rm = T),
    # Errori = per-species standard error (CI?)
    Errori = sd_samples/(sqrt(n_samples)),
    #cum_norm_dist = cum_norm_dist_function(x = Li, M = sp_mean, S = sd_samples) # Måste kanske köra en for-loop???
  ) %>% 
  group_by(CAS.Number) %>% 
  summarise(
    # n records per substance
    n_data_recs = sum(n_samples),
    # n species per substance
    n_species = n(),
    # n extrapolated ECx values
    sum_extrapolations = sum(n_extrapolations),
    #n taxonomic groups per substance
    n_Taxonomy.Group = n_distinct(Taxonomy.Group[!is.nan(sp_mean)]),
    # the µ per substance
    mean_EC10eq = mean(sp_mean, na.rm = T),
    # the sigma, standard deviation, per substance
    Sd_EC10eq = sd(sp_mean, na.rm = T),
    #least_squares = least_square_min_res_function(Ei, Mi, n_species, cum_norm_dist)  # Måste kanske köra en for-loop??????
  )  %>% 
  mutate(HC20EC10eq = mean_EC10eq +(-0.842 * Sd_EC10eq),
         CRF = (20/100)/(10^HC20EC10eq)) 















meanEC10eq = mean(EC10eq),
medianEC10eq = median(EC10eq),
n_data_points = n(), # n data points in each species group
n_acute = sum(!is.na(EC10eq_Acute)), # makes a logical operation to evaluate if NA or not NA, sum(all that is true)
n_chronic = sum(!is.na(EC10eq_Chronic)), # makes a logical operation to evaluate if NA or not NA, sum(all that is true)
log10EC10eq = mean(log10(EC10eq)), # Arithmetic mean of log10-transformed EC10eq-data based on acute AND Chronic effect data per species group

log10EC10eq_acute = mean(log10(EC10eq_Acute)), # Arithmetic mean of log10-transformed EC10eq-data based on acute effect data
log10EC10eq_chronic = mean(log10(EC10eq_Chronic)), # Arithmetic mean of log10-transformed EC10eq-data based on acute effect data
n_extrapolations = sum(No.Extrapolations),
sd_all = sd(EC10eq, na.rm = TRUE), # Standard deviation of all data points per species group
sd_acute = sd(EC10eq_Acute, na.rm = TRUE), # Standard deviation of all acute data points per species group
sd_chronic = sd(EC10eq_Chronic, na.rm = TRUE), # Standard deviation of all chronic data points per species group
CoV_all = goeveg::cv(EC10eq, na.rm = TRUE), # Coefficient of variation for all data points per species group  
CoV_acute = goeveg::cv(EC10eq_Acute, na.rm = TRUE), # Coefficient of variation for all acute data points per species group
CoV_chronic = goeveg::cv(EC10eq_Chronic, na.rm = TRUE) # Coefficient of variation for all chronic data points per species group
#EC10eq_range = case_when(n_data_points == 1 | max(EC10eq) - min(EC10eq) == 0 ~ as.numeric(NA),
#                         TRUE ~ max(EC10eq) - min(EC10eq)),
#EC10eq_acute_range = case_when(n_acute %in% c(0, 1) | max(EC10eq_acute) - min(EC10eq_acute) == 0 ~ as.numeric(NA),
#                         TRUE ~ max(EC10eq_acute) - min(EC10eq_acute)),
#EC10eq_chronic_range = case_when(n_chronic %in% c(0, 1) | max(EC10eq_chronic) - min(EC10eq_chronic) ==0 ~ as.numeric(NA),
#                         TRUE ~ max(EC10eq_chronic) - min(EC10eq_chronic))
#time_hours_min = min(Time.Hours, na.rm = TRUE), # Minimum test time data is based on
#time_hours_max = max(Time.Hours, na.rm = TRUE) # Maximum test time data is based on
) %>% 
  mutate(across(starts_with("log"), ~ case_when(
    is.nan(.x) ~ as.numeric(NA), TRUE ~ .x))
  )

HESTIA_HC20_dataset <- HESTIA_envirotox_EC10eq_first_avg %>% 
  # Summarizing a second time to get averages across all taxa as well as counts of n = acute data points, n = chronic data points, n = USEtox taxa per acute/chronic factor
  group_by(CAS.Number) %>% 
  summarise(
    mean_sp_EC10eq = mean(meanEC10eq),
    median_sp_EC10eq = median(medianEC10eq),
    meanlog10EC10eq = mean(log10EC10eq), # arithmetic mean of (chronic) logEC10eq for all species per chemical
    meanlog10EC10eq_chronic = mean(log10EC10eq_chronic, na.rm = TRUE), # arithmetic mean of (chronic) logEC10eq for all species per chemical
    sd_log = sd(log10EC10eq, na.rm = TRUE), # Standard deviation of all log-transformed data for each substance
    sd_log_chronic = sd(log10EC10eq_chronic, na.rm = TRUE), # Standard deviation of all log-transformed chronic data for each substance
    logHC20EC10eq_all = meanlog10EC10eq + (sd_log * -0.842), # Calculating the log HC20EC10eq which uses mean, sd and z0.2; the inverse of the standard normal
    logHC20EC10eq_chronic = meanlog10EC10eq_chronic + (sd_log_chronic * -0.842), # Calculating the log HC20EC10eq which uses mean, sd and z0.2; the inverse of the standard normal distribution at the 0.2 probability level (-0.842) (Owsianiak et al., 2022)
    logHC20EC10eq_best = case_when(
      sum(!is.na(log10EC10eq_chronic)) <5 & sum(!is.na(log10EC10eq)) >=5 ~ meanlog10EC10eq + (sd_log * -0.842),
      TRUE ~ logHC20EC10eq_chronic),# Using Owsianiak's principle to calculate EC10eq from chronic data, unless fewer than 5 species-specific chronic EC10eq data points exist. Then calculate SSDss using both acute & chronic data
    # median_range_all_EC10eq = median(EC10eq_range),
    # Concentration-response slope factor calculation represents an incremental change in the potentially affected fraction of species due to an incremental exposure to the bioavailable concentration of a chemical at the HC20 level  
    CRF_all = 0.2/10^logHC20EC10eq_all, # Concentration-response slope factor calculation based on all available data
    
    CRF_best = 0.2/10^logHC20EC10eq_best,# Concentration-response slope factor calculation for >5 chronic datapoints or if not sufficient data, acute data is used. According to Owsianiak et al. 2022 
    CoV_log_all = goeveg::cv(log10EC10eq, na.rm = TRUE), # Coefficient of variation for all data points per species group  
    CoV_log_acute = goeveg::cv(log10EC10eq_acute, na.rm = TRUE), # Coefficient of variation for all data points per species group  
    CoV_log_chronic = goeveg::cv(log10EC10eq_chronic, na.rm = TRUE), # Coefficient of variation for all data points per species group  
    n_data_tot = sum(n_data_points), # Number of data points available for this substance
    sum_extrapolations = sum(n_extrapolations),
    n_acute_sum = sum(n_acute), # Number of acute data points available for this substance
    n_chronic_sum = sum(n_chronic), # Number of chronic data points available for this substance
    n_best = case_when(
      sum(!is.na(log10EC10eq_chronic)) <5 & sum(!is.na(log10EC10eq)) >=5 ~ (sum(n_chronic) + sum(n_acute)),
      TRUE ~ n_chronic_sum), # Number of data points available for this substance when considering only >5 species' chronic data or, if not sufficient, all data points.
    n_species_all = sum(!is.na(log10EC10eq)),
    n_species_chronic = sum(!is.na(log10EC10eq_chronic)),
    n_CoV_all = sum(!is.na(CoV_all)), # makes a logical operation to evaluate if NA or not NA, sum(all that is true)
    n_CoV_chronic = sum(!is.na(CoV_chronic)), # makes a logical operation to evaluate if NA or not NA, sum(all that is true)
    CoV_sp_mean_EC10eq = goeveg::cv(0.2/10^meanEC10eq, na.rm = TRUE), # CoV of species mean values
    n_Taxonomy.Group_all = n_distinct(Taxonomy.Group[!is.nan(log10EC10eq)]),
    n_Taxonomy.Group_chronic = n_distinct(Taxonomy.Group[!is.nan(log10EC10eq_chronic)]),
    CRF_bad = case_when(
      n_species_all == 1 ~ 0.2/10^(meanlog10EC10eq + (0.8710236 * -0.842)), # Median sd_log is the median sd of the best and intermediate quality CRFs sd of ONLY chronic data: a generic value of 0.8710236. THIS IS AN AD HOC OPERATION!!! using `summarize(median = median(sd_log))`
      TRUE ~ as.numeric(NA)),
    QSe = log(n_species_all)*log(n_Taxonomy.Group_all)*(1/(1+sum_extrapolations)^0.1),
    QSe_desc = case_when(
      QSe > 1.77 ~ "High Quality",
      QSe > 1.48 & QSe <= 1.77 ~ "Intermediate Quality",
      QSe >= 0 & QSe <= 1.48 ~ "Low Quality")
  ) %>%
  ungroup() %>% 
  mutate(across(starts_with(c("mean", "log", "CRF")), ~ case_when(
    is.nan(.x) ~ as.numeric(NA), # making sure empty cells are specified as "NA".
    TRUE ~ .x))) %>% 
  # Joining the HC50_EC50 data to have "backwards compatibility" with the database
  left_join(
    x = .,
    y = HESTIA_HC50_calc %>% 
      select(
        CAS.Number, avlog10_HC50_EC50, HC50_AcuteChronic, 
        n_data_points_HC50, n_species_acute_HC50, 
        n_species_Chronic_HC50, USEtox_HC50_recommend),
    by = "CAS.Number")

# Description of the generic SD
#HESTIA_HC20_dataset %>% 
#    filter(QSe_desc %in% c("High Quality", "Intermediate Quality")) %>% 
#  summarize(median = median(sd_log_chronic, na.rm=T))
