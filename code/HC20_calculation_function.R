
# function to calculate HC20EC10eq from EC10eq, as well as the CI of the EC10eq values
HCx_calculator <- function(df, CAS = CAS.Number, Tax = Taxonomy.Group, sp = Species, EC10 = EC10eq, HCx = 20){
#library(dplyr)
  df_out <- df %>% 
    mutate(Li = log10({{EC10}})) %>% 
    group_by({{CAS}}, {{Tax}}, {{sp}}) %>% 
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
    group_by({{CAS}}) %>% 
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
           CRF = ({{HCx}}/100)/(10^HC20EC10eq)) 
  colnames(df_out)[colnames(df_out) == "CRF"] = paste0("HC", deparse(substitute(HCx)))
  return(df_out)
}
