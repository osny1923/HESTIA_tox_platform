library(dplyr)


# Define the dataset and calculate Species-specific averages:
# MC_df <- HESTIA_BASE_EnviroTox_FILL %>% 
#   filter(CAS.Number %in% H_CAS[3]) %>% 
#   select(CAS.Number, Species, Taxonomy.Group, EC10eq) %>% 
#   mutate(ln_EC10eq = log(EC10eq),
#          log_EC10eq = log10(EC10eq)) %>% 
#   group_by(CAS.Number,Taxonomy.Group, Species) %>% 
#   summarize(n = n(),
#             mean_log_EC10eq = mean(log_EC10eq),
#             mean_ln_EC10eq = mean(ln_EC10eq),
#             sd_ln_EC10eq = sd(ln_EC10eq),
#             CV_ln_EC10eq = sd_ln_EC10eq/mean_ln_EC10eq) %>% 
#   ungroup()

Uncertain_norm_dataset_func <- function(x) {
  # Make a CAS.Number vector to loop over
  H_CAS <- x %>% 
    distinct(CAS.Number) %>% pull()
  
  for (i in 1:length(H_CAS)) {
    MC_df[i,] <- x %>% 
      filter(CAS.Number %in% H_CAS[i]) %>% 
      select(CAS.Number, Species, Taxonomy.Group, EC10eq) %>% 
      # mutate(ln_EC10eq = log(EC10eq),
      #        #log_EC10eq = log10(EC10eq)
      #        ) %>% 
      group_by(CAS.Number) %>% 
      summarize(n = n(),
                #mean_log_EC10eq = mean(log_EC10eq),
                Q0.05 = quantile(EC10eq, probs = 0.05),
                Q0.95 = quantile(EC10eq, probs = 0.95),
                mean_EC10eq = mean(EC10eq),
                median_EC10eq = median(EC10eq),
                sd_EC10eq = sd(EC10eq),
                CV_EC10eq = sd_EC10eq/mean_EC10eq) %>% 
      ungroup()
    
  }
  return(MC_df)  
}

system.time(uncert_data <- Uncertain_norm_dataset_func(HESTIA_BASE_EnviroTox_FILL))


# Calculating HC20EC10eq-value
HC20EC10eq_dat <- MC_df %>% 
  group_by(CAS.Number) %>% 
  summarize(
    n = sum(n),
    all_mean = mean(mean_ln_EC10eq, na.rm = T),
    all_sd = sd(sd_ln_EC10eq, na.rm = T)
    ) %>% 
  ungroup() %>% 
  mutate(HC10EC10eq = all_mean + (-0.842 * all_sd),
         CRF_HC20 = 0.2/(10^HC10EC10eq))

## Using the ecdf function (empirical cumulative distribution function) to derive the 20th percentile (HC20EC10eq) point of 
ggplot(MC_df, aes())
plot(ecdf(MC_df$mean_ln_EC10eq))
line(ecdf(MC_df$mean_ln_EC10eq))
     
g <- ecdf(MC_df$mean_ln_EC10eq)
plot(g)
g <- ecdf(MC_df$mean_log_EC10eq)
quantile(g, probs = 0.2, type = 2)


# Calculate mu and sigma for the specified CAS.number level
# i = 1 
# for (i in 1:length(H_CAS)) {
cas_number <- H_CAS[3]  # Specify the CAS.number level
subset_df <- subset(HESTIA_BASE_EnviroTox_FILL, CAS.Number == cas_number)

species_stats <- aggregate(EC10eq ~ Species, subset_df, function(x) {
  if (length(x) >= 2) c(mu = mean(x), sigma = sd(x))
  else c(mu = NA, sigma = NA)
})

# Perform Monte Carlo simulations
set.seed(42)  # Set a seed for reproducibility
n_samples <- 1000
simulated_data <- do.call(rbind, lapply(species_stats$Species, function(species) {
  mu <- species_stats[species_stats$Species == species, "value"]["mu"]
  sigma <- species_stats[species_stats$Species == species, "value"]["sigma"]
  
  if (!is.na(mu) && !is.na(sigma)) {
    replicate(n_samples, rnorm(1, mean = mu, sd = sigma))
  } else {
    rep(NA, n_samples)
  }
}))

#}
# Print the simulated data
print(simulated_data)
