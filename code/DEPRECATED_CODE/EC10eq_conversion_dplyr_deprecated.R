# Deprecated `dplyr` approach for EC10eq extrapolation factors.
  # This code is messy and long, but very efficient since it's leaning on dplyrs' mutate & case_when functions
library(dplyr)


Q_dat <- HESTIA_HC20_DB_duplicated_filter %>% 
  mutate(DB = "HESTIA",
         version = "HESTIA 1.3") %>%
  
  mutate(
    # Applying the extrapolation factors available in Aurisano et al., 2019, to each Endpoint_conv to both account for acute/chronic conversion as well as endpoint conversion into EC10eq
    ##################
    EC10eq = as.numeric(case_when(
      Endpoint_conv == "EC50" & AcuteChronic == "Chronic" & Taxonomy.Group == "Fish" ~ Value.mg_l/1.55,
      Endpoint_conv == "EC50" & AcuteChronic == "Chronic" & Taxonomy.Group == "Crustacean" ~ Value.mg_l/1.94,
      Endpoint_conv == "EC50" & AcuteChronic == "Chronic" & Taxonomy.Group == "Algae" ~ Value.mg_l/2.24,
      Endpoint_conv == "NOEC" & AcuteChronic == "Chronic" & Taxonomy.Group == "Crustacean" ~ Value.mg_l/0.95,
      Endpoint_conv == "NOEC" & AcuteChronic == "Chronic" & Taxonomy.Group == "Algae" ~ Value.mg_l/0.44,
      Endpoint_conv == "NOEC" & AcuteChronic == "Chronic" & !Taxonomy.Group %in% c("Crustacean", "Algae") ~ Value.mg_l/0.6, # Generalized conversion factor To EC10eq-chronic from NOECchronic
      Endpoint_conv == "EC50" & AcuteChronic == "Chronic" & !Taxonomy.Group %in% c("Fish","Crustacean", "Algae") ~ Value.mg_l/2, # Generalized conversion factor To EC10eq-chronic from EC50chronic
      Endpoint_conv == "EC10" & AcuteChronic == "Chronic" ~ Value.mg_l,
      Endpoint_conv == "EC50" & AcuteChronic == "Acute" & Taxonomy.Group == "Fish" ~ Value.mg_l/7.44,
      Endpoint_conv == "EC50" & AcuteChronic == "Acute" & Taxonomy.Group == "Crustacean" ~ Value.mg_l/3.38,
      Endpoint_conv == "NOEC" & AcuteChronic == "Acute" & Taxonomy.Group == "Fish" ~ Value.mg_l/3.97,
      Endpoint_conv == "NOEC" & AcuteChronic == "Acute" & Taxonomy.Group == "Crustacean" ~ Value.mg_l/1.55,
      Endpoint_conv == "NOEC" & AcuteChronic == "Acute" & !Taxonomy.Group %in% c("Crustacean", "Fish") ~ Value.mg_l/1.8, # Generalized conversion factor To EC10eq-chronic from NOECacute 
      Endpoint_conv == "EC50" & AcuteChronic == "Acute" & !Taxonomy.Group %in% c("Crustacean", "Fish") ~ Value.mg_l/4, # Generalized conversion factor To EC10eq-chronic from EC50acute
      Endpoint_conv == "EC10" & AcuteChronic == "Acute"~ Value.mg_l, 
      TRUE ~ as.numeric(NA))),
    No.Extrapolations = as.numeric(case_when(Endpoint_conv == "EC10" ~ 0,
                                             TRUE ~ 1)),
    EC10eq_chronic = as.numeric(case_when(
      Endpoint_conv == "EC50" & AcuteChronic == "Chronic" & Taxonomy.Group == "Fish" ~ Value.mg_l/1.55,
      Endpoint_conv == "EC50" & AcuteChronic == "Chronic" & Taxonomy.Group == "Crustacean" ~ Value.mg_l/1.94,
      Endpoint_conv == "EC50" & AcuteChronic == "Chronic" & Taxonomy.Group == "Algae" ~ Value.mg_l/2.24,
      Endpoint_conv == "NOEC" & AcuteChronic == "Chronic" & Taxonomy.Group == "Crustacean" ~ Value.mg_l/0.95,
      Endpoint_conv == "NOEC" & AcuteChronic == "Chronic" & Taxonomy.Group == "Algae" ~ Value.mg_l/0.44,
      Endpoint_conv == "NOEC" & AcuteChronic == "Chronic" & !Taxonomy.Group %in% c("Crustacean", "Algae") ~ Value.mg_l/0.6, # Generalized conversion factor To EC10eq-chronic from NOECchronic
      Endpoint_conv == "EC50" & AcuteChronic == "Chronic" & !Taxonomy.Group %in% c("Fish","Crustacean", "Algae") ~ Value.mg_l/2, # Generalized conversion factor To EC10eq-chronic from EC50chronic
      Endpoint_conv == "EC10" & AcuteChronic == "Chronic" ~ Value.mg_l,
      TRUE ~ as.numeric(NA))),# <- all the rest, very many data points. 
    EC10eq_acute = as.numeric(case_when( 
      Endpoint_conv == "EC50" & AcuteChronic == "Acute" & Taxonomy.Group == "Fish" ~ Value.mg_l/7.44,
      Endpoint_conv == "EC50" & AcuteChronic == "Acute" & Taxonomy.Group == "Crustacean" ~ Value.mg_l/3.38,
      Endpoint_conv == "NOEC" & AcuteChronic == "Acute" & Taxonomy.Group == "Fish" ~ Value.mg_l/3.97,
      Endpoint_conv == "NOEC" & AcuteChronic == "Acute" & Taxonomy.Group == "Crustacean" ~ Value.mg_l/1.55,
      Endpoint_conv == "NOEC" & AcuteChronic == "Acute" & !Taxonomy.Group %in% c("Crustacean", "Fish") ~ Value.mg_l/1.8, # Generalized conversion factor To EC10eq-chronic from NOECacute 
      Endpoint_conv == "EC50" & AcuteChronic == "Acute" & !Taxonomy.Group %in% c("Crustacean", "Fish") ~ Value.mg_l/4, # Generalized conversion factor To EC10eq-chronic from EC50acute
      Endpoint_conv == "EC10" & AcuteChronic == "Acute" ~ Value.mg_l,
      TRUE ~ as.numeric(NA))),# <- all the rest, very many data points. 
    
    EC10eq_low = as.numeric(case_when(
      Endpoint_conv == "EC50" & AcuteChronic == "Chronic" & Taxonomy.Group == "Fish" ~ Value.mg_l/0.67, # low-CI
      Endpoint_conv == "EC50" & AcuteChronic == "Chronic" & Taxonomy.Group == "Crustacean" ~ Value.mg_l/1.56, # low-CI
      Endpoint_conv == "EC50" & AcuteChronic == "Chronic" & Taxonomy.Group == "Algae" ~ Value.mg_l/1.9, # low-CI
      Endpoint_conv == "NOEC" & AcuteChronic == "Chronic" & Taxonomy.Group == "Crustacean" ~ Value.mg_l/0.77, # low-CI
      Endpoint_conv == "NOEC" & AcuteChronic == "Chronic" & Taxonomy.Group == "Algae" ~ Value.mg_l/0.39, # low-CI
      Endpoint_conv == "EC50" & AcuteChronic == "Acute" & Taxonomy.Group == "Fish" ~ Value.mg_l/2.92, # low-CI
      Endpoint_conv == "EC50" & AcuteChronic == "Acute" & Taxonomy.Group == "Crustacean" ~ Value.mg_l/2.14, # low-CI
      Endpoint_conv == "NOEC" & AcuteChronic == "Acute" & Taxonomy.Group == "Fish" ~ Value.mg_l/0.9, # low-CI
      Endpoint_conv == "NOEC" & AcuteChronic == "Acute" & Taxonomy.Group == "Crustacean" ~ Value.mg_l/0.91, # low-CI
      Endpoint_conv == "NOEC" & AcuteChronic == "Acute" & !Taxonomy.Group %in% c("Crustacean", "Fish") ~ Value.mg_l/1.0, # low-CI
      Endpoint_conv == "NOEC" & AcuteChronic == "Chronic" & !Taxonomy.Group %in% c("Crustacean", "Algae") ~ Value.mg_l/0.5, # low-CI
      Endpoint_conv == "EC50" & AcuteChronic == "Acute" & !Taxonomy.Group %in% c("Crustacean", "Fish") ~ Value.mg_l/2.6, # low-CI
      Endpoint_conv == "EC50" & AcuteChronic == "Chronic" & !Taxonomy.Group %in% c("Fish","Crustacean", "Algae") ~ Value.mg_l/1.8, # low-CI
      Endpoint_conv == "EC10" ~ Value.mg_l,
      TRUE ~ as.numeric(NA))),
    
    EC10eq_high = as.numeric(case_when(
      Endpoint_conv == "EC50" & AcuteChronic == "Chronic" & Taxonomy.Group == "Fish" ~ Value.mg_l/3.66, # High-CI
      Endpoint_conv == "EC50" & AcuteChronic == "Chronic" & Taxonomy.Group == "Crustacean" ~ Value.mg_l/2.41, # High-CI
      Endpoint_conv == "EC50" & AcuteChronic == "Chronic" & Taxonomy.Group == "Algae" ~ Value.mg_l/2.65, # High-CI
      Endpoint_conv == "NOEC" & AcuteChronic == "Chronic" & Taxonomy.Group == "Crustacean" ~ Value.mg_l/1.16, # High-CI
      Endpoint_conv == "NOEC" & AcuteChronic == "Chronic" & Taxonomy.Group == "Algae" ~ Value.mg_l/0.49, # High-CI
      Endpoint_conv == "EC50" & AcuteChronic == "Acute" & Taxonomy.Group == "Fish" ~ Value.mg_l/18.95, # Notice this incredibly large spread in the data. 3.97 (0.90–17.39)
      Endpoint_conv == "EC50" & AcuteChronic == "Acute" & Taxonomy.Group == "Crustacean" ~ Value.mg_l/5.34, # High-CI
      Endpoint_conv == "NOEC" & AcuteChronic == "Acute" & Taxonomy.Group == "Fish" ~ Value.mg_l/17.39, # Notice this incredibly large spread in the data. 3.97 (0.90–17.39)
      Endpoint_conv == "NOEC" & AcuteChronic == "Acute" & Taxonomy.Group == "Crustacean" ~ Value.mg_l/2.64, # High-CI
      Endpoint_conv == "NOEC" & AcuteChronic == "Acute" & !Taxonomy.Group %in% c("Crustacean", "Fish") ~ Value.mg_l/2.7, # High-CI
      Endpoint_conv == "NOEC" & AcuteChronic == "Chronic" & !Taxonomy.Group %in% c("Crustacean", "Algae") ~ Value.mg_l/0.7, # High-CI
      Endpoint_conv == "EC50" & AcuteChronic == "Acute" & !Taxonomy.Group %in% c("Crustacean", "Fish") ~ Value.mg_l/6.1, # High-CI
      Endpoint_conv == "EC50" & AcuteChronic == "Chronic" & !Taxonomy.Group %in% c("Fish","Crustacean", "Algae") ~ Value.mg_l/2.5, # High-CI 
      Endpoint_conv == "EC10" ~ Value.mg_l,
      TRUE ~ as.numeric(NA)))
  ) %>% 
  # Making sure that the Taxonomy.Group is a factor      
  mutate(Taxonomy.Group = as.factor(Taxonomy.Group))
```