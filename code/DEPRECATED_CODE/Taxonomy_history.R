# Binding together as a final taxonomy dataset
Final_taxonomy <- rbind(
  NCBI_Taxonomy %>% filter(!is.na(Source)),
  GBIF_Taxonomy %>% filter(!is.na(Source)),
  no_match_GBIF_fixed) %>%
  # fixing a few strange class-descriptions
  mutate(Phylum = case_when(
    Class == "Dinophyceae" ~ "Myzozoa", # algae
    Class == "Cryptophyceae" ~ "Cryptophyta", # algae
    Class == "Chrysophyceae" ~ "Ochrophyta",
    Class == "Choanoflagellata" ~ "Choanozoa",
    Class == "Raphidophyceae" ~ "Ochrophyta",
    Class == "Synurophyceae" ~ "Ochrophyta",
    Class == "Eustigmatophyceae" ~ "Ochrophyta",
    Class == "Phaeophyceae" ~ "Ochrophyta",
    TRUE ~ Phylum),
  ) %>%
  mutate(Taxonomy.Group = case_when(
    Phylum == "Chordata" & Class %in% c("Amphibia", "Amphibian") ~ "Amphibian",
    Phylum == "Chordata" & Class %in% c("Actinopterygii", "Actinopteri") ~ "Fish",
    Phylum == "Arthropoda" & Class %in% c("Insecta") ~ "Insect",
    Phylum == "Arthropoda" & !Class %in% c("Insecta", "Arachnida") ~ "Crustacean",
    Phylum == "Mollusca" ~ "Mollusca",
    Phylum == "Streptophyta" ~ "Plant",
    Phylum == "Rotifera" ~ "Rotifera",
    Phylum == "Annelida" ~ "Annellidae",
    Phylum %in% c(
      "Bacillariophyta",
      "Cyanobacteria",
      "Chlorophyta",
      "Ochrophyta",
      "Rhdophyta",
      "Miozoa",
      "Cryptophyta",
      "Haptophyta",
      "Charophyta",
      "Euglenophyta") ~ "Algae",
    is.na(Phylum) ~ as.character(NA),
    TRUE ~ "Others"
  )  )

