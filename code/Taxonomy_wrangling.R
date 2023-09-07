library(tidyverse)
library(taxize)


Species_list <- read.csv("../data/RAW_DB/HESTIA_HC20_DB_raw_toxdata.csv") %>% 
  distinct(Test.organisms..species.) 

#----------Binomial species definitions query at NCBI----------------#

# Separating the Species_list into columns with binomial species names, and "non-binomial" annotations    
main_taxa_df <- Species_list %>% 
  separate(col = Test.organisms..species., into = "species", sep = " sp\\.| ssp\\.| var\\.| x ", remove = F, convert = F) %>% 
  mutate(non_binom = case_when(!grepl("\\s", species) ~ species),
         species = case_when(!grepl("\\s", species) ~ as.character(NA), TRUE ~ species))

# All binomial species
main_species_df <- main_taxa_df %>% filter(!is.na(species))

# Checking for names status and their correct taxonomic name
correct_names <- tol_resolve(main_species_df$species)

# Joining the correct names to main species list
main_species_df_corr <- left_join(
  main_species_df %>% mutate(species = tolower(species)),
  correct_names %>% rename(species = search_string),
  by = "species")

# Removing everything following a parenthesis & selecting correctly named species
Species_taxonomy_search_df <- main_species_df_corr %>%
  mutate(unique_name = sub("\\(.*", "", unique_name)) %>%
  filter(!is.na(unique_name))

# Removing punctuation & selecting species missing correct names for later processing
Species_taxonomy_leftover_names <- main_species_df_corr %>%
  mutate(unique_name = sub("\\(.*", "", unique_name)) %>%
  filter(is.na(unique_name))

#------------Query the NCBI database--------------------#

# Run the NCBI query and hang on the results to the main_species_df
Species_taxonomy_search_df$uid <- get_uid(Species_taxonomy_search_df$unique_name, rank_query = "species")

# Querying NCBI for the main_taxa_df species list, and wrangling this output into a df.
NCBI_Species_query <- rbind(
  classification(
    Species_taxonomy_search_df$uid,
    db = 'ncbi',
    rank_query = "species", # <- for the df with species defined!
    rows = 1,
    verbose = TRUE)
    ) %>%
  select(-id) %>%
  distinct(rank, query, .keep_all = T) %>%
  pivot_wider(names_from = rank, values_from = name) %>%
  select(query, phylum, class, order, genus, species) %>%
  # The NCBI Query returns the "sister"-class "Actinopteri" for bony fishes, 
  # not the proper class "Actinopterygii"!! Manually changing this.
  mutate(class = case_when(class == "Actinopteri" ~ "Actinopterygii", TRUE ~ class))
# Writing the output to file.
write.csv(NCBI_Species_query, "../data/Taxonomy/NCBI_Species_query.csv", row.names = F)

# Creating a subset with correct species annotations
Species_list_NCBI <- left_join(
x = Species_taxonomy_search_df %>%
  mutate(uid = as.character(uid)) %>%
  rename(query = uid),
y = NCBI_Species_query %>%
  filter(!is.na(species)),
by = "query") %>%
  mutate(source = case_when(!is.na(query) ~"NCBI"))

# Export the NCBI-query part to a .csv
write.csv(Species_list_NCBI, "../data/Taxonomy/Species_list_NCBI.csv", row.names = F)

#-------------Missing species from NCBI---------------------------#

## Defining taxa that are missing information  
missing_taxa <- main_taxa_df %>% 
  left_join(
    x = ., 
    y = Species_list_NCBI %>% select(-c("species.x", "non_binom")), 
    by = "Test.organisms..species.") %>% 
  filter(is.na(source)) %>% 
  mutate(unique_name = gsub("\\[|\\]", "", unique_name),
         species = case_when(!is.na(unique_name) & species != unique_name ~ unique_name, TRUE ~ species))

# Separating the missing taxonomy into three datasets: 
  
  # species not described in the NCBI query, 
missing_species <- missing_taxa[,1:9] %>% 
  filter(!is.na(species) & !is.na(unique_name)) %>% 
  mutate(unique_name = sub("Chlorella fusca", "Desmodesmus abundans", unique_name))

  # Common names  
missing_common_names <- missing_taxa[,1:4] %>% 
  filter(!is.na(species) & is.na(unique_name)) %>% 
  mutate(species = sub("\\(.*", "", species))

  # non-binomial definitions
missing_non_binom <- missing_taxa[,c(1,3)] %>% 
  filter(!is.na(non_binom))

#----------------------GBIF query--------------------#

# Querying the GBIF database for additional taxonomic information.
# Code for the query
missing_species$query <- get_gbifid(missing_species$unique_name)

# saving the slow query to file.
write.csv(missing_species, "../data/Taxonomy/missing_species_GBIF.csv", row.names = F)

# With GBIF_id, look up taxonomic information
GBIF_non_binom_query <- rbind(
  classification(
   missing_species %>% pull(query),
    db = 'gbif',
    rank_query = "species", # <- for the targeted taxonomic level!
    rows = 1,
    verbose = TRUE)
    ) %>%
  select(-id) %>%
  distinct(rank, query, .keep_all = T) %>%
  pivot_wider(names_from = rank, values_from = name) %>%
  mutate(query = as.character(query)) %>%
  select(query, phylum, class, order, genus, species)

# Saving slow query to file
write.csv(GBIF_non_binom_query, "../data/Taxonomy/GBIF_non_binom_query.csv", row.names = F)

# Compiling a subset with correct taxonomic identification
GBIF_non_binom_Taxonomy <- left_join(
  x = missing_species %>% 
    select(-species) %>% 
    mutate(query = as.integer(query)) %>%
    rename(species = unique_name),
  y = GBIF_non_binom_query, 
  by = "query") %>% 
  mutate(source = case_when(!is.na(query) ~ "GBIF")) %>% 
  select(-species.x) %>% 
  rename(species = species.y)

#-------------------------Common names-----------------------------#

# Using the object created above to look for identiiers using common names
missing_common_names$uid <- get_uid(missing_common_names$species, modifier = "Common Name")
# write slow query to file
write.csv(missing_common_names, "../data/Taxonomy/missing_common_names_UID.csv", row.names = F)

NCBI_missing_common_query <- rbind(
  classification(
    missing_common_names %>%
      filter(!is.na(uid)) %>%
               pull(uid),
    db = 'ncbi',
    rank_query = "species", # <- for the targeted taxonomic level!
    rows = 1,
    verbose = TRUE)
    ) %>%
  select(-id) %>%
  distinct(rank, query, .keep_all = T) %>%
  pivot_wider(names_from = rank, values_from = name) %>%
  select(query, phylum, class, order, genus, species)

# Write slow query to file
write.csv(NCBI_missing_common_query, "../data/Taxonomy/NCBI_missing_common_query.csv", row.names = F)

# Create a subset with correct taxonomic information based on common names
common_names_taxonomy <- left_join(
  x = missing_common_names %>%
    mutate(uid = as.integer(uid)) %>% 
    rename(query = uid, 
           common_name = species),
  y = NCBI_missing_common_query,
  by = "query") %>% 
  mutate(source = case_when(!is.na(query) ~ "NCBI"))

#----------------non-binomial names "hotfix"-------------------------#
NCBI_Taxonomy_unique <- Species_list_NCBI %>% select(phylum:genus) %>% distinct(genus, .keep_all = T)

non_binom_taxonomy <- missing_non_binom %>% 
  mutate(higher_taxonomy = case_when(
    !grepl(" sp\\.", Test.organisms..species.) ~ non_binom, 
    TRUE ~ as.character(NA)),
    non_binom = case_when(
      !is.na(higher_taxonomy) ~ as.character(NA), 
      TRUE ~ non_binom)) %>% 
  rename(genus = non_binom) 

NCBI_non_binom_query <- rbind(
  classification(
    get_uid(non_binom_taxonomy %>%
              filter(!is.na(genus)) %>%
              pull(genus)),
    db = 'ncbi',
    rank_query = "genus", # <- for the targeted taxonomic level!
    rows = 1,
    verbose = TRUE)
    ) %>%
  select(-id) %>%
  distinct(rank, query, .keep_all = T) %>%
  pivot_wider(names_from = rank, values_from = name) %>%
  select(query, phylum, class, order, genus)
# Write slow query to file
write.csv(NCBI_non_binom_query, "../data/Taxonomy/NCBI_non_binom_query.csv", row.names = F)

NCBI_non_binom_genus_query <- left_join(
  x = non_binom_taxonomy, 
  y = NCBI_non_binom_query,
  by = "genus") %>% 
  mutate(species = as.character(NA),
         source = case_when(!is.na(query) ~"NCBI"))

#-----------------Merging taxonomic information----------------#

# Define relevant columns to include in the merge
selection_list <- c("Test.organisms..species.", "phylum", "class", "order", "genus", "species", "source", "query")

# merge the various queries into one comprehensive object
Species_taxonomy <- rbind(
  Species_list_NCBI %>% 
    filter(!is.na(query)) %>% 
    rename(species = species.y) %>% 
    select(all_of(selection_list)),
  GBIF_non_binom_Taxonomy %>% 
    select(all_of(selection_list)),
  common_names_taxonomy %>% 
    select(all_of(selection_list)),
  NCBI_non_binom_genus_query %>% 
    select(all_of(selection_list))
) %>% 
  distinct(Test.organisms..species., .keep_all = T) %>% 
  # Manual adjustments to phylum for some cryptic algal species
  mutate(phylum = case_when(
    is.na(phylum) & class == "Phaeophyceae" ~ "Sar", 
    is.na(phylum) & class == "Chrysophyceae" ~ "Sar", 
    is.na(phylum) & class == "Eustigmatophyceae" ~ "Sar", 
    is.na(phylum) & class == "Pelagophyceae" ~ "Sar", 
    is.na(phylum) & class == "Raphidophyceae" ~ "Sar", 
    is.na(phylum) & class == "Synurophyceae" ~ "Sar", 
    TRUE ~ phylum)) %>% 
  # Categorization of phyla into Taxonomic groups
  mutate(Taxonomy.Group = case_when( 
    phylum == "Chordata" & class %in% c("Amphibia", "Amphibian") ~ "Amphibian",
    phylum == "Chordata" & class %in% c("Actinopterygii", "Actinopteri") ~ "Fish", 
    phylum == "Arthropoda" & class %in% c("Insecta") ~ "Insect",
    phylum == "Arthropoda" & !class %in% c("Insecta", "Arachnida") ~ "Crustacean",
    phylum == "Mollusca" ~ "Mollusca",
    phylum == "Streptophyta" ~ "Plant",
    phylum == "Rotifera" ~ "Rotifera",
    phylum == "Annelida" ~ "Annellidae",
    phylum %in% c(
      "Sar",
      "Rhodophyta",
      "Bacillariophyta", 
      "Cyanobacteriota", 
      "Chlorophyta", 
      "Ochrophyta",
      "Rhdophyta",
      "Miozoa",
      "Cryptophyta",	
      "Haptophyta",
      "Charophyta",
      "Euglenophyta") ~ "Algae",
    TRUE ~ "Others"),
    Taxonomy.Group = case_when(is.na(query) ~ as.character(NA), TRUE ~ Taxonomy.Group)
  )%>%
  rename_with(str_to_title) %>% 
  rename(Taxonomy.Group = Taxonomy.group,
         Test.organisms..species. = Test.organisms..Species.)

# Write results to file
write.csv(Species_taxonomy, "../data/Taxonomy/Species_taxonomy.csv", row.names = F)




# Exporting a list of updated species names for the Envirotox database names update!
read.csv("../data/Taxonomy/Species_list_NCBI.csv") %>%
  mutate(species.x = stringr::str_to_sentence(species.x)) %>%
  filter(species.x != unique_name) %>%
  select(species.x, unique_name) %>%
  rename(old_name = species.x,
         current_name = unique_name) %>%
  distinct(old_name, .keep_all = T) %>%
  write.csv(., "../data/Taxonomy/names_update.csv", row.names = F)
