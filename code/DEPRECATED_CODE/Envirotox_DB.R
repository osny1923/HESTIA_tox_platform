## ----setup------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(verbose = TRUE)
library(rmarkdown)
library(xlsx)
library(tidyverse)


## ---------------------------------------------------------------------------------------------------------------------------------------
# Managing the Endpoints, included are the ones selected in the file created in the operation above.
# Selections and convertions of endpoints are made according to Aurisano et al., 2019
Endpoints_filter <- read.xlsx("../data/excel_references/Endpoints_selector.xlsx", sheetIndex = 1, header = T)%>% 
                      filter(Group != "N") %>% 
                      select(Endpoint, Group)


## ---- warning=FALSE---------------------------------------------------------------------------------------------------------------------

# Taxonomy wrangling. Adapting the Envirotox Species list to the HESTIA species list.
EnviroTox_taxaDB <- readxl::read_xlsx(path = "../data/Envirotox_2023-01-10/envirotox_DB_20230110081208.xlsx",sheet = "taxonomy", col_names = T) %>% 
              rename(Species = `Latin name`,
                     Class = `Taxonomic class`) %>%
              select(-`Taxonomic superclass`) %>% 
  # Harmonizing some Class annotations and correcting classifications
              mutate(
                Class = gsub("Phaeophycaea", "Phaeophyceae", Class),
                Class = gsub("Cyanophyceae", "Cyanobacteriia", Class),
                Class = gsub("Euglenophyceae", "Euglenoidea", Class),
                Class = gsub("Monogonta", "Monogononta", Class),
                Class = gsub("Cephalaspidomorphi", "Hyperoartia", Class),
                Class = gsub("Maxillipoda", "Maxillopoda", Class),
                Class = gsub("Branchiopoda\\?", "Branchiopoda", Class),
                Class = case_when(Species == "Tubifex tubifex" ~ "Clitellata", 
                                    TRUE ~ Class)) %>% 
              left_join(x = .,
                        y = read.xlsx("../data/Final_Taxonomy_dataset.xlsx", sheetName = "Data1", header = T) %>% 
                            distinct(Class, .keep_all = T) %>% 
                            select(Class, Taxonomy.Group) %>% 
                            filter(!is.na(Class)),
              by = "Class") %>% 
  # Managing the last few NAs manually
    mutate(Taxonomy.Group = case_when(Species == "Lepidodermella squamatum" ~ "Others", 
                                   TRUE ~ Taxonomy.Group),
           Taxonomy.Group = case_when(Species == "Moina dubia" ~ "Crustacean", 
                                   TRUE ~ Taxonomy.Group),
           Taxonomy.Group = case_when(is.na(Taxonomy.Group) & `Taxonomic phylum or division` == "Annelida" ~ "Annellidae", 
                                   TRUE ~ Taxonomy.Group),
           Taxonomy.Group = case_when(is.na(Taxonomy.Group) & `Taxonomic phylum or division` == "Rotifera" ~ "Rotifera", 
                                   TRUE ~ Taxonomy.Group),
           Taxonomy.Group = case_when(is.na(Taxonomy.Group) & `Taxonomic phylum or division` == "Nemata" ~ "Others", 
                                   TRUE ~ Taxonomy.Group),
           Taxonomy.Group = case_when(is.na(Taxonomy.Group) & `Taxonomic phylum or division` == "Ciliophora" ~ "Others", 
                                   TRUE ~ Taxonomy.Group),
           Taxonomy.Group = case_when(is.na(Taxonomy.Group) & `Taxonomic phylum or division` == "Chlorophyta" ~ "Algae", 
                                   TRUE ~ Taxonomy.Group),
           Taxonomy.Group = case_when(is.na(Taxonomy.Group) & `Taxonomic phylum or division` == "Rhodophycota" ~ "Algae", 
                                   TRUE ~ Taxonomy.Group)
           ) %>% 
  # Also removing data where Species annotations are missing by:
    # 1. Removing any sub-species information to allow for more homogenus species averaging
        # Should do somethign about the var. sub-specification as well!!
              separate(col = Species, sep = " ssp ", into = "Species", remove = F, convert = T) %>% 
              separate(col = Species, sep = " var ", into = "Species", remove = F, convert = T) %>% 
              separate(col = Species, sep = " ", into = "Genus", remove = F, convert = T) #%>% 
    # 2. Removing rows that contain "sp." and lack a "space" character, but keeping "spp" annotations.
              #filter(!grepl(" sp\\.", Species))
  
EnviroTox_toxDB <- readxl::read_xlsx(
  path = "../data/Envirotox_2023-01-10/envirotox_DB_20230110081208.xlsx", 
  sheet = "test", col_names = T, 
  col_types = c("skip", rep("text", times = 4), "guess", rep("text", times = 4), 
                "guess", "guess", "logical", rep("text", times = 3), "numeric"
                )
  ) %>% 
  rename(
   Endpoint = `Test statistic`,
   AcuteChronic= `Test type`,
   Species = `Latin name`,
   Value.mg_l = `Effect value`,
   Time.Hours = `Duration (hours)`) %>% 
  mutate(
    DB = "EnviroTox",
    CAS.Number = as.cas(`original CAS`),
    AcuteChronic = case_when(
      AcuteChronic == "A" ~ "Acute",
      TRUE ~ "Chronic")) %>% 
  left_join(x = ., # Adding an Endpoint conversion column - a manually curated .xlsx file with annotations for endpoints conversion IDs
              y = Endpoints_filter, 
              by = "Endpoint") %>% 
          rename(Endpoint_conv = Group) %>% 
  left_join(x = ., # Adding Taxonomic Information, wrangled prior to this operation.
            y = EnviroTox_taxaDB %>% 
                select(Species, Genus, Medium, Taxonomy.Group),
            by = "Species") %>%
  select(CAS.Number, 1:22) %>% 
  select(-c("original CAS")) %>% 
  filter(!is.na(Endpoint_conv)) %>% 
  # Adding a dot after the "sp." annotation of species defined only down to genus.
  mutate(Species = sub("sp$", "sp.", Species))



## ---------------------------------------------------------------------------------------------------------------------------------------
# load EC10eq converson function
source("EC10eq_conversion_functions.R")

EnviroTox_toxDB_EC10eq <- EnviroTox_toxDB %>% 
  filter(
    Medium == "Freshwater", # including freshwater species only
    Time.Hours != "NA"
  ) %>% # removing time duration defined as "NA"
  mutate(
    EC10eq = mapply(ec10eq_extrapolation_function, Value.mg_l, Endpoint, AcuteChronic, Taxonomy.Group, "extpl"),
    EC10eq_high = mapply(ec10eq_extrapolation_function, Value.mg_l, Endpoint, AcuteChronic, Taxonomy.Group, "high_CI"),
    EC10eq_low = mapply(ec10eq_extrapolation_function, Value.mg_l, Endpoint, AcuteChronic, Taxonomy.Group, "low_CI"),
    EC10eq_Acute = case_when(AcuteChronic == "Acute" ~ EC10eq, TRUE ~ as.numeric(NA)),
    EC10eq_Chronic = case_when(AcuteChronic == "Chronic" ~ EC10eq, TRUE ~ as.numeric(NA)),
    # Making sure to annotate whether an effect value is extrapolated or not for downstream applications.
    No.Extrapolations = as.numeric(case_when(Endpoint_conv == "EC10" ~ 0, TRUE ~ 1)), 
    Database = "EnviroTox"
  ) %>%  
  # selecting relevant columns and merging some substance data
  select(
    CAS.Number, Species, Genus, Effect, Value.mg_l, AcuteChronic, Endpoint, 
    Endpoint_conv, Time.Hours, Source, version, DB, Medium, Taxonomy.Group, 
    EC10eq, No.Extrapolations, EC10eq_Chronic, EC10eq_Acute, Database
  )  

 write.csv(EnviroTox_toxDB_EC10eq, "../data/EnviroTox_toxDB_EC10eq.csv", row.names = F)


