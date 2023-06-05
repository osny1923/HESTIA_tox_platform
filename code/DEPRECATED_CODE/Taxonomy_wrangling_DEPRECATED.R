## DEPRECATED Taxonomy wrangling code


## Making a Taxonomic classification list for all species

```{r, include = F, eval=FALSE}
# NCBI_Taxon_query <- rbind(classification(get_uid(Species_list$Test.organisms..species.), db = 'ncbi', rows = 1, verbose = TRUE))
# saveRDS(NCBI_Taxon_query, file = "../data/Taxonomy/NCBI_Taxon_query.txt")
NCBI_Taxon_query <- readRDS("../data/Taxonomy/NCBI_Taxon_query.txt")
# NCBI_Taxon_query_matrix <- NCBI_Taxon_query %>%
#      select(-id) %>% 
#    pivot_wider(names_from = rank, values_from = name)

# Extracting each desired taxonomic level
NCBI_phylum <- NCBI_Taxon_query %>%
  filter(rank == "phylum") 
NCBI_genus <- NCBI_Taxon_query %>%
  filter(rank == "genus") 
NCBI_order <- NCBI_Taxon_query %>%
  filter(rank == "order") 
NCBI_class <- NCBI_Taxon_query %>%
  filter(rank == "class") 
NCBI_species <- NCBI_Taxon_query %>%
  filter(rank == "species") 

# Binding together the query into one comprehensive df based on the species query number.
NCBI_taxonomy_df <- NCBI_phylum %>% 
  left_join(
    ., 
    NCBI_genus, 
    by = "query") %>% 
  rename(Genus = name.y, Phylum = name.x) %>% 
  select(c(4, 5,1)) %>% 
  distinct(Genus, .keep_all = TRUE) %>% 
  left_join(
    ., 
    NCBI_order %>% 
      select(-c(id, rank)) %>% 
      rename(Order = name) %>% 
      distinct(query, .keep_all = T),
    by = "query") %>% 
  left_join(
    ., 
    NCBI_class %>% 
      select(-c(id, rank)) %>% 
      rename(Class = name) %>% 
      distinct(query, .keep_all = T),
    by = "query") %>% 
  left_join(
    ., 
    NCBI_species %>% 
      select(-c(id, rank)) %>% 
      rename(Species = name) %>% 
      distinct(query, .keep_all = T),
    by = "query") %>% 
  select(c(Species, Genus, Order, Class, Phylum)) %>% 
  mutate(Source = "NCBI")

NCBI_DB_merge <- rbind(
  NCBI_taxonomy_df,
  Species_list %>% 
    filter(!Test.organisms..species. %in% NCBI_taxonomy_df$Species,
           grepl(" ", Test.organisms..species.)) %>% 
    separate(Test.organisms..species., sep = " ", remove = F, into = "Genus") %>% 
    left_join(x = .,
              y = NCBI_taxonomy_df %>% 
                select(-Species),
              by = "Genus") %>% 
    rename(Species = Test.organisms..species.) %>% 
    filter(!is.na(Source))
) %>% 
  filter(!is.na(Species))



# For the missing taxonomic descriptors, where names are old or common or just bad, a second search is performed.
# Extracting the recods without a match
missing_names <- Species_list %>% 
  filter(!Test.organisms..species. %in% NCBI_DB_merge$Species,
         # this grepl filters out non-binomial taxonomic names, e.g "Gammarus pulex" is kept, but "Gammarus" is filtered out, since it is missing the " ".
         grepl(" ", Test.organisms..species.))

# Taxonomy _identification number for NCBI_ Query (UID - NCBI database identifier)
# At two points, a question to select data is asked. "1" is selected.
good_names <- get_uid(missing_names$Test.organisms..species.[1:5], rank_query = "Species")
# write_rds(good_names, "../data/Taxonomy/Taxonomy_UID.txt")
good_names <- readRDS("../data/Taxonomy/Taxonomy_UID.txt")

good_names
missing_names$uid <- good_names

# NCBI_missing_names <- classification(id = good_names, db = "ncbi")
# write_rds(NCBI_missing_names, "Taxonomy_classification_UID.txt")
NCBI_missing_names <- readRDS("../data/Taxonomy/Taxonomy_classification_UID.txt")

NCBI_missing_names_rbind <- rbind(NCBI_missing_names)
pivot_wider(NCBI_missing_names_rbind, names_from = c(id, rank), values_from = name)
# Binding rows with each taxonomic rank and then compiling into one df
NCBI_phylum <- NCBI_missing_names_rbind %>%
  filter(rank == "phylum") 
NCBI_genus <- NCBI_missing_names_rbind %>%
  filter(rank == "genus") 
NCBI_order <- NCBI_missing_names_rbind %>%
  filter(rank == "order") 
NCBI_class <- NCBI_missing_names_rbind %>%
  filter(rank == "class") 
NCBI_species <- NCBI_missing_names_rbind %>%
  filter(rank == "species") 

missing_names_second <- NCBI_phylum %>% 
  left_join(NCBI_genus, by = "query") %>% 
  rename(Genus = name.y, Phylum = name.x) %>% 
  select(c(4, 5,1)) %>% 
  distinct(Genus, .keep_all = TRUE) %>% 
  left_join(., 
            NCBI_order %>% 
              select(-c(id, rank)) %>% 
              rename(Order = name) %>% 
              distinct(query, .keep_all = T),
            by = "query") %>% 
  left_join(., 
            NCBI_class %>% 
              select(-c(id, rank)) %>% 
              rename(Class = name) %>% 
              distinct(query, .keep_all = T),
            by = "query") %>% 
  left_join(., 
            NCBI_species %>% 
              select(-c(id, rank)) %>% 
              rename(Species = name) %>% 
              distinct(query, .keep_all = T),
            by = "query") %>% 
  select(c(Species, Genus, Order, Class, Phylum, query)) %>% 
  mutate(Source = "NCBI")


third_taxonomy <- rbind(NCBI_DB_merge %>% 
                          mutate(Species_update = NA_character_),
                        
                        left_join(x = missing_names %>% 
                                    rename(Species = Test.organisms..species.,
                                           query = uid) %>% 
                                    mutate(query = as.character(query)), 
                                  y = missing_names_second %>% 
                                    rename(Species_update = Species),
                                  by = "query") %>% 
                          select(-query)
)

baddest_names <- Species_list %>% 
  filter(!Test.organisms..species. %in% third_taxonomy$Species)

# Manualy curated list
# Last_taxonomy_list <- read.xlsx("missing_phyla_manual.xlsx", sheetName = "Data1", header = T, colIndex = c(1:6))

fourth_taxonomy <- rbind(third_taxonomy,
                         left_join(x = baddest_names %>% 
                                     rename(Species = Test.organisms..species.),
                                   y = Last_taxonomy_list,
                                   by = "Species") %>% 
                           mutate(Species_update = NA_character_) 
)

# completing the Taxonomic descriptions and merging with the big dataset.
HESTIA_HC20_DB_taxonomy <- HESTIA_HC20_DB_prefil %>% 
  left_join(x = .,
            fourth_taxonomy %>% 
              rename(Test.organisms..species. = Species),
            by = "Test.organisms..species.") %>% 
  separate(col = Phylum.x, sep = " ", into = "Phylum_internal", remove = F, convert = T) %>% 
  separate(col = Division, sep = " ", into = "Division_internal", remove = F, convert = T) %>% 
  separate(col = Class.x, sep = " ", into = "Class_internal", remove = F, convert = T) %>% 
  separate(col = Order.x, sep = " ", into = "Order_internal", remove = F, convert = T) %>% 
  mutate(Phylum.y = case_when(is.na(Phylum.y) ~ Phylum_internal, 
                              TRUE ~ Phylum.y),
         Phylum.y = case_when(is.na(Phylum.y) & is.na(Phylum_internal) ~ Division_internal, 
                              TRUE ~ Phylum.y),
         Class.y = case_when(is.na(Class.y) ~ Class_internal, 
                             TRUE ~ Class.y),
         Order.y = case_when(is.na(Order.y) ~ Order_internal, 
                             TRUE ~ Order.y),
  ) %>% 
  select(-c("Superphylum", "Superdivision", "Phylum.x", "Division", "Subphylum", "Subdivision", "Infraphylum", "Superclass", "Class.x", "Subclass", "Infraclass", "Superorder", "Order.x", "Suborder", "Infraorder", "Superfamily", "Family", "Subfamily", "Tribe", "Genus.x", "Subgenus", "Section", "Subsection", "Phylum_internal", "Division_internal", "Class_internal","Order_internal")) %>% 
  rename(Genus = Genus.y,
         Order = Order.y,
         Class = Class.y,
         Phylum = Phylum.y)


Missing_genus_list <- HESTIA_HC20_DB_taxonomy %>% 
  filter(is.na(Source)) %>% 
  distinct(Test.organisms..species.) %>% 
  separate(col = Test.organisms..species., sep = " ", into = "Genus", remove = F, convert = T) %>% 
  filter(!Genus %in% Final_taxa_list$Genus,
         !grepl("not", Genus)) %>% 
  distinct(Genus,.keep_all = T) 



#Miss_class <- classification(get_uid(Missing_genus_list$Genus, rank_query = "genus"))
# write_rds(Miss_class, "Taxonomy_missing_classes.txt")
Miss_class <- readRDS("Taxonomy_missing_classes.txt")

Miss_class_nona <- Miss_class[!is.na(Miss_class)]

Miss_class_rbind <- rbind(Miss_class) # Flatten the nested lists into one comprehensive data frame

miss_class_df <- Miss_class_rbind %>% 
  filter(rank %in% c("clade", "class", "order", "genus", "phylum")) %>% 
  pivot_wider(names_from = rank, values_from = name, id_cols = query)


merge_names_missing <- left_join(x = Missing_genus_list,
                                 y = miss_class_df %>% 
                                   rename(Genus = genus,
                                          Order = order) %>%
                                   mutate(across(everything(), ~ as.character(.))) %>%
                                   select(-c(query, clade)), 
                                 by = "Genus")

fifth_taxonomy <- rbind(
  fourth_taxonomy %>% 
    filter(
      !Species %in% merge_names_missing$Test.organisms..species.),
  merge_names_missing %>% 
    mutate(Source = "NCBI",
           phylum = gsub("NULL", NA_character_, phylum),
           Order = gsub("NULL", NA_character_, Order)) %>% 
    # Filling in the blank Phylum-annotations by using the "fourth taxonomy" list.
    left_join(
      x = ., 
      y = fourth_taxonomy %>% 
        distinct(Order, .keep_all = T) %>% 
        select(c(Order, Phylum)) %>% 
        filter(!is.na(Order),
               !Order == ""),
      by = "Order") %>% 
    mutate(phylum = case_when(is.na(phylum) & !is.na(Phylum) ~ Phylum,
                              TRUE ~ phylum)) %>% 
    select(-Phylum) %>% 
    rename(Class = class,
           Phylum = phylum,
           Species = Test.organisms..species.) %>% 
    mutate(Species_update = NA_character_) 
) %>% 
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
  )

GBIF_class <- fifth_taxonomy %>% 
  filter(is.na(Phylum)) %>% 
  select(Species)

#GBIF_class$gbif_id <- get_gbifid(Wiki_class$Species, match = "first", verbose = T)
# write_rds(GBIF_class$gbif_id, "Taxonomy_GBIF_id_query.txt")
GBIF_class$gbif_id <- readRDS("Taxonomy_GBIF_id_query.txt")


GBIF_query <- GBIF_class %>% 
  filter(!is.na(gbif_id))
GBIF_query <- GBIF_query$gbif_id
#GBIF_class_list <- classification(GBIF_query, db = "gbif")
#write_rds(GBIF_class_list, "Taxonomy_GBIF_query.txt")
GBIF_class_list <- readRDS("Taxonomy_GBIF_query.txt")

GBIF_class_rbind <- rbind(GBIF_class_list)


GBIF_tax_list <- left_join(x = GBIF_class %>% 
                             rename(query = gbif_id) %>% 
                             mutate(query = as.character(query)),
                           y = GBIF_class_rbind %>% 
                             filter(rank %in% c("clade", "class", "order", "genus", "phylum")) %>% 
                             pivot_wider(names_from = rank, values_from = name, id_cols = query)) %>% 
  mutate(across(everything(), ~ as.character(.)),
         across(everything(), ~ case_when(.x == "NULL" ~ NA_character_,
                                          TRUE ~ .x))) %>% 
  select(-query)

# write.xlsx(GBIF_tax_list, "GBIF_class_list.xlsx", sheetName = "Data1", col.names = T, row.names = F, append = T, showNA = F)
GBIF_tax_list_cur <- read.xlsx("GBIF_class_list.xlsx", sheetName = "Data1", header = T)

#Combining the last part and cleaning up the taxonomy list from errors
sixth_taxonomy <-  rbind(GBIF_tax_list_cur %>% 
                           mutate(Source = "GBIF",
                                  Species_update = NA_character_) %>% 
                           rename(Phylum = phylum, 
                                  Class = class,
                                  Order = order,
                                  Genus = genus), #%>% 
                         #select(Species, Genus, Order, Class, Phylum, Source, Species_update),
                         fifth_taxonomy %>% 
                           filter(!is.na(Phylum))) %>% 
  mutate(Species_update = case_when(grepl("Selenastrum capricornutum", Species_update) ~ "Raphidocelis subcapitata",
                                    TRUE ~ Species_update),
         Class = case_when(grepl("Anura", Species) ~ "Amphibian",
                           TRUE ~ Class),
         across(everything(), ~ case_when(grepl("Algae", Species) ~ NA_character_,
                                          TRUE ~ .x)),
         Source = case_when(is.na(Phylum) ~ NA_character_,
                            TRUE ~ Source),
         Species_update = case_when(Species == "Common hornwort (live plants)" ~ "Ceratophyllum demersum",
                                    Species == "Green cabomba (live plants)"~ "Cabomba caroliniana",
                                    Species == "Harlequin fly" ~ "Chironomus riparius",
                                    Species == "Mysid shrimp" ~"Mysida",
                                    Species == "Orthocladiinae" ~ "Orthocladiinae",
                                    Species == "Zygoptera" ~ "Zygoptera",
                                    Species == "Copepoda" ~ "Copepoda",
                                    Species == "Tanypodinae" ~ "Tanypodinae",
                                    Species == "Common duckweed" ~ "Lemna minor",
                                    Species == "Parrot feather" ~ "Myriophyllum aquaticum",
                                    Species == "Green alga" ~ "Chlorophyta",
                                    Species == "Brown trout" ~ "Salmo trutta",
                                    Species == "Eurasian watermilfoil" ~ "Myriophyllum spicatum",
                                    Species == "Tanytarsini" ~ "Tanytarsini",
                                    Species == "Daggerblade grass Shrimp" ~ "Palaemonetes pugio",
                                    Species == "Ectyphinae" ~ "Ectyphinae",
                                    Species == "Basommatophora" ~ "Gastropoda",
                                    Species == "Conchostraca" ~ "Conchostraca",
                                    Species == "Floating sweet-grass" ~ "Glyceria fluitans",
                                    Species == "Cryptophycophyta" ~ "",
                                    Species == "Wild celery (water celery) (live plants)" ~ "Vallisneria americana",
                                    Species == "Louisiana crayfish" ~ "Procambarus clarkii",
                                    Species == "Pyrrophycophyta" ~ "",
                                    Species == "American waterweed (Pondweed) (live plants)" ~ "Elodea canadensis"),
         Source = case_when(!is.na(Species_update) & is.na(Source) ~ "Wiki",
                            TRUE ~ Source)
  ) %>%
  filter(!is.na(Species))  

# write.xlsx(sixth_taxonomy, "Final_Taxonomy_dataset.xlsx", sheetName = "Data1", append = T, col.names = T, row.names = F, showNA = F)

# I wiggled a little with this final list.
# It is done and i think i need a new document to merge it all into!

# also, need to match the Taxonomic.Group thingy with Saouter!!


seventh_taxonomy <- read.xlsx("Final_Taxonomy_dataset.xlsx", sheetName = "Data1", header = T)  

Final_taxonomy <- seventh_taxonomy %>% 
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
    TRUE ~ "Others")
  )

# write.xlsx(Final_taxonomy, "Final_Taxonomy_dataset.xlsx", sheetName = "Data1", append = T, col.names = T, row.names = F, showNA = F)

```
