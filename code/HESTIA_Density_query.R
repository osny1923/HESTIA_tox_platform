# Adding a request from Joseph on 2023-03-23 - density query 
library(tidyverse)
library(webchem)


PC_density_inventory <- read.csv("data/HESTIA_additional_data/PubChem_compound_list_density_2023-03-23.csv") %>% 
  rename(PubChemID = ï..cid) %>% 
  select(PubChemID, cmpdname, isosmiles, iupacname, inchikey)

# Dependency on `Pesticide_annotations.Rmd`
source(knitr::purl("code/Pesticide_annotations.Rmd"))

Density_query <- PC_density_inventory %>% 
  filter(PubChemID %in% HESTIA_Comp_info_7$PubChemId)

# Saving the query which takes half an hour-ish.
# Density_HESTIA <- pc_sect(Density_query$PubChemID, section = "Density", domain = "compound", verbose = T)
# write.csv(Density_HESTIA, "data/HESTIA_additional_data/Density_HESTIA.csv", row.names = F)
Density_HESTIA <- read.csv("data/HESTIA_additional_data/Density_HESTIA.csv")
Density_HESTIA_deliverable <- Density_HESTIA %>% 
  rename(PubChemId = CID,
         PC_name = Name) %>% 
  mutate(PubChemId = as.integer(PubChemId)) %>% 
  left_join(x = ., 
            y = HESTIA_Comp_info_7 %>% 
              filter(PubChemId %in% Density_HESTIA$CID) %>% 
              select(CAS.Number, PubChemId, PesticideAI_name),
            by = "PubChemId") %>% 
  select(CAS.Number, PesticideAI_name, PubChemId, 2:5)

write.csv(Density_HESTIA_deliverable, "data/HESTIA_additional_data/Density_HESTIA_long_format.csv", row.names = F)

Density_HESTIA_deliverable %>% 
  filter(PesticideAI_name == "Lindane")

Density_HESTIA_wide_format <- pivot_wider(Density_HESTIA_deliverable, id_cols = c("CAS.Number", "PesticideAI_name", "PubChemId", "PC_name"), names_from = SourceName, values_from = Result) %>% 
  # unnest(`CAMEO Chemicals`, .drop = F) %>% 
  # unnest(`Hazardous Substances Data Bank (HSDB)`, .drop = F) %>% 
  # unnest(`ILO International Chemical Safety Cards (ICSC)`, .drop = F) %>% 
  # unnest(`Joint FAO/WHO Expert Committee on Food Additives (JECFA)`, .drop = F) %>% 
  # unnest(`Occupational Safety and Health Administration (OSHA)`, .drop = F) %>% 
  # unnest(`The National Institute for Occupational Safety and Health (NIOSH)`, .drop = F) %>% 
  mutate(across(5:10, ~ as.character(.x))) %>% 
  arrange(PubChemId)

write.csv(Density_HESTIA_wide_format, "data/HESTIA_additional_data/Density_HESTIA_wide_format.csv", row.names = F)
# UNNESTING the wider format df is tricky!!
# here is a function that does it, but not great. (https://stackoverflow.com/questions/56372140/how-to-unnest-multiple-list-columns-of-a-dataframe-in-one-go-with-dplyr-pipe)
# unnest_cross <- function(data, cols, ...) {
#     .df_out <- data
#     .cols <- tidyselect::eval_select(rlang::enquo(cols), data)
#     purrr::walk(
#         .cols,
#         function(col) {
#             .df_out <<- unnest(.df_out, {{ col }}, ...)
#         }
#     )
#     .df_out
# }
# 
# Unnested_density_HESTIA_wide_format <- unnest_cross(Density_HESTIA_wide_format, tidyselect::ends_with("\\)"))
