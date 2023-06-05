# Reading and initial filtering of the raw toxicological data output from QSAR Toolbox.

# Dependencies:
library(tidyverse)

# Making a vector containing the columns that I want to include, based on the vector Full_name_include where many redundant columns are excluded:
columns_to_incl <- c('CAS.Number', 'SMILES', 'Molecular.formula', 'Predefined.substance.type', 'Additional.Ids', 'Identity', 'Database', 'Effect', 'Endpoint', 'Comments',
                     'Endpoint.comments', 'Title', 'Author', 'Year', 'Publication.year', 'Reference.source', 'URL', 'DOI', "Substance.Test.material.equality", 'Test.method', 'Test.method.comments', 
                     'Test.guideline', 'Data.quality', 'Qualifier', 'Measurement', 'Control.type', 'Organism.lifestage', 'Experimental.design', 
                     'Test.type', 'Media.type', 'Water.type', 'Subhabitat','Test.organisms..species.', 'Duration.MeanValue', 'Duration.Qualifier',
                     'Duration.Unit', 'Exposure.duration.MeanValue', 'Exposure.duration.Qualifier', 'Exposure.duration.Unit', 'Value.MeanValue', 'Value.Qualifier', 'Value.Unit', 'Value.Scale',
                     'Value.MinValue', 'Value.MinQualifier', 'Value.MaxValue', 'Value.MaxQualifier', 
                     'Superphylum', 'Superdivision', 'Phylum', 'Division', 'Subphylum', 'Subdivision', 'Infraphylum', 'Superclass', 'Class', 'Subclass', 'Infraclass', 'Superorder', 'Order', 
                     'Suborder', 'Infraorder', 'Superfamily', 'Family', 'Subfamily', 'Tribe', 'Genus', 'Subgenus', 'Section', 'Subsection'
                     )

# Defining which endpoints that should be included. missing or other endpoints are to be left out of the analysis
endpoints_to_incl <- c(
  "NOEC", "NOEL", "LC0", "EC0", "NOER", "NOAEC", # NOEC equivalent endpoints
  "LC50", "EC50", "IC50", "LD50", "ER50", "EL50", # EC50 equivalent endpoints
  "LOEC","EC10","LC10","IC10","LD10","ER10","EL10" # EC10 equivalent endpoints
  )

# Read each raw data and rbind into one complete df
HESTIA_HC20_DB_raw <- rbind(
  read.csv("data/OECD_2023-03-28/QSAR_Aquatic_ALL_TOX_2023-03-28_with_metadata[1-4k].csv",header = T, sep = "\t", na.strings = "", fileEncoding = "UTF-16LE",  stringsAsFactors = FALSE) %>% 
    select(all_of(columns_to_incl)),
  read.csv("data/OECD_2023-03-28/QSAR_Aquatic_ALL_TOX_2023-03-28_with_metadata[4k-8k].csv",header = T, sep = "\t", na.strings = "", fileEncoding = "UTF-16LE",  stringsAsFactors = FALSE) %>% 
    select(all_of(columns_to_incl)),
  read.csv("data/OECD_2023-03-28/QSAR_Aquatic_ALL_TOX_2023-03-28_with_metadata[8k-12k].csv",header = T, sep = "\t", na.strings = "", fileEncoding = "UTF-16LE",  stringsAsFactors = FALSE) %>% 
    select(all_of(columns_to_incl)),
  read.csv("data/OECD_2023-03-28/QSAR_Aquatic_ALL_TOX_2023-03-28_with_metadata[12k-16k].csv",header = T, sep = "\t", na.strings = "", fileEncoding = "UTF-16LE",  stringsAsFactors = FALSE) %>% 
    select(all_of(columns_to_incl)),
  read.csv("data/OECD_2023-03-28/QSAR_Aquatic_ALL_TOX_2023-03-28_with_metadata[16k-end].csv",header = T, sep = "\t", na.strings = "", fileEncoding = "UTF-16LE",  stringsAsFactors = FALSE) %>% 
    select(all_of(columns_to_incl))
) %>% 
  # Removing empty rows of data (If 'Database' is empty, no data exists in other columns)
  filter(!is.na(Database)) %>% 
  # Selecting only relevant endpoints
  filter(Endpoint %in% endpoints_to_incl)%>% 
  distinct() # Instantly removing duplicated rows

#######
#OUTPUT
#######

# The main dataset, NAs removed and ready to be dealt with.
write.csv(HESTIA_HC20_DB_raw, "data/RAW_DB/HESTIA_DB_pre_filtered.csv", row.names = FALSE)

# Creating a more condensed data file for the toxicological data wrangle, to reduce runtime
write.csv(HESTIA_HC20_DB_raw[, 1:48], "data/RAW_DB/HESTIA_HC20_DB_raw_toxdata.csv", row.names = FALSE)
