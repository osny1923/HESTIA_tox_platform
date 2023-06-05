library(tidyverse)
library(webchem)
# Read the list of chemicals in ´PesticideAI.csv` and fetch chemical information on these from CIR
CAS.list <- read.csv("data/excel_references/pesticideAI.csv", header = T, sep = ",") %>%
  separate(term.id, c("cas.termstz","CAS.id"), sep = "CAS-", remove = FALSE) %>%
  select(CAS.id) %>%
  rename(CAS.number = CAS.id) 

# Perform the CIR_query to fetch information on chemicals.           

#####################################################################
## WARNING - THIS OPERATION TAKES HOURS!!                          ##
## cir_query has a hard coded timeout of 1 sec between queries.    ##
## 17 000 sec equals 4.7 hours                                     ##
#####################################################################

# CIR_CAS_SMILES <- cir_query(CAS.list$CAS.number, "SMILES", match = "first", verbose = TRUE)
 # Output is exported as .csv file to
# write.csv(CIR_CAS_SMILES, file = "data/excel_references/CIR_CAS_SMILES.csv") 
CIR_CAS_SMILES <- read.csv("data/excel_references/CIR_CAS_SMILES.csv") 
### I should add any available SMILES from other repositories. ###

# looking for SMILES in ChemIdPlus
p <- CIR_CAS_SMILES %>% 
  filter(is.na(SMILES)) %>% 
  pull(CAS.Number)

#  ChemIDPlus query
# b <- webchem::ci_query(p, from = "rn", verbose = T)
# Saving the ChemIDPlus query as a .txt object
# saveRDS(b, file = "data/excel_references/ChemIDPlus_query.txt")
ChemIDPlus_query <- readRDS("data/excel_references/ChemIDPlus_query.txt")
ChemIDPlus_query_dropna <- ChemIDPlus_query[!is.na(ChemIDPlus_query)]
ChemIDPlus_query_rbind <- as.data.frame(do.call(rbind, ChemIDPlus_query_dropna)) %>% 
  select(cas, smiles) %>% 
  filter(!is.na(smiles),
         smiles != "")
rownames(ChemIDPlus_query_rbind) = NULL
ChemIDPlus_query_rbind$cas <- as.character(ChemIDPlus_query_rbind$cas)
ChemIDPlus_query_rbind$smiles <- as.character(ChemIDPlus_query_rbind$smiles)

ChemIDPlus_query_rbind_filtered <- ChemIDPlus_query_rbind %>% 
  rename(CAS.Number = cas) %>% 
  filter(!grepl(c("['Rgp']"), smiles),
         !grepl(c("[*]"), smiles))

# adding the 167 SMILES configurations from ChemIDPlus query to QSAR_TOX df.
CAS_SMILES <- CIR_CAS_SMILES %>% 
  left_join(., 
            ChemIDPlus_query_rbind_filtered,
            by = "CAS.Number") %>% 
  mutate(SMILES = coalesce(SMILES, smiles)) %>% 
  select(-smiles)

# Ok, There are still 299 SMILES configurations missing. 


# Querying the remaining SMILES in the Pubchem database
outside_ci <- CAS_SMILES %>% 
  distinct(CAS.Number, .keep_all = T) %>% 
  filter(is.na(SMILES)) %>% 
  pull(CAS.Number)

# pubchem query
#pubchem_id <- get_cid(outside_ci, from = "xref/rn", match = "all", verbose = T)
# Saving the ChemIDPlus query as a .csv object
# write.csv(pubchem_id, "data/excel_references/pubchem_id.csv")

# load the query from file:
pubchem_id <- read.csv("data/excel_references/pubchem_id.csv", header = T) %>% 
    distinct(query, .keep_all = T)
 
# pubchem chemical properties (SMILES) query
#pubchem_prop <- pc_prop(pubchem_id$cid, "IsomericSMILES", verbose = T)

# Saving the ChemIDPlus query as a .csv object
# write.csv(pubchem_prop, "data/excel_references/pubchem_prop.csv")
# load the query from file
pubchem_prop <- read.csv("data/excel_references/pubchem_prop.csv", header = T)

pubchem_smiles <- data.frame(CAS.Number = pubchem_id$query, 
                             IsomericSMILES = pubchem_prop$IsomericSMILES)
                              
# Query returned another 134 SMILES!
# Adding these to the dataset
FULL_CAS_SMILES <- CAS_SMILES %>% 
  left_join(., 
            pubchem_smiles,
            by = "CAS.Number") %>% 
  mutate(SMILES = coalesce(SMILES, IsomericSMILES)) %>% 
  select(CAS.Number, SMILES)

write.csv(FULL_CAS_SMILES, file = "data/excel_references/Full_CAS_SMILES.csv", row.names = F) 

# Exporting the full dataset from QSAR toolbox (16797 CAS.Numbers with duplicates and all tox data = 600k+ data points. Server threw an exception and i had to split the data into 5 parts.)
# After this, i can just `rbind` those data sets along with a filter function (include the endpoints i want to work with)
# What follows is the documentation of CAS-list and method to query SMILES-configurations from CAS numbers. 
# I use the CAS numbers supplied in the ´PesticideAI.csv` file from the HESTIA home page to query for matching SMILES configurations. These two vectors are then used as input into the QSAR toolbox application.
# Exporting the CAS_SMILES object as 4k row subsets due to issues when wxporting too large datasets from OECD QSAR Toolbox.
cas_smiles_list <- read.csv("data/excel_references/Full_CAS_SMILES.csv")

cas_smiles_list[1:4000,] %>% 
  write.table("data/cas_smiles_list[1-4k].txt", quote = FALSE, col.names = F, row.names = F, sep = "\t")

cas_smiles_list[4001:8000,] %>% 
  write.table("data/cas_smiles_list[4k-8k].txt", quote = FALSE,  col.names = F, row.names = F, sep = "\t")

cas_smiles_list[8001:12000,] %>% 
  write.table("data/cas_smiles_list[8k-12k].txt", quote = FALSE, col.names = F, row.names = F, sep = "\t")

cas_smiles_list[12001:16000,] %>% 
  write.table("data/cas_smiles_list[12-16k].txt",  quote = FALSE, col.names = F, row.names = F, sep = "\t")

cas_smiles_list[16001:nrow(cas_smiles_list),] %>% 
  write.table("data/cas_smiles_list[16k-end].txt",  quote = FALSE, col.names = F, row.names = F, sep = "\t")

