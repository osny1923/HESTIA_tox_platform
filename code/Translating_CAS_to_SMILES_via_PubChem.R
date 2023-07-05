library(tidyverse)
library(webchem)
# Purpose: collect correct SMILES configurations for as many CAS.Numbers as possible in the pesticideAI.csv Hestia inventory:

# Read the HESTIA list of chemials and fetch chemical identification-information on these from PubChem
  # Pubchem is selected since the queries are fast and have much other relevant information if a PubChemId is matched.
# Importing the HESTIA db and create a list of the CAS.numbers
CAS.list_HESTIA <- read.csv("data/excel_references/pesticideAI.csv", header = T, sep = ",") %>% 
  separate(term.id, c("cas.termstz","CAS.id"), sep = "CAS-", remove = FALSE) %>% 
  select(CAS.id, term.name) %>% 
  rename(CAS.number = CAS.id)

# Retrieving PubChemId, based on CAS RN.
  # CAS_CID_list_fromname <- get_cid(CAS.list_HESTIA$PesticideAI_name, from = "name", domain = "compound", match = "all", verbose = T)
# Saving the returned data as file
  # write_rds(CAS_CID_list_fromname, file = "data/excel_references/PubChem_CID_fromname.txt")

# Read returned data as .txt
PC_id_fromname <- read_rds("data/excel_references/PubChem_CID_fromname.txt") %>% 
  rename(term.name = query) 

# Or i could join PesticideAI_names & cid with the CAS.Numbers and thereafter look at which CAS.Numbers match a query together with SMILES
PC_id_fromname_CAS_name_smiles <- PC_id_fromname %>% 
  left_join(x = ., 
            y = CAS.list_HESTIA, 
            by = "term.name")

# And now do a query for SMILES and Title where Title should match PesticideAI_name"
  # PC_id_fromname_CAS_name_smiles[,4:6] <- pc_prop(PC_id_fromname_CAS_name_smiles$cid, properties = c("CanonicalSMILES", "Title"), verbose = T)
# Saving the returned data as file
  # write.csv(PC_id_fromname_CAS_name_smiles, "data/excel_references/PC_id_fromname_CAS_name_smiles.csv", row.names = F)
# Read returned data as .txt
PC_id_fromname_CAS_name_smiles <- read.csv("data/excel_references/PC_id_fromname_CAS_name_smiles.csv") %>% 
# Select relevant information
  select(CAS.number, term.name, Title, cid, CID, CanonicalSMILES) %>% 
# Better descriptions of columns
  rename(PC_title = Title, 
         PesticideAI_name = term.name, 
         PC_SMILES = CanonicalSMILES)

# Cleaning up the duplicates: (The unique matches are used further down)
duplicate_matches <- PC_id_fromname_CAS_name_smiles %>% 
  # Selecting all duplicated substances
  group_by(CAS.number) %>% filter(n() > 1) %>% ungroup() %>% 
# This gets matches for many substances except where letter-cases are heterogenous and in 72 cases no matches are found.
# Making all letters lower case 
  mutate(PesticideAI_name = tolower(PesticideAI_name),
         PC_title = tolower(PC_title))

# Adding a "MATCH" annotation when the correct name is matched between PesticideAI and PubChem query
duplicate_matches_matches <- duplicate_matches %>% 
  mutate(MATCH = case_when(PesticideAI_name == PC_title ~ "MATCH", TRUE ~ as.character(NA))) 

# Generating a data frame of unique substances with matching names and a SMILES that is correct
Matched_items_list <- duplicate_matches_matches %>% 
  group_by(CAS.number) %>% 
  filter(MATCH %in% "MATCH")

# These are the missing cases. By using the distinct function across SMILES, if there is only one CAS.number left, all other had identical smiles and should be valid, despite no matching PestcideAI and PubChem names
filtered_noMatch_unique_smiles <- duplicate_matches_matches %>% 
  filter(!CAS.number %in% Matched_items_list$CAS.number) %>% 
  group_by(CAS.number) %>% 
  distinct(PC_SMILES, .keep_all = T) %>% 
  group_by(CAS.number) %>% filter(n() == 1) %>% ungroup() %>% 
  mutate(MATCH = "MATCH")

# lastly, a reduced subset for where i use distinct(), as i cannot select substance in a reproducible way, nor have a way to distinguish between CAS/PesticideAI matches to PC_ID/PC_Title/PC_SMILES 
manually_curated_smiles <- duplicate_matches_matches %>% 
  filter(!CAS.number %in% Matched_items_list$CAS.number) %>% 
  group_by(CAS.number) %>% 
  distinct(PC_SMILES, .keep_all = T) %>% 
  group_by(CAS.number) %>% filter(n() > 1) %>% ungroup()

# Selecting the first match here
manually_curated_smiles_distinct <- manually_curated_smiles %>%   
  distinct(CAS.number, .keep_all = T) %>% 
  mutate(MATCH = "non-MATCH")

# Binding together a final list of substances
CAS_CID_list_first <- rbind(
  PC_id_fromname_CAS_name_smiles %>% 
    # Selecting all duplicated substances
    group_by(CAS.number) %>% filter(n() == 1) %>% ungroup() %>% 
    mutate(MATCH = "insta-MATCH"),
  Matched_items_list,
  filtered_noMatch_unique_smiles,
  manually_curated_smiles_distinct
) %>% 
  select(-c("MATCH", "cid")) %>% 
  rename(CAS.Number = CAS.number, 
         CanonicalSMILES = PC_SMILES)

##########################
# Second phase of looking up curated SMILES-CAS Matches
##########################
# From the original query, there were 3911 missing PC_IDs
missing_3911_PC_ID <- PC_id_fromname_CAS_name_smiles %>% 
  filter(is.na(x = cid))

# Looking into PC_ID using CAS RN
# missing_3911_PC_ID_CAS_query <- get_cid(missing_3911_PC_ID$CAS.number, from = "xref/RN", domain = "compound", match = "all", verbose = T)
#  write.csv(missing_3911_PC_ID_CAS_query, "data/excel_references/missing_3911_PC_ID_CAS_query.csv", row.names = F)
missing_3911_PC_ID_CAS_query <- read.csv("data/excel_references/missing_3911_PC_ID_CAS_query.csv")

# joining in the names for substances to investigate matches 
missing_3911_PC_ID_CAS_query_w_names <- missing_3911_PC_ID_CAS_query %>% 
  rename(CAS.number = query) %>% 
  left_join(
    x = ., 
    y = PC_id_fromname_CAS_name_smiles %>% 
      select(CAS.number, PesticideAI_name),
    by = "CAS.number") %>% 
  filter(!is.na(cid))

# A lot of information is missing due to inconsistent spelling when querying "name". commonly sulfate // sulphate spelling is causing this.
# To not have to go in and gsub a bunch of "sulf" --> "sulph" and maybe fix 10-30% of data, i do a query back to CASRN, so: 
# CASRN -> PC-ID (multiple hits per CAS!) -> CanonicalSMILES -> CASRN. to see which SMILES are correct.
missing_3911_PC_ID_CAS_query_w_names[,4:5] <- pc_prop(missing_3911_PC_ID_CAS_query_w_names$cid,properties = c("CanonicalSMILES", "Title"), verbose = T)[,2:3]


missing_3911_PC_ID_CAS_query_w_names_updated <- missing_3911_PC_ID_CAS_query_w_names %>% 
  mutate(Title = tolower(gsub("sulf", "sulph", Title, ignore.case = T)),
         PesticideAI_name = tolower(PesticideAI_name),
         MATCH = case_when(PesticideAI_name == Title ~ "MATCH", TRUE ~ as.character(NA)))

# Creating a df with the true matches, no double matches.
missing_3911_PC_ID_CAS_query_w_names_matched <- missing_3911_PC_ID_CAS_query_w_names_updated %>% 
  filter(MATCH == "MATCH")

# Removing all matched CAS from the search
missing_386_fewer_CAS <- missing_3911_PC_ID_CAS_query_w_names_updated %>% 
  filter(!CAS.number %in% missing_3911_PC_ID_CAS_query_w_names_matched$CAS.number)

missing_386_fewer_CAS_more_IDs <- missing_386_fewer_CAS %>% 
  distinct(CAS.number, CanonicalSMILES, .keep_all = T) %>% 
  group_by(CAS.number) %>% filter(n() == 1) %>% ungroup() %>%
  mutate(MATCH = "MATCH") 

# Merging two datasets with controlled SMILES match
Filling_final_SMILES <- rbind(missing_3911_PC_ID_CAS_query_w_names_matched,
                              missing_386_fewer_CAS_more_IDs)

last_remaining_CAS_SMILES <- missing_386_fewer_CAS %>% 
  filter(!CAS.number %in% Filling_final_SMILES$CAS.number)

### Using the pc_sect() function 
# cid_to_cas_query <- pc_sect(last_remaining_CAS_SMILES$cid, section = "CAS", domain = "compound", verbose = T)
# write.csv(cid_to_cas_query, "data/excel_references/cid_to_cas_query.csv", row.names = F)
# Reading the data from pc_sect query
cid_to_cas_query <- read.csv("data/excel_references/cid_to_cas_query.csv")

cid_to_cas_query_CAS_repo <- cid_to_cas_query %>% 
  filter(SourceName == "CAS Common Chemistry") 

joined_missing_df <- left_join(
    x = cid_to_cas_query_CAS_repo %>% 
      rename(cid = CID) %>% 
      mutate(cid = as.integer(cid)), 
    y = last_remaining_CAS_SMILES, 
    by = "cid") %>% 
  mutate(MATCH = case_when(
    CAS.number == Result ~ "MATCH",
      TRUE ~ as.character(NA))) %>% 
  filter(MATCH == "MATCH") %>% 
  distinct(CAS.number, CanonicalSMILES, .keep_all = T)

Unique_last_missing <- joined_missing_df %>% 
  group_by(CAS.number) %>% filter(n() == 1) %>% ungroup() %>% 
  select(6,1,7:10)

# For the last missing substances:
# Simply perform a `distinct()`operation, 
last_ambiguous_SMILES <- joined_missing_df %>% 
  group_by(CAS.number) %>% filter(n() > 1) %>% ungroup() %>% 
  distinct(CAS.number, .keep_all = T) %>% 
  select(6,1,7:10)

# Merging the third and the fourth dataset with SMILES MATCH
CAS_SMILES_filler <- rbind(
  Filling_final_SMILES, 
  Unique_last_missing, 
  last_ambiguous_SMILES)

write.csv(CAS_SMILES_filler, "data/excel_references/CAS_SMILES_filler.csv", row.names = F)

## Reading in the the data to not have loose dependencies above
last_missing_SMILES <- read.csv("data/excel_references/CAS_SMILES_filler.csv") %>% 
  rename(PC_title = Title,
         CID = cid,
         CAS.Number = CAS.number) %>% 
  select(names(CAS_CID_list_first))

# Third lookup
######################

# Last effort to secure the missing 2259 CAS numbers without SMILES
CAS_CID_Second_Step <- left_join(x = CAS_CID_list_first, y = last_missing_SMILES, by = "CAS.Number") %>% 
  mutate(PC_title = coalesce(PC_title.x, PC_title.y),
         CID = coalesce(CID.x, CID.y),
         CanonicalSMILES = coalesce(CanonicalSMILES.x, CanonicalSMILES.y)) %>% 
  select(-c("PesticideAI_name.y", "PC_title.x", "PC_title.y", "CID.x", "CID.y", "CanonicalSMILES.x", "CanonicalSMILES.y")) %>% 
  rename(PesticideAI_name = PesticideAI_name.x) %>% 
  mutate(CID = as.integer(CID))


####################
# select chemicals missing SMILES-configuration 
Missing_list_1 <- CAS_CID_Second_Step %>% 
  filter(is.na(CanonicalSMILES))

# Missing_list_1_cid <- get_cid(Missing_list_1$CAS.Number, from = "xref/RN", domain = "compound", match = "all", verbose = T)
# write.csv(Missing_list_1_cid, "data/excel_references/Missing_list_1_cid.csv", row.names = F)
Missing_list_1_cid <- read.csv("data/excel_references/Missing_list_1_cid.csv")
Missing_list_1_cid[,3] <- pc_prop(Missing_list_1_cid$cid, properties = "CanonicalSMILES",verbose = T)[,2]

Missing_list_2_cid <- Missing_list_1_cid %>% 
  filter(!is.na(cid)) %>% 
  rename(CanonicalSMILES = ...3,
         CAS.Number = query) 

Missing_list_3_cid <- left_join(
  x = Missing_list_2_cid, 
  y = Missing_list_1,
  by = "CAS.Number") %>% 
  select(-c("PC_title", "CID", "CanonicalSMILES.y")) %>% 
  distinct(CAS.Number, CanonicalSMILES.x, .keep_all = T)

# Missing_list_4_cid <- pc_sect(Missing_list_3_cid$cid, section = "CAS", domain = "compound", verbose = T)
# write.csv(Missing_list_4_cid, "data/excel_references/Missing_list_4_cid.csv", row.names = F)

Missing_list_4_cid <- read.csv("data/excel_references/Missing_list_4_cid.csv")

Missing_list_5 <- left_join(
  x = Missing_list_4_cid,
  y = Missing_list_3_cid %>% 
    rename(CID = cid),
  by = "CID")

Missing_list_5_MATCH <- Missing_list_5 %>% 
  filter(CAS.Number == Result) %>% 
  distinct(CAS.Number, CanonicalSMILES.x, .keep_all = T) %>% 
  distinct(CAS.Number, .keep_all = T) %>% 
  mutate(MATCH = "MATCH", 
         PC_title = "Not Queried") %>% 
  rename(CanonicalSMILES = CanonicalSMILES.x) %>% 
  select(names(CAS_CID_Second_Step))

write.csv(Missing_list_5_MATCH, "data/excel_references/Missing_list_5_MATCH.csv", row.names = F)

Missing_list_5_NOT_MATCH <- Missing_list_5 %>% 
  filter(CAS.Number != Result)

Missing_list_5_refined_NOT_MATCH <- Missing_list_5_NOT_MATCH %>% 
  filter(!CAS.Number %in% Missing_list_5_MATCH$CAS.Number) %>% 
  distinct(CAS.Number, CanonicalSMILES.x, .keep_all = T) %>% 
  mutate(MATCh = "NOT_MATCH") %>% 
  rename(CanonicalSMILES = CanonicalSMILES.x)

Missing_list_5_MATCH <- read.csv("data/excel_references/Missing_list_5_MATCH.csv")

# compiling all identified SMILES-configurations for the third operation 
last_missing_SMILES <- rbind(last_missing_SMILES, Missing_list_5_MATCH)

# nrow(last_missing_SMILES)

# compiling a list of SMILES-information for the operations.
CAS_CID_list_final <- left_join(x = CAS_CID_list_first, y = last_missing_SMILES, by = "CAS.Number") %>% 
  mutate(PC_title = coalesce(PC_title.x, PC_title.y),
         CID = coalesce(CID.x, CID.y),
         CanonicalSMILES = coalesce(CanonicalSMILES.x, CanonicalSMILES.y)) %>% 
  select(-c("PesticideAI_name.y", "PC_title.x", "PC_title.y", "CID.x", "CID.y", "CanonicalSMILES.x", "CanonicalSMILES.y")) %>% 
  rename(PesticideAI_name = PesticideAI_name.x) %>% 
  mutate(CID = as.integer(CID))

write.csv(CAS_CID_list_final, "data/excel_references/CAS_CID_list_final.csv", row.names = F)
#DONE

#####
# Printing a .txt file that is used to query Physchem data in OECD QSAR Toolbox
write.table(CAS_CID_list_final[,c("CAS.Number", "CanonicalSMILES")], "data/excel_references/CAS_CID_list_final.txt", col.names = F, row.names = F, quote = FALSE, sep = "\t", na = "")

# Export subsets of CAS_SMILES matched .txt documents for input into OECD QSAR Toolbox for the Toxicological query
CAS_CID_list_final[1:4000,c("CAS.Number", "CanonicalSMILES")] %>% 
  write.table("data/excel_references/cas_smiles_list[1-4k].txt", quote = FALSE, col.names = F, row.names = F, sep = "\t")

CAS_CID_list_final[4001:8000,c("CAS.Number", "CanonicalSMILES")] %>% 
  write.table("data/excel_references/cas_smiles_list[4k-8k].txt", quote = FALSE,  col.names = F, row.names = F, sep = "\t")

CAS_CID_list_final[8001:12000,c("CAS.Number", "CanonicalSMILES")] %>% 
  write.table("data/excel_references/cas_smiles_list[8k-12k].txt", quote = FALSE, col.names = F, row.names = F, sep = "\t")

CAS_CID_list_final[12001:16000,c("CAS.Number", "CanonicalSMILES")] %>% 
  write.table("data/excel_references/cas_smiles_list[12-16k].txt",  quote = FALSE, col.names = F, row.names = F, sep = "\t")

CAS_CID_list_final[16001:nrow(cas_smiles_list),c("CAS.Number", "CanonicalSMILES")] %>% 
  write.table("data/excel_references/cas_smiles_list[16k-end].txt",  quote = FALSE, col.names = F, row.names = F, sep = "\t")

