library(tidyverse)
library(webchem)

latest_list <- read.csv("results/HESTIA_chem_list_slim.csv") %>%  as.tibble()
early_list <- read.csv("data/HESTIA_chem_list_slim_deprecated.csv") %>%  as.tibble()

Look_into_list <- latest_list %>% 
  left_join(x = ., 
            y = early_list %>% 
              select(CAS.Number, Group, CanonicalSMILES) %>% 
              rename(Group_OLD = Group),
            by = "CAS.Number") %>% 
  mutate(MATCH = case_when(Group == Group_OLD ~ "MATCH",
                           TRUE ~ as.character(NA)))
Checkup_list <- Look_into_list %>% 
  filter(is.na(MATCH)) %>% 
  select(CAS.Number, CanonicalSMILES.x, Group, Group_OLD, MATCH, CanonicalSMILES.y, PesticideAI_name)

# write.csv(Checkup_list %>% 
#   count(Group, Group_OLD), "data/excel_references/Group_comparison.csv", row.names = F)

# dependency on `Pesticide_anotations.Rmd`
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

# Creating a df with the true matches! No double matches! Great!
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

cid_to_cas_query <- pc_sect(last_remaining_CAS_SMILES$cid, section = "CAS", domain = "compound", verbose = T)

cid_to_cas_query_CAS_repo <- cid_to_cas_query %>% 
filter(SourceName == "CAS Common Chemistry") 


joined_missing_df <- left_join(x = cid_to_cas_query_CAS_repo %>% 
              rename(cid = CID) %>% 
              mutate(cid = as.integer(cid)), 
            y = last_remaining_CAS_SMILES, 
            by = "cid") %>% 
  mutate(MATCH = case_when(CAS.number == Result ~ "MATCH",
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


# Last effort to secure the missing 2259 CAS numbers without SMILES
#Dependency on `Pesticide_annotations.Rmd`
Missing_list_1 <- CAS_CID_list_final %>% 
  filter(is.na(CanonicalSMILES))

# Missing_list_1_cid <- get_cid(Missing_list_1$CAS.Number, from = "xref/RN", domain = "compound", match = "all", verbose = T)
# write.csv(Missing_list_1_cid, "data/excel_references/Missing_list_1_cid.csv", row.names = F)
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

Missing_list_4_cid <- pc_sect(Missing_list_3_cid$cid, section = "CAS", domain = "compound", verbose = T)

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
  select(names(CAS_CID_list))

write.csv(Missing_list_5_MATCH, "data/excel_references/Missing_list_5_MATCH.csv", row.names = F)


Missing_list_5_NOT_MATCH <- Missing_list_5 %>% 
  filter(CAS.Number != Result)

Missing_list_5_refined_NOT_MATCH <- Missing_list_5_NOT_MATCH %>% 
  filter(!CAS.Number %in% Missing_list_5_MATCH$CAS.Number) %>% 
  distinct(CAS.Number, CanonicalSMILES.x, .keep_all = T) %>% 
  mutate(MATCh = "NOT_MATCH") %>% 
  rename(CanonicalSMILES = CanonicalSMILES.x)
  

