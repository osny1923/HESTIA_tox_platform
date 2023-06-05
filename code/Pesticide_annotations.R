## ----setup--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "C:/PhD Work/R data work/HESTIA_Project")
knitr::opts_chunk$set(warning = FALSE)

   #install.packages("webchem")
    library(rmarkdown)
    library(xlsx)
    library(readr)
    library(tidyverse)
    library(webchem)
   



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Make a vector of file names that are to be read in upcoming function
USEPA_filenames <- list.files(path = "data/USEPA - Comptox_categories", pattern = "Chemical List", full.names = TRUE)

# Function loading all the USEPA_COMPTOX substance information annotation
# need to load 16 files with a lot of specific stuff to read in. 
USEPA_CHEMLIST <- do.call(rbind, lapply(USEPA_filenames, function(x){
  read.csv(x, header = T, sep = ",", colClasses = c("NULL", rep("character", 7), rep("numeric", 7), "character", rep("numeric", 2))) %>%
    select(CASRN, PREFERRED.NAME, IUPAC.NAME, SMILES,  MOLECULAR.FORMULA) %>%
    mutate(Source = "USEPA_COMPTOX",
           Sub_source = paste(as.character(            # Inserting part of substance information in column "Group", based on file name
             gsub("data/USEPA - Comptox_categories/Chemical List ", "",
              gsub("-2023-01-20.csv", "",
                gsub("WIKI", "Wikipedia", 
                  gsub(",.*", "", x)))))),
          Group = paste(as.character(                  # Inserting part of substance information in column "Group", based on file name
             gsub("data/USEPA - Comptox_categories/Chemical List ", "",
              gsub("-2023-01-20.csv", "",
                gsub("WIKI", "Wikipedia", 
                  sub(".*?,", "", x)))))))
            }
    )
  ) %>% 
  pivot_wider(id_cols = c(CASRN, PREFERRED.NAME, IUPAC.NAME, SMILES,  MOLECULAR.FORMULA), names_from = Sub_source, values_from = Group) %>% 
  unnest(cols = everything()) %>% 
  filter(!grepl("NOCAS", CASRN)) 

#write.xlsx(USEPA_CHEMLIST, "data/excel_references/USEPA_CHEMLIST_2023-01_20.xlsx", sheetName = "Data1", col.names = T, append = T, showNA = F)



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Read the HESTIA list of chemials and fetch chemical information on these from CIR
    # Importing the HESTIA db and create a list of the CAS.numbers
CAS.list <- read.csv("data/excel_references/pesticideAI.csv", header = T, sep = ",") %>% 
  separate(term.id, c("cas.termstz","CAS.id"), sep = "CAS-", remove = FALSE) %>% 
  select(CAS.id, term.name) %>% 
  rename(CAS.number = CAS.id)

# Retrieving substance information, based on CAS from  the PubChem database.
  # CAS_CID_list <- get_cid(CAS.list$CAS.number, from = "xref/RN", domain = "compound", match = "first", verbose = TRUE)
  # write_rds(CAS_CID_list, file = "data/excel_references/PubChem_CID.txt")

PC_id <- read_rds("data/excel_references/PubChem_CID.txt") %>% 
        rename(CAS.Number = query) 


CAS_CID_list <- CAS.list %>% 
                rename(CAS.Number = CAS.number) %>% 
                left_join(x = .,
                          y = PC_id,
                          by = "CAS.Number")


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# US EPA Chemical Type Annotations
USEPA_Ecotox_Chem_DB <- read.delim("..\\..\\Backup HESTIA\\USEPA-ECOTOX DB\\ecotox_ascii_03_10_2022\\validation\\chemicals.txt", header = T, sep = "|") %>% 
                    mutate(CAS.Number = as.cas(cas_number))

# Merging data with HESTIA CAS numbers list
HESTIA_Comp_info <- CAS_CID_list %>% 
    left_join(x = . ,
              y = USEPA_Ecotox_Chem_DB %>% 
                  select(-cas_number),
              by = "CAS.Number") %>% 
    rename(PubChemId = cid)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Query the BCPC Compendium of Pesticide Common Names https://pesticidecompendium.bcpc.org
  # Query and code has been made as: 

# Pest_query <- bcpc_query(CAS_CID_list$CAS.Number, from = "cas", verbose = T)
# save.rds(Pest_query, "Pest_query.txt")
Pest_query <- readRDS("data/excel_references/Pest_query.txt")
#Pest_drop.na <- Pest_query[!is.na(Pest_query)] # Removing NAs list objects. not really needed here.

Pest_rbind <- as.data.frame(do.call(rbind, Pest_query)) # Flatten the nested lists into one comprehensive data frame
Pest_rbind$Query_CAS <- rownames(Pest_rbind)

rownames(Pest_rbind) = NULL # Removing rownames

Pest_rbind[1:12] <- lapply(Pest_rbind[1:12], as.character) # Coercing all data into character. Previously each cell was defined as a list.

Pest_info <- Pest_rbind %>% # selecting the cols i want to add in the final output
                rename(CAS.Number = Query_CAS) %>% 
                select(CAS.Number, formula, cname, iupac_name, activity, subactivity, source_url) 


# Merging on to the HESTIA_Comp_info
HESTIA_Comp_info_2 <- HESTIA_Comp_info %>% 
    left_join(., 
              Pest_info, 
              by = "CAS.Number") %>% 
    rename(BCPC_formula = formula,
           BCPC_cname = cname,
           BCPC_iupac_name = iupac_name,
           BCPC_activity = activity,
           BCPC_subactivity = subactivity,
           BCPC_source_url = source_url)


# Using the PubChem id to query for canonical SMILES configuration for all substances.
  # Pubchem_CanonicalSMILES_query <- webchem::pc_prop(HESTIA_Comp_info$PubChemId, properties = "CanonicalSMILES")
  # write_rds(Pubchem_CanonicalSMILES_query, "data/excel_references/Pubchem_SMILES_query.txt")
PubChem_SMILES_Query <- read_rds("data/excel_references/Pubchem_SMILES_query.txt")
  HESTIA_Comp_info_3 <- HESTIA_Comp_info_2 %>%
     mutate(PubChemId = as.integer(PubChemId)) %>%
   left_join(x = .,
             y = PubChem_SMILES_Query %>%
                 filter(!is.na(CID)) %>%
                 mutate(CID = as.integer(CID)) %>%
                 rename(PubChemId = CID),
               by = "PubChemId") %>%
     distinct(CAS.Number, .keep_all = T) 
#write_rds(HESTIA_Comp_info_3, "data/excel_references/Pesticide_info_plus_PubChem.txt")


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Checking out the new ChEMBL service in latest Webchem update (as of early 2023 ;)
  # ATC_classes <- chembl_atc_classes(verbose = T)
  # write_rds(ATC_classes, "data/excel_references/ATC_classes.txt")
ATC_classes <- read_rds("data/excel_references/ATC_classes.txt")

ChEMBL_agents <- ATC_classes %>% 
    filter(level1 %in% c("P", "J"))

# defining names to be used in the CIR query
names <- ChEMBL_agents %>% 
    pull(who_name)
  # ChEMBL_CAS_CIR <- cir_query(names, representation = "cas", resolver = "name_by_cir", verbose = T)
  # write_rds(ChEMBL_CAS_CIR, "data/excel_references/ChEMBL_ATC_classes_CAS.txt")
ChEMBL_CAS_CIR <- read_rds("data/excel_references/ChEMBL_ATC_classes_CAS.txt")

# Merging with the ATC_CLASS df
ChEMBL_data <- ChEMBL_CAS_CIR %>% 
    rename(who_name = query) %>% 
    left_join(x = .,
              y = ChEMBL_agents,
              by = "who_name") %>% 
  rename(CAS.Number = cas)

# Not really sure what this code is doing...
        #   #ChEMBL_ATC_cid_names <- get_cid(names, from = "name", domain = "substance", match = "all", verbose = T)
        #   #write_rds(ChEMBL_ATC_cid_names, "ChEMBL_ATC_cid_names.txt")
        # ChEMBL_ATC_cid_names <- read_rds("data/excel_references/ChEMBL_ATC_cid_names.txt")
        # ChEMBL_ATC_cid_names_2 <- ChEMBL_ATC_cid_names %>%
        #     filter(!is.na(cid)) %>% 
        #   distinct(cid, .keep_all = T) %>% 
        #     pull(cid)
        # 
        # #ChEMBL_ATC_cid_CAS <- pc_sect(ChEMBL_ATC_cid_names_2[1:10], section = "CAS", domain = "substance", verbose = T)
        # # write_rds(ChEMBL_ATC_cid_CAS, "ChEMBL_ATC_cid_CAS.txt")
        # ChEMBL_ATC_cid_CAS <- read_rds("data/excel_references/ChEMBL_ATC_cid_CAS.txt")
        # 
        # #If the code is working, then i can just load this file below!
        # ChEMBL_ATC_cid_CAS <- read_rds("data/excel_references/ChEMBL_ATC_cid_CAS.txt") # This is the result. All available CAS data from the CID (PubChem identifiers) from the list of WHO ATC codes. 



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# dependency:`HESTIA_Comp_info_3` object created in 

#HESTIA_ChEBI_name <- get_chebiid(HESTIA_Comp_info_3$CAS.Number, from = "registry numbers", match = "best", stars = "all", verbose = T) # OK, Done
#write_rds(HESTIA_ChEBI_name, "HESTIA_Chem_Chebi_IDs.txt")
HESTIA_ChEBI_name <- read_rds("data/excel_references/HESTIA_Chem_Chebi_IDs.txt")

#HESTIA_ChEMBL_all_list <- chebi_comp_entity(HESTIA_ChEBI_name$chebiid, verbose = T) #OK, Done
#write_rds(HESTIA_ChEMBL_all_list, "HESTIA_Chem_all_ChEBI_data.txt")
HESTIA_ChEMBL_all_list <- read_rds("data/excel_references/HESTIA_Chem_all_ChEBI_data.txt")

# Starting with the wrangling of lists:
HESTIA_ChEMBL_all_list_drop.na <- HESTIA_ChEMBL_all_list[!is.na(HESTIA_ChEMBL_all_list)] # Removing NAs list objects
#HESTIA_ChEMBL_all_list__flat <- as.data.frame(do.call(rbind, HESTIA_ChEMBL_all_list_drop.na))
HESTIA_ChEMBL_all_list_parents <- do.call("rbind", lapply(HESTIA_ChEMBL_all_list_drop.na, "[[", 10) ) # Flattening the "Parents" data frame with substance information
HESTIA_ChEMBL_all_list_parents$ChEBId <- rownames(HESTIA_ChEMBL_all_list_parents) # Creating a substance name column
rownames(HESTIA_ChEMBL_all_list_parents) = NULL # Removing rownames

#########
# Investigating if there is a logic to the naming of IDs for "has role" stuff
# list of agrochemicals;
HESTIA_ChEMBL_all_list_parents %>% 
    #filter(type == "has role") %>% 
   separate(col = ChEBId, into = c("CheBId", "row.nr"), sep = "\\.", remove = T)

# list of pesticides;
HESTIA_ChEMBL_definitions <- HESTIA_ChEMBL_all_list_parents %>% 
    #filter(type == "has role") %>% 
    separate(col = ChEBId, into = c("CheBId", "row.nr"), sep = "\\.", remove = T) %>% 
    mutate(type = case_when(chebiId %in% c("CHEBI:38805", "CHEBI:22153", "CHEBI:39215", "CHEBI:22583", "CHEBI:51076", "CHEBI:86328", "CHEBI:33289",
                                           "CHEBI:23092", "CHEBI:39276", "CHEBI:24527", "CHEBI:24852", "CHEBI:33904", "CHEBI:25491", "CHEBI:25943",
                                           "CHEBI:167183", "CHEBI:136643", "CHEBI:33288", "CHEBI:25944", "CHEBI:38601", "CHEBI:38656", "CHEBI:9442",
                                           "CHEBI:73333", "CHEBI:38956"
                                           ) ~ "Pesticide",
                            
                            chebiId %in%  c("CHEBI:33286", "CHEBI:71692"
                                            ) ~ "agrochemical",
                            
                            chebiId %in%  c("CHEBI:33281", "CHEBI:9337", "CHEBI:87230", "CHEBI:87228", "CHEBI:86324", "CHEBI:86322", "CHEBI:26179",  
                                            "CHEBI:26177", "CHEBI:25807", "CHEBI:25105", "CHEBI:22507", "CHEBI:156449", "CHEBI:15369", "CHEBI:23066",
                                            "CHEBI:23765", "CHEBI:35358", "CHEBI:23007", "CHEBI:25558", "CHEBI:25605", "CHEBI:25903", "CHEBI:27933",
                                            "CHEBI:39217", "CHEBI:49318", "CHEBI:49322", "CHEBI:86478", "CHEBI:87211", "CHEBI:88225", "CHEBI:36043"
                                            ) ~ "antimicrobial agent",
                            
                            chebiId %in%  c("CHEBI:33282", "CHEBI:36047") ~ "antibacterial agent",
                            chebiId %in%  c("CHEBI:22587", "CHEBI:36044") ~ "antiviral agent",
                            chebiId %in%  c("CHEBI:39094", "CHEBI:35444"
                                            ) ~ "nematocide",
                            
                            chebiId %in%  c("CHEBI:33287", "CHEBI:16199"
                                            ) ~ "fertilizer",
                            
                            chebiId %in%  c("CHEBI:82599", "CHEBI:25540", "CHEBI:25705", "CHEBI:26409", "CHEBI:38455", "CHEBI:38611", "CHEBI:38804",
                                            "CHEBI:39116", "CHEBI:39190", "CHEBI:39213", "CHEBI:136644", "CHEBI:22917", "CHEBI:23457", "CHEBI:25708",
                                            "CHEBI:25715", "CHEBI:38060", "CHEBI:38461", "CHEBI:38488", "CHEBI:38494", "CHEBI:38577", "CHEBI:39090",
                                            "CHEBI:39167", "CHEBI:39191", "CHEBI:39208", "CHEBI:39210", "CHEBI:39277", "CHEBI:39295", "CHEBI:39351",
                                            "CHEBI:39398", "CHEBI:39415"
                                            ) ~ "insecticide",
                            
                            chebiId %in%  c("CHEBI:60575", "CHEBI:60575", "CHEBI:133673", "CHEBI:15930", "CHEBI:137504"
                                            ) ~ "herbicide",
                            
                            chebiId %in%  c("CHEBI:36053", "CHEBI:38489", "CHEBI:38602", "CHEBI:38612", "CHEBI:38657", "CHEBI:38806", "CHEBI:38820",
                                            "CHEBI:39219", "CHEBI:39259", "CHEBI:39292", "CHEBI:39296", "CHEBI:39298", "CHEBI:39301", "CHEBI:39318",
                                            "CHEBI:39363", "CHEBI:39369", "CHEBI:39412"
                                            ) ~ "acaricide",
                            
                            chebiId %in%  c("CHEBI:24127", "CHEBI:87208", "CHEBI:87198", "CHEBI:87197", "CHEBI:87195", "CHEBI:87114", "CHEBI:87100",
                                            "CHEBI:87068", "CHEBI:87067", "CHEBI:87066", "CHEBI:87064", "CHEBI:87061", "CHEBI:87039", "CHEBI:87038",
                                            "CHEBI:87071", "CHEBI:87069", "CHEBI:87035", "CHEBI:87034", "CHEBI:87019", "CHEBI:87018", "CHEBI:87015",
                                            "CHEBI:60600", "CHEBI:87113", "CHEBI:86327", "CHEBI:86494", "CHEBI:86488", "CHEBI:86487", "CHEBI:86484",
                                            "CHEBI:86482", "CHEBI:35718", "CHEBI:38819", "CHEBI:86417", "CHEBI:86440", "CHEBI:86478", "CHEBI:86485",
                                            "CHEBI:87013", "CHEBI:87021", "CHEBI:87022", "CHEBI:87025", "CHEBI:87026", "CHEBI:87027", "CHEBI:87029",
                                            "CHEBI:87036", "CHEBI:87037", "CHEBI:87101", "CHEBI:87110", "CHEBI:87127", "CHEBI:87128", "CHEBI:87134",
                                            "CHEBI:87135", "CHEBI:87207", "CHEBI:136645" 
                                            ) ~ "fungicide",
                            
                           TRUE ~ as.character(NA))) 

#HESTIA_ChEMBL_definitions %>% 
#    #filter(chebiId ==  "CHEBI:33282") %>% 
#    #separate(col = ChEBId, into = c("CheBId", "row.nr"), sep = "\\.", remove = T) %>% 
#    filter(#row.nr == "4",
#           !grepl("metabolite", chebiName),
#           grepl("anti", chebiName),
#           is.na(type)
#           ) %>% 
#    distinct(chebiId, .keep_all = T) %>% 
#    arrange(chebiId) #%>% 
#    pull(chebiId)

ChEBI_annot <-  HESTIA_ChEMBL_definitions %>% 
                distinct(CheBId, .keep_all = T) %>% 
                select(chebiName, CheBId) %>%
                left_join(x =.,
                          y = HESTIA_ChEMBL_definitions %>% 
                                filter(!grepl(paste(c("metabolite", "inhibitor"), collapse = "|"), chebiName),
                                       !is.na(type)) %>% 
                                select(-c(row.nr, status, cyclicRelationship, chebiName)) %>% 
                                pivot_wider(values_from = type, names_from = chebiId) %>% 
                                unnest(cols = everything() ) %>% 
                                unite("Use_as", 2:160, sep = ";", na.rm = T),
                          by = "CheBId") %>% 
                rename(chebiid = CheBId) %>% 
# a separate-and-clean exercise of the CheBI_annot object
    separate(Use_as, into = c("Use_as", "Use_as.1", "Use_as.2", "Use_as.3", "Use_as.4", "Use_as.5", "Use_as.6"), sep = ";", remove = TRUE) %>% 
    mutate(Use_as.1 = case_when(Use_as.1 == Use_as ~ as.character(NA),
              TRUE ~ Use_as.1),
           Use_as.2 = case_when(Use_as.2 == Use_as | 
                                  Use_as.2 == Use_as.1 ~ as.character(NA),
              TRUE ~ Use_as.2),
           Use_as.3 = case_when(Use_as.3 == Use_as | 
                                Use_as.3 == Use_as.1 |
                                Use_as.3 == Use_as.2 ~ as.character(NA),
              TRUE ~ Use_as.3),
           Use_as.4 = case_when(Use_as.4 == Use_as | 
                                Use_as.4 == Use_as.1 |
                                Use_as.4 == Use_as.2 |
                                Use_as.4 == Use_as.3 ~ as.character(NA),
              TRUE ~ Use_as.4),
           Use_as.5 = case_when(Use_as.5 == Use_as | 
                                Use_as.5 == Use_as.1 |
                                Use_as.5 == Use_as.2 |
                                Use_as.5 == Use_as.3 |
                                Use_as.5 == Use_as.4 ~ as.character(NA),
              TRUE ~ Use_as.5),
           Use_as.6 = case_when(Use_as.6 == Use_as | 
                                Use_as.6 == Use_as.1 |
                                Use_as.6 == Use_as.2 |
                                Use_as.6 == Use_as.3 |
                                Use_as.6 == Use_as.4 |
                                Use_as.6 == Use_as.5 ~ as.character(NA),
              TRUE ~ Use_as.6)
           ) %>% 
  unite("Use_as", 3:9, sep = ";", na.rm = T)

ChEBI_annot$Use_as <- sapply(strsplit(as.character(ChEBI_annot$Use_as), ';\\s*'), 
                     function(x) toString(sort(x, decreasing = T)))
      
# Flattening the "regnumbers" data frame with Registry numbers
# But it's "too" flat! I have lost the ChEBI ID nr.
HESTIA_ChEMBL_all_list_regnumbers <- do.call("bind_rows", lapply(HESTIA_ChEMBL_all_list_drop.na, "[[", 7) ) %>%  # The dplyr "bind_rows" function inserts NAs to force together dfs with uneven amount of columns!! 
                                      filter(type == "CAS Registry Number") %>% 
                                      distinct(data, .keep_all = T)

as.data.frame(do.call(rbind, HESTIA_ChEMBL_all_list_drop.na))

reg_nr <- as.data.frame(do.call(rbind, HESTIA_ChEMBL_all_list_drop.na)) %>% 
    select(regnumbers)


reg_nr <- do.call("rbind", lapply(reg_nr$regnumbers, "[[", 1) )
reg_nr <- as.data.frame(reg_nr)
ChEBI_ID_CAS <- reg_nr %>% 
        mutate(across(everything(), ~ case_when(is.cas(.x) ~ .x, 
                                                   TRUE ~ as.character(NA)))) %>% 
        rename(CAS.Number = V1) %>% 
        mutate(CAS.Number = case_when(is.na(CAS.Number) ~ V2, 
                                      TRUE ~ CAS.Number),
               CAS.Number = case_when(is.na(CAS.Number) ~ V3, 
                                      TRUE ~ CAS.Number),
               CAS.Number = case_when(is.na(CAS.Number) ~ V4, 
                                      TRUE ~ CAS.Number),
               CAS.Number = case_when(is.na(CAS.Number) ~ V5, 
                                      TRUE ~ CAS.Number),
               CAS.Number = case_when(is.na(CAS.Number) ~ V6, 
                                      TRUE ~ CAS.Number),
               CAS.Number = case_when(is.na(CAS.Number) ~ V7, 
                                      TRUE ~ CAS.Number)
               )
  
ChEBI_ID_CAS$ChebId <- rownames(ChEBI_ID_CAS) 
rownames(ChEBI_ID_CAS) = NULL # Removing rownames

# Adding lot's of inchl and stuff... is it important?? NO. Not now..
p <- as.data.frame(do.call(rbind, HESTIA_ChEMBL_all_list_drop.na)) %>% 
    select(properties)
p$ChEBId <- rownames(p)
rownames(p) = NULL # Removing rownames

pq <- do.call(rbind, p$properties) %>% 
      select(chebiid, chebiasciiname, definition, smiles, inchi, inchikey, mass)
rownames(pq) = NULL # Removing rownames


#####
# MERGING ALL THE ChEBI DATA#
#####
ChEBI_annot_merge <- ChEBI_annot %>% 
  left_join(x = ., 
            y = ChEBI_ID_CAS %>% 
                select(ChebId, CAS.Number) %>%
                mutate(ChebId = gsub("\\.", ":", ChebId)) %>% 
                rename(chebiid = ChebId), 
            by = "chebiid") %>% 
  left_join(x = ., 
            y = pq, 
            by = "chebiid")
  



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

HESTIA_Comp_info_4 <- HESTIA_Comp_info_3 %>% 
    left_join(x = ., 
              y = ChEBI_annot_merge, 
              by = "CAS.Number") %>% 
    left_join(x = ., 
              y = ChEMBL_data %>% 
                  filter(!is.na(CAS.Number)),
              by = "CAS.Number") %>% 
    mutate(dtxsid = gsub("None", as.character(NA), dtxsid),
           ecotox_group = gsub("None", as.character(NA), ecotox_group),
           across(starts_with("BCPC"), ~ case_when(.x == "NA" ~ as.character(NA),
                                                 TRUE ~ .x))) %>% 
    rename_at(vars(matches("level")), ~ paste0("ATC_", .)) %>% 
    rename(WHO_ATC_name = who_name,
           USEPA_ecotox_group = ecotox_group) %>% 
    select(c(
"CAS.Number", "CanonicalSMILES", "smiles", "dtxsid", # Chem_id's
"term.name", "chebiasciiname", "chemical_name", "BCPC_cname", "WHO_ATC_name", "BCPC_iupac_name", # "chebiName",  # Names
"mass", "BCPC_formula", # Chem_props
"ATC_level1_description", "USEPA_ecotox_group",  "BCPC_activity", "BCPC_subactivity", "Use_as", # Pest IINFO
"PubChemId", "chebiid",  "definition", "inchi", "inchikey", # ID_codes
"BCPC_source_url",  "ATC_level1",  "ATC_level2", "ATC_level2_description", "ATC_level3", "ATC_level3_description", "ATC_level4", "ATC_level4_description", "ATC_level5" 
    )) %>% 
    distinct(CAS.Number, .keep_all = T) 



## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Halogenated: contains c(F, Cl, Br, I)

# Investigating smiles that contain C in the smiles (Cu, Cr, Co),yet excluding such elements with a density <5g/cm3 and Carbon (C). !
# Creatng a vector of substances that I can filter against.
# This vector is based on 3 operations: 
  # 1. Looking up all SMILES containing elements with the letter "C" and filtering out all possible carbon components
  # 2. Looking up all elements not containing carbon whatsoever
  # 3. Merging both lists as CAS.numbers to filter out.
# 332 substances are identified. Although substances without SMILES configuration will be remain in the data set!

Inorganic_vector_x <- rbind(HESTIA_Comp_info_4 %>% 
                              distinct(CAS.Number, .keep_all = T) %>% # this filter detects all SMILES 
                              filter(str_detect(CanonicalSMILES, paste(c("Cd", "Cs", "Co", "Cs", "Cr", "Cu", "Cl", "Ca"), collapse = "|")), # Making sure to select all substances with elements containing a "C" or "c" in SMILES
                                     !str_detect(CanonicalSMILES, paste(c("CC", "c","C+[0-9]", "C+\\(", "C+\\)", "C+\\[", "C+\\]", "C+\\#", "C+\\="), collapse="|"))),
                          HESTIA_Comp_info_4 %>% # This filter detects all SMILES without Carbon.
                              distinct(CAS.Number, .keep_all = T) %>% 
                              filter(!str_detect(CanonicalSMILES, paste(c("c","C"), collapse = "|")))) %>% 
                          pull(CAS.Number)

HM_vector_x <- HESTIA_Comp_info_4 %>% # Just as in the operation above with an inorganic annotation, This creates a vector of CAS.Numbers to annotate Heavy metals presence in the smiles configuration
            distinct(CAS.Number, .keep_all = T) %>% 
            filter(str_detect(CanonicalSMILES, paste(read.xlsx("data/excel_references/Elements_list.xlsx", sheetName = "Sheet1", header = T) %>% # Using a list of elements (From Wikipedia) with annotations of Elements and their density
                                                      mutate(Density = as.numeric(Density)) %>% 
                                                      filter(Density >=5) %>% 
                                                      pull(Abbr), 
                                                      collapse = "|"))) %>% 
            pull(CAS.Number)

Halogenated_vector_x <- HESTIA_Comp_info_4 %>% # Just as in the operation above with an inorganic annotation, This creates a vector of CAS.Numbers to annotate Halogens presence in the SMILES configuration
                        distinct(CAS.Number, .keep_all = T) %>% 
                        filter(str_detect(CanonicalSMILES, paste(c("F", "Cl", "Br", "I"), collapse = "|"))) %>% 
                        pull(CAS.Number)


## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
HESTIA_Comp_info_5 <-  HESTIA_Comp_info_4 %>% 
                                mutate(Substance_type = case_when(CAS.Number %in% Inorganic_vector_x ~ "Inorganic",
                                                                  BCPC_subactivity == "inorganic" ~ "Inorganic",
                                                                  grepl("mixt.", term.name) ~ "Chemical mixture",
                                                                  is.na(CanonicalSMILES) ~ "Unknown",
                                                                TRUE ~ "Organic"),
                                       `Heavy Metals` = case_when(CAS.Number %in% HM_vector_x ~ "1", 
                                                                  TRUE ~ "0"),
                                       Halogenated = case_when(CAS.Number %in% Halogenated_vector_x ~ "1", 
                                                                  TRUE ~ "0" ))

HESTIA_Comp_info_6 <- HESTIA_Comp_info_5 %>% 
    mutate_all(~na_if(., '')) %>% 
    select(CAS.Number, CanonicalSMILES, term.name, Use_as, ATC_level1_description, USEPA_ecotox_group, BCPC_activity, BCPC_subactivity, ATC_level1, ATC_level1_description, ATC_level4_description, Substance_type, `Heavy Metals`, Halogenated, definition) %>% 
    left_join(x = ., 
              y = USEPA_CHEMLIST %>% 
                  rename(CAS.Number = CASRN) %>% 
                  select(-c("PREFERRED.NAME", "IUPAC.NAME", "SMILES", "MOLECULAR.FORMULA")),
              by = "CAS.Number") %>% 
    unite("Subgroup", c("Use_as", "ATC_level1_description", "USEPA_ecotox_group", 
                        "BCPC_activity", "BCPC_subactivity", "DRUGBANK", "EPAPCS", 
                        "HEALTHY_BUILDING_NETWORK", "NORMAN_ITN", "OPPIN", "PPDB", 
                        "USEPA", "Wikipedia"), sep = "; ", na.rm = T, remove = F) %>%
    mutate(Group = case_when(
                             grepl(paste(c("Antimicrobial", "antibacterial", "antibiotic"), collapse = "|"), Subgroup, fixed = F) ~ "Antibiotic",
                             grepl("Pesticide", Subgroup) ~ "Pesticide",
                             !is.na(BCPC_activity) ~ "Pesticide",
                             grepl("PPCP", USEPA_ecotox_group, fixed = F) ~ "PPCP",
                             grepl("Pharmaceutical", Subgroup) ~ "Pharmaceutical",
                             ATC_level1 == "P" ~ "Pesticide",
                             ATC_level1 == "J" ~ "Antibiotic",
                             TRUE ~ case_when(Substance_type == "Organic" ~ "Other organic chemicals",
                                              is.na(CanonicalSMILES) ~ "Unknown",
                                              TRUE ~ "Other inorganic chemicals")),
           Use_as = gsub("\\,", ";", Use_as)) %>%
      rowwise() %>% # Removing the Group category from Subgroup.
      mutate(Subgroup = case_when(grepl(paste((Group)), Subgroup) ~ gsub(paste(Group, ";", sep=""), "", Subgroup),
                                  TRUE ~ Subgroup)) %>% 
      mutate(Subgroup = case_when(grepl(paste((Group)), Subgroup) ~ gsub(paste(Group), "", Subgroup),
                                  TRUE ~ Subgroup),
             Subgroup = sub(" $", "", Subgroup),
             Subgroup = sub(";$", "", Subgroup),
             Subgroup = gsub(",", "", Subgroup),
             Subgroup = sub(" \\(s\\)", "", Subgroup),
             Subgroup = sub("^ ", "", Subgroup)
             ) %>% 
      ungroup() %>% 
      mutate_all(~na_if(., '')) %>% 
      rename(ChEBI_DB = Use_as) %>% 
      select(CAS.Number, term.name, CanonicalSMILES, Group, Subgroup, Substance_type, `Heavy Metals`, Halogenated, definition, ChEBI_DB, ATC_level1_description, 
             USEPA_ecotox_group, BCPC_activity, BCPC_subactivity, ATC_level1, ATC_level1_description, ATC_level4_description, DRUGBANK, EPAPCS, 
                        HEALTHY_BUILDING_NETWORK, NORMAN_ITN, OPPIN, PPDB, USEPA, Wikipedia) %>% 
      mutate(Subgroup = tolower(Subgroup)) %>% 
      distinct(CAS.Number, .keep_all = T)
                      
#write.xlsx(as.data.frame(HESTIA_Comp_info_6), "results/HESTIA_chem_prop_list_full.xlsx", sheetName = "Full_DB", col.names = T, row.names = F, append = T, showNA = F)
#write.csv(HESTIA_Comp_info_6, "results/HESTIA_chem_prop_list_full.csv")


HESTIA_Comp_info_7 <- HESTIA_Comp_info_6 %>% 
    select(CAS.Number, term.name, CanonicalSMILES, Group, Subgroup, Substance_type, `Heavy Metals`, Halogenated, definition)
    
#xlsx::write.xlsx(as.data.frame(HESTIA_Comp_info_7), "results/HESTIA_chem_prop_list_full.xlsx", sheetName = "Slim_DB", col.names = T, row.names = F, append = T, showNA = F)
#write.csv(HESTIA_Comp_info_7, "results/HESTIA_chem_list_slim.csv")


