
## Substance Properties Annotations
Creating an annotation on substances being 1. Inorganic/Organic, 2. containing Heavy metals, and 3. being Halogenated.
```{r}
# Halogenated: contains c(F, Cl, Br, I)
# Investigating smiles that contain C in the smiles (Cu, Cr, Co) yet excluding such elements with a density <5g/cm3 and Carbon (C). !
# Creatng a vector of substances that I can filter against.
# This vector is based on 3 operations: 
# 1. Looking up all SMILES containing elements with the letter "C" and filtering out all possible carbon components
# 2. Looking up all elements not containing carbon whatsoever
# 3. Merging both lists as CAS.numbers to filter out.
# 332 substances are identified. Although substances without SMILES configuration will be remain in the data set!

Inorganic_vector <- rbind(NEW_PHYSCHEM %>% 
                            filter(str_detect(SMILES, paste(c("Cd", "Cs", "Co", "Cs", "Cr", "Cu", "Cl", "Ca"), collapse = "|")), 
                                   !str_detect(SMILES, paste(c("CC", "c","C+[0-9]", "C+\\(", "C+\\)", "C+\\[", "C+\\]", "C+\\#", "C+\\="), collapse="|"))),
                          NEW_PHYSCHEM %>% 
                            filter(!str_detect(SMILES, paste(c("c","C"), collapse = "|")))) %>% 
  pull(CAS.Number)

HM_vector <- NEW_PHYSCHEM %>% # Just as in the operation above with an inorganic annotation, This creates a vector of CAS.Numbers to annotate Heavy metals presence in the SMILES configuration
  filter(str_detect(SMILES, paste(read.xlsx("../data/excel_references/Elements_list.xlsx", sheetName = "Sheet1", header = T) %>% # Using a list of elements (From Wikipedia) with annotations of Elements and their density
                                    # only looking at elements >5 g/cm3 (Ca and Cl are not counted here)
                                    mutate(Density = as.numeric(Density)) %>% 
                                    filter(Density >=5) %>% 
                                    pull(Abbr), 
                                  collapse = "|"))) %>% 
  pull(CAS.Number)

Halogenated_vector <- NEW_PHYSCHEM %>% # Just as in the operation above with an inorganic annotation, This creates a vector of CAS.Numbers to annotate Halogens presence in the SMILES configuration
  filter(str_detect(SMILES, paste(c("F", "Cl", "Br", "I"), collapse = "|"))) %>% 
  pull(CAS.Number)

NEW_PHYSCHEM <-  NEW_PHYSCHEM %>% 
  mutate(Substance_type = case_when(CAS.Number %in% Inorganic_vector ~ "Inorganic",
                                    TRUE ~ "Organic"),
         `Heavy Metals` = case_when(CAS.Number %in% HM_vector ~ "1", 
                                    TRUE ~ "0"),
         Halogenated = case_when(CAS.Number %in% Halogenated_vector ~ "1", 
                                 TRUE ~ "0" ))
rm(list = c("Inorganic_vector", "HM_vector", "Halogenated_vector"))
```

### USEPA Comptox categorized lists

Copied Chunk from Pesticide_annotations.Rmd - "Compiling lists of chemicals from the USEPA Comptox categorized lists"
Reading here to decrease dependency on the above mentioned file, which clutters the work space
This dataset comes from [USEPA Comptox](https://comptox.epa.gov/dashboard/chemical-lists) where free search is possible.
Search terms are:  
  * Antibiotic
* Pesticide
* PPCP
* Pharmaceutical
Each match is downloaded, imported and the use-type is ascribed depending on list name (e.g. WIKI_Antibiotics = Antibiotic)
CAS.Numbers are later joined to the 

```{r}
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

```


### Pesticide information
collected data on potential pesticides from [BCPC](https://pesticidecompendium.bcpc.org) using `webchem::bcpc_query()` and substance CAS.Number

```{r}
# Query the BCPC Compendium of Pesticide Common Names https://pesticidecompendium.bcpc.org
# Query and code has been made as: 
# Pest_query <- bcpc_query(CAS_CID_list$CAS.Number, from = "cas", verbose = T)
# save.rds(Pest_query, "Pest_query.txt")
Pest_query <- readRDS("data/excel_references/Pest_query.txt")
Pest_drop.na <- Pest_query[!is.na(Pest_query)] # Removing NAs list objects

Pest_rbind <- as.data.frame(do.call(rbind, Pest_drop.na)) # Flatten the nested lists into one comprehensive data frame
Pest_rbind$Query_CAS <- rownames(Pest_rbind)

rownames(Pest_rbind) = NULL # Removing rownames

Pest_rbind[1:12] <- lapply(Pest_rbind[1:12], as.character) # Coercing all data into character. Previously each cell was defined as a list.

Pest_info <- Pest_rbind %>% # selecting the cols i want to add in the final output
  rename(CAS.Number = Query_CAS) %>% 
  select(CAS.Number, cname, iupac_name, activity, subactivity) 

#joining the pesticide activity columns with the data frame
NEW_PHYSCHEM <- NEW_PHYSCHEM %>% 
  left_join(., 
            Pest_info,
            by = "CAS.Number")

rm(list = c("Pest_query", "Pest_drop.na", "Pest_rbind", "Pest_info"))

```

### Substance names and identifiers for both HESTIA and Envirotox
Merging all of the various data queries into one(two) outputs
```{r}
# Selected identifiers for the Envirotox database
Envirotox_subst <- readxl::read_xlsx(path = "..\\..\\Backup HESTIA\\EnviroTox DB\\envirotox_DB_20230110081208.xlsx", sheet = "substance", col_names = T, col_types = c("skip", "guess", rep("text", times = 2), rep("skip", times = 19))) %>% 
  mutate(CAS.Number = as.cas(`original CAS`)) %>% 
  select(CAS.Number, `Chemical name`, `Canonical SMILES`) 

# Making a list of CAS.Numbers to select whether a substance is Organic or inorganic
Inorganic_vector_Enviro <- rbind(Envirotox_subst %>% 
                                   filter(str_detect(`Canonical SMILES`, paste(c("Cd", "Cs", "Co", "Cs", "Cr", "Cu", "Cl", "Ca"), collapse = "|")), 
                                          !str_detect(`Canonical SMILES`, paste(c("CC", "c","C+[0-9]", "C+\\(", "C+\\)", "C+\\[", "C+\\]", "C+\\#", "C+\\="), collapse = "|"))),
                                 Envirotox_subst %>% 
                                   filter(!str_detect(`Canonical SMILES`, paste(c("c","C"), collapse = "|")))) %>% 
  pull(CAS.Number)

HESTIA_Enviro_subst_type_annot <- rbind(
  # loading the HESTIA substance type annotations
  readxl::read_xlsx(path = "data/HESTIA_chem_prop_list_full.xlsx", sheet = "Slim_DB", col_names = T, col_types = c("guess", rep("text", times = 4), rep("skip", times = 3), "text")), 
  # loading and wrangling substance type annotations for the envirotox substances that are missing from the HESTIA list
  # Keep in mind that this whole operations is to have physchem data annotations for the envirotox-unique substances. it's messy and repetitive, but works!
  Envirotox_subst %>% 
    left_join(x = .,
              y = USEPA_CHEMLIST %>% 
                rename(CAS.Number = CASRN),
              by = "CAS.Number") %>% 
    # Using the Above created "inorganic_vector_Enviro" to define inorganic/organic substance types
    mutate(Substance_type = case_when(CAS.Number %in% Inorganic_vector_Enviro ~ "Inorganic",
                                      TRUE ~ "Organic")) %>% 
    unite("Subgroup", c("DRUGBANK", "EPAPCS", "HEALTHY_BUILDING_NETWORK", "NORMAN_ITN", "OPPIN", "PPDB", "USEPA", "Wikipedia"), sep = "; ", na.rm = T, remove = T)  %>% 
    mutate(Group = case_when(
      grepl(paste(c("Antimicrobial", "antibacterial", "antibiotic", "Antibiotics"), collapse = "|"), Subgroup, fixed = F) ~ "Antibiotic",
      grepl("Pesticide", Subgroup) ~ "Pesticide",
      grepl("Pharmaceutical", Subgroup) ~ "Pharmaceutical",
      TRUE ~ case_when(Substance_type == "Organic" ~ "Other organic chemicals",
                       is.na(`Canonical SMILES`) ~ "Unknown",
                       TRUE ~ "Other inorganic chemicals"))
    ) %>%
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
    # Warning comes in here, since many term.names are removed. Only the first is selected in the `separate` function.
    separate(col = "Chemical name", into = "term.name", sep = ";") %>% 
    mutate(term.name = case_when(!is.na(PREFERRED.NAME) ~ PREFERRED.NAME,
                                 TRUE ~ term.name),
           definition = as.character(NA)) %>% 
    rename(CanonicalSMILES = `Canonical SMILES`) %>% 
    select(CAS.Number, term.name, CanonicalSMILES, Group, Subgroup, definition)
) %>% 
  distinct(CAS.Number, .keep_all = T)

NEW_PHYSCHEM <- NEW_PHYSCHEM %>% 
  left_join(x = ., 
            y = HESTIA_Enviro_subst_type_annot %>% 
              select(-CanonicalSMILES),
            by = "CAS.Number")
```

