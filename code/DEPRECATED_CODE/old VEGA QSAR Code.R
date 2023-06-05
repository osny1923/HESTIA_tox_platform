
#### DEPREATED CODE #####
Wrangling the VEGA QSAR data

VEGA models are available for Chronic NOEC data for three taxonomic groups, Fish (vertebrate), D. magna (Invertebrate/Crustacean) and Algae (well.... Algae.)
Following the methodology above, by Aurisano et al., 2019, the coefficient to transform Chronic NOEC data into EC10eq is "0.6" (Effect concentration/0.6). 
The amount of species data to construct SSD models is not sufficient (n species >5), but the demand of >=3 taxonomic groups is fulfilled. 

```{r VEGA_wrangle}
ID_index_HESTIA_dataset <- read.csv("../data/2022-11-07/ID_index_HESTIA_dataset.csv")
# I have to load the original QSAR_toolbox query, as this dataset makes the basis for VEGA QSAR SMILES-configurations.
# By importing CAS.Numbers and SMILES-configurations and adding a unique ID nr to each substance, i can trace and match the ambiguous SMILES configurations exported by VEGA.
#names(VEGA_QSARS)
# A list of possible QSARS created using VEGA software.
# VEGA has, however, changed the SMILES configuration during the operation. Adding a "unique" identifier to each row to avoid having copies of rows. (Both list follow the same chronological order)
VEGA_QSARS <- 
  bind_rows(
    read.delim("../data/VEGA_REPORT/2022-11-28 1-4k/report_summary.txt", sep = "\t", skip = 6, na.strings = c("[Error]-", "", "NA")),
    read.delim("../data/VEGA_REPORT/2022-11-28 4-8k/report_summary.txt", sep = "\t", skip = 6, na.strings = c("[Error]-", "", "NA")),
    read.delim("../data/VEGA_REPORT/2022-11-28 8-12k/report_summary.txt", sep = "\t", skip = 6, na.strings = c("[Error]-", "", "NA")),
    read.delim("../data/VEGA_REPORT/2022-11-28 12-16k/report_summary.txt", sep = "\t", skip = 6, na.strings = c("[Error]-", "", "NA")) 
  ) %>% 
  mutate(ID = seq(1, nrow(.), by = 1)) %>% # <- Adding a numerical index to the dataset. Enabling  merge with VEGA, since SMILES is reconfigured in the VEGA output
  select(ID, everything(), -c(No., tId)) %>% 
  rename(SMILES_VEGA = SMILES) %>% 
  filter(!is.na(Fish.Chronic..NOEC..Toxicity.model..IRFMN..assessment), # removing all rows lacking estimations (NA's) 
         !grepl("\\?", Fish.Chronic..NOEC..Toxicity.model..IRFMN..assessment), # removing all rows with outrageous estimations (marked with an ?)
         !grepl("\\?", Daphnia.Magna.Chronic..NOEC..Toxicity.model..IRFMN..assessment), # removing all rows with outrageous estimations (marked with an ?)
         !grepl("\\?", Algae.Chronic..NOEC..Toxicity.model..IRFMN..assessment), # removing all rows with outrageous estimations (marked with an ?)
  ) %>% 
  distinct(CAS.Number, .keep_all = T)

# I now have a list with VEGA_QSARs attached to CAS.Numbers, both "QSAR toolbox-smiles" and VEGA-SMILES (slightly different configuration).  

# Tidying up in all cells and across the data frame
VEGA <- VEGA_QSARS %>% 
  select(!contains("prediction")) %>%  #<- remove all predictions columns and keep the assessment columns
  # Transforming the assessment columns into numeric effect concentrations and keeping the quality parameter in separate column 
  separate(Fish.Chronic..NOEC..Toxicity.model..IRFMN..assessment, c("Fish_Chrnic_NOEC_IRFMN_mgL", "Unit_Fish_IRFMN", "Fish_IRFMN_eval"), sep = (" "), convert = TRUE, extra = "merge") %>% 
  separate(Daphnia.Magna.Chronic..NOEC..Toxicity.model..IRFMN..assessment, c("D.magna_Chronic_NOEC_IRFMN_mgL", "Unit_D.magna_IRFMN", "D.magna_IRFMN_eval"), convert = TRUE, sep = (" "), extra = "merge") %>% 
  separate(Algae.Chronic..NOEC..Toxicity.model..IRFMN..assessment, c("Algae_Chronic_NOEC_IRFMN_mgL", "Unit_Algae_IRFMN", "Algae_IRFMN_eval"), sep = (" "), convert = TRUE, extra = "merge") %>% 
  # removing the "unit" column. Unit is specified in name of column ("mg/L")    
  select(!contains("Unit")) %>% 
  # Assigning columns as numeric type
  mutate(Fish_Chrnic_NOEC_IRFMN_mgL = as.numeric(Fish_Chrnic_NOEC_IRFMN_mgL),
         D.magna_Chronic_NOEC_IRFMN_mgL = as.numeric(D.magna_Chronic_NOEC_IRFMN_mgL),
         Algae_Chronic_NOEC_IRFMN_mgL = as.numeric(Algae_Chronic_NOEC_IRFMN_mgL),
         # removing redundant text in quality parameter column
         across(contains("_eval"), ~ as.character(gsub("\\(", "", .x))),
         across(contains("_eval"), ~ as.character(gsub(" reliability)", "", .x))),
         across(contains("_eval"), ~ as.character(gsub(" value)", "", .x)))
  ) 

#VEGA %>%   
#write.csv("Data_Output\\VEGA_QSAR_Wrangle_2022-06-10.csv", sep = ",", col.names = T)

```

VEGA HC20EC10eq calculations

Dealing with toxicity averages so that I can compare the USEtox database with in silico data.

```{r}

# Create a subset of only relevant data to create HC20 averages
VEGA_CF <- VEGA %>%
  #  filter(!CAS.Number %in% Inorganic_vector_VEGA) %>% # removing inorganic chemicals
  filter(Fish_Chrnic_NOEC_IRFMN_mgL != 0) %>%  # Removing data reported as 0 mg/L which further down creates infinite avlogHC50 values. these are obviously misrepresenting the toxicity.
  mutate(Fish_EC10eq_mgL = Fish_Chrnic_NOEC_IRFMN_mgL/0.6, # transforming To EC10eq-chronic from NOECchronic = 0.6 (0.4–0.7)
         D.magna_EC10eq_mgL = D.magna_Chronic_NOEC_IRFMN_mgL/0.6,
         Algae_EC10eq_mgL = Algae_Chronic_NOEC_IRFMN_mgL/0.6,
         # Adding an annotation of data quality 
         QSAR_Quality_fish_crust_algae = as.character(paste(Fish_IRFMN_eval, D.magna_IRFMN_eval, Algae_IRFMN_eval))) %>%
  # Calculating EC20EC10eq for VEGA estimations       
  rowwise() %>% 
  mutate(VEGA_log_Mean = mean(log10(c(Fish_EC10eq_mgL, D.magna_EC10eq_mgL, Algae_EC10eq_mgL)), na.rm = TRUE),
         sd_VEGA_log_Mean = sd(log10(c(Fish_EC10eq_mgL, D.magna_EC10eq_mgL, Algae_EC10eq_mgL)), na.rm = TRUE), # STRUL SD across three columns??
         CoV_VEGA_Mean = goeveg::cv(c(Fish_EC10eq_mgL, D.magna_EC10eq_mgL, Algae_EC10eq_mgL), na.rm = TRUE) # Coefficient of variation for all data points per species group
  ) %>% 
  ungroup() %>% 
  mutate(logHC20EC10eq_VEGA = VEGA_log_Mean + (sd_VEGA_log_Mean * -0.842),
         CRF_VEGA = 0.2/10^logHC20EC10eq_VEGA) %>%  # Concentration-response slope factor calculation based on all available data)
  distinct(CAS.Number, .keep_all = T) # Removing duplicate CAS.Numbers from the initial VEGA-data frame
```

Step 7. PHYSICOCHEMICAL Properties
Code loaded from file "Physicochemical_properties.Rmd".
```{r}
# Loading Physicochemial properties object from file.
NEW_PHYSCHEM <- read.csv("../data/PHYSCHEM_INFO.csv")
```

Merging the QSAR dataset with Physchem data

Issue: 
  The Physchem dataset has been revised to fit the HESTIA_EnviroTox database. Several physchem calculations and estimations have been changed along with it. 
There are 4k missing datapoints for the merge below.
However, there is no pressing need to add Physchem data to the QSAR dataset, since we will refrain from characterizing QSAR data at the endpoint level. Also, comparisons are based on the CRF slope factor, not at the midpoint CF-level.
```{r VEGA_Tox_wrangle}
# Making sure to get all the essential Physicochemical parameters into the VEGA_QSAR dataset
VEGA_CF_full <- left_join(
  x = VEGA_CF,
  y = NEW_PHYSCHEM %>% 
    select(
      CAS.Number, 
      MW.g.mol,
      pKaChemClass,
      pKa.gain, 
      pKa.loss, 
      Kow_L.L, 
      Koc_L.kg_MCI,
      kH25C_Pa.m3.mol,
      Vapor.Pressure_Pa,
      Sol_mg.L, 
      KdegA, 
      KdegW, 
      KdegSd, 
      KdegSl, 
      BAF_L.Kg,
      contains("Source")
    ),
  by = "CAS.Number") %>% 
  select(
    ID, 
    CAS.Number, 
    SMILES, 
    SMILES_VEGA,
    MW.g.mol,
    pKaChemClass,
    pKa.gain, 
    pKa.loss, 
    Kow_L.L, 
    Koc_L.kg_MCI,
    kH25C_Pa.m3.mol,
    Vapor.Pressure_Pa,
    Sol_mg.L, 
    KdegA, 
    KdegW, 
    KdegSd, 
    KdegSl, 
    CRF_VEGA, 
    BAF_L.Kg,
    logHC20EC10eq_VEGA,
    sd_VEGA_log_Mean,
    CoV_VEGA_Mean,
    QSAR_Quality_fish_crust_algae,
    contains("Source")
  ) 

# Export the data frame: 
write.csv(VEGA_CF_full, "../data/VEGA_HC20.csv", row.names = F)

# summary table


#write.xlsx(x= 
# cbind(data.frame(
#   VEGA %>% count(Fish_IRFMN_eval, sort = T),
#   VEGA %>% count(D.magna_IRFMN_eval, sort = T),
#   VEGA %>% count(D.magna_IRFMN.combase_eval, sort = T),
#   VEGA %>% count(Algae_IRFMN_eval, sort = T),
#   VEGA %>% count(Algae_ProtoQSAR.combase_eval, sort = T)
#   )),
#         file = "Data_output\\VEGA_Quality_table.xlsx",
#         sheetName = "VEGA_quality",
#         append = T
#         )

```
