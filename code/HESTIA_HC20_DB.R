## ----setup, echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(message = F, warning = F, echo = F)
# These libraries are used for analysis
   #install.packages("webchem")
   #install.packages("taxize")   # <- Installing the "taxize" library http://dx.doi.org/10.5281/zenodo.7097
   #install.packages("networkD3") <- For the Sankey flow chart visualization
   #install.packages("goeveg") # <- for simple coefficient of Variation calculation at summary of data
    library(rmarkdown)
    library(tidyverse)
    library(webchem)
    library(xlsx)
    library(readr)
    library(kableExtra)
    library(taxize)
    library(ggpubr)
    library(networkD3)
    library(htmlwidgets)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# This step takes place in the `code/Translating_CAS_to_SMILES_via_PubChem.R`-file
# Output from this operation creates three files: a .csv doc "data/excel_references/CAS_CID_list_final.csv", a .txt file containing CASRN-SMILES matches: "data/excel_references/CAS_CID_list_final.txt", as well as five subsets from the previous .txt file named "data/excel_references/cas_smiles_list[XXXX-XXXXk].txt"


## ----Physchem info, echo=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Physchem wrangle function
source("../code/Physchem_read_wrangle_function.R")
# Processing and wrangling the raw data through function loaded above
Physchem_HESTIA <- Physchem_read_wrangle_function(read.csv("../data/QSAR_Toolbox_physchem_data_2023-03-24_RAW.csv", header = T, sep = "\t", na.strings = "", fileEncoding = "UTF-16LE",  stringsAsFactors = FALSE))

# Loading the slim version of the Pesticide annotations data for substance use properties and merge this here.
HESTIA_chem_list_slim <- read.csv("../results/HESTIA_chem_list_slim.csv")

NEW_PHYSCHEM <- left_join(
  x = Physchem_HESTIA, 
  y = HESTIA_chem_list_slim, 
  by = "CAS.Number"
) %>% 
  select(CAS.Number, CanonicalSMILES, PesticideAI_name, 2:31)



## ----reading TOX data-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Reading in the data input ready for wrangling. "Pre-filtered" implying that the data has been 
# Dependency -> "../data/raw_data_read_and_wrangle.R"
HESTIA_HC20_DB <- read.csv("../data/RAW_DB/HESTIA_HC20_DB_raw_toxdata.csv") %>% 
  select(-SMILES) %>% 
  mutate(Year = coalesce(Year, Publication.year)) %>% 
  select(-Publication.year)



## ----Endpoint selection---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## Applying conversions for Endpoints
source("../code/EC10eq_conversion_functions.R")

HESTIA_HC20_DB_endpoint_conversions <- HESTIA_HC20_DB %>%
  mutate(Endpoint_conv = mapply(endpoint_conv_function, Endpoint)) 



## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
HESTIA_HC20_DB_qualifier_filter <- HESTIA_HC20_DB_endpoint_conversions %>% 
  # Qualifiers, step one. in Saouter et al., 2018 "When there is a lower range value with the descriptors ‘>=, ca., or empty’, the lowest value is selected. If, within this group, a test has also a higher value, this higher value is ignored. "case_when(is.na(Value.MeanValue)
  mutate(
    Value.MeanValue = case_when(is.na(Value.MeanValue) ~
      case_when(
        # If a qualifyer is stated as "ca." or "larger than" 
        Value.MinQualifier %in% c("ca.", ">=") ~ Value.MinValue,
        # When there is a lower range value with the descriptors "empty", the lowest value is selected
        !is.na(Value.MeanValue) & is.na(Value.MinQualifier) & !is.na(Value.MinValue) ~ Value.MinValue, 
        # Mean value is NOT empty, but the qualifier is empty => The original value is selected (Most data!)
        !is.na(Value.MeanValue) & is.na(Value.MinQualifier) ~ Value.MeanValue,
        # Lower range qualifier is empty => lowest value is selected
        is.na(Value.MeanValue) & is.na(Value.MinQualifier) ~ Value.MinValue, 
          # In case of NOEC > than, the value was kept since it is still representing a concentration with no observed effect."  
          # Qualifiers, step Two. "All lower range values described as ‘>’ are ignored, unless the higher value is described as ‘=<’. 
        Value.MinQualifier == ">" ~ case_when(
            Endpoint_conv == "NOEC" ~ Value.MinValue,
            Value.MaxQualifier == "<="  ~ Value.MinValue,
            TRUE ~ "ignore"),
        # Qualifiers, step 3. All higher values described as ‘< than’ are ignored, unless the lower range value is described as ‘>=’. Then the lower value is used.
        Value.MaxQualifier == "<" ~ case_when(
          Value.MinQualifier == ">=" ~ Value.MinValue, TRUE ~ "ignore"),
        # Qualifiers, step 4. When a lower range value is missing (0 or blank) and a higher value is available described as ‘<=’, the higher value is used. 
          Value.MinValue = 0 | is.na(Value.MinValue) & Value.MaxQualifier == "<=" ~ Value.MaxValue,
        # Qualifiers, step 5. When a lower value is described as >= and the higher value is described as <=, the lowest value is used. 
          Value.MinQualifier == ">=" & Value.MaxQualifier == "<="  ~ Value.MinValue,
        # Qualifiers, step 5. Values expressed as ‘<’ are excluded.
          Value.MinQualifier == "<" ~ "ignore", TRUE ~ Value.MeanValue),
      TRUE ~ Value.MeanValue
      )
    ) %>% 
    filter(is.na(Value.MeanValue)|!Value.MeanValue == "ignore") # Removing all data with where qualifier filter is set to "ignore"

# Creating a summary table for value.qualifiers reported.
HESTIA_HC20_DB_endpoint_conversions %>%
  group_by(Value.Qualifier) %>% 
  summarise(n_val.qualifier = n()) %>% 
    left_join(
      x = .,
      y = HESTIA_HC20_DB_endpoint_conversions %>%
            group_by(Value.MinQualifier) %>% 
            summarise(n_min.qualifier = n()) %>% 
            rename(Value.Qualifier = Value.MinQualifier),
      by = "Value.Qualifier"
        ) %>% 
    left_join(
      x = .,
      y = HESTIA_HC20_DB_endpoint_conversions %>%
            group_by(Value.MaxQualifier) %>% 
            summarise(n_max.qualifier = n()) %>% 
            rename(Value.Qualifier = Value.MaxQualifier),
      by = "Value.Qualifier"
        ) %>% 
  write.csv(., "../data/excel_references/summary table for value.qualifiers.csv", row.names = F)


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
HESTIA_HC20_DB_Control.Type_filter <- HESTIA_HC20_DB_qualifier_filter %>% 
  filter(
    # Removing data based on insufficient/unsatisfactory controls
    !Control.type %in% c("Insufficient", "Unsatisfactory"),
    # Removing QSAR in title
    !grepl("QSAR", x = Title, ignore.case = T),
    # Removing bioassay in title
    !grepl("bioassay", x = Title, ignore.case = T),
    # Removing Quantitative in title
    !grepl("Quantitative", x = Title, ignore.case = T),
    # Removing QSAR in experimental design annotation
    !grepl("QSAR", x = Experimental.design, ignore.case = T),
    # Removing bioassay in experimental design annotation
    !grepl("bioassay", x = Experimental.design, ignore.case = T),
    # Removing QSAR study
    !grepl("Assessment of Aquatic Experimental Versus Predicted and Extrapolated Chronic Toxicity Data of Four Structural Analogues", Title),
    # Removing QSAR study
    !grepl("Effects of Chlorpyrifos, Carbendazim, and Linuron on the Ecology of a Small Indoor Aquatic Microcosm", Title),
    # Removing QSAR study
    !grepl("Altenburger,R., H. Walter, and M. Grote", Title),
    # Removing QSAR study
    !grepl("PREDICTING MODES OF TOXIC ACTION FROM CHEMICAL STRUCTURE: ACUTE TOXICITY IN THE FATHEAD MINNOW", Title),
    # Data entered incorrectly, mixing cardiovascular injections with soaking solutions, also, data reported as minimum non-lethal concentration is entered as EC10. Remove all of these.
    !grepl("Human Cardiotoxic Drugs Delivered by Soaking and Microinjection Induce Cardiovascular Toxicity in Zebrafish", Title),
    # Data on species-specific toxicity entered incorrectly, attributing D magna toxicity to S. capricirnutum.
    !grepl("A Multi-Battery Toxicity Investigation on Fungicides", Title),
    # Bioassay test-results
    !grepl("Rainbow Trout Larvae Compared with Daphnia", Title),
    # No source available, yet substance 107-21-1 from this study causes an extreme outlier
    !grepl("GERISH", Title),
    # 95% CI has one report of negative conventrations, which the qualifyer filter above selects. removing all data
    !grepl("Acute Effects of Binary Mixtures of Imidacloprid and Tebuconazole on 4 Freshwater Invertebrates", Title),
    # 95% CI has one report of negative conventrations, which the qualifyer filter above selects. removing all data
    !grepl("Relative Chronic Sensitivity of Neonicotinoid Insecticides to Ceriodaphnia dubia and Daphnia magna", Title),
    #Presents ridiculously low effect concentrations from micrplastic study on D. rerio embryos concentrations reported are used as an effect test result! incorrect entry into database
    !grepl("Multi-Laboratory Hazard Assessment of Contaminated Microplastic Particles by Means of Enhanced Fish Embryo Test with the Zebrafish", Title),
    # Study runs toxicological effect testing across a range of temperatures, leading to skewed data output.
    !grepl("Water Toxicology and Radioecology. Acute Toxicity of Heavy Metals to Aquatic Invertebrates at Different Temperatures", Title),
    # Duplicate data exists from the same author, also recorded in "Little, L.W., J.C. Lamb III, M.A. Chillingworth, and W.B. Durkin, Acute Toxicity of Selected Commercial Dyes to the Fathead Minnow and Evaluation of Biological Treatment for Reduction of Toxicity"
    !grepl("Acute Toxicity of 46 Selected Dyes to the Fathead Minnow, Pimephales promelas", Title),
    # The Japan MoE Ecotoxicological tests on chemical Hexabromobenzene, CASRN=87-82-1, is reported with data 6 orders of magnitude from all other data, and the specific test report is not to traceable.
    Database != "Aquatic Japan MoE" | CAS.Number != "87-82-1"
  )



## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
HESTIA_HC20_DB_water_filter <- HESTIA_HC20_DB_Control.Type_filter %>% 
  mutate(
    Media.type = gsub(" ", "", Media.type),
    Media.type = coalesce(
      Media.type, Water.type),
    Media.type = case_when(
      is.na(Media.type) ~ "Freshwater", 
      TRUE ~ Media.type)
    ) %>% 
  # Filtering out all "saltwater" and "no substrate" tests. 
  filter(!Media.type %in% c("Saltwater", "Nosubstrate")
         ) 
  


## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
HESTIA_HC20_DB_effect_crit_filter <- HESTIA_HC20_DB_water_filter %>% 
filter(Effect %in% c(
  "Behavior", "Growth", "Intoxication", "Population", "Reproduction", 
  "Acute", "Cell(s)", "Growth Rate", "Mortality", "Feeding Behavior",
  "Biomass", "Body Weight", "Chronic", "Frond Number", "Development", 
  "Mobility","Seedling Emergence", "Immobilisation", "Behaviour" 
  ) )

HESTIA_HC20_DB_effect_crit_filter %>% 
  group_by(Effect) %>% 
  summarise(n = n()) %>% 
  arrange(desc(n)) %>% 
  write.csv(., "../data/excel_references/summary_table_effect_crit.csv", row.names = F)



## ---- warning=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Removing poor taxonomic descriptions and Species data annotated at higher taxonomic level. 
HESTIA_HC20_DB_taxonomy_filter <-  HESTIA_HC20_DB_effect_crit_filter %>% 
  left_join(x = .,
  # Joining with the (in part) manually curated taxonomy list for the HESTIA DB
    read.csv("../data/Taxonomy/Species_taxonomy.csv"), 
    by = "Test.organisms..species.") %>%
  filter(!is.na(Species),
         !is.na(Taxonomy.Group))

# Overview of the filtered out species, or genus and higher taxonomy rather...
# HESTIA_HC20_DB_effect_crit_filter %>% 
#   left_join(x = .,
#   # Joining with the (in part) manually curated taxonomy list for the HESTIA DB
#     read.csv("../data/Taxonomy/Species_taxonomy.csv"), 
#     by = "Test.organisms..species.") %>%
#   filter(is.na(Taxonomy.Group) | is.na(Species)) %>% 
#   count(Test.organisms..species., sort = T)



## ----dur_and_val harmonization, warning = FALSE---------------------------------------------------------------------------------------------------------------------------------------------------------------------

HESTIA_HC20_DB_unit_and_value_conv <- HESTIA_HC20_DB_taxonomy_filter %>% 
  left_join(
    x = .,
    y = NEW_PHYSCHEM %>% 
      select(CAS.Number, MW.g.mol), 
    by = "CAS.Number") %>%
  # Some mean value cells are empty, but can be calculated from the min and max value
  mutate(
  # replacing excel-format commas for dots.
    Value.MeanValue = gsub(pattern = ",", ".", Value.MeanValue), 
    # Converting effect concentrations into numeric data
    Value.MeanValue = as.numeric(Value.MeanValue), 
    # selecting values of Exposure duration as the determined duration, but when missing, test duration is used.
    Duration.MeanValue = case_when(
      is.na(Exposure.duration.MeanValue) & !is.na(Duration.MeanValue) ~ Duration.MeanValue, 
      TRUE ~ Exposure.duration.MeanValue),
    # selecting Units of Exposure.duration when "original" duration is missing
    Duration.Unit = case_when(
      is.na(Exposure.duration.Unit) & !is.na(Duration.Unit) ~ Duration.Unit, 
      TRUE ~ Exposure.duration.Unit),               
    # Converting Test duration values into numeric data
    Duration.MeanValue = as.numeric(Duration.MeanValue),
    # Converting test-duration to a coherent unit (hour)
    Time.Hours = case_when(
      Duration.Unit == "min" ~ Duration.MeanValue/60,
      Duration.Unit == "Second(s)" ~ Duration.MeanValue/3600,
      Duration.Unit == "h" ~ Duration.MeanValue,
      Duration.Unit == "d" ~ Duration.MeanValue*24,
      Duration.Unit == "wk" ~ Duration.MeanValue*(24*7),
      Duration.Unit == "mo" ~ Duration.MeanValue*(24*30),
      Duration.Unit == "yr" ~ Duration.MeanValue*(24*365),
      TRUE ~ as.numeric(NA)
    ),
  # Converting all eligible values into a coherent unit (mg/L)
  Value.Unit = gsub(" ", "", as.character(Value.Unit)), 
  Value.mg_l = case_when(
    Value.Unit %in% c("µg/L", "ng/mL", "µg/dm³") ~ Value.MeanValue/1E3,
    Value.Unit == "ppb" ~ Value.MeanValue/1E3,
    Value.Unit %in% c("mg/L", "µg/mL", "g/m³", "mg/dm³")	 ~ Value.MeanValue,
    Value.Unit == "µg/3.5L" ~ Value.MeanValue/0.000285714286,
    Value.Unit %in% c("ppm", "µg/cm³") ~ Value.MeanValue,
    Value.Unit %in% c("ng/L","pg/mL", "µg/µL") ~ Value.MeanValue/1E6,
    Value.Unit %in% c("g/L", "g/dm³", "mg/mL", "µg/mm³") ~ Value.MeanValue*1E3,
    Value.Unit == "µg/10L" ~ Value.MeanValue/1E4,
    Value.Unit == "pg/L" ~ Value.MeanValue/1E9,
    Value.Unit == "g/mL" ~ Value.MeanValue*1E6,
    Value.Unit == "µg/100mL" ~ Value.MeanValue/100,
    Value.Unit == "mg/100cm³" ~ Value.MeanValue*10,
    Value.Unit == "µg/5mL" ~ Value.MeanValue/0.2,
    Value.Unit == "mg/200mL" ~ Value.MeanValue*5,
    Value.Unit %in% c("mol/L", "M", "mol","mol/dm³") ~ (Value.MeanValue*MW.g.mol)*1E3,
    Value.Unit %in% c("mmol/L", "mM", "mmol/dm³", "mol/m³") ~ (Value.MeanValue*MW.g.mol),
    Value.Unit %in% c("µmol/L", "µM", "mmol","µm", "µmol/dm³", "uM/L", "µM/L", "mmol/m³", "nmol/mL") ~ (Value.MeanValue*MW.g.mol)/1E3,
    Value.Unit %in% c("nmol/L", "nmol", "nM/L", "nM")~ (Value.MeanValue*MW.g.mol)/1E6,
    Value.Unit == c("pM")~ (Value.MeanValue*MW.g.mol)/1E9,
    TRUE ~ as.numeric(NA)
  ),
  # Making sure that all units that are convertable into mg/L are actually converted to the same unit
  # including values reported as molar units, where 1 M = 1 mol/L (https://en.wikipedia.org/wiki/Molar_concentration)
  Value.unit.mg_l = as.factor(case_when(
    Value.Unit %in% c(
      "µg/L", "µg/cm³","ng/L","ng/mL","ppm", "ppb", "pg/mL", "mg/mL", "g/L", "g/dm³", 
      "µM", "µm", "mM", "nM", "M", "pM", "mol", "nmol", "mmol",
      "mg/dm³", "g/m³", "µg/dm³", "µg/100mL", "µg/µL", "µg/5mL", "µg/10L", "mg/L",
      "µg/mm³", "pg/L", "g/mL", "mg/200mL", "mol/dm³", "µg/mL", "nmol/L", "nM/L", "mg/100cm³",
      "µmol/L", "µmol/dm³", "uM/L", "µM/L", "mmol/L", "mmol/dm³", 
      "mol/L", "mmol/m³", "µg/3.5L", "mol/m³", "nmol/mL" ) ~ "mg/L",
      TRUE ~ as.character(NA)))
  )


## ----dur_and_val_filtering, warning = FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# FILTER STEP! 
# Removing effect concentrations that was reported as NA or "0", as well as Value.mg_l == NA

HESTIA_HC20_DB_effect_conc_filter <- HESTIA_HC20_DB_unit_and_value_conv %>% 
  filter(
    # Removing records not convertible into mg/L due to inapplicable units 
    !is.na(Value.mg_l),
    # Removing records missing mol weights.
    !is.na(MW.g.mol), 
    # Removing effect concentration "0"
    Value.MeanValue != 0
  ) 



## -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
HESTIA_HC20_DB_val.unit_filter <- HESTIA_HC20_DB_effect_conc_filter %>% 
  mutate(Time.Hours = as.numeric(Time.Hours)) %>%
  # Defining time chronic/Acute exposure time units for USEtox classification
  mutate(AcuteChronic = as.factor(
    case_when(
      Taxonomy.Group %in% c("Fish", "Plant", "Insect", "Mollusca", "Annellidae", "Amphibian") ~ case_when(Time.Hours < 168 ~ "Acute",
                                                                                                          Time.Hours >= 168 ~ "Chronic"),
      Taxonomy.Group == "Crustacean" ~ case_when(Time.Hours < 96 ~ "Acute",
                                    Time.Hours >= 96 ~ "Chronic"),
      Taxonomy.Group %in% c("Algae", "Rotifera") ~ case_when(Time.Hours < 24 ~ "Acute",
                                            Time.Hours >= 24 ~ "Chronic"),
      Taxonomy.Group == "Others" & !Phylum %in% c("Chordata", "Arthropoda", "Cnidaria") ~ case_when(Time.Hours < 24 ~ "Acute",
                             Time.Hours >= 24 ~ "Chronic"),
      Taxonomy.Group == "Others" & Phylum == "Arthropoda" ~ case_when(Time.Hours < 96 ~ "Acute",
                             Time.Hours >= 96 ~ "Chronic"),
                             TRUE ~ case_when(Time.Hours < 168 ~ "Acute",
                                              Time.Hours >= 168 ~ "Chronic") ) 
    ))



## ----dur_unit_filter------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
HESTIA_HC20_DB_time.unit_filter <- HESTIA_HC20_DB_val.unit_filter %>% 
  mutate(Time.Hours = as.numeric(Time.Hours)) %>% 
  filter(Time.Hours != 0,
         Time.Hours >= 24, 
         !is.na(Time.Hours))


## ----EC10eq_conversions---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# load EC10eq converson function
source("../code/EC10eq_conversion_functions.R")

system.time(Q_dat <- HESTIA_HC20_DB_time.unit_filter %>% 
  mutate(
    EC10eq = mapply(ec10eq_extrapolation_function, Value.mg_l, Endpoint, AcuteChronic, Taxonomy.Group, "extpl"),
    EC10eq_high = mapply(ec10eq_extrapolation_function, Value.mg_l, Endpoint, AcuteChronic, Taxonomy.Group, "high_CI"),
    EC10eq_low = mapply(ec10eq_extrapolation_function, Value.mg_l, Endpoint, AcuteChronic, Taxonomy.Group, "low_CI"),
    EC10eq_Acute = case_when(AcuteChronic == "Acute" ~ EC10eq, TRUE ~ as.numeric(NA)),
    EC10eq_Chronic = case_when(AcuteChronic == "Chronic" ~ EC10eq, TRUE ~ as.numeric(NA)),
         No.Extrapolations = as.numeric(case_when(Endpoint_conv == "EC10" ~ 0, TRUE ~ 1)),# Making sure to annotate whether an effect value is extrapolated or not for downstream applications.
    DB = "HESTIA",
    version = "HESTIA 1.3"
    )  )

# write.csv(expf_df, "../data/excel_references/EC10_extrapolation_factor_summary_table.csv", row.names = F)

write.csv(Q_dat, "../results/HESTIA_EC10eq_DB.csv", row.names = F)




## ---- echo=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

nodes = data.frame("name" = c(
   "OECD_QSAR_Toolbox", 
   "ECOTOX", "Aquatic ECETOC", "Food TOX Hazard EFSA", "Aquatic Japan MoE", "Aquatic Oasis",
   "1. Validated data", 
   "2. Data in range-filter", "3. Controls insufficient or unsatisfactory", "4. Non-Fresh Water data", "5. Effect criterions irrelevant", "6. Poor taxonomic descriptions",
   "7. Effect unit/value '0'or'NA'", "8. Duration 'NA' or <24h",
   "EC50-acute", "EC50-chronic", "EC10-acute", "EC10-chronic", "NOEC-acute", "NOEC-chronic",
   "Chronic EC10eq"),
   "Node_Group" = c("Toolbox", rep("DBs", 5), "Validated_data", rep("Filter1", 7),  rep(c("EndpointA", "EndpointC"), 3), rep("EC10eq",1))
               )

links = as.data.frame(matrix(c(
    1,0,nrow(HESTIA_HC20_DB %>% filter(Database == "ECOTOX")),
    2,0,nrow(HESTIA_HC20_DB %>% filter(Database == "Aquatic ECETOC")),
    3,0,nrow(HESTIA_HC20_DB %>% filter(Database == "Food TOX Hazard EFSA")),
    4,0,nrow(HESTIA_HC20_DB %>% filter(Database == "Aquatic Japan MoE")),
    5,0,nrow(HESTIA_HC20_DB %>% filter(Database == "Aquatic OASIS")),
    0,6,nrow(Q_dat), # Stuff kept from filter operations
    0,7,nrow(HESTIA_HC20_DB) - nrow(HESTIA_HC20_DB_qualifier_filter), # Selecting effect data within a range 
    0,8,nrow(HESTIA_HC20_DB_qualifier_filter) - nrow(HESTIA_HC20_DB_Control.Type_filter),  # Test.Control == "insufficient" or "unsatisfactory"
    0,9,nrow(HESTIA_HC20_DB_Control.Type_filter) - nrow(HESTIA_HC20_DB_water_filter), # Non-FW media type
    0,10,nrow(HESTIA_HC20_DB_water_filter) - nrow(HESTIA_HC20_DB_effect_crit_filter), # Irrelevant effect criterions
    0,11,nrow(HESTIA_HC20_DB_effect_crit_filter) - nrow(HESTIA_HC20_DB_taxonomy_filter), # Taxonomy
    0,12,nrow(HESTIA_HC20_DB_unit_and_value_conv) - nrow(HESTIA_HC20_DB_effect_conc_filter), # Effect Unit/Value "na" or "0"
    #0,13,nrow(HESTIA_HC20_DB_effect_conc_filter) - nrow(HESTIA_HC20_DB_val.unit_filter), # Incorrect units to describe effect
    0,13,nrow(HESTIA_HC20_DB_val.unit_filter) - nrow(HESTIA_HC20_DB_time.unit_filter), # Incorrect time unit
    
    6,14,nrow(Q_dat %>% # How many EC50 do we have? 
                    filter(Endpoint_conv == "EC50",
                           AcuteChronic == "Acute")),
    6,15,nrow(Q_dat %>% # How many LOEC do we have? 
                    filter(Endpoint_conv == "EC50",
                           AcuteChronic == "Chronic")),
    6,16,nrow(Q_dat %>% # How many NOEC do we have? 
                    filter(Endpoint_conv == "EC10",
                           AcuteChronic == "Acute")),
    6,17,nrow(Q_dat %>% # How many EC50 do we have? 
                    filter(Endpoint_conv == "EC10",
                           AcuteChronic == "Chronic")),
    6,18,nrow(Q_dat %>% # How many LOEC do we have? 
                    filter(Endpoint_conv == "NOEC",
                           AcuteChronic == "Acute")),
    6,19,nrow(Q_dat %>% # How many NOEC do we have? 
                    filter(Endpoint_conv == "NOEC",
                           AcuteChronic == "Chronic")),
    
    14,20,nrow(Q_dat %>% # How many Acute EC50 do we have? 
                    filter(Endpoint_conv == "EC50",
                           AcuteChronic == "Acute")),
    15,20,nrow(Q_dat %>% # How many Chronic EC50 do we have? 
                    filter(Endpoint_conv == "EC50",
                           AcuteChronic == "Chronic")),
    16,20,nrow(Q_dat %>% # How many Acute EC10 do we have? 
                    filter(Endpoint_conv == "EC10",
                           AcuteChronic == "Acute")),
    17,20,nrow(Q_dat %>% # How many Chronic do we have? 
                    filter(Endpoint_conv == "EC10",
                           AcuteChronic == "Chronic")),
    18,20,nrow(Q_dat %>% # How many Acute NOEC do we have? 
                    filter(Endpoint_conv == "NOEC",
                           AcuteChronic == "Acute")),
    19,20,nrow(Q_dat %>% # How many Chronic NOEC do we have? 
                    filter(Endpoint_conv == "NOEC",
                           AcuteChronic == "Chronic"))
                                        ),
  byrow = TRUE, ncol = 3))

names(links) = c("source", "target", "value")
links$groups = c(rep("A", 5), "B", rep("C", 7), rep("D", 6), "A_EC10", "C_EC10","A_EC10", "C_EC10","A_EC10", "C_EC10")

Wrangle_plot <- sankeyNetwork(Links = links, Nodes = nodes,
 Source = "source", Target = "target",
 Value = "value", NodeID = "name", 
 LinkGroup = "groups", NodeGroup = "Node_Group",
 units = "n", fontFamily = "sans-serif",
 fontSize = 10, nodeWidth = 30, sinksRight = FALSE)

Wrangle_plot

saveNetwork(Wrangle_plot, "Wrangle_plot.html", selfcontained = TRUE)
saveWidget(Wrangle_plot, file = "Wrangle_plot_widget.html", selfcontained = TRUE, title = "Overview of data curation steps")


