# Making the physicochemical data wrangling as a function to just insert into runnig code. 
# The Physchem data can be read from file and through this function become ready to present in USEtox-friendly format.

# REMEMBER! The output file from OECD QSAR Toolbox is formatted in an UFT-16 format and need some specifications within the read.csv() operation:
  #read.csv("filepath.csv", header = T, sep = "\t", na.strings = "", fileEncoding = "UTF-16LE",  stringsAsFactors = FALSE)


# Libraries needed
library(dplyr)

Physchem_read_wrangle_function <- function(raw_physchem_input){
  # Apply function to the read.csv("OECD QSAR Toolbox output.csv") directly. 
  # Some wrangling and sub() expressions are needed to read the data as a useful format from the raw output. But this enables reading of an unprocessed raw output file from OECD QSAR Toolbox
physchem_output  <- {{raw_physchem_input}} %>% 
    select(-X., -contains("value", ignore.case = T)) %>% 
    # Multi constituent substances are are filtered out. I have no good way to deal with these substances
    filter(Predefined.substance.type != "Multi constituent") %>% 
    mutate(
      across(8:length(names(.)), ~ sub("No value", "", .)),
      across(8:length(names(.)), ~ sub(",", ".", .)),
      across(8:length(names(.)), ~ sub("\\s.*", "", .)), # <- removes anything after a blankspace which works great for removing all unit annotations!
      across(8:length(names(.)), ~ as.numeric(.))
    ) %>% 
    # If substances have been described with invalid CASRN, the OECD QSAR Toolbox will tag these as "Invalid Registry Number: XXXX-XX-X" 
    filter(!grepl("Invalid", CAS.Number)) %>% 
    # Removing fugacity models of biodegradation half-lives. not needed.
    select(-starts_with("FM")) %>% 
    # Unit conversions to USEtox-friendly format
    mutate(
      log.Kow = 10^log.Kow, # Converting Kow to L/L, removing logarithmic value
      Koc..MCI. = 10^Koc..MCI., # Converting Koc to L/kg, removing logarithmic value
      Koc..Log.Kow. = 10^Koc..Log.Kow., # Converting Koc to L/kg, removing logarithmic value
      Exp.Vapor.Pressure = 133.3*Exp.Vapor.Pressure, # Converting from mmHg to Pascal
      Selected.Vapor.Pressure = 133.3*Selected.Vapor.Pressure, # Converting from mmHg to Pascal
      Exp.Henrys.Law.Constant = 101325*Exp.Henrys.Law.Constant, # Converting from atm to Pa
      BAF..upper.trophic. = 10^BAF..upper.trophic., # Removing logarithmm
      Exp.Log.P = 10^Exp.Log.P  # Removing logarithm
    ) %>% 
    # Renaming columns to better fit USEtox input format data frame
    rename(
      Est.Kow_L.L = log.Kow,
      Exp.Kow_L.L = Exp.Log.P,
      Koc_L.kg_MCI = Koc..MCI.,
      Koc_L.kg_Kow = Koc..Log.Kow.,
      Biodeg_BIOWIN3 = Ultimate.biodeg..Biowin.3.,
      Exp.Water.Solubility_mg.L = Exp.Water.Solubility,
      Exp.Vapor.Pressure_Pa = Exp.Vapor.Pressure,
      BAF_L.Kg = BAF..upper.trophic.,
      pKa.gain = Basic.pKa..OASIS.Regression.,
      pKa.loss = Acidic.pKa..OASIS.Regression.,
      MW.g.mol = Molecular.Weight
   ) %>% 
   # Making a prioritization strategy, prioritizing experimental over estimated data into a new column, with an annotation column.
   mutate(
    # Biodegradation rates needs conversion from probability
    Biodeg_1s = case_when(
      Biodeg_BIOWIN3 > 4.75 & Biodeg_BIOWIN3 < 5 ~ 4.7E-5,
      Biodeg_BIOWIN3 > 4.25 & Biodeg_BIOWIN3 <= 4.75 ~ 6.4E-6,
      Biodeg_BIOWIN3 > 3.75 & Biodeg_BIOWIN3 <= 4.25 ~ 3.4E-5,
      Biodeg_BIOWIN3 > 3.25 & Biodeg_BIOWIN3 <= 3.75 ~ 9.3E-7,
      Biodeg_BIOWIN3 > 2.75 & Biodeg_BIOWIN3 <= 3.25 ~ 5.3E-7,
      Biodeg_BIOWIN3 > 2.25 & Biodeg_BIOWIN3 <= 2.75 ~ 2.1E-7,
      Biodeg_BIOWIN3 > 1.75 & Biodeg_BIOWIN3 <= 2.25 ~ 1.3E-7,
      Biodeg_BIOWIN3 <= 1.75 ~ 4.7E-5
      ),
      # applying the division factors for extrapolating W, Sl & Sd biodegradation rates with 1:2:9 respectively, according to Usetox manual.
      KdegW = Biodeg_1s,
      KdegSl = Biodeg_1s/2,
      KdegSd = Biodeg_1s/9,
      # Degradation rates in air is different. Prioritization of OH rate constant method, taking Episuite's SUMMARY (AOP v1.92): HYDROXYL RADICALS (25 deg C) 
      KdegA = case_when(
        is.na(OVERALL.OH.rate.constant) ~ as.numeric(NA),
          TRUE ~ (OVERALL.OH.rate.constant * 1.5E6)/2),
      KdegA = case_when(
        KdegA == 0 ~ as.numeric(NA),
          TRUE ~ KdegA),
      Kow_L.L = case_when(
        is.na(Exp.Kow_L.L) ~ Est.Kow_L.L,
          TRUE ~ Exp.Kow_L.L),
      Koc_L.kg = case_when(
        is.na(Koc_L.kg_MCI) ~ Koc_L.kg_Kow,
        #Koc_L.kg_MCI == 1 ~ Koc_L.kg_Kow,
          TRUE ~ Koc_L.kg_MCI),
      pKaChemClass = case_when(
        !is.na(pKa.loss) & is.na(pKa.gain) ~ "acid",
        is.na(pKa.loss) & !is.na(pKa.gain) ~ "base",
        !is.na(pKa.loss) & !is.na(pKa.gain) ~ "amphoter",
        Predefined.substance.type == "Multi constituent" ~ "undefined",
        is.na(pKa.loss) & is.na(pKa.gain) ~ "neutral"),
      Source_Kow = case_when(
        is.na(Est.Kow_L.L) & is.na(Exp.Kow_L.L) ~ "",
        is.na(Exp.Kow_L.L) ~ "Estimated",
          TRUE ~ "Experimental"),
      Source_Koc = case_when(
        is.na(Koc_L.kg_MCI) & is.na(Koc_L.kg_Kow) ~ "",
        is.na(Koc_L.kg_MCI) | Koc_L.kg_MCI == 1 ~ "Biowin_LogKow",
          TRUE ~ "Biowin_MCI"),
      Vapor.Pressure_Pa = case_when(
        is.na(Exp.Vapor.Pressure_Pa) ~ Selected.Vapor.Pressure,
          TRUE ~ Exp.Vapor.Pressure_Pa),
      Source_Pvap = case_when(
        is.na(Selected.Vapor.Pressure) & is.na(Exp.Vapor.Pressure_Pa) ~ as.character(NA),
        is.na(Exp.Vapor.Pressure_Pa) ~ "Estimated",
          TRUE ~ "Experimental"),
      Sol_mg.L = case_when(
        is.na(Exp.Water.Solubility_mg.L) ~ Water.Solubility,
          TRUE ~ Exp.Water.Solubility_mg.L),
      Source_Sol = case_when(
        is.na(Exp.Water.Solubility_mg.L) & is.na(Water.Solubility) ~ as.character(NA),
        is.na(Exp.Water.Solubility_mg.L) ~ "Estimated",
          TRUE ~ "Experimental"),
      kH25C_Pa.m3.mol = case_when(
        is.na(Exp.Henrys.Law.Constant) ~ (Vapor.Pressure_Pa*MW.g.mol)/Sol_mg.L, 
          TRUE ~ Exp.Henrys.Law.Constant),
      Source_KH25C = case_when(
        is.na(kH25C_Pa.m3.mol)  ~ as.character(NA),
        is.na(Exp.Henrys.Law.Constant) ~ "Calculated",
          TRUE ~ "Experimental"),
      # According to the USEtox manual, The estimation procedures (regression equations) are only suitable under certain conditions. "The regressions used in USEtox for calculating the Koc for the electrolytes are suited for acids within the pKa range 0–12 and with a log Kow between -2.18 and 8.50. For bases the pKa needs to be above 2 and log Kow is between -1.66 and 7.03 (Franco & Trapp 2008)."
    Koc_L.kg = case_when(
        pKaChemClass == "acid" | pKaChemClass == "amphoter" & pKa.loss <12 & pKa.loss > 0 & log(Kow_L.L) >= -2.18 & log(Kow_L.L) <= 8.5  ~ Koc_L.kg,
        pKaChemClass == "base" | pKaChemClass == "amphoter" & pKa.gain >2 & log(Kow_L.L) >= -1.66 & log(Kow_L.L) <= 7.03 ~ Koc_L.kg,
        pKaChemClass == "neutral" ~ Koc_L.kg)
    ) %>% 
    select(CAS.Number, Molecular.formula, Predefined.substance.type, MW.g.mol, pKaChemClass, pKa.gain, pKa.loss, Kow_L.L, Koc_L.kg, kH25C_Pa.m3.mol, Vapor.Pressure_Pa, Sol_mg.L, KdegA, KdegW, KdegSl, KdegSd, BAF_L.Kg, starts_with("Source"))
  return(physchem_output)
}
