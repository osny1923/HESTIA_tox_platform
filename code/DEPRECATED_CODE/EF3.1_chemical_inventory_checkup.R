# read libraries
library(tidyverse)
library(readxl)
library(xlsx)

# Read the EF3.1 CF's inventory
Owsianiak <- readxl::read_xlsx("data/EF3.1/EF3.1CFs.xlsx", sheet = 2, col_names = TRUE, na = "n/a")
HESTIA_inventory <- readxl::read_xlsx("data/EF3.1/pesticideAI-HESTIA summary.xlsx", sheet = "modified chems", col_names = TRUE) %>% 
  mutate(CAS.Number = gsub("CAS-", "", CAS.Number))

# Looking for chemicals in the EF3.1 Toxicological database
not_in_EF3.1 <- HESTIA_inventory %>% 
  filter(!CAS.Number %in% Owsianiak$CAS.Number)

present_in_EF3.1 <- HESTIA_inventory %>% 
  filter(CAS.Number %in% Owsianiak$CAS.Number)

# Looking for chemicals in the HESTIA Toxicological database
not_in_HEST_TOX <- HESTIA_inventory %>% 
  filter(!CAS.Number %in% HESTIA_HC20_dataset$CAS.Number)

present_in_HEST_TOX <- HESTIA_inventory %>% 
  filter(CAS.Number %in% HESTIA_HC20_dataset$CAS.Number)

# Looking for chemicals in the HESTIA chemical inventory
not_in_HEST_INVENTORY <- HESTIA_inventory %>% 
  filter(!CAS.Number %in% CAS_SMILES_list$CAS.Number)

present_in_HEST_INVENTORY <- HESTIA_inventory %>% 
  filter(CAS.Number %in% CAS_SMILES_list$CAS.Number)

library(openxlsx)

write_to_excel <- function(object_list, file_path) {
  wb <- createWorkbook()
  
  for (name in names(object_list)) {
    addWorksheet(wb, sheetName = name)
    writeData(wb, sheet = name, object_list[[name]])
  }
  
  saveWorkbook(wb, file = file_path, overwrite = TRUE)
}


object_list <- list(not_in_EF3.1 = not_in_EF3.1,
                    present_in_EF3.1 = present_in_EF3.1,
                    not_in_HEST_TOX = not_in_HEST_TOX,
                    present_in_HEST_TOX = present_in_HEST_TOX,
                    not_in_HEST_INVENTORY = not_in_HEST_INVENTORY,
                    present_in_HEST_INVENTORY = present_in_HEST_INVENTORY)

write_to_excel(object_list, "data/EF3.1/Inventory_data_availability.xlsx")

