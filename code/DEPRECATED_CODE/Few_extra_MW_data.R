library(tidyverse)
library(webchem)


missing_CAS_again <- HESTIA_HC20_DB_unit_and_value_conv %>% 
  filter(is.na(MW.g.mol)) %>% 
  distinct(CAS.Number) %>% 
  pull(CAS.Number)

CAS_CID_list <- webchem::get_cid(missing_CAS_again, from = "xref/rn", domain = "compound", match = "first", verbose = T)

CAS_CID_list$MW.g_mol <- webchem::pc_prop(CAS_CID_list$cid, properties = "MolecularWeight", verbose = T)[[2]]

str(read.csv("data/excel_references/CAS_CID_MW_list.csv") %>% 
      select(query, MW.g_mol) %>% 
    filter(!is.na(MW.g_mol)))

write.csv(CAS_CID_list, "data/excel_references/CAS_CID_MW_list.csv", row.names = F)
