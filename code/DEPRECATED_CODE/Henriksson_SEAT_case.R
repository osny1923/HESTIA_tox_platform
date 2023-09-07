Henriksson_inventory <- readxl::read_xlsx('data/HESTIA_additional_data/Henriksson_et_al_2015_inventory.xlsx', sheet = 1, col_names = TRUE)

Henriksson_inventory_inputs <- Henriksson_inventory[-1,] %>% 
  select(cycle.description, contains("cycle.inputs")) 

Henriksson_inventory_inputs_clean <- rename_with(Henriksson_inventory_inputs, ~ gsub("cycle.inputs.", "", .x), .cols = contains("cycle.inputs"))

Henriksson_inventory_long <- Henriksson_inventory_inputs_clean %>% 
  select(contains("@id")) %>% 
  pivot_longer(everything()) %>% 
  filter(value != "-") %>% 
  arrange(name) %>% 
  separate(value, into = c("name", "CASRN"), sep = "CAS-", remove = FALSE, convert = TRUE) %>% 
  mutate(CASRN = case_when(value == "oxytetracycline" ~ "79-57-2", TRUE ~ CASRN)) %>% 
  mutate(CASRN = case_when(value == "aminoglycosides" ~ "1403-66-3", TRUE ~ CASRN)) %>% 
  mutate(name = case_when(value == "aminoglycosides" ~ "gentamicin", TRUE ~ name)) %>% 
  filter(!is.na(CASRN)) %>% 
  distinct(CASRN, .keep_all = T)
  
Henriksson_inventory_w_names<- Henriksson_inventory_long %>% 
  rename(CAS.Number = CASRN) %>% 
  left_join(., 
            cas_smiles,
            by = "CAS.Number") %>% 
  mutate(term.name = coalesce(PesticideAI_name, name)) %>% 
  select(-name) %>% 
  left_join(., 
            HESTIA_HC20_dataset %>% 
              select(CAS.Number, HC20EC10eq, HC20), 
            by = "CAS.Number") %>% 
  rename(unweighted_HC20EC10eq = HC20EC10eq)

Henriksson_inventory_uncertainty <- Henriksson_inventory_w_names %>% 
  left_join(., 
            nls_output_df_OK,
            by = "CAS.Number") %>% 
  mutate(quantile_range = (Q2.5-Q97.5)/2, # <- average quantile range from the central value
         uncertainty_ratio = as.numeric(abs(quantile_range)/abs(log_HC20EC10eq)))

write.csv(Henriksson_inventory_uncertainty, "results/Henriksson_inventory_uncertainty.csv", row.names = F)

#########################################################################
# Trying again with the SEAT uncertainty inventory
library(tidyverse)
sheet_names <- readxl::excel_sheets("data/HESTIA_additional_data/Comparison of Asian Aquaculture Products by Use of Statistically Supported Life Cycle Assessment HENRIKSSON 2015.xls")[2:73]
final_sheet <- data.frame()
for (i in 1:length(sheet_names)) {
df_inventory <- readxl::read_xls("data/HESTIA_additional_data/Comparison of Asian Aquaculture Products by Use of Statistically Supported Life Cycle Assessment HENRIKSSON 2015.xls", sheet = sheet_names[i], col_names = TRUE)[,3]
names(df_inventory) <- sheet_names[i]
  
start_row <- na.omit((1:nrow(df_inventory))[df_inventory[,1] == "Environmental outputs"])
end_row <- nrow(df_inventory)

df_inventory <- df_inventory[start_row:end_row,]

removals <- df_inventory[3:10, 1] %>% pull(sheet_names[i])

df_inventory_chems <- df_inventory %>% 
  distinct() 

df_inventory_chems <- na.omit(df_inventory_chems[-c(1,3:10), 1]) %>% pull(sheet_names[i])

df_inventory_chems <- df_inventory_chems[!grepl(", ", df_inventory_chems)]

tmp_sheet <- list(inventory = NULL, process = NULL)
tmp_sheet$inventory <- df_inventory_chems 
tmp_sheet$process <- c(rep(sheet_names[i],length(df_inventory_chems)))

final_sheet <- rbind(final_sheet, tmp_sheet)

i = i + 1

}

SEAT_chems <- left_join(
  x = final_sheet %>% 
    rename(term.name = inventory),
  y = CAS.list_HESTIA,
  by = "term.name")

SEAT_chems_correct <- SEAT_chems %>% 
  mutate(
    CAS.number = case_when(
      term.name == "Oxytetracycline" ~ "79-57-2", 
      term.name == "Copper sulfate" ~ "17599-81-4",
      term.name == "2,4-D" ~ "94-75-7",
      term.name == "Fluazifop-P-buty" ~ "79241-46-6",
      term.name == "Ametryne" ~ "834-12-8",
      term.name == "Abamectin" ~ "71751-41-2",
      term.name == "Bensulfuron methyl ester" ~ "83055-99-6",
      term.name == "Chloropyrifos" ~ "2921-88-2",
      term.name == "Thiamethoxam" ~ "153719-23-4",
      term.name == "Sulfadimethoxine" ~ "122-11-2",
      term.name == "Trifloxystrobin" ~ "141517-21-7",
      term.name == "Prothioconazol" ~ "178928-70-6",
      term.name == "Imidacloprid" ~ "138261-41-3",
      term.name == "Metaconazole" ~ "125116-23-6",
      term.name == "Sulfometuron" ~ "74223-56-6",
      term.name == "Nicosulforon" ~ "111991-09-4",
      term.name == "Mancozeb" ~ "8018-01-7",
      term.name == "Chlorothaloni" ~ "1897-45-6",
      term.name == "Sulfamethoxazole" ~ "723-46-6",
      term.name == "Enrofoxacin" ~ "93106-60-6",
      term.name == "Avermectin" ~ "71751-41-2",
      term.name == "Dinitrigen monoxide" ~ "10024-97-2",
      term.name == "Dinitrogen monoxide" ~ "10024-97-2",
      term.name == "Dinitrogen monoixde" ~ "10024-97-2",
      term.name == "Dinitrogen monooxide" ~ "10024-97-2",
      term.name == "Colistin" ~ "1264-72-8",
      term.name == "Levofloxacin hydrate" ~ "100986-85-4",
      term.name == "Florfenical" ~ "73231-34-2",
      term.name == "Gentamycin sulfate" ~ "1405-41-0",
      term.name == "Zinc sulfate" ~ "7733-02-0",
      term.name == "DBDMH" ~ "77-48-5",
      term.name == "Benzalkonium bromide" ~ "7281-04-1",
      term.name == "BCDMH" ~ "16079-88-2",
      term.name == "Trichorfon" ~ "52-68-6",
      term.name == "Norlfoxacin" ~ "70458-96-7",
      term.name == "Amoxicillin" ~ "26787-78-0",
      term.name == "Chlorotetracycline" ~ "57-62-5",
      term.name == "Benzalkonium chloride" ~ "8001-54-5",
      term.name == "Benzalkonium chloride" ~ "8001-54-5",
      TRUE ~ CAS.number)) %>% 
  filter(!term.name %in% c("Wastewater", "Waste water", "waste water", "Chemical oxygen demand (COD)", "Ammonia to air"),
         !grepl("Inherent uncertainty", term.name)) %>% 
  rename(CAS.Number = CAS.number)

SEAT_uncertainty <- left_join( 
  x = SEAT_chems_correct,
  y = nls_output_df_OK, 
  by = "CAS.Number")


SEAT_uncertainty %>% 
filter(is.na(CRF)) %>% 
  distinct(term.name, .keep_all = T)

SEAT_summary <- SEAT_uncertainty %>% 
  mutate(count_na = case_when(is.na(Q2.5) ~ 0, TRUE ~ 1 )) %>% 
  group_by(process) %>% 
  summarise(n = n(),
            n_chems = sum(count_na),
            pct_covered = round((n_chems/n)*100, 0)) %>% 
  mutate(all_chems_covered = case_when(n == n_chems ~ "YES", TRUE ~ "NO")) %>% 
  ungroup()

SEAT_priority_chems <- left_join(x = SEAT_uncertainty, 
          y = SEAT_summary, 
          by = "process") %>% 
  arrange(desc(pct_covered)) 

SEAT_sum <- SEAT_priority_chems[, c(2, 15:16)] %>% distinct(process, .keep_all = T) %>% 
  arrange(desc(pct_covered)) 

write.csv(SEAT_sum, "results/SEAT_sum.csv", row.names = F)
write.csv(SEAT_priority_chems, "results/SEAT_priority_chems.csv", row.names = F)
write.csv(SEAT_uncertainty, "results/SEAT_uncertainty.csv", row.names = F)

