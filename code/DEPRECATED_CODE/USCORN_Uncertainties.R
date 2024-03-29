library(xlsx)
library(tidyverse)

Process_data <- read.xlsx("data/HESTIA_additional_data/HESTIA_data_distribution_summary.xlsx", sheetName = "Sheet1") %>% 
  filter(process == "CORN US") %>% 
  mutate(CAS.Number = case_when(term.name == "Mancozeb" ~ "8018-01-7", TRUE ~ CAS.Number))

# reading the USEtox input dataframe
US_CORN_USEtox <- read.csv("results/HESTIA_HC20_USEtox_format.csv") %>% 
  filter(CASRN %in% Process_data$CAS.Number)

# CV & Phi calculated from the Monte Carlo simulation vector of 100k iterations across a given mu and sigma generated by fitting a nonlinear least squares model to each chemical's data
US_CORN_phi <- read.csv("results/nls_output_df.csv") %>% 
    filter(CAS.Number %in% Process_data$CAS.Number) %>%
  select(CAS.Number, n_recs, log_HC20EC10eq, CV_HCx, Phi_HCx) %>% 
  left_join(., 
            US_CORN_USEtox %>% 
              select(CASRN, Name) %>% 
              rename(CAS.Number = CASRN), 
            by = "CAS.Number") %>% 
  select(1, 6, 2:5)

write.csv(US_CORN_phi, "results/US_CORN_phi.csv", row.names = F)


write.xlsx(US_CORN_USEtox, "results/US_CORN_inventory_USEtox_HC20.xlsx", sheetName = "Substance data", col.names = T, row.names = F, showNA = FALSE)


# Looking up missing datapoints
Process_data %>% 
  filter(!CAS.Number %in% HESTIA_BASE_dat$CAS.Number)
# Chemicals not present in the HESTIA_BASE_dat dataframe!!
# Azoxystrobin CORN US 131860-33-8
# Dicrotophos CORN US    141-66-2
# These are ONLY present in the EnviroTox database :(

# Old and bad CV/Phi calculations based on plain EC10eq data
# US_CORN_phi <- HESTIA_BASE_dat %>% 
#   filter(CAS.Number %in% Process_data$CAS.Number) %>% 
#   select(CAS.Number, EC10eq) %>% 
#   group_by(CAS.Number) %>% 
#   summarize(
#     n_data = n(), 
#     sd_EC10eq = sd(EC10eq),
#     cv_EC10eq = cv(EC10eq),
#     phi_EC10eq = sqrt(log((cv_EC10eq^2)+1))
#     ) %>% 
#   ungroup()


function(x, y = 20) {
  df <- data.frame()
  df <- df %>% mutate(!!paste("HC", y, sep = "_") := 2 * 2)
  return(df)
}
