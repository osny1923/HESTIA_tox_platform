library(dplyr)

##### Filter decided on ####
HESTIA_HC20_DB_qualifier_filter <- HESTIA_HC20_DB_endpoint_conversions %>% 
  # Qualifiers, step one. in Saouter et al., 2018 "When there is a lower range value with the descriptors ‘>=, ca., or empty’, the lowest value is selected. If, within this group, a test has also a higher value, this higher value is ignored.
  mutate(
    Value.MeanValue = case_when(Qualifier == "Data in range" ~ case_when(
      # If a qualifyer is stated as "ca." or "larger than" 
      Value.MinQualifier %in% c("ca.", "=")| is.na(Value.MinQualifier) ~ Value.MinValue,
      # In case of NOEC > than, the value was kept since it is still representing a concentration with no observed effect. 
      Value.MinQualifier == ">" ~ case_when(
        Endpoint_conv == "NOEC" ~ Value.MinValue, 
        TRUE ~ "ignore"),
      Value.MinQualifier == "<" ~ "ignore"),
        TRUE ~ Value.MeanValue)
    )

##### First filter, that is dealing with values where a man is already appointed. (not neccesary).
HESTIA_HC20_DB_qualifier_filter <- HESTIA_HC20_DB_endpoint_conversions %>% 
  # Qualifiers, step one. in Saouter et al., 2018 "When there is a lower range value with the descriptors ‘>=, ca., or empty’, the lowest value is selected. If, within this group, a test has also a higher value, this higher value is ignored.
  mutate(
    Value.MeanValue = case_when(is.na(Value.MeanValue) ~
                                  case_when(
                                    # If a qualifyer is stated as "ca." or "larger than" 
                                    Value.MinQualifier == "ca." ~ Value.MinValue,
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
  # Removing all data with where qualifier filter is set to "ignore"
  filter(is.na(Value.MeanValue)|!Value.MeanValue == "ignore") 

nrow(HESTIA_HC20_DB_endpoint_conversions) - nrow(HESTIA_HC20_DB_qualifier_filter)

### Summary tables ##########
# Creating a overview counts table for value.qualifiers reported.
left_join(
  x = HESTIA_HC20_DB_endpoint_conversions %>% 
    filter(is.na(Value.MeanValue)) %>% 
    count(Value.MinQualifier) %>% 
    rename(n_min_qual = n, 
           Qualifier= Value.MinQualifier),
  y = HESTIA_HC20_DB_endpoint_conversions %>% 
    filter(is.na(Value.MeanValue)) %>% 
    count(Value.MaxQualifier) %>% 
    rename(n_max_qual = n, 
           Qualifier= Value.MaxQualifier),
  by = "Qualifier"
) %>% 
  mutate(Qualifier = case_when(is.na(Qualifier) ~ "No qualifier", TRUE ~ Qualifier)) 

# HESTIA_HC20_DB_endpoint_conversions %>%
#   filter(is.na(Value.MeanValue)) %>% 
#   group_by(Value.Qualifier) %>%
#   summarise(n_val.qualifier = n()) %>% 
#     left_join(
#       x = .,
#       y = HESTIA_HC20_DB_endpoint_conversions %>%
#             group_by(Value.MinQualifier) %>% 
#             summarise(n_min.qualifier = n()) %>% 
#             rename(Value.Qualifier = Value.MinQualifier),
#       by = "Value.Qualifier"
#         ) %>% 
#     left_join(
#       x = .,
#       y = HESTIA_HC20_DB_endpoint_conversions %>%
#             group_by(Value.MaxQualifier) %>% 
#             summarise(n_max.qualifier = n()) %>% 
#             rename(Value.Qualifier = Value.MaxQualifier),
#       by = "Value.Qualifier"
#         ) %>% 
#   mutate(Value.Qualifier = case_when(is.na(Value.Qualifier) ~ "No qualifier", TRUE ~ Value.Qualifier)) %>% 
#   arrange(-n_val.qualifier) %>% 
#   write.csv(., "../data/excel_references/summary table for value.qualifiers.csv", row.names = F)
   

HESTIA_HC20_DB_endpoint_conversions %>% 
  filter(Qualifier == "Data in range") %>% 
  count(Value.MinQualifier, Value.MaxQualifier)

HESTIA_HC20_DB_endpoint_conversions %>% 
  filter(is.na(Value.MeanValue)) %>%
  count(Qualifier)

HESTIA_HC20_DB_qualifier_filter %>% 
  #filter(is.na(Value.MeanValue)) %>%
  filter(grepl("ignore", Value.MeanValue)) %>% 
  count(Value.MeanValue)


nrow(HESTIA_HC20_DB_endpoint_conversions) - nrow(HESTIA_HC20_DB_qualifier_filter)
