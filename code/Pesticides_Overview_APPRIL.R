library(tidyverse)
library(xlsx)

APPRIL_DB <- read.csv("data/APPRIL_Database_2023-07-09.csv") %>% 
  mutate(
    Status.Date = as.Date(Status.Date, format = "%m/%d/%Y"),
    Date.First.Registered = as.Date(Date.First.Registered, format = "%m/%d/%Y"),
    Latest.Label.Date = as.Date(Latest.Label.Date, format = "%m/%d/%Y"),
    Sort.by.Reg.. = as.factor(Sort.by.Reg..)
    )


## Table of counts of Registered pesticides
APP_use_cat_tab <- APPRIL_DB %>% 
  filter(Status == "Registered") %>% 
  group_by(Pesticide.Category) %>% 
  #filter(!grepl("\\,", Pesticide.Category)) %>% 
  summarize(n = n()) %>% 
  mutate(Pesticide.Category = case_when(Pesticide.Category == "" ~ "No Category", TRUE ~ Pesticide.Category)) %>% 
  arrange(-n)

write.csv(APP_use_cat_tab, "figures/use_cat_list.csv", row.names = FALSE)

APPRIL_FULL_DB %>% 
  count(Status)

nrow(APPRIL_FULL_DB)
nrow(APPRIL_DB %>% 
  filter(Status != "Registered"))

APPRIL_FULL_DB <- read.csv("data/APPRIL_Database_FULL_2023-07-09.csv") %>% 
  mutate(
    Status.Date = as.Date(Status.Date, format = "%m/%d/%Y"),
    Date.First.Registered = as.Date(Date.First.Registered, format = "%m/%d/%Y"),
    Latest.Label.Date = as.Date(Latest.Label.Date, format = "%m/%d/%Y"),
    Sort.by.Reg.. = as.factor(Sort.by.Reg..)
  ) %>% 
  filter(Date.First.Registered != "0002-02-22")

# Summarize registration data per year
{
## Pesticides - Year registered
Y_reg <- APPRIL_FULL_DB %>% 
  group_by(Date.First.Registered) %>% 
  summarize(n = n()) %>% 
  #filter(Date.First.Registered != 0002-02-22) %>% 
  mutate(cum_n = cumsum(n),
         Status = "Registered") %>% 
  ungroup() %>% 
  rename(Year = Date.First.Registered)

## Pesticides - Year Conditionally Registered"
Y_conReg <- APPRIL_FULL_DB %>% 
  filter(Status == "Conditionally Registered") %>% 
  group_by(Status.Date) %>% 
  summarize(n = n()) %>% 
  #filter(Date.First.Registered != 2) %>% 
  mutate(cum_n = cumsum(n),
         Status = "Conditionally Registered") %>% 
  ungroup() %>% 
  rename(Year = Status.Date)

## Pesticides - Year canceled
Y_cancel <- APPRIL_FULL_DB %>% 
  filter(Status == "Canceled") %>% 
  group_by(Status.Date) %>% 
  summarize(n = n()) %>% 
  #filter(Date.First.Registered != 2) %>% 
  mutate(cum_n = cumsum(n),
         Status = "Cancelled") %>% 
  ungroup() %>% 
  rename(Year = Status.Date)

## Pesticides - Year Reregistered
Y_reReg <- APPRIL_FULL_DB %>% 
  filter(Status == "Reregistered") %>% 
  group_by(Status.Date) %>% 
  summarize(n = n()) %>% 
  #filter(Date.First.Registered != 2) %>% 
  mutate(cum_n = cumsum(n),
         Status = "Reregistered") %>% 
  ungroup() %>% 
  rename(Year = Status.Date)

## Pesticides - Year Conditionally Reregistered
Y_conReReg <- APPRIL_FULL_DB %>% 
  filter(Status == "Conditionally Reregistered") %>% 
  group_by(Status.Date) %>% 
  summarize(n = n()) %>% 
  #filter(Date.First.Registered != 2) %>% 
  mutate(cum_n = cumsum(n),
         Status = "Conditionally Reregistered") %>% 
  ungroup() %>% 
  rename(Year = Status.Date)

## Pesticides - Year Conditionally Reinstated
Y_reInst <- APPRIL_FULL_DB %>% 
  filter(Status == "Reinstated") %>% 
  group_by(Status.Date) %>% 
  summarize(n = n()) %>% 
  #filter(Date.First.Registered != 2) %>% 
  mutate(cum_n = cumsum(n),
         Status = "Reinstated") %>% 
  ungroup() %>% 
  rename(Year = Status.Date)
}

Pest_Reg_cumulative <- rbind(Y_reg,# Y_conReg, 
                             Y_cancel)
Pest_Reg_cumulative$Status <- factor(Pest_Reg_cumulative$Status, levels = c("Registered", 
                                                                            "Cancelled"#, 
                                                                            #"Conditionally Registered"
                                                                            ))

Pest_Reg_cumulative %>%   
ggplot(aes(x = Year, y = cum_n, color = Status)) +
  #geom_point(size = 1, color = "black", alpha=0.3) + 
  geom_line(size = 1.5) +
  #geom_area(aes(fill=Action)) +
  theme_bw() +
  
  theme(
        axis.text.x = element_text(size = 14, angle = 90), 
        axis.text.y = element_text(size = 14), 
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)
        ) +
  scale_x_date(date_labels = "%Y", date_breaks = "5 years") +
  scale_y_continuous(breaks = c(seq(0, 3.5e+5, 50000)), expand = c(0, 0), limits = c(0, 3.7e+5)) +
  ylab("Number of chemicals") #+
  #ggtitle("Cumulative number of registered vs. cancelled pesticides by the USEPA 1947-2023", subtitle = "Data obtained from the APPRIL database on July 9, 2023")
