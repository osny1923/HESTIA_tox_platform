library(dplyr)

use_prio_function <- function(x) {
  col_to_count <- c("ChEBI_DB", "ATC_Type", "USEPA_ecotox_group", 
                    "BCPC_pesticide", "DRUGBANK", "EPAPCS", 
                    "HEALTHY_BUILDING_NETWORK", "NORMAN_ITN", "OPPIN", "PPDB", 
                    "USEPA", "Wikipedia") 
  
prio_vec <- matrix(nrow = nrow(x), ncol = 2)

prio_vec <- t(sapply(1:nrow(x), function(i) {
  count_df <- x[i,] %>%
    select(all_of(col_to_count)) %>%
    mutate(across(all_of(col_to_count), ~ na_if(., ''))) %>%
    mutate(across(all_of(col_to_count), ~ na_if(., 'NA'))) %>%
    mutate(across(all_of(col_to_count), ~ tolower(.))) %>%
    mutate(across(all_of(col_to_count), ~ gsub(".*cide.*", "pesticide", .))) %>%
    mutate(across(all_of(col_to_count), ~ gsub(".*antib.*", "antibiotic", .))) %>%
    pivot_longer(cols = all_of(col_to_count)) %>%
    group_by(value) %>%
    na.omit(value) %>%
    summarize(n = n())%>% 
    arrange(desc(n), value) %>% 
    rename(Select_use = value)
  
  # Fetch the top row using count_df[1, c(1:2)]
  top_row <- count_df[1, c(1:2)]
  
  # Check if there is a tie in the top rows of column "n"
  if (sum(count_df$n == top_row$n) > 1) {
    # Check if "pesticide" or "antibiotic" is among the top tied strings
    if ("pesticide" %in% count_df$Select_use[which(count_df$n == top_row$n)]) {
      # Set top_row$Select_use to "pesticide"
      count_df$Select_use <- ifelse(count_df$Select_use == "pesticide", "pesticide", "pesticide tie")
    } else if ("antibiotic" %in% count_df$Select_use[which(count_df$n == top_row$n)]) {
      # Set top_row$Select_use to "antibiotic"
      count_df$Select_use <- ifelse(count_df$Select_use == "antibiotic", "antibiotic", "antibiotic tie")
    } else {
      # Add "tie" to the fetched cell (count_df[1,1])
      count_df$Select_use <- ifelse(count_df$Select_use == top_row$Select_use, paste(top_row$Select_use, "tie", sep = " "), count_df$Select_use)
    }
  }
  
  data.frame(Select_use = as.character(count_df[1, "Select_use"]), n = as.numeric(count_df[1, "n"]), stringsAsFactors = FALSE)
}))

# Replace "NA" or "NANA" strings with empty cells
prio_vec[prio_vec == "NA" | prio_vec == "NANA"] <- ""

prio_vec <- as.data.frame(prio_vec, stringsAsFactors = FALSE)

return(prio_vec)
}
