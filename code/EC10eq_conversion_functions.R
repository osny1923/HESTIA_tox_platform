# Creating the endpoint conversion function

# I have a data frame describing which endpoints to harmonize and which to exclude. 
# I just need to figure out how to pair an endpoint to the harmonized endpoint within another data frame. 
# The dplyr way would be to do a `mutate(harmonized_endpoint = case_when(Endpoint == Endpoints_filter$Endpoint ~ Endpoints_filter$Group), TRUE ~ as.integer(NA))` by using the indexed parameters from the Endpoints_filter df.
# Another dplyr way is to read the conversion table and simply use left_join() - probably a much faster way of doing it, but these methods relies on reading in a dataframe for conversions. 

# Single value functions, not minding the whole vector of comparisons

# function to convert input ECx endpoints into EC50/EC10/NOEC
endpoint_conv_function <- function(e_p){
  # Defining Endpoints-filter data frame here so the function becomes portable.
  filter_df <- data.frame(Endpoint = c("NOEC","LC50","LOEC","EC50","NOEL","IC50","EC10","LC10","LD50","LC0","EC0","IC10","NOER","ER50","LD0","LD10","ER10","EL50","NOAEC","EL10"),
                          Group = c("NOEC","EC50","EC10","EC50","NOEC","EC50","EC10","EC10","EC50","NOEC","NOEC","EC10","NOEC","EC50","NOEC","EC10","EC10","EC50","NOEC","EC10")
                          ) 
  converted_endpoint <- as.character(NA) # defining output as an integer vector
  l_ep <- length(filter_df$Group) # length of the vector to check conversions against
  i = 1 # Ticks for the filter_df length
  for (i in 1:l_ep) {
    if ((is.na(e_p)) | (!e_p %in% filter_df$Endpoint) ) {
      converted_endpoint = "Endpoint is NA or invalid"
      i = 1
    } else if (e_p == filter_df$Endpoint[i]) {
      converted_endpoint <- filter_df$Group[i] # if statement is true, then populate the vector "result" with the value from "endpoints_filter$group".
      i = 1
      } else if (i <= l_ep) {
        i <- i + 1
        } 
  }
return(converted_endpoint)
}


# Taxonomic conversions function that applies to the extrapolation factors used in Aurisano et al., 2019
taxa_conversion_function <- function(tax){
  if (tax %in% c("Fish", "Crustacean", "Algae")) { # if the input corresponds to the vector
    other_vector = tax # then use that value as output
} else{
  other_vector = "Others" # otherwise, call the output value as "Others"
  } 
  return(other_vector)
}


# Function to Extrapolating effect concentrations of all ECx data into EC10eq
# `input_type` needs to be defined as either "extpl", "high_CI" or "low_CI"!
ec10eq_extrapolation_function <- function(x, endp, ac, tax, input_type ) { # This function will rely on the effect concentration value and information within a data frame to apply extrapolation factors
  # Here is the table of conversions defined in Aurisano et al., 2019
  expf_df <- data.frame(q = c(rep("EC50", 8), rep("EC10", 8), rep("NOEC", 8)), # Endpoints (q)
                        a = c(rep(c(rep("Acute", 4), rep("Chronic", 4)), 3)), # Acute or chronic exposure (a)
                        t = c(rep(c("Fish", "Crustacean", "Algae", "Others"), 6)), # Taxonomy group (t)
                        g = as.numeric(c(7.44, 3.38, 4, 4, 1.55, 1.94, 2.24, 2, rep(1, 8), 3.97, 1.55, 1.8, 1.8, 0.6, 0.95, 0.44, 0.6)), # extrapolation factor (g)
                        h = as.numeric(c(18.95, 5.34, 6.1, 6.1, 3.66, 2.41, 2.65, 2.5, rep(1, 8), 17.39, 2.64, 2.7, 2.7, 0.7, 1.16, 0.49, 0.7)), # High CI (95%) extrapolation factor (h)
                        l = as.numeric(c(2.92, 2.14, 2.6, 2.6, 0.67,  1.56, 1.9, 1.8, rep(1, 8), 0.9, 0.91, 1, 1, 0.4, 0.77, 0.39, 0.4)) # Low CI (95%) extrapolation factor (l)  
                        )
  
  if (input_type == "high_CI") expf_df$g <- expf_df$h 
  if (input_type == "low_CI") expf_df$g <- expf_df$l 
  if (!input_type %in% c("extpl", "high_CI", "low_CI")) warning("Input type is neither 'extpl', 'high' or 'low'. n/ Need to define if input should be converted as extrapolation factor or high or low CI") 
  
  extpl_endpoint <- as.numeric(NA) # pre-allocating an output-vector
  l_expf_df <- length(expf_df$q) # length of the conversion df to delimit the lookup
  p = 1 # ticks for expl_df length
  ep <- endpoint_conv_function(endp) # Applying the endpoint conversion function to input endpoints -> output either EC50/EC10/NOEC
  Tax <- taxa_conversion_function(tax) # Defining Taxonomy.Group according to the Taxonomy in Aurisano et al., 2019.
  for (p in 1:l_expf_df) {
    if (is.na(ep) | is.na(ac) | is.na(Tax)) {
        extpl_endpoint = as.numeric(NA)
        p = 1
    } else if ((ep == expf_df$q[p]) & (ac == expf_df$a[p]) & (Tax == expf_df$t[p])) {
      extpl_endpoint = x / expf_df$g[p]
      p = 1
      } else if (p <= l_expf_df) {
        p = p + 1
        } else {
          extpl_endpoint = "NA' last else - broken"
          p = 1
          }
  }
return(extpl_endpoint)
}

