#install.packages("pracma")
#install.packages("proto")
#install.packages("nlraa")
library(dplyr)
library(pracma)
library(rlang)
library(nlraa) # <- contains the "predict2_nls" function
require(EnvStats) # <- to calculate Geometric St.Dev as a comparative unitless measure across all data

nls_across_shiny <- function(dataset, CAS, HCx = 20, MC_n = 10000) {
  options(dplyr.summarise.inform = FALSE)
  ### Loading functions to base nls on ###
  # Cumulative normal distribution function:
  cum_norm_dist_function <- function(x, mu, sig){
    0.5 + (0.5*erf((x - mu)/(sig*sqrt(2))))
  }
  
  ### Selfstart function ###
  cnormSS <- function(formula, data, start=NULL, weights=NULL){
    if(is.null(start)){
      x <- data[[all.vars(formula)[1]]]
      y <- data[[all.vars(formula)[2]]]
      mu <- mean(x, na.rm = TRUE)
      sig <- mean(y, na.rm = TRUE)
      if(is.na(mu)) mu <- 0
      if(is.na(sig)) sig <- 1
      start <- list(mu=mu, sig=sig)
    }
    return(start)
  }
  
  ### selecting data with n >=5 distinct Species and >= 3 taxonomic groups ###
  error_df <- dataset %>%
    filter(CAS.Number == {{CAS}}) %>% 
    count(CAS.Number, Taxonomy.Group, Species) %>%
    group_by(CAS.Number) %>% 
    summarise(n_sp = length(unique(Species)), 
              n_tax.grp = length(unique(Taxonomy.Group)), 
              n_recs = sum(n))
  

  output_list <- list(CAS.Number = NULL,
                    log_HCxEC10 = NULL,
                    non_W_HCxEC10 = NULL,
                    MC_logHC20ec10eq = NULL, 
                    Resp_lvl = {{HCx}},
                    CRF = NULL,
                    n_sp = NULL,
                    n_tax.grp = NULL,
                    n_recs = NULL,
                    mu = NULL,
                    sigma = NULL,
                    mean_HCx = NULL,
                    sd_HCx  = NULL,
                    # CV_HCx = NULL, 
                    # Phi_HCx = NULL, 
                    Q2.5 = NULL,
                    Q97.5 = NULL,
                    Geo_St.Dev = NULL,
                    status = NULL,
                    nls_results = NULL,
                    HCx_vec = NULL,
                    MC_mu = NULL,
                    MC_sig = NULL,
                    responseSequence = NULL,
                    predicted_values = NULL,
                    prediction_data = NULL
                    )

  
      ### Using the input dataset to get counts, means and Sd.
      d.frame_1 <- dataset %>%
         filter(CAS.Number == {{CAS}}) %>% 
        # Perform first averaging of species-specific tox data to get species mean and standard error
        mutate(Li = log10(EC10eq)) %>%
        group_by(Taxonomy.Group, Species) %>%
        # Log10-transforming effect concentrations, then calculating the arithmetic mean per species' log10 effect concentration
        summarise(
          n_samples = n(),
          sp_mean = mean(Li, na.rm = TRUE),
          sd_Li = sd(Li, na.rm = TRUE),
        ) %>%
        ungroup() 
      # Special dataset for non-weighted data
      d.frame_2 <- d.frame_1 %>%
        # Denoting the order of appearance in cumulative form:
        mutate(y_rank = (rank(sp_mean, na.last = NA, ties.method = "random") - 0.5)/length(unique(Species))) %>% 
        arrange(y_rank) %>%
        mutate(Taxonomy.Group = as.factor(Taxonomy.Group))
      
      # Duplicating d.frame_1 to allow non-weighted method to pull from all available data, as in HC20_calc.
      d.frame <- d.frame_1 %>% 
               #Need for removing cases where sigma (Sd_Li) is missing
                filter(!is.na(sd_Li), 
                        sd_Li != 0) %>% 
        mutate(y_rank = (rank(sp_mean, na.last = NA, ties.method = "random") - 0.5)/length(unique(Species))) %>% 
        arrange(y_rank) %>%
        mutate(Taxonomy.Group = as.factor(Taxonomy.Group))
      
      z_HCx <- sqrt(2)*erfinv(2*({{HCx}}/100)-1)
      output_list$non_W_HCxEC10 <- mean(d.frame_2$sp_mean, na.rm = T) + (z_HCx * sd(d.frame_2$sp_mean, na.rm = T))
      
      
      # adding a tryCatch to identify status occurring due to few records, resulting in non-convergence.
      output_list$status <- tryCatch({
  
      # The nonlinear least squares (nls()) function assuming a cumulative normal distribution
      nlc <- nls.control(maxiter = 1000, warnOnly = TRUE)
      nls_out <- nls(y_rank ~ cum_norm_dist_function(sp_mean, mu, sig), 
                     control = nlc, 
                     data=d.frame, 
                     weights = 1/(sd_Li^2),
                     start = cnormSS(sp_mean ~ sd_Li, 
                                     data = d.frame)
                     )
      if (nrow(d.frame)<5 | nrow(d.frame %>% distinct(Taxonomy.Group))<3) {
        "Convergence; Data insufficient"
      } else {
        "Convergence"
      }
      },
      error = function(e) {
        if (nrow(d.frame)<5 | nrow(d.frame %>% distinct(Taxonomy.Group))<3) {
          "Fail.Init; data insufficient"
        } else {
          "Fail.Init; data present"
        }
      },
      warning = function(w){
        if (nrow(d.frame) < 5 | nrow(d.frame %>% distinct(Taxonomy.Group)) <3) {
          "Fail.Conv; data insufficient"
        } else {
          "Fail.Conv; data present"
        }
      }
      )
      
      # Generate predictions from the nls model
      output_list$responseSequence <- data.frame(sp_mean = seq(min(d.frame$sp_mean), max(d.frame$sp_mean), length.out = 10000))
      output_list$predicted_values <- predict(nls_out, newdata = output_list$responseSequence)
      
      # Create a data frame with predictions
      output_list$prediction_data <- data.frame(sp_mean = output_list$responseSequence[[1]], predicted_values = predicted_values)
      
      # If the model fails, print number of records and terminate
      if (grepl("Fail", output_list$status)) {
        output_list$CAS.Number <- {{CAS}}
        output_list$n_sp <- as.numeric(error_df[i,"n_sp"])
        output_list$n_tax.grp <- as.numeric(error_df[i,"n_tax.grp"])
        output_list$n_recs <- as.numeric(error_df[i,"n_recs"])
        output_list$nls_results <- as.numeric(NA)
        
      } else { # Run Monte Carlo analysis on the distribution of HC20EC10eq data to assess uncertainty!
        output_list$nls_results <- list(nls_out)
        # assign CAS.Number to output df
        output_list$CAS.Number <- {{CAS}}
        
        # Defining the number of MC runs (MC_n) is done in function input.
        # Monte Carlo of the HC20EC10eq distribution for the mu using the Standard ERROR for sigma
        output_list$MC_mu <- list(rnorm({{MC_n}}, mean = summary(output_list$nls_results[[1]])$parameters[1,1], sd = summary(output_list$nls_results[[1]])$parameters[1,2]))
        # Monte carlo of the HC20EC10eq distribution for the sigma using the Standard ERROR for sigma
        output_list$MC_sig <- list(rnorm({{MC_n}}, mean = summary(output_list$nls_results[[1]])$parameters[2,1], sd = summary(output_list$nls_results[[1]])$parameters[2,2]))
        # these two vectors is used to construct HC20 using the qnorm()
        output_list$HCx_vec <- list(qnorm(({{HCx}}/100), mean = output_list$MC_mu[[1]], sd = output_list$MC_sig[[1]]))
        
        # assigning all the lists with respectively generated data
        # Automating extraction of the mu & sig nls results
        output_list$log_HCxEC10 <- qnorm(({{HCx}}/100), mean = coef(nls_out)[1], sd =  coef(nls_out)[2]) # gives me a point value of the corresponding HC20 (data fetched from console output)
        output_list$CRF <- ({{HCx}}/100)/(10^output_list$log_HCxEC10)
        output_list$mu <- coef(output_list$nls_results[[1]])[1]
        output_list$sigma <- coef(output_list$nls_results[[1]])[2]
        # Mean HC_x based on MC vector
        output_list$mean_HCx <- mean(output_list$HCx_vec[[1]], na.rm=T)
        output_list$sd_HCx <- sd(output_list$HCx_vec[[1]], na.rm=T)
        # output_df$CV_HCx <- (output_df$sigma/output_df$mu)*100
        # output_df$Phi_HCx <- sqrt(log((output_df$CV_HCx^2)+1))
        output_list$Q2.5 <- quantile(output_list$HCx_vec[[1]], 0.025, na.rm=T)
        output_list$Q97.5 <- quantile(output_list$HCx_vec[[1]], 0.975, na.rm=T)
        output_list$MC_logHC20ec10eq = output_list$mean_HCx + (-0.842 * output_list$sd_HCx)
        output_list$n_sp <- as.numeric(error_df[,"n_sp"])
        output_list$n_tax.grp <- as.numeric(error_df[,"n_tax.grp"])
        output_list$n_recs <- as.numeric(error_df[,"n_recs"])
        # Comparing variation across chemicals using the geoSD and first back-log-transform all the log-data
        output_list$Geo_St.Dev <- geoSD(10^output_list$HCx_vec[[1]], na.rm = TRUE)
        #output_list$HCx_vec <- output_list$HCx_vec[[1]] # for plotting a histogram over the HC20EC10eq Monte Carlo data distribution
      
        # Making all numbers to be presented as maxmum 4 decimals
        output_list[c(2:4, 6, 10:16)] <- lapply(output_list[c(2:4, 6, 10:16)], round, digits = 4)
        
  # Plotting output!
  if (class(output_list$nls_results[[1]]) == "nls") {
      d.frame <- cbind(d.frame, predict2_nls(nls_out, interval = "conf"))
      # Conditionally plot the confidence intervals at HC20 working point, based on if the confint() gives an output or throws an error due to "Conf.int infinity produced" (Starting values are out of bounds)
      c_int_on <- !is.na(output_list$Q2.5) && !is.na(output_list$Q97.5)
      if(c_int_on == FALSE) output_list$status <- "Conf.int contains 'NA'"
      if(nrow(d.frame)<5 | nrow(d.frame %>% distinct(Taxonomy.Group))<3) output_list$status <- "Convergence; Data insufficient"
      
     ## Catching statement if ggplot did output ssd curves
     # output_list$status <-  tryCatch({
      
      # Adding a conditional hjust-value for geom_textprinting text on either side of the geom_point.
      d.frame <- d.frame %>%
        mutate(hjust_value = ifelse(y_rank >= 0.5, 1.5, -0.5))
      
     # Generate a plot object that can be put in a list at the end of the function.
     plt <- ggplot(d.frame, aes(x = sp_mean, y = y_rank)) +
        geom_point(aes(fill = Taxonomy.Group), shape = 21, size = 3, alpha = 0.7) +
        geom_text(aes(label = Species, hjust = hjust_value, fontface=3)) +
     # Conditionally plotting percentiles if there are values available
          {if(c_int_on)geom_point(aes(x = output_list$Q2.5,
                         y = ({{HCx}}/100), color = "2.5 percentile", shape = "2.5 percentile"),
                     inherit.aes = F)} +
          {if(c_int_on)geom_point(aes(x = output_list$Q97.5,
                         y = ({{HCx}}/100), color = "97.5 percentile", shape = "97.5 percentile"),
                     inherit.aes = F)} +
           geom_smooth(
             data = output_list$prediction_data, 
             aes(x = sp_mean, y = output_list$predicted_values), 
             color = "red", 
             method = "auto", 
             se = FALSE) +
       # plotting the nls function using the same arguments as the nls model above. (DEPRECATED, This plots the curve of the nonweighted nls model)
          # geom_smooth(method = "nls",
          #             formula = y ~ cum_norm_dist_function(x, mu, sig),
          #             method.args = list(start = cnormSS(sp_mean ~ sd_Li, data = d.frame)),
          #             se =  FALSE) + # this is important
           #geom_ribbon(aes(ymin=Q2.5, ymax=Q97.5), alpha=0.2, na.rm = TRUE) +
       
          geom_point(aes(x = output_list$log_HCxEC10, y = ({{HCx}}/100), color = "logHC20EC10eq", shape = "logHC20EC10eq"), inherit.aes = F) +
          geom_point(aes(x = output_list$MC_logHC20ec10eq, y = ({{HCx}}/100), color = "MC_logHC20ec10eq", shape = "MC_logHC20ec10eq"), inherit.aes = F) +
          geom_point(aes(x = output_list$non_W_HCxEC10, y = ({{HCx}}/100), color = "non_W_HCxEC10", shape = "non_W_HCxEC10"), inherit.aes = F) +
          scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0", "20", "40", "60", "80", "100")) +
          scale_color_manual(
           name = element_blank(),
           aesthetics = c("color"),
           labels = if(c_int_on) {
             c("2.5 percentile", "97.5 percentile", "logHC20EC10eq", "MC_logHC20ec10eq", "non_W_HCxEC10") 
             } else {
               c("logHC20EC10eq", "MC_logHC20ec10eq", "non_W_HCxEC10")},
           values = if(c_int_on) {
             c("2.5 percentile" = "green", "97.5 percentile" = "blue", "logHC20EC10eq" = "red", "MC_logHC20ec10eq" = "purple", "non_W_HCxEC10" = "black") 
             } else {
               c("logHC20EC10eq" = "red", "MC_logHC20ec10eq" = "purple", "non_W_HCxEC10" = "black")}
          ) +
          scale_shape_manual(
           name = element_blank(),
           labels = if(c_int_on) {
             c("2.5 percentile", "97.5 percentile", "logHC20EC10eq", "MC_logHC20ec10eq", "non_W_HCxEC10") 
           } else {
             c("logHC20EC10eq", "MC_logHC20ec10eq", "non_W_HCxEC10")},
           values = if(c_int_on) {
             c("2.5 percentile" = 25, "97.5 percentile" = 24, "logHC20EC10eq" = 15, "MC_logHC20ec10eq" = 19, "non_W_HCxEC10" = 13) 
             } else {
               c("logHC20EC10eq" = 15, "MC_logHC20ec10eq" = 19, "non_W_HCxEC10" = 13)}
          ) +
          labs(title = paste("CAS", output_list$CAS.Number, sep = " "),
               subtitle = if(c_int_on) {
                  paste("log10HC", {{HCx}}, " = ", round(output_list$log_HCxEC10, digits = 4), sep = "")
               } else {
                 paste("log10HC", {{HCx}}, " = ", round(output_list$log_HCxEC10, digits = 4), ", DATA INSUFFICIENT. No CI available", sep = "")
                 },
               x = "mean(log) EC10eq (mg L-1)",
               y = "Response level (%)") +
          theme_linedraw() +
          theme(
            axis.text.y = element_text(),
            legend.background = element_blank(),
            legend.position = c(0.1, 0.75),
            legend.box.background = element_rect(fill = "white"),
            legend.box.margin = margin(-3,-3,-3,-3),
            legend.spacing.y = unit(1, "mm"),
            legend.key.size = unit(4, 'mm'),
            legend.title = element_text(size=10),
            legend.text = element_text(size=8)
          )
     
        #"OK"
      # }, error = function(Err) {
      #   "Error: Blank Plot. step factor reduced below 'minFactor'"
      # }, warning = function(Warn){
      #   "Warning: too many missing values"
      # } )
    
    } 
    #     else {
    #   warning(print(paste("No nls available")))
    # }
  }
#}
  
return(list(output_list, plt))
}
