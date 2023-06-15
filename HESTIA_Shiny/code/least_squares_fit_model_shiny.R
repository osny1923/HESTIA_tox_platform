#install.packages("pracma")
#install.packages("proto")
#install.packages("nlraa")
library(dplyr)
library(pracma)
library(rlang)
library(nlraa) # <- contains the "predict2_nls" function

# # #### TEST_bullshit
#  dataset <- HESTIA_BASE_EnviroTox_FILL %>% filter(CAS.Number == "100-00-5")
# 
# nls_data_test <- nls_across_shiny(dataset = dataset, "100-00-5", HCx = 20)
# nls_t <- nls_data_test[[1]]
# test <- data.frame(do.call(cbind, nls_t[1:11])) %>%
#   mutate(across(c(2:10), ~ as.numeric(.x)))

                             
nls_across_shiny <- function(dataset, CAS, HCx = 20) {
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
  
  ### create a template list-object to fill ###
  output_list <- list(CAS.Number = NULL, log_HC20EC10eq = NULL,
                      Effect_level = {{HCx}}, 
                    CRF  = NULL, mu  = NULL, sigma  = NULL,
                    Q2.5 = NULL, Q97.5 = NULL,
                    n_sp = NULL, n_tax.grp = NULL,
                    n_recs = NULL,
                    status = NULL, nls_results  = NULL)

    # If insufficient data (number of species <5 and/or number of taxonomic groups <3), print statement "not enough data" and move to the next substance. 
    if (error_df[,"n_sp"] <5 | error_df[,"n_tax.grp"] <3) {
      output_list$CAS.Number <- {{CAS}}
      output_list$log_HC20EC10eq <- as.numeric(NA)
      output_list$CRF <- as.numeric(NA)
      output_list$mu <- as.numeric(NA)
      output_list$sigma <- as.numeric(NA)
      output_list$Q2.5 <- as.numeric(NA)
      output_list$Q97.5 <- as.numeric(NA)
      output_list$status <- "not enough data"
      output_list$n_sp <- as.numeric(error_df[,"n_sp"])
      output_list$n_tax.grp <- as.numeric(error_df[,"n_tax.grp"])
      output_list$n_recs <- as.numeric(error_df[,"n_recs"])
      output_list$nls_results <- as.numeric(NA)
      # Create a blank plot
      plt <- ggplot() +
        theme_void() +
        labs(title = "WARNING: Too few data points to produce an SSD curve", 
             subtitle = NULL, 
             x = NULL, 
             y = NULL)
      
        
    }
    
    else {
      ### Using the input dataset to get counts, means and Sd.
      d.frame <- dataset %>%
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
        ungroup() %>%
        # Denoting the order of appearance in cumulative form:
        mutate(y_rank = (rank(sp_mean, na.last = NA) - 0.5)/length(unique(Species)),
               # Some occurrences where sd == 0 cause error in the nls model. making these into NA's
               sd_Li = case_when(sd_Li == 0 ~ as.numeric(NA), TRUE ~ sd_Li)) %>%
        arrange(y_rank) %>%
        mutate(Taxonomy.Group = as.factor(Taxonomy.Group))
     
      # adding a tryCatch to identify status occurring due to few records, resulting in non-convergence.
      output_list$status <- tryCatch({
        
      # The nonlinear least squares (nls()) function assuming a cumulative normal distribution
      nlc <- nls.control(maxiter = 1000, warnOnly = TRUE)
      nls_out <- nls(y_rank ~ cum_norm_dist_function(sp_mean, mu, sig), 
                     control = nlc, 
                     data=d.frame, 
                     start = cnormSS(sp_mean ~ sd_Li, 
                                     data = d.frame)
                     )
      "Convergence" # <- returning an OK in the status-column
      },
      error = function(e) {
        "singular gradient matrix at initial parameter estimates"
      },
      warning = function(w){
        "warning: step factor reduced below 'minFactor'"
      })
      if (output_list$status %in% c("warning: step factor reduced below 'minFactor'",
                                     "singular gradient matrix at initial parameter estimates")) {
        output_list$CAS.Number <- {{CAS}}
        output_list$n_sp <- as.numeric(error_df[i,"n_sp"])
        output_list$n_tax.grp <- as.numeric(error_df[i,"n_tax.grp"])
        output_list$n_recs <- as.numeric(error_df[i,"n_recs"])
        output_list$nls_results <- as.numeric(NA)
        
      } else {
        # assigning all the lists with respectively generated data
        output_list$nls_results <- list(nls_out)
        # assign CAS.Number to output df
        output_list$CAS.Number <- {{CAS}}
        # Automating extraction of the mu & sig nls results
        output_list$log_HC20EC10eq <- qnorm(({{HCx}}/100), mean = coef(nls_out)[1], sd =  coef(nls_out)[2]) # gives me a point value of the corresponding HC20 (data fetched from console output)
        output_list$CRF <- ({{HCx}}/100)/(10^output_list$log_HC20EC10eq)
        output_list$mu <- coef(nls_out)[1]
        output_list$sigma <- coef(nls_out)[2]
        output_list$n_sp <- as.numeric(error_df[,"n_sp"])
        output_list$n_tax.grp <- as.numeric(error_df[,"n_tax.grp"])
        output_list$n_recs <- as.numeric(error_df[,"n_recs"])
        
      catch_warning <- list(warn = NULL)
      catch_warning <- tryCatch({
        # Calculating the confidence intervals for mu and sigma at HC20 working point
        confint_res <- confint(output_list$nls_results[[1]], parm = c("mu", "sig"), level = 0.95)
        output_list$Q2.5 <- qnorm(({{HCx}}/100), mean = confint_res[1,1], sd = confint_res[2,1])
        output_list$Q97.5 <- qnorm(({{HCx}}/100), mean = confint_res[1,2], sd = confint_res[2,2], lower.tail = FALSE)
        # Adding "nothing" to the tryCatch output, since i don't want a warning inside the plot
        ""
      }, error = function(e) {
        "Conf.int infinity produced"
      },
      warning = function(w){
        "Conf.ints. warning"
      }
      )
  # Plotting output!
  if (class(output_list$nls_results[[1]]) == "nls") {
      d.frame <- cbind(d.frame, predict2_nls(nls_out, interval = "conf"))
      # Conditionally plot the confidence intervals at HC20 working point, based on if the confint() gives an output or throws an error due to "Conf.int infinity produced" (Starting values are out of bounds)
      c_int_on <- catch_warning == "" && !is.na(output_list$Q2.5) && !is.na(output_list$Q97.5)
      if(c_int_on == FALSE) catch_warning <- "Conf.int contains 'NA'"
      
     ## Catching statement if ggplot did output ssd curves
     output_list$status <-  tryCatch({
      
     # Generate a plot object that can be put in a list at the end of the function.
     plt <- ggplot(d.frame, aes(x = sp_mean, y = y_rank)) +
        geom_point(aes(fill = Taxonomy.Group), shape = 21, size = 1.5, alpha = 0.7) +
       
     # Conditionally plotting percentiles if there are values available
          {if(c_int_on)geom_point(aes(x = output_list$Q2.5,
                         y = ({{HCx}}/100), color = "low_CI", shape = "low_CI"),
                     inherit.aes = F)} +
          {if(c_int_on)geom_point(aes(x = output_list$Q97.5,
                         y = ({{HCx}}/100), color = "high_CI", shape = "high_CI"),
                     inherit.aes = F)} +
          
       # plotting the nls function usng the same arguments as the nls model above.
          geom_smooth(method = "nls",
                      formula = y ~ cum_norm_dist_function(x, mu, sig),
                      method.args = list(start = cnormSS(sp_mean ~ sd_Li, data = d.frame)),
                      se =  FALSE) + # this is important
          #geom_ribbon(aes(ymin=Q2.5, ymax=Q97.5), alpha=0.2, na.rm = TRUE) +
          geom_point(aes(x = output_list$log_HC20EC10eq, y = ({{HCx}}/100), color = "HC20EC10eq", shape = "HC20EC10eq"), inherit.aes = F) +
          scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0", "20", "40", "60", "80", "100")) +
          scale_color_manual(name = element_blank(),
                             aesthetics = c("color"),
                             labels = c("logHC20EC10eq", "97.5 percentile", "2.5 percentile"),
                             values = c("HC20EC10eq"="red", "low_CI"="green", "high_CI" = "blue")
          ) +
          scale_shape_manual(name = element_blank(),
                             labels = c("logHC20EC10eq", "97.5 percentile", "2.5 percentile"),
                             values = c(15, 17, 19)
          ) +
          labs(title = paste("CAS", output_list$CAS.Number, sep = " "),
               subtitle = paste("log10HC", {{HCx}}, "EC10eq"," = ", round(output_list$log_HC20EC10eq, digits = 4), " ", catch_warning, sep = ""),
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
            legend.title = element_text(size=8),
            legend.text = element_text(size=6)
          )
     
        #ggsave(filename = paste("SSD", output_list$CAS.Number,".png", sep = "_"), plot = last_plot(), device = "png", path = paste(Plot_destination, "/", sep = ""))
     
        "OK"
      }, error = function(Err) {
        "Error: Blank Plot. step factor reduced below 'minFactor'"
      }, warning = function(Warn){
        "Warning: too many missing values"
      } )
    
    } else {
      warning(print(paste("No nls available")))
    }
  }
}
  
return(list(output_list, plt))
}
