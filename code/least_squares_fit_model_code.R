#install.packages("pracma")
#install.packages("proto")
#install.packages("nlraa")
library(tidyverse)
library(pracma)
library(rlang)
library(nlraa) # <- contains the "predict2_nls" function

nls_across_all <- function(dataset, CAS = CAS.Number, Tax = Taxonomy.Group, EC = EC10eq, Sp = Species, HCx = 20, Plot_output = c("YES", "NO"), Plot_destination = "folder") {
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
  
  ### making a CAS.Number vector to loop over ###
  CAS_list <- dataset %>%
    distinct({{CAS}}) %>%
    arrange(CAS.Number) %>%
    pull({{CAS}})
  
  ### selecting data with n >=5 distinct Species and >= 3 taxonomic groups ###
  error_df <- dataset %>%
    count({{CAS}}, {{Tax}}, {{Sp}}) %>%
    group_by({{CAS}}) %>%
    summarise(n_sp = length(unique(Species)), n_tax.grp = length(unique(Taxonomy.Group)), n_recs = sum(n)) %>%
    arrange(CAS.Number)
  
  ### create a template list-object to fill ###
  output_df <- list(CAS.Number = NULL, log_HC20EC10eq = NULL,
                    CRF  = NULL, mu  = NULL, sigma  = NULL,
                    Q2.5 = NULL, Q97.5 = NULL,
                    n_sp = NULL, n_tax.grp = NULL,
                    n_recs = NULL,
                    status = NULL, nls_results  = NULL)
  # Create a dummy integer "i"
  i = 1
  
  ### Start looping each substance ###
  # Producing a list object where data per substance is gathered, along with the nls models in last list-column.
  while (i <= length(CAS_list)) {
    # If insufficient data (number of species <5 and/or numbre of taxonomic groups <3), print statement "not enough data" and move to the next substance. 
    if (error_df[i,"n_sp"] <5 | error_df[i,"n_tax.grp"] <3) {
      output_df$CAS.Number[i] <- CAS_list[i]
      output_df$log_HC20EC10eq[i] <- as.numeric(NA)
      output_df$CRF[i] <- as.numeric(NA)
      output_df$mu[i] <- as.numeric(NA)
      output_df$sigma[i] <- as.numeric(NA)
      output_df$Q2.5[i] <- as.numeric(NA)
      output_df$Q97.5[i] <- as.numeric(NA)
      output_df$status[i] <- "not enough data"
      output_df$n_sp[i] <- as.numeric(error_df[i,"n_sp"])
      output_df$n_tax.grp[i] <- as.numeric(error_df[i,"n_tax.grp"])
      output_df$n_recs[i] <- as.numeric(error_df[i,"n_recs"])
      output_df$nls_results[i] <- as.numeric(NA)
      i = i+1
    }
    
    else {
      ### Using the input dataset to get counts, means and Sd.
      d.frame <- dataset %>%
        filter({{CAS}} == CAS_list[i]) %>%
        # Perform first averaging of species-specific tox data to get species mean and standard error
        mutate(Li = log10({{EC}})) %>%
        group_by({{Tax}}, {{Sp}}) %>%
        # Log10-transforming effect concentrations, then calculating the arithmetic mean per species' log10 effect concentration
        summarise(
          n_samples = n(),
          sp_mean = mean(Li, na.rm = TRUE),
          sd_Li = sd(Li, na.rm = TRUE),
        ) %>%
        ungroup() %>%
        # Denoting the order of appearance in cumulative form:
        mutate(y_rank = (rank(sp_mean, na.last = NA) - 0.5)/length(unique({{Sp}})),
               # Some occurrences where sd == 0 cause error in the nls model. making these into NA's
               sd_Li = case_when(sd_Li == 0 ~ as.numeric(NA), TRUE ~ sd_Li)) %>%
        arrange(y_rank) %>%
        mutate(Taxonomy.Group = as.factor(Taxonomy.Group))
     
      # adding a tryCatch to identify status occurring due to few records, resulting in non-convergence.
      output_df$status[i] <- tryCatch({
        
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
      if (output_df$status[i] %in% c("warning: step factor reduced below 'minFactor'",
                                     "singular gradient matrix at initial parameter estimates")) {
        output_df$CAS.Number[i] <- CAS_list[i]
        output_df$n_sp[i] <- as.numeric(error_df[i,"n_sp"])
        output_df$n_tax.grp[i] <- as.numeric(error_df[i,"n_tax.grp"])
        output_df$n_recs[i] <- as.numeric(error_df[i,"n_recs"])
        output_df$nls_results[i] <- as.numeric(NA)
        i = i+1
      } else {
        # assigning all the lists with respectively generated data
        output_df$nls_results[i] <- list(nls_out)
        # assign CAS.Number to output df
        output_df$CAS.Number[i] <- CAS_list[i]
        # Automating extraction of the mu & sig nls results
        output_df$log_HC20EC10eq[i] <- qnorm(({{HCx}}/100), mean = coef(nls_out)[1], sd =  coef(nls_out)[2]) # gives me a point value of the corresponding HC20 (data fetched from model output)
        output_df$CRF[i] <- ({{HCx}}/100)/(10^output_df$log_HC20EC10eq[i])
        output_df$mu[i] <- coef(nls_out)[1]
        output_df$sigma[i] <- coef(nls_out)[2]
        output_df$n_sp[i] <- as.numeric(error_df[i,"n_sp"])
        output_df$n_tax.grp[i] <- as.numeric(error_df[i,"n_tax.grp"])
        output_df$n_recs[i] <- as.numeric(error_df[i,"n_recs"])
        
      catch_warning <- list(warn = NULL)
      catch_warning[i] <- tryCatch({
        # Calculating the confidence intervals for mu and sigma at HC20 working point
        confint_res <- confint(output_df$nls_results[[i]], parm = c("mu", "sig"), level = 0.95)
        output_df$Q2.5[i] <- qnorm(({{HCx}}/100), mean = confint_res[1,1], sd = confint_res[2,1])
        output_df$Q97.5[i] <- qnorm(({{HCx}}/100), mean = confint_res[1,2], sd = confint_res[2,2], lower.tail = FALSE)
        # Adding "nothing" to the tryCatch output, since i don't want a warning inside the plot
        ""
      }, error = function(e) {
        "Conf.int infinity produced"
      },
      warning = function(w){
        "Conf.ints. warning"
      }
      )
        # Do i want plot output or not?
        # If i don't want plots, i skip nls model fit calculations and ggplot operation
        if (Plot_output == "YES") {
          if (class(output_df$nls_results[[i]]) == "nls") {
            d.frame <- cbind(d.frame, predict2_nls(nls_out, interval = "conf"))
            # Conditionally plot the confidence intervals at HC20 working point, based on if the confint() gives an output or throws an error due to "Conf.int infinity produced" (Starting values are out of bounds)
            c_int_on <- catch_warning[i] == "" && !is.na(output_df$Q2.5[i]) && !is.na(output_df$Q97.5[i])
            if(c_int_on == FALSE) catch_warning[i] <- "Conf.int contains 'NA'"
            ## ggplot output ssd curves
            output_df$status[i] <-  tryCatch({
              ggplot(d.frame, aes(x = sp_mean, y = y_rank)) +
                geom_point(aes(fill = Taxonomy.Group), shape = 21, size = 1.5, alpha = 0.7) +
            # conditionally plotting precentiles if there are values available
                {if(c_int_on)geom_point(aes(x = output_df$Q2.5[i],
                               y = ({{HCx}}/100), color = "low_CI", shape = "low_CI"),
                           inherit.aes = F)} +
                {if(c_int_on)geom_point(aes(x = output_df$Q97.5[i],
                               y = ({{HCx}}/100), color = "high_CI", shape = "high_CI"),
                           inherit.aes = F)} +
                # plotting the nls function usng the same arguments as the nls model above.
                geom_smooth(method = "nls",
                            formula = y ~ cum_norm_dist_function(x, mu, sig),
                            method.args = list(start = cnormSS(sp_mean ~ sd_Li, data = d.frame)),
                            se =  FALSE) + # this is important
                geom_ribbon(aes(ymin=Q2.5, ymax=Q97.5), alpha=0.2) +
                geom_point(aes(x = output_df$log_HC20EC10eq[i], y = ({{HCx}}/100), color = "HC20EC10eq", shape = "HC20EC10eq"), inherit.aes = F) +
                scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0", "20", "40", "60", "80", "100")) +
                scale_color_manual(name = element_blank(),
                                   aesthetics = c("color"),
                                   labels = c("logHC20EC10eq", "2.5 quantile", "97.5 quantile"),
                                   values = c("HC20EC10eq"="red", "low_CI"="blue", "high_CI" = "blue")
                ) +
                scale_shape_manual(name = element_blank(),
                                   labels = c("logHC20EC10eq", "2.5 quantile", "97.5 quantile"),
                                   values = c(15, 8, 8)
                ) +
                labs(title = paste("CAS", output_df$CAS.Number[i], sep = " "),
                     subtitle = paste("log10HC", {{HCx}}, "EC10eq"," = ", round(output_df$log_HC20EC10eq[i], digits = 4), " ", catch_warning[i], sep = ""),
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
              ggsave(filename = paste("SSD", output_df$CAS.Number[i],".png", sep = "_"), plot = last_plot(), device = "png", path = paste(Plot_destination, "/", sep = ""))
              "OK"
            }, error = function(Err) {
              "Error: Blank Plot. step factor reduced below 'minFactor'"
            }, warning = function(Warn){
              "Warning: too many missing values"
            } )
          }else {
            warning(print(paste("No nls available")))
          }
        } # If plot_output == "YES"
        # Increase the tick by one
        i = i+1
      }
    }
  }
  return(output_df)
}
