#install.packages("pracma")
#install.packages("proto")
#install.packages("nlraa")
library(tidyverse)
library(pracma)
library(rlang)
library(nlraa) # <- contains the "predict2_nls" function
library(EnvStats)

nls_across_all <- function(dataset, CAS = CAS.Number, Tax = Taxonomy.Group, EC = EC10eq, Sp = Species, HCx = 20, MC_n = 100000, rm_singles = c("YES", "NO"), Plot_output = c("YES", "NO"), Plot_destination = "folder") {
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
      #if(is.na(mu)) mu <- 0
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
  
  ### counting data to later sort data conditionally, with n >=5 distinct Species and >= 3 taxonomic groups ###
  error_df <- dataset %>%
    count({{CAS}}, {{Tax}}, {{Sp}}) %>%
    group_by({{CAS}}) %>%
    summarise(n_sp = n_distinct({{Sp}}), n_tax.grp = n_distinct({{Tax}}), n_recs = sum(n)) %>%
    arrange(CAS.Number)
  
  output_df <- tibble(
    CAS.Number = CAS_list,
    log_HC20EC10eq = rep(NA_real_, length(CAS_list)),
    CRF = rep(NA_real_, length(CAS_list)),
    mu = rep(NA_real_, length(CAS_list)),
    sigma = rep(NA_real_, length(CAS_list)),
    mean_HCx = rep(NA_real_, length(CAS_list)),
    sd_HCx = rep(NA_real_, length(CAS_list)),
    CV_HCx = rep(NA_real_, length(CAS_list)),
    Phi_HCx = rep(NA_real_, length(CAS_list)),
    Q2.5 = rep(NA_real_, length(CAS_list)),
    Q97.5 = rep(NA_real_, length(CAS_list)),
    MC_CRFQ2.5 = rep(NA_real_, length(CAS_list)),
    MC_CRFQ97.5 = rep(NA_real_, length(CAS_list)),
    n_sp = error_df$n_sp,
    n_tax.grp = error_df$n_tax.grp,
    n_recs = error_df$n_recs,
    GStDev = rep(NA_real_, length(CAS_list)),
    Iterations.to.Convergence = rep(NA_real_, length(CAS_list)),
    Achieved.convergence.tolerance = rep(NA_real_, length(CAS_list)),
    status = rep("", length(CAS_list)),
    nls_results = vector("list", length = length(CAS_list)),
    MC_mu = vector("list", length = length(CAS_list)),
    MC_sig = vector("list", length = length(CAS_list)),
    log_HCx_vec = vector("list", length = length(CAS_list)),
    MC_CRF = vector("list", length = length(CAS_list)),
    n_sigma = rep(NA_real_, length(CAS_list)),
    n_taxa_sigma = rep(NA_real_, length(CAS_list))
  )
  
  ### Start looping each substance ###
  # Producing a list object where data per substance is gathered, along with the nls models in last list-column.
  for (i in seq_along(CAS_list)) {
    # If insufficient data (n = 1), print statement "not enough data" and move to the next substance.
    if (error_df[i,"n_sp"] == 1 | error_df[i,"n_tax.grp"] == 1) {
      output_df$status[i] <- "Not enough data"
      next
    } else {
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
        mutate(y_rank = (rank(sp_mean, na.last = NA) - 0.5)/length(unique({{Sp}})
        )
        ) %>% 
        arrange(y_rank) %>%
        mutate(Taxonomy.Group = as.factor(Taxonomy.Group)) 
      # Remove all data with missing sigma or sigma = 0 
      if (rm_singles == "YES") {
        d.frame <- d.frame %>% 
          filter(!is.na(sd_Li), 
                 sd_Li != 0)
        output_df$n_sigma[i] <- nrow(d.frame)
        output_df$n_taxa_sigma[i] <- nrow(d.frame %>% distinct(Taxonomy.Group))
      }
      
      # Check if sigma == 0 or NA, if the subset is reduced to 1 or fewer rows 
      if (nrow(d.frame) <= 1) {
        output_df$status[i] <- "One sigma or fewer"
        next
      }
      
      # adding a tryCatch to identify status occurring due to few records, resulting in non-convergence.
      output_df$status[i] <- tryCatch({
        
        # The nonlinear least squares (nls()) function assuming a cumulative normal distribution
        nlc <- nls.control(maxiter = 1000, warnOnly = TRUE)
        nls_out <- nls(y_rank ~ cum_norm_dist_function(sp_mean, mu, sig), 
                       control = nlc, 
                       weights = 1/(sd_Li),
                       data=d.frame, 
                       start = cnormSS(sp_mean ~ sd_Li, 
                                       data = d.frame)
        )
        if ((error_df[i,"n_sp"] >= 5 | error_df[i,"n_tax.grp"] >= 3) & (nrow(d.frame)<5|nrow(d.frame %>% distinct(Taxonomy.Group))<3)) {
          "Converge:<5 sigma|<3 tax.grp"
        } else if (error_df[i,"n_sp"] < 5 | error_df[i,"n_tax.grp"] < 3) {
          "Data insufficient"
        } else {
          "Convergence"
        }
      },
      error = function(e) {
        if ((error_df[i,"n_sp"] >= 5 | error_df[i,"n_tax.grp"] >= 3) & (nrow(d.frame)<5 | nrow(d.frame %>% distinct(Taxonomy.Group))<3)) {
          "Fail.Init.<5 sigma|<3 tax.grp"
        } else if (error_df[i,"n_sp"] >= 5 | error_df[i,"n_tax.grp"] >= 3) {
          "Fail.Init.data insufficient"
        } else {
          "Fail.Init.full.data.present"
        }
      },
      warning = function(w){
        if ((error_df[i,"n_sp"] >= 5 | error_df[i,"n_tax.grp"] >= 3) & (nrow(d.frame) < 5 | nrow(d.frame %>% distinct(Taxonomy.Group)) <3)) {
          "Fail.Conv.<5 sigma|<3 tax.grp"
        } else if (error_df[i,"n_sp"] >= 5 | error_df[i,"n_tax.grp"] >= 3) {
          "Fail.Conv.data insufficient"
        } else {
          "Fail.Conv.full.data.present"
        }
      }
      )
      if (output_df$status[i] %in% c("Fail.Conv.full.data.present",
                                     "Fail.Init.data insufficient",
                                     "Fail.init.<5 sigma|<3 tax.grp",
                                     "Fail.Conv.<5 sigma|<3 tax.grp")) {
        output_df$CAS.Number[i] <- CAS_list[i]
        output_df$n_sp[i] <- as.numeric(error_df[i,"n_sp"])
        output_df$n_tax.grp[i] <- as.numeric(error_df[i,"n_tax.grp"])
        output_df$n_recs[i] <- as.numeric(error_df[i,"n_recs"])
        
        next
      } else {
        
        # assigning all the lists with respectively generated data
        output_df$nls_results[i] <- list(nls_out)
        # assign CAS.Number to output df
        output_df$CAS.Number[i] <- CAS_list[i]
        # Automating extraction of the mu & sig nls results
        output_df$log_HC20EC10eq[i] <- qnorm(({{HCx}}/100), mean = coef(output_df$nls_results[[i]])[1], sd =  coef(output_df$nls_results[[i]])[2]) # gives me a point value of the log10(HC20EC10eq) (data fetched from model output)
        output_df$CRF[i] <- ({{HCx}}/100)/(10^output_df$log_HC20EC10eq[i])
        output_df$mu[i] <- coef(nls_out)[1]
        output_df$sigma[i] <- coef(nls_out)[2]
        output_df$CV_HCx[i] <- (output_df$sigma[i]/output_df$mu[i])*100
        output_df$Phi_HCx[i] <- sqrt(log((output_df$CV_HCx[i]^2)+1))
        output_df$n_sp[i] <- as.numeric(error_df[i,"n_sp"])
        output_df$n_tax.grp[i] <- as.numeric(error_df[i,"n_tax.grp"])
        output_df$n_recs[i] <- as.numeric(error_df[i,"n_recs"])
        
        ##############
        # Monte Carlo of the HC20EC10eq distribution for the mu using the Standard ERROR for sigma
        ###############
        # Defining the number of runs is done in the function using parameter (MC_n)
        output_df$MC_mu[i] <- list(rnorm({{MC_n}}, mean = summary(output_df$nls_results[[i]])$parameters[1,1], sd = summary(output_df$nls_results[[i]])$parameters[1,2]))
        # Monte carlo of the HC20EC10eq distribution for the sigma using the Standard ERROR for sigma
        output_df$MC_sig[i] <- list(rnorm({{MC_n}}, mean = summary(output_df$nls_results[[i]])$parameters[2,1], sd = summary(output_df$nls_results[[i]])$parameters[2,2]))
        # these two vectors is used to construt HC20 using the qnorm()
        output_df$log_HCx_vec[i] <- list(qnorm(({{HCx}}/100), mean = output_df$MC_mu[[i]], sd = output_df$MC_sig[[i]]))
        output_df$mean_HCx[i] <- mean(output_df$log_HCx_vec[[i]], na.rm = T)
        output_df$sd_HCx[i] <- sd(output_df$log_HCx_vec[[i]], na.rm = T)
        output_df$GStDev[i] <- geoSD(10^output_df$log_HCx_vec[[i]], na.rm = T)
        output_df$Q2.5[i] <- quantile(output_df$log_HCx_vec[[i]], 0.025, na.rm = T)
        output_df$Q97.5[i] <- quantile(output_df$log_HCx_vec[[i]], 0.975, na.rm = T)
        output_df$Iterations.to.Convergence[i] <- summary(output_df$nls_results[[i]])[7]$convInfo$finIter
        output_df$Achieved.convergence.tolerance[i] <-summary(output_df$nls_results[[i]])[7]$convInfo$finTol
        output_df$MC_CRF[i] <- list(0.2/10^output_df$log_HCx_vec[[i]])
        output_df$MC_CRFQ2.5[i] <- quantile(output_df$MC_CRF[[i]], 0.025, na.rm = T)
        output_df$MC_CRFQ97.5[i] <- quantile(output_df$MC_CRF[[i]], 0.975, na.rm = T)
        gc(verbose = F)
        # 
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
              "Convergence"
            }, error = function(Err) {
              "Error: Blank Plot. step factor reduced below 'minFactor'"
            }, warning = function(Warn){
              "Warning: too many missing values"
            } )
          } else {
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
