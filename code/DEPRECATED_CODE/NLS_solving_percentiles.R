#install.packages("pracma")
#install.packages("proto")
#install.packages("nlraa")
library(tidyverse)
library(pracma)
library(rlang)
library(nlraa) # <- contains the "predict2_nls" function
require(EnvStats) # <- to calculate Geometric St.Dev as a comparative unitless measure across all data
# "100-25-4"
test_df <- read.csv("results/FINAL_HESTIA.csv")  %>% filter(CAS.Number == "27554-26-3")

#Loading all prerequisites
{
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
      #if(is.na(sig)) sig <- 1
      start <- list(mu=mu, sig=sig)
    }
    return(start)
  }
  
  ### making a CAS.Number vector to loop over ###
  CAS_list <- test_df[1,1]
  
  ### selecting data with n >=5 distinct Species and >= 3 taxonomic groups ###
  error_df <- test_df %>%
    count(CAS.Number, Taxonomy.Group, Species) %>%
    group_by(CAS.Number) %>%
    summarise(n_sp = n_distinct(Species), n_tax.grp = n_distinct(Taxonomy.Group), n_recs = sum(n)) %>%
    arrange(CAS.Number)
 
  
  ### create a template list-object to fill ###
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
}

   {
     i = 1
      ### Using the input dataset to get counts, means and Sd.
      d.frame <- test_df %>%
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
        mutate(y_rank = (rank(sp_mean, na.last = NA) - 0.5)/length(unique(Species))
               # Some occurrences where sd == 0 or missing, to ensure weighting to work, transforming sd == 1. 
               # sd_Li = case_when(sd_Li == 0 ~ 1,
               #                   is.na(sd_Li) ~ 1,
               #                   TRUE ~ sd_Li)
               ) %>%
        arrange(y_rank) %>%
        mutate(Taxonomy.Group = as.factor(Taxonomy.Group)) %>% 
        # Removing all single species from the each chemical-subset
        filter(!is.na(sd_Li), 
               sd_Li != 0)
      
      output_df$n_sigma[i] <- nrow(d.frame)
      output_df$n_taxa_sigma[i] <- nrow(d.frame %>% distinct(Taxonomy.Group))
      
        
      
    
         
        # defining the inverse of the standard normal distribution at the 0.2 probability level (z_HCx)
        z_HCx <- sqrt(2)*erfinv(2*(20/100)-1)
        output_df$non_W_HCxEC10 <- mean(d.frame$sp_mean, na.rm = T) + (z_HCx * sd(d.frame$sp_mean, na.rm = T))
      
        # The nonlinear least squares (nls()) function assuming a cumulative normal distribution
        nlc <- nls.control(maxiter = 1000, warnOnly = TRUE)
        nls_out <- nls(y_rank ~ cum_norm_dist_function(sp_mean, mu, sig), 
                       control = nlc, 
                       data = d.frame, 
                       weights = 1/(sd_Li^2), 
                       start = cnormSS(sp_mean ~ sd_Li, 
                                       data = d.frame))
        
        # Generate predictions from the nls model
        #new_data <- data.frame(sp_mean = seq(min(d.frame$sp_mean), max(d.frame$sp_mean), length.out = 1000))
        #predicted_values <- predict(nls_out, newdata = new_data)
        
        # Create a data frame with predictions
        #prediction_data <- data.frame(sp_mean = new_data$sp_mean, predicted_values = predicted_values)
        
        {
        # assigning all the lists with respectively generated data
        output_df$nls_results[i] <- list(nls_out)
        # assign CAS.Number to output df
        output_df$CAS.Number[i] <- CAS_list[i]
        # Automating extraction of the mu & sig nls results
        output_df$log_HC20EC10eq[i] <- qnorm((20/100), mean = coef(nls_out)[1], sd = coef(nls_out)[2]) # gives me a point value of the corresponding HC20 (data fetched from console output)
        output_df$CRF[i] <- (20/100)/(10^output_df$log_HC20EC10eq[i])
        output_df$mu[i] <- coef(nls_out)[1]
        output_df$sigma[i] <- coef(nls_out)[2]
        
        output_df$n_sp[i] <- as.numeric(error_df[i,"n_sp"])
        output_df$n_tax.grp[i] <- as.numeric(error_df[i,"n_tax.grp"])
        output_df$n_recs[i] <- as.numeric(error_df[i,"n_recs"])
        
        catch_warning <- list(warn = "")
    
        summary(nls_out)$parameters
        d.frame
        #summary(nls_out)$parameters[1,2]
        # summary(nls_out)$parameters[2,2]
        #coef(nls_out)[1]
        # Define the number of runs
        MC_n = 10000  
        # Monte carlo of the HC20EC10eq distribution for the mu using the Standard ERROR for sigma
        MC_mu <- rnorm(MC_n, mean = summary(nls_out)$parameters[1,1], sd = summary(nls_out)$parameters[1,2])
        # Monte carlo of the HC20EC10eq distribution for the sigma using the Standard ERROR for sigma
        MC_sig <- rnorm(MC_n, mean = summary(nls_out)$parameters[2,1], sd = summary(nls_out)$parameters[2,2])
        # these two vectors is used to construt HC20 using the qnorm()
        HC20_vec <- qnorm((20/100), mean = MC_mu, sd = MC_sig)
        # hist(HC20_vec)
        output_df$HC20_vec[i] <- list(qnorm((20/100), mean = MC_mu, sd = MC_sig))
        mean_HC20_vec <- mean(output_df$HC20_vec[[i]], na.rm = T)
        sd_HC20_vec <- sd(HC20_vec, na.rm = T)
        output_df$CV[i] <- (sd_HC20_vec/mean_HC20_vec)*100
        output_df$Phi[i] <- sqrt(log((output_df$CV[i]^2)+1))
        output_df$MC_logHC20ec10eq[i] = mean_HC20_vec + (-0.842 * sd_HC20_vec) 
        output_df$Q2.5[i] <- quantile(HC20_vec, 0.025, na.rm = T)
        output_df$Q97.5[i] <- quantile(HC20_vec, 0.975, na.rm = T)
        
         output_df$MC_CRF[i] <- list(0.2/10^output_df$HC20_vec[[i]])
         output_df$MC_CRFQ2.5 <- quantile(output_df$MC_CRF[[i]], 0.025, na.rm = T)
         output_df$MC_CRFQ97.5 <- quantile(output_df$MC_CRF[[i]], 0.975, na.rm = T)

        # Presenting the uncertainties across the whole dataset
        # Comparing variation across chemicals using the geoSD and first back-log-transform all the log-data
        # Back log-transform data to be represented in the geo.sd.dev
        Geometric_St.Dev <- geoSD(10^HC20_vec)
        
          # # Calculating the confidence intervals for mu and sigma at HC20 working point
          # confint_res <- confint(output_df$nls_results[[i]], parm = c("mu", "sig"), level = 0.95)
          # output_df$Q2.5[i] <- qnorm(0.2, mean = coef(nls_out)[1], sd = coef(nls_out)[2])
          # output_df$Q97.5[i] <- qnorm(0.2, mean = coef(nls_out)[1], sd = coef(nls_out)[2], lower.tail = FALSE)
          # # Shoul I use the following input instead of a point parameter?
          # # predict(nls_out)
      
            d.frame <- cbind(d.frame, predict2_nls(nls_out, interval = "confidence"))
          
            ## ggplot output ssd curves
            
              ggplot(d.frame, aes(x = sp_mean, y = y_rank)) +
                geom_point(aes(fill = Taxonomy.Group), shape = 21, size = 1.5, alpha = 0.7) +
                # conditionally plotting precentiles if there are values available
                geom_point(aes(x = output_df$Q2.5[i],
                                            y = (20/100), color = "low_CI", shape = "low_CI"),
                                        inherit.aes = F) +
                geom_point(aes(x = output_df$Q97.5[i],
                                            y = (20/100), color = "high_CI", shape = "high_CI"),
                                        inherit.aes = F) +
                geom_smooth(data = prediction_data, 
                            aes(x = sp_mean, y = predicted_values), 
                            color = "red", 
                            method = "auto", 
                            se = FALSE) +
              
                # plotting the nls function using the same arguments as the nls model above.
                geom_smooth(method = "nls",
                            formula = y ~ cum_norm_dist_function(x, mu, sig),
                            method.args = list(start = cnormSS(sp_mean ~ sd_Li, data = d.frame)),
                            color = "blue",
                            se =  FALSE) + # this is important

                #geom_ribbon(aes(xmin=sp_mean - Q2.5, xmax= sp_mean + Q97.5), alpha=0.2) +
                geom_point(aes(x = output_df$log_HC20EC10eq[i], y = (20/100), color = "HC20EC10eq", shape = "HC20EC10eq"), inherit.aes = F) +
                #geom_point(aes(x = output_df$MC_logHC20ec10eq[i], y = (20/100), color = "MC_logHC20ec10eq", shape = "MC_logHC20ec10eq"), inherit.aes = F) +
                #geom_point(aes(x = output_df$non_W_HCxEC10[i], y = (20/100), color = "non_W_HCxEC10", shape = "non_W_HCxEC10"), inherit.aes = F) +
                scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0", "20", "40", "60", "80", "100")) +
                scale_color_manual(name = element_blank(),
                                   aesthetics = c("color"),
                                   labels = c("logHC20EC10eq", "97.5 percentile", "2.5 percentile", "MC_logHC20ec10eq", "non_W_HCxEC10"),
                                   values = c("HC20EC10eq"="red", "low_CI"="green", "high_CI" = "blue", "MC_logHC20ec10eq" = "purple", "non_W_HCxEC10" = "black")
                ) +
                scale_shape_manual(name = element_blank(),
                                   labels = c("logHC20EC10eq", "97.5 percentile", "2.5 percentile", "MC_logHC20ec10eq", "non_W_HCxEC10"),
                                   values = c(15, 17, 19, 8, 14)
                ) +
                labs(title = paste("CAS", output_df$CAS.Number[i], sep = " "),
                     subtitle = paste("log10HC", 20, " = ", round(output_df$log_HC20EC10eq[i], digits = 4), " ", catch_warning[i], "| BLUE = Unweighted, RED = Weighted", sep = ""),
                     x = "mean(log) EC10eq (mg L-1)",
                     y = "Response level (%)") +
                theme_linedraw() +
                theme(
                  axis.text.y = element_text(),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  legend.background = element_blank(),
                  legend.position = c(0.1, 0.75),
                  legend.box.background = element_rect(fill = "white"),
                  legend.box.margin = margin(-3,-3,-3,-3),
                  legend.spacing.y = unit(1, "mm"),
                  legend.key.size = unit(4, 'mm'),
                  legend.title = element_text(size=8),
                  legend.text = element_text(size=6)
                )
              #ggsave(filename = paste("SSD", output_df$CAS.Number[i],".png", sep = "_"), plot = last_plot(), device = "png", path = paste(Plot_destination, "/", sep = ""))
              
          
      }
    }
  
# sum <- summary(nls_out)
#summary(nls_out)[7]$convInfo$finIter # Number of iterations to convergence 
#summary(nls_out)[7]$convInfo$finTol # Achieved convergence tolerance 
# sum[7]$convInfo$finIter # get data on Number of iterations to convergence 
# sum[7]$convInfo$finTol # get data on Achieved convergence tolerance