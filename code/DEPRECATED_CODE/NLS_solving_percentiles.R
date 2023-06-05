#install.packages("pracma")
#install.packages("proto")
#install.packages("nlraa")
library(tidyverse)
library(pracma)
library(rlang)
library(nlraa) # <- contains the "predict2_nls" function

test_df <- read.csv("results/FINAL_HESTIA_BASE_EnviroTox_FILL.csv") %>% filter(CAS.Number == "100-61-8")
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
      if(is.na(mu)) mu <- 0
      if(is.na(sig)) sig <- 1
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
    summarise(n_sp = length(unique(Species)), n_tax.grp = length(unique(Taxonomy.Group)), n_recs = sum(n)) %>%
    arrange(CAS.Number)
  
  ### create a template list-object to fill ###
  output_df <- list(CAS.Number = NULL, log_HC20EC10eq = NULL,
                    CRF  = NULL, mu  = NULL, sigma  = NULL,
                    Q2.5 = NULL, Q97.5 = NULL,
                    n_sp = NULL, n_tax.grp = NULL,
                    n_recs = NULL,
                    status = NULL, nls_results  = NULL)
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
        mutate(y_rank = (rank(sp_mean, na.last = NA) - 0.5)/length(unique(Species)),
               # Some occurrences where sd == 0 cause error in the nls model. making these into NA's
               sd_Li = case_when(sd_Li == 0 ~ as.numeric(NA), TRUE ~ sd_Li)) %>%
        arrange(y_rank) %>%
        mutate(Taxonomy.Group = as.factor(Taxonomy.Group))
      
     
        # The nonlinear least squares (nls()) function assuming a cumulative normal distribution
        nlc <- nls.control(maxiter = 1000, warnOnly = TRUE)
        nls_out <- nls(y_rank ~ cum_norm_dist_function(sp_mean, mu, sig), 
                       control = nlc, 
                       data=d.frame, 
                       start = cnormSS(sp_mean ~ sd_Li, 
                                       data = d.frame))
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
        
          # Calculating the confidence intervals for mu and sigma at HC20 working point
          confint_res <- confint(output_df$nls_results[[i]], parm = c("mu", "sig"), level = 0.95)
          output_df$Q2.5[i] <- qnorm((20/100), mean = confint_res[1,1], sd = confint_res[2,1])
          output_df$Q97.5[i] <- qnorm((20/100), mean = confint_res[1,2], sd = confint_res[2,2])
          # Shoul I use the following input instead of a point parameter?
          # predict(nls_out)
      
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
                
                #plotting the nls function usng the same arguments as the nls model above.
                geom_smooth(method = "nls",
                            formula = y ~ cum_norm_dist_function(x, mu, sig),
                            method.args = list(start = cnormSS(sp_mean ~ sd_Li, data = d.frame)),
                            se =  FALSE) + # this is important
                #geom_ribbon(aes(xmin=sp_mean - Q2.5, xmax= sp_mean + Q97.5), alpha=0.2) +
                geom_point(aes(x = output_df$log_HC20EC10eq[i], y = (20/100), color = "HC20EC10eq", shape = "HC20EC10eq"), inherit.aes = F) +
                scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0", "20", "40", "60", "80", "100")) +
                scale_color_manual(name = element_blank(),
                                   aesthetics = c("color"),
                                   labels = c("logHC20EC10eq", "2.5 percentile", "97.5 percentile"),
                                   values = c("HC20EC10eq"="red", "low_CI"="green", "high_CI" = "blue")
                ) +
                scale_shape_manual(name = element_blank(),
                                   labels = c("logHC20EC10eq", "2.5 percentile", "97.5 percentile"),
                                   values = c(15, 8, 8)
                ) +
                labs(title = paste("CAS", output_df$CAS.Number[i], sep = " "),
                     subtitle = paste("log10HC", 20, "EC10eq"," = ", round(output_df$log_HC20EC10eq[i], digits = 4), " ", catch_warning[i], sep = ""),
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
# sum[7]$convInfo$finIter # get data on Number of iterations to convergence 
# sum[7]$convInfo$finTol # get data on Achieved convergence tolerance
