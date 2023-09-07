library(tidyverse)
library(pracma)
library(nlraa)
library(EnvStats)

HESTIA_nls_function <- function(dataset, HCx = 20, MC_n = 100000, rm_singles = "YES") {
  
  options(dplyr.summarise.inform = FALSE)
  
  cum_norm_dist_function <- function(x, mu, sig) {
    0.5 + 0.5 * erf((x - mu) / (sig * sqrt(2)))
  }
  
  cnormSS <- function(formula, data, start = NULL, weights = NULL) {
    if (is.null(start)) {
      x <- data[[all.vars(formula)[1]]]
      y <- data[[all.vars(formula)[2]]]
      mu <- mean(x, na.rm = TRUE)
      sig <- mean(y, na.rm = TRUE)
      start <- list(mu = mu, sig = sig)
    }
    return(start)
  }
  
  CAS_list <- dataset %>%
    distinct(CAS.Number) %>%
    arrange(CAS.Number) %>%
    pull(CAS.Number)
  
  error_df <- dataset %>%
    count(CAS.Number, Taxonomy.Group, Species) %>%
    group_by(CAS.Number) %>%
    summarise(n_sp = n_distinct(Species), n_tax.grp = n_distinct(Taxonomy.Group), n_recs = sum(n)) %>%
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
    status = rep("", length(CAS_list))
  )
  
  for (i in seq_along(CAS_list)) {
    if (output_df$n_sp[i] == 5 || output_df$n_tax.grp[i] == 3) {
      output_df$status[i] <- "Not enough data"
      next
    }
    
    d.frame <- dataset %>%
      filter(CAS.Number == CAS_list[i]) %>%
      mutate(Li = log10(EC10eq)) %>%
      group_by(Taxonomy.Group, Species) %>%
      summarise(
        n_samples = n(),
        sp_mean = mean(Li, na.rm = TRUE),
        sd_Li = sd(Li, na.rm = TRUE)
      ) %>%
      ungroup() %>%
      mutate(
        y_rank = (rank(sp_mean, na.last = NA) - 0.5) / length(unique(Species)),
        sd_Li = ifelse(sd_Li == 0 | is.na(sd_Li), 1, sd_Li)
      ) %>%
      arrange(y_rank) %>%
      mutate(Taxonomy.Group = as.factor(Taxonomy.Group))
    
    if (rm_singles == "YES") {
      d.frame <- d.frame %>%
        filter(sd_Li != 1)
    }
    
    if (nrow(d.frame) == 1) {
      output_df$status[i] <- "Sigma available for only one species"
      next
    }
    
    nlc <- nls.control(maxiter = 1000, warnOnly = TRUE)
    nls_out <- tryCatch(
      nls(
        y_rank ~ cum_norm_dist_function(sp_mean, mu, sig),
        control = nlc,
        weights = 1 / (sd_Li^2),
        data = d.frame,
        start = cnormSS(sp_mean ~ sd_Li, data = d.frame)
      ),
      error = function(e) {
        output_df$status[i] <- "Failure to initialize"
        return(NULL)
      },
      warning = function(w){
        output_df$status[i] <- "Failure to converge"
        return(NULL)
      }
    )
    
    if (!is.null(nls_out)) {
      output_df$status[i] <- "Convergence"
      output_df$log_HC20EC10eq[i] <- qnorm((HCx / 100), mean = coef(nls_out)[1], sd = coef(nls_out)[2])
      output_df$CRF[i] <- (HCx / 100) / (10^output_df$log_HC20EC10eq[i])
      output_df$mu[i] <- coef(nls_out)[1]
      output_df$sigma[i] <- coef(nls_out)[2]
      output_df$CV_HCx[i] <- (output_df$sigma[i] / output_df$mu[i]) * 100
      output_df$Phi_HCx[i] <- sqrt(log((output_df$CV_HCx[i]^2) + 1))
      output_df$Iterations.to.Convergence[i] <- summary(nls_out)[7]$convInfo$finIter
      output_df$Achieved.convergence.tolerance[i] <- summary(nls_out)[7]$convInfo$finTol
      
      output_df$MC_mu[i] <- rnorm(MC_n, mean = summary(nls_out)$parameters[1, 1], sd = summary(nls_out)$parameters[1, 2])
      output_df$MC_sig[i] <- rnorm(MC_n, mean = summary(nls_out)$parameters[2, 1], sd = summary(nls_out)$parameters[2, 2])
      output_df$log_HCx_vec[i] <- qnorm((HCx / 100), mean = output_df$MC_mu[[i]], sd = output_df$MC_sig[[i]])
      output_df$mean_HCx[i] <- mean(output_df$log_HCx_vec[[i]], na.rm = TRUE)
      output_df$sd_HCx[i] <- sd(output_df$log_HCx_vec[[i]], na.rm = TRUE)
      output_df$GStDev[i] <- geoSD(10^output_df$log_HCx_vec[[i]], na.rm = TRUE)
      output_df$Q2.5[i] <- quantile(output_df$log_HCx_vec[[i]], 0.025, na.rm = TRUE)
      output_df$Q97.5[i] <- quantile(output_df$log_HCx_vec[[i]], 0.975, na.rm = TRUE)
      output_df$MC_CRFQ2.5[i] <- quantile(0.2 / 10^output_df$log_HCx_vec[[i]], 0.025, na.rm = TRUE)
      output_df$MC_CRFQ97.5[i] <- quantile(0.2 / 10^output_df$log_HCx_vec[[i]], 0.975, na.rm = TRUE)
    }
  }
  
  return(output_df)
}
