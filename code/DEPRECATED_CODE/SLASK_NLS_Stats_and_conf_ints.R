### testing of the least_squares_model script
system.time(nls_output <- nls_across_all(dataset = test_df, Plot_output = "YES", Plot_destination = "figures/SSD_plots/test"))
test_df <- read.csv("results/FINAL_HESTIA_BASE_EnviroTox_FILL.csv") %>% filter(CAS.Number == "104-40-5")
#test_df <- read.csv("results/FINAL_HESTIA_BASE_EnviroTox_FILL.csv") %>% filter(CAS.Number == "500008-45-7")
nls_output_df <- as.data.frame(do.call(cbind, nls_output[1:13])) %>%
  mutate(across(c(2:10), ~ as.numeric(.x)),
         CAS.Number = as.factor(CAS.Number))



test_df <- read.csv("results/FINAL_HESTIA_BASE_EnviroTox_FILL.csv") %>% filter(CAS.Number == "104-40-5")
d.frame <- test_df %>% 
  # filter(CAS.Number == CAS_list[i]) %>% # in testing redundant, since we have only one substance
  # Perform first averaging of species-specific tox data to get species mean and standard error
  mutate(Li = log10(EC10eq)) %>% 
  group_by(Taxonomy.Group, Species) %>% 
  # Log10-transforming effect concentrations, then calculating the arithmetic mean per species' log10 effect concentration
  summarise(
    n_samples = n(),
    sp_mean = mean(Li),
    sd_Li = sd(Li, na.rm = TRUE),
  ) %>% 
  ungroup() %>% 
  # Denoting the order of appearance in cumulative form:
  mutate(y_rank = (rank(sp_mean, na.last = NA) - 0.5)/length(unique(Species)),
         # Some occurrences where sd == 0 cause error in the nls model. making these into NA's
         sd_Li = case_when(sd_Li == 0 ~ as.numeric(NA), TRUE ~ sd_Li)) %>% 
  arrange(y_rank) %>% 
  mutate(Taxonomy.Group = as.factor(Taxonomy.Group))

### Functions to apply ###
# Cumulative normal distribution function:
cum_norm_dist_function <- function(x, mu, sig){
  0.5 + (0.5*erf((x - mu)/(sig*sqrt(2))))
}

## Selfstart function
cnormSS <- function(formula, data, start=NULL, weights=NULL){
  if(is.null(start)){
    x <- data[[all.vars(formula)[1]]]
    y <- data[[all.vars(formula)[2]]]
    mu <- mean(x, na.rm = TRUE)
    sig <- mean(y, na.rm = TRUE)
    if(is.na(mu)) mu <- 1
    if(is.na(sig)) sig <- 1
    start <- list(mu=mu, sig=sig)
  }
  return(start)
}

# cnormSS <- function(formula, data, start=NULL, weights=NULL){
#   if(is.null(start)){
#     x <- data[[all.vars(formula)[1]]]
#     y <- data[[all.vars(formula)[2]]]
#     mu <- mean(x, na.rm = TRUE)
#     sig <- mean(y, na.rm = TRUE)
#     if(all(x < 0)) mu <- abs(mu)
#     if(!is.finite(mu)) mu <- 0
#     if(!is.finite(sig)) sig <- 1
#     start <- list(mu=mu, sig=sig)
#   }
#   return(start)
# }
# 
# cnormSS <- function(formula, data, start=NULL, weights=NULL) {
#   if(is.null(start)) {
#     x <- data[[all.vars(formula)[1]]]
#     y <- data[[all.vars(formula)[2]]]
#     
#     # Compute initial values for mu and sigma
#     mu0 <- median(x)
#     q25 <- quantile(x, 0.25)
#     q75 <- quantile(x, 0.75)
#     sig0 <- (q75 - q25) / 1.34  # Estimate sigma using IQR / 1.34
#     
#     # Return starting values as a list
#     start <- list(mu = mu0, sig = sig0)
#   }
#   return(start)
# }
# 
# normSS_custom <- function(x, y, ...) {
#   if (is.null(start)) {
#     mu <- mean(x)
#     sig <- sd(x)
#     if (sig == 0) sig <- mean(x)
#     start <- list(mu = mu, sig = sig)
#   }
#   else start <- as.list(start)
#   start$mu <- max(-10, min(10, start$mu))
#   start$sig <- max(0.01, min(10, start$sig))
#   nls(y ~ pnorm(x, start$mu, start$sig, lower.tail = TRUE), 
#       start = start, ...)
#}

cnormSS_custom <- function(formula, data,start=NULL) {
  if(is.null(start)) {
    x <- data[[all.vars(formula)[1]]]
    y <- data[[all.vars(formula)[2]]]
  mu <- sum(x*y)/sum(y)
  sig <- sqrt(sum(y*(x-mu)^2)/sum(y))
  start <- list(mu = mu, sig = sig)
  }
  return(start)
}
# nls model
nls_out <- nls(y_rank ~ cum_norm_dist_function(sp_mean, mu, sig), data=d.frame, start = cnormSS(sp_mean ~ sd_Li, data = d.frame), trace = TRUE)

confint(nls_out, parm = c("mu", "sig"), level = 0.95)
c_int <- confint(nls_out, method = "profile")
# Extracting information form the nls
# coef() extracts the estimated parameter values from the nls-model. 
coef(nls_out)
# Each of the below operations extract estimated mu and sigma form the model fit.
est_mu <- coef(nls_out)[1]
est_sig <- coef(nls_out)[2]
# i use these estimated parameters from the fitted model to calculate the value of the 0.2/20% point on the model curve as such
logHC20EC10eq <- qnorm(0.2, est_mu, est_sig)
# But this only gives me a point value for the logHC20EC10eq, I want to find the 95% probability where this value could be found based on the input data.
# here i calculate the 2.5 and 97.5 quantile (95% confidence interval) based on the nls model output, which delimits the range of values that is likely to contain the true value of a population parameter with 95% probability.
confint_res <- confint(profile(nls_out), parm = c("mu", "sig"), level = 0.95) 

# Obtaining estimates for the parameters mu & sig of a normal distribution. These are to be used when predicting the respective low and high value of my independent variable
# Below are the calculated lower and upper range for the 95% confidence interval of the logHC20EC10eq value (at the 0.2/20% working level). 
Q2.5 <- qnorm(0.2, mean = confint_res[1,1], sd = confint_res[2,1]) # <- the 2.5% quantile 
Q97.5 <- qnorm(0.2, mean = confint_res[1,2], sd = confint_res[2,2]) # <- the 97.5% quantile

# Now i have the estimated 95% confidence interval for the 20% working point generated from the nls model fit.
# Then, when it comes to plotting the area of 95% confidence interval 
# here i extract the Estimated points for the curve, standard error and predictions for 2.5 and 97.5 quartiles for each estimated point of the nls model. 
d.frame <- cbind(d.frame, predict2_nls(nls_out, interval = "confidence", level = 0.95))

sum_nls_out <- summary(nls_out)
sum_nls_out$coefficients
sum_nls_out$parameters[1:2, 2]
sum(resid(nls_out)^2)
AIC(nls_out)
BIC(nls_out)


d.frame <- cbind(d.frame, data.frame(pred_fit = predict(nls_out, interval = "confidence")))
d.frame$low <- qnorm(d.frame$y_rank, mean = confint_res[1,1], sd = confint_res[2,1])  # upper mean, upper sig quantile of the normal distribution at 20%
d.frame$high <- qnorm(d.frame$y_rank, mean = confint_res[1,2], sd = confint_res[2,2])  # upper mean, upper sig quantile of the normal distribution at 20%
predict(as.lm.nls(nls_out), interval = "confidence")
predict2_nls(nls_out, interval = "confidence", level = 0.95)


testsiffra <- 0.2
## Plot data and bands
ggplot(data = d.frame, aes(x = sp_mean, y = y_rank)) + 
  geom_point() + 
  geom_line(aes(y = fitted(nls_out))) + 
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = "purple", alpha = 0.2) + 
  geom_point(aes(x = qnorm(0.2, mean = confint_res[1,1], sd = confint_res[2,1]), 
                 y = 0.2, 
                 color = "low_CI", 
                 shape = "low_CI"), 
             inherit.aes = F) +
  geom_point(aes(x = qnorm(0.2, 
                           mean = confint_res[1,2], 
                           sd = confint_res[2,2]), 
                 y = 0.2, 
                 color = "high_CI", 
                 shape = "high_CI"), 
             inherit.aes = F) +
  ggtitle(paste("Test of ChatGPTs", testsiffra, "self-start function", sep = " "))




########################################################################################################

# # Setup for selecting substances across a vector using indexes
# CAS_lookup <- read.csv("results/FINAL_HESTIA_BASE_EnviroTox_FILL.csv") %>%
#   distinct(CAS.Number) %>% pull(CAS.Number)
# # # Subsetting the toxicity dataset
# Tox_data_one_CAS_test <- read.csv("results/FINAL_HESTIA_BASE_EnviroTox_FILL.csv") %>%
#   #filter(CAS.Number == "114-07-8")
#   filter(CAS.Number %in% CAS_lookup[1:25])

### Happened to plot the taxonomic group responses!
# # Cool way to show the different responses. possibility to estimate HC20 response per taxonomic group!
# ggplot(Tox_data_one_CAS %>% 
#          group_by(Taxonomy.Group, Species) %>%        
#          summarise(sp_mean = mean(EC10eq)) %>% 
#          arrange(sp_mean) %>% 
#          mutate(y_rank = (rank(sp_mean, na.last = NA) - 0.5)/length(unique(Species))),
#        aes(x = sp_mean, y = y_rank)) +
#   geom_point(aes(color = Taxonomy.Group), size = 1.5, alpha = 0.7)+
#   scale_x_log10()
# 
# 
# # plot of actual data distribution # Works poorly with all the data points right now.
# # Skip this part for now!!!
# plot_df <- Tox_data_one_CAS %>% 
#   left_join(
#     x = ., 
#     y = Tox_data_one_CAS %>% 
#       group_by(Species) %>% 
#       summarise(sp_mean = mean(EC10eq)) %>% 
#       arrange(sp_mean) %>% 
#       mutate(y_rank = (rank(sp_mean, na.last = NA) - 0.5)/length(unique(Species))),
#   by = "Species") %>% 
#   left_join(x = ., 
#             y = EC10eq_df %>% 
#                 select(Species, model_lower_fit, model_upper_fit),
#             by = "Species")
# 
# ggplot(plot_df, aes(x = EC10eq, y = y_rank)) +
#   geom_point(aes(color = Taxonomy.Group), size = 1.5, alpha = 0.7)+
#   scale_x_log10()+
#   geom_line(stat = "smooth", method = "nls", 
#             # x is mapped to EC10eq
#             # y is mapped to y_rank
#             formula = y ~ cum_norm_dist_function(x, mu, sig), 
#             se =  FALSE, # this is important 
#             method.args = list(start = c(mu = full_df_mu, sig = full_df_sig)),
#             size = 0.5) +
#   geom_ribbon(aes(x=EC10eq, ymin=model_lower_fit, ymax=model_upper_fit), alpha=0.5 )
