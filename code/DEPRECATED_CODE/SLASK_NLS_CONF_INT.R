# To improve the runtime of the nls_across_all function, here are some suggestions:
#   
#  - Instead of loading the entire tidyverse package, load only the necessary packages, such as dplyr and ggplot2.
#  - Remove the require(pracma) statement inside the cum_norm_dist_function function since the package is already required outside of the function.
#  - Use data.table package instead of dplyr for faster data manipulation.
#  - Preallocate output_df as a data.frame instead of a list.
#  - Use lapply or purrr::map instead of a while loop for looping through each substance.
#  - Replace the if...else statement inside the loop with a vectorized version using if_else from dplyr.
#  - Simplify the as.lm.nls function by removing the warning since the warning already appears when using the function.
#  - Remove the options(dplyr.summarise.inform = FALSE) statement since it only affects the printing of dplyr summaries.

?qnorm

test_df <- read.csv("results/FINAL_HESTIA_BASE_EnviroTox_FILL.csv") %>% filter(CAS.Number == "50-28-2")

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

mu_guess <- mean(d.frame$sp_mean, na.rm = TRUE)
sig_guess <- mean(d.frame$sd_Li, na.rm = TRUE)

cum_norm_dist_function <- function(x, mu, sig){
  0.5 + (0.5*erf((x - mu)/(sig*sqrt(2))))
}
#### NLS MODEL ####
nlc <- nls.control(maxiter = 1000, warnOnly = TRUE)
nls_out <- nls(
  y_rank ~ cum_norm_dist_function(sp_mean, mu, sig),
  data = d.frame,
  control = nlc,
  start = list(mu = mu_guess, sig = sig_guess))
#### NLS MODEL ####
fitted(nls_out)
# as.lm.nls function from "the internet".
fit_as.lm <- predict(as.lm.nls(nls_out), interval = "confidence")
fit <- predict(nls_out, interval = "confidence")
summary_nls <- summary(nls_out, interval = "confidence")
nls_out_dm <- predict2_nls(nls_out, interval = "conf")
plot(d.frame$sp_mean, d.frame$y_rank)
line()
# adding a y-axis ranking to the `fit` matrix
fit <- cbind(fit, vector(mode = "numeric", length = length(fit[1])))
d.frame$model_lower_fit <- fit[,2]
d.frame$model_upper_fit <- fit[,3]
qnorm(nls_out)
# I can add profile(output... ) here to supress warnnings and the "Waiting for profiling to be done..." note in output
confint_res <- confint(nls_out, parm = c("mu", "sig"), level = 0.95)
# quantile function = uncertainty of the HC20 value
CI_HC20s <- as.vector(rbind(
  qnorm(0.2, mean = confint_res[1,1], sd = confint_res[2,1]), # lower mean, lower sig quantile of the normal distribution at 20%
  #qnorm(0.2, mean = confint_res[1,2], sd = confint_res[2,1]), # upper mean, lower sig
  #qnorm(0.2, mean = confint_res[1,1], sd = confint_res[2,2]), # lower mean, upper sig
  qnorm(0.2, mean = confint_res[1,2], sd = confint_res[2,2])  # upper mean, upper sig quantile of the normal distribution at 20%
))

### Functions to apply ###
# Cumulative normal distribution function:
cum_norm_dist_function <- function(x, mu, sig){
  0.5 + (0.5*erf((x - mu)/(sig*sqrt(2))))
}


cnormSS <- function(formula, data, start=NULL, weights=NULL){
  if(is.null(start)){
    x <- data[[all.vars(formula)[2]]]
    mu <- mean(x)
    sig <- sd(x)
    start <- list(mu=mu, sig=sig)
  }
  return(start)
}

fit <- nls(y_rank ~ cum_norm_dist_function(sp_mean, mu, sig), data=d.frame, start=cnormSS(y_rank ~ sp_mean, data = d.frame), trace = TRUE)
summary(fit)
fit_pred <- predict(fit, interval = "confidence")
qnorm(0.2, coef(fit)[1], coef(fit)[2])
fit_out_dm <- predict2_nls(fit, interval = "conf")

qnorm()
#####################
cnormSS <- function(formula, data, start=NULL, weights=NULL){
  if(is.null(start)){
    x <- data[[all.vars(formula)[2]]]
    mu <- mean(x)
    sig <- sd(x)
    start <- list(mu=mu, sig=sig)
  }
  return(start)
}
set.seed(123)
x <- seq(-5, 5, length.out=100)
y <- cum_norm_dist_function(x, mu=1, sig=2) + rnorm(length(x), 0, 0.1)

data <- data.frame(x=x, y=y)
plot(data$x, data$y)
# Fit the nls model with starting values and the port algorithm
fit <- nls(y ~ cum_norm_dist_function(x, mu, sig), data=data, start=list(mu=0, sig=1), algorithm="port", trace=TRUE)
summary(fit)
######################
?selfStart
## The "initializer" (finds initial values for parameters from data):
initLogis <- function(mCall, data, LHS) {
  xy <- data.frame(sortedXyData(mCall[["input"]], LHS, data))
  if(nrow(xy) < 4)
    stop("too few distinct input values to fit a logistic model")
  z <- xy[["y"]]
  ## transform to proportion, i.e. in (0,1) :
  rng <- range(z); dz <- diff(rng)
  z <- (z - rng[1L] + 0.05 * dz)/(1.1 * dz)
  xy[["z"]] <- log(z/(1 - z))		# logit transformation
  aux <- coef(lm(x ~ z, xy))
  pars <- coef(nls(y ~ 1/(1 + exp((xmid - x)/scal)),
                   data = xy,
                   start = list(xmid = aux[[1L]], scal = aux[[2L]]),
                   algorithm = "plinear"))
  setNames(pars[c(".lin", "xmid", "scal")], nm = mCall[c("sp_mean", "mu", "sig")])
}

pchisq(2*x^2,1)*sign(x)

SScum_norm_dist <- selfStart( ~ 0.5 + (0.5*((x - mu)/(sig*sqrt(2)))), initial = initLogis, parameters = c("x", "mu", "sig"))

nls_out_ss <- nls(
  y_rank ~ SScum_norm_dist(sp_mean, mu, sig),
  data = d.frame,
  control = nlc)

confint_res <- confint(nls_out, parm = c("mu", "sig"), level = 0.95)
fit <- data.frame(predict(nls_out, interval = "confidence"))
fit$low <- qnorm(predict(nls_out, interval = "confidence"), mean = confint_res[1,1], sd = confint_res[2,1])  # upper mean, upper sig quantile of the normal distribution at 20%
fit$high <- qnorm(predict(nls_out, interval = "confidence"), mean = confint_res[1,2], sd = confint_res[2,2])  # upper mean, upper sig quantile of the normal distribution at 20%

## Confidence bands using the bootstrap method
nls_out_bt <- boot_nls(nls_out)
pairs(fm1.P.bt$t, labels = c("mu", "sig"))
## Bootstrapped confidence bands
nls_out_bt_ft <- boot_nls(nls_out, fitted) ## This takes about 5s

nls_out_bt_ft_prd <- summary_simulate(t(nls_out_bt_ft$t))
d.frame_boot <- cbind(d.frame, nls_out_bt_ft_prd)
## Plot data and bands
ggplot(data = d.frame, aes(x = sp_mean, y = y_rank)) + 
  geom_point() + 
  geom_line(aes(y = fitted(fit))) + 
  geom_ribbon(aes(ymin = fit$low, ymax = fit$high), fill = "purple", alpha = 0.2) + 
  geom_point(aes(x = qnorm(0.2, mean = confint_res[1,1], sd = confint_res[2,1]), 
                 y = 0.2, 
                 color = "low_CI", 
                 shape = "low_CI"), 
             inherit.aes = F) +
  #geom_point(aes(x = CI_HC20s[2], y = ({{HCx}}/100), shape=8, color = "blue", show.legend = F, inherit.aes = F)) +
  #geom_point(aes(x = CI_HC20s[3], y = ({{HCx}}/100), shape=8, color = "blue", show.legend = F, inherit.aes = F)) +
  geom_point(aes(x = qnorm(0.2, 
                           mean = confint_res[1,2], 
                           sd = confint_res[2,2]), 
                 y = 0.2, 
                 color = "high_CI", 
                 shape = "high_CI"), 
             inherit.aes = F) +
  ggtitle("95% Bootstrap Confidence Bands")

## Predictions at 0.4
prd_fun <- function(x) predict(x, newdata = data.frame(conc = 0.2))
prd_fun(nls_out) ## It also provides the gradient

nls_out_x_0.2 <- boot_nls(nls_out, prd_fun) ## This takes about 6s

boot::boot.ci(nls_out_x_0.2, type = "perc") 

# Monte Carlo method
ndat <- data.frame(sp_mean = seq(min(d.frame$sp_mean), max(d.frame$sp_mean), length.out = 50))
Pprd <- predict_nls(nls_out, interval = "conf",
                    newdata = ndat)
Pprdd <- data.frame(ndat, Pprd)
ggplot() + 
  geom_point(data = d.frame, aes(x = sp_mean, y = y_rank)) + 
  geom_line(data = Pprdd, aes(x = sp_mean, y = Estimate)) + 
  geom_ribbon(data = Pprdd, aes(x = sp_mean, ymin = Q2.5, ymax = Q97.5), 
              fill = "purple", alpha = 0.4) + 
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
  ggtitle("Monte Carlo 95% Confidence Bands")

fit <- data.frame(fitten = predict(nls_out, interval = "confidence"))
fit$low <- qnorm(fit$fitten, mean = confint_res[1,1], sd = confint_res[2,1])
fit$high <- qnorm(fit$fitten, mean = confint_res[1,2], sd = confint_res[2,2])
lm(cum_norm_dist_function(16, 2, 0.5))
# Checking out the as.ml.nls function
##########################################
gradient_test<- nls_out$m$gradient()
colnames(gradient) <- names(nls_out$m$getPars())
as.character(formula(nls_out)[[2]])
response.name <- as.character(formula(nls_out)[[2]])
lhs <- nls_out$m$lhs()
L <- data.frame(lhs, gradient_test)
names(L)[1] <- response.name
fo <- sprintf("%s ~ %s - 1", response.name, 
              paste(colnames(gradient), collapse = "+"))
fo <- as.formula(fo, env = proto:::as.proto.list(L))
do.call("lm", list(fo, offset = substitute(fitted(nls_out))))

###########################################
# The as.lm.nls function that enables plotting linear models with confidence intervals in ggplot
as.lm.nls <- function(object, ...) {
  if (!inherits(object, "nls")) {
    w <- paste("expected object of class nls but got object of class:", 
               paste(class(object), collapse = " "))
    warning(w)
  }
  gradient <- object$m$gradient()
  if (is.null(colnames(gradient))) {
    colnames(gradient) <- names(object$m$getPars())
  }
  response.name <- if (length(formula(object)) == 2) "0" else 
    as.character(formula(object)[[2]])
  lhs <- object$m$lhs()
  L <- data.frame(lhs, gradient)
  names(L)[1] <- response.name
  fo <- sprintf("%s ~ %s - 1", response.name, 
                paste(colnames(gradient), collapse = "+"))
  fo <- as.formula(fo, env = proto:::as.proto.list(L))
  do.call("lm", list(fo, offset = substitute(fitted(object))))
}

# So this function is fitting a linear model to the predicted nls() output. 