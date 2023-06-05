
  
confint(row_1$nls_results, parm = c("mu", "sig"), level = 0.95)
class(row_1$nls_results)
predict(as.lm.nls(row_1$nls_results), interval = "confidence")
coef(row_1$nls_results)[1]
coef(row_1$nls_results)[2]
  
predict(row_1$nls_results)
predict(test_nls_results)
#the single nls output 
?predict
class(test_nls_results)
?confint


nls_results <- nls(
  y_rank ~ cum_norm_dist_function(sp_mean, mu, sig),
  data = EC10eq_df,
  control = nlc,
  start = list(mu = mean(EC10eq_df$sp_mean, na.rm = T), sig = 0.01))

summary(nls_results)
predict(nls_results)

# skapa nls_results list-object med samma namn som 

# Single substance test code

# scripting plot output for each 
# Need to get $nls_result[i] (model results) applied to the 
# my thought is that the apply function looks like:

##plot_generator <- lapply(nls_output, plot_nls_function())

plot_SSDs <- function(x){
  
  if (class(x[i]) == "nls") {
    # paste the as.lm.nls function here
    # predicting curve-fit and generating a df with fit, lower & upper predictions 
    fit <- predict(as.lm.nls(x), interval = "confidence") 
    # adding a y-axis ranking to the `fit` matrix
    fit <- cbind(fit, vector(mode = "numeric", length = length(fit[,1])))
    fit[,4] <- (rank(fit[,1])-0.5)/length(fit[,1])
    # gives me the confidence intervals
    confint_res <- confint(profile(x), parm = c("mu", "sig"), level = 0.95)
    # quantile function = uncertainty of the HC20 value
    CI_HC20s <- as.vector(rbind( 
      qnorm(0.20, mean = confint_res[1,1], sd = confint_res[2,1]),
      qnorm(0.20, mean = confint_res[1,2], sd = confint_res[2,1]),
      qnorm(0.20, mean = confint_res[1,1], sd = confint_res[2,2]),
      qnorm(0.20, mean = confint_res[1,2], sd = confint_res[2,2]) 
    ))
    
  }else {
    warning(print(paste("No nls available")))
  }
  # I can add the plot-operation here
  png(filename = paste("figures/SSD_plots/SSD_curve_CAS", output_df$CAS.Number[i],".png"), width = 480, height = 480)
  plot(x = d.frame$sp_mean, y = d.frame$y_rank, col = "grey", log = "x", main = paste("SSD for CAS", output_df$CAS.Number[i], sep = " "),
       xlab="mean EC10eq", ylab="Ranking per species (0-1)")
  points(x = output_df$log_HC20EC10eq[i], y = 0.2, col = "black", pch = 4)
  points(x = CI_HC20s[1:4], y = c(rep(0.2, 4)), col= "green", pch = 2)
  lines(x = fit[,1], y = predict(x), col = "red")
  text(x = mean(d.frame$sp_mean), y = 0.9, paste("CRF = HC",{{HCx}}, ", value = ", output_df$CRF[i], sep = ""))
  dev.off()
}


# Selecting only first ro to inspect if function is working
row_1 <- lapply(nls_output, '[[', 3)
row_1 <- sapply(nls_output, '[', 3)
plot_nls_function(row_1)


for (i in nls_output) {
  row_i <- lapply(nls_output, '[[', i)
  plot_nls_function(row_i)
}



#################################################################################################################
f = 1
y_axis <- vector(mode = "numeric", length = length(fit1[,1]))
for (f in fit1) {
  while (f <= length(fit1)) {
    y_axis <- (rank(fit1[,f])-0.5)/length(y_axis)
    f = f+1
  }
}

#################################################################################################################
#SLASK
?AIC
?predict
#################################################################################################################
# Need make a summary output for the nls modelling operation
nls_results_df %>% 
  filter(warnings == "OK") %>% 
  summarise(min_sp = min(n_sp),
            min_tax = min(n_tax.grp))

#################################################################################################################
# Setup for selecting substances across a vector using indexes
CAS_lookup <- read.csv("results/FINAL_HESTIA_BASE_EnviroTox_FILL.csv") %>% 
  distinct(CAS.Number) %>% pull(CAS.Number)

CAS_list <- read.csv("results/FINAL_HESTIA_BASE_EnviroTox_FILL.csv") %>% 
  count(CAS.Number, Taxonomy.Group) %>%
  group_by(CAS.Number) %>%
  summarise(n_samples = sum(n), n_tax = n())

#################################################################################################################
# SINGLE substance operation (for debugging)
nls_results <- list(a = NULL, name = NULL, c = NULL)
# Perform first averaging of species-specific tox data to get species mean and standard error
Test_single <- read.csv("results/FINAL_HESTIA_BASE_EnviroTox_FILL.csv") %>% 
  filter(CAS.Number == "63-25-2")
EC10eq_df <- Test_single %>% 
  mutate(Li = log10(EC10eq)) %>% 
  group_by(CAS.Number, Taxonomy.Group, Species) %>% 
  # Log10-transforming effect concentrations, then calculating the arithmetic mean per species' log10 effect concentration
  summarise(
    n_samples = n(),
    sp_mean = mean(Li, na.rm = T),
    sd_Li = sd(Li, na.rm = T)
  ) %>% 
  ungroup() %>% 
  # Denoting the order of appearance in cumulative form:
  mutate(y_rank = (rank(sp_mean, na.last = NA) - 0.5)/length(unique(Species))) %>% 
  arrange(y_rank)

### Functions to apply ###
# Cumulative normal distribution function:
cum_norm_dist_function <- function(x, mu, sig){
  require(pracma)
  0.5 + (0.5*erf((x - mu)/(sig*sqrt(2))))
}
full_df_mu <- mean(EC10eq_df$sp_mean, na.rm = TRUE)
full_df_sig <- mean(EC10eq_df$sd_Li, na.rm = TRUE)
# The nonlinear least squares (nls()) function assuming a cumulative normal distribution
nls_ctrl <- nls.control(maxiter = 1000)
nls_results <-  nls(y_rank ~ cum_norm_dist_function(sp_mean, mu, sig), 
                    data = EC10eq_df, 
                    control = nls_ctrl,
                    start = list(mu = full_df_mu, sig = full_df_sig)
)

fit <- predict(as.lm.nls(nls_results), interval = "confidence")
# adding a y-axis ranking to the `fit` matrix
fit <- cbind(fit, vector(mode = "numeric", length = length(fit[,1])))
fit[,2]
fit[,3]
#fit[,4] <- (rank(fit[,1])-0.5)/length(fit[,1]) # <- kasta skiten!! SKRÄP!

#gives me the confidence intervals
# I can add profile(output... ) here to supress warnnings and the "Waiting for profiling to be done..." note in output
confint_res <- confint(nls_results, parm = c("mu", "sig"), level = 0.95)
# make sure to warn if there are bad models for quantiles

  # quantile function = uncertainty of the HC20 value
  CI_HC20s <- as.vector(rbind(
    qnorm(0.20, mean = confint_res[1,1], sd = confint_res[2,1]),
    qnorm(0.20, mean = confint_res[1,2], sd = confint_res[2,1]),
    qnorm(0.20, mean = confint_res[1,1], sd = confint_res[2,2]),
    qnorm(0.20, mean = confint_res[1,2], sd = confint_res[2,2])
  ))
  
class(nls_results)
plot(EC10eq_df$sp_mean, EC10eq_df$y_rank)
lines(EC10eq_df$sp_mean,  )

##################################################
#testing out a fit_nls_model
nls_lm_out_df <- data.frame("sp_mean" = EC10eq_df$sp_mean, "nls_predict" = predict(nls_results))
fit_nls_lm_out_df <- lm(nls_lm_out_df)
AIC(fit_nls_lm_out_df)
plot(nls_lm_out_df$sp_mean, nls_lm_out_df$nls_predict)
lines(nls_lm_out_df$sp_mean, nls_lm_out_df$nls_predict)
rsquared(fit_nls_lm_out_df)

# Get the R2-values for the model compared to measured results using this function:
summary(lm(obs ~ mod, data=df))$r.squared 
#I can use a summary for linear models
rsqr <- summary(lm(EC10eq_df$sp_mean ~ predict(nls_results)))$r.squared 





# Function to extract metadata from nls() to generate 95% CI smooth line, 
# found at [$http://www.leg.ufpr.br/~walmes/cursoR/ciaeear/as.lm.R$]
as.lm.nls <- function(object, ...) {
  require(proto)
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

# Automating extraction of the mu & sig nls results
coef(test_nls_results)
test_list$b$m$getAllPars()[1] # [1] # <- possible way to extract single value
predict(test_nls_results)

# These parameters are only to enable plotting the confidence intervals in GGplot!!
fit1 <- predict(as.lm.nls(nls_results), interval = "confidence") # using the function as.lm.nls from the internet
EC10eq_df$model_lower_fit <- fit1[,2]
EC10eq_df$model_upper_fit <- fit1[,3]
# We can extract the dispersion of the HC20 values
point_HC20 <- qnorm(0.20, mean = test_list$nls$testa$m$getAllPars()[1], sd = coef(test_list$m)) # gives me a point value of the corresponding HC20 (data fetched from console output)

# gives me the confidence intervals
confint_res <- confint(test_nls_results, parm = c("mu", "sig"), level = 0.95)

# 95% Confidence intervals for the 
CI_HC20s <- as.vector(rbind( 
  qnorm(0.20, mean = confint_res[1,1], sd = confint_res[2,1]),
  qnorm(0.20, mean = confint_res[1,2], sd = confint_res[2,1]),
  qnorm(0.20, mean = confint_res[1,1], sd = confint_res[2,2]),
  qnorm(0.20, mean = confint_res[1,2], sd = confint_res[2,2])
))

# plot the predictions
plot(x = EC10eq_df$sp_mean, y = EC10eq_df$y_rank, col = "grey")
points(x = CI_HC20s[2:5], y = c(rep(0.2, 4)), col= "orange")
points(x = qnorm(0.20, mean = coef(test_nls_results)[1], sd = coef(test_nls_results)[2]), y = 0.2, col = "black", pch = 4)
lines(x = EC10eq_df$sp_mean, y = predict(test_nls_results), col = "red")



###
# GGPLOT code to plot curves for each substance!
###
ggplot(EC10eq_df, aes(x = sp_mean, y = y_rank)) +
  geom_point(aes(color = Taxonomy.Group), shape=1, size = 1.5, alpha = 0.7) +
  geom_point(aes(x = point_HC20, y = 0.2), shape=15, color = "red", show.legend = FALSE, inherit.aes = FALSE) +
  geom_point(aes(x = CI_HC20s[1], y = 0.2), shape=18, color = "orange", show.legend = FALSE, inherit.aes = FALSE) +
  geom_point(aes(x = CI_HC20s[2], y = 0.2), shape=8, color = "orange", show.legend = FALSE, inherit.aes = FALSE) +
  geom_point(aes(x = CI_HC20s[3], y = 0.2), shape=8, color = "orange", show.legend = FALSE, inherit.aes = FALSE) +
  geom_point(aes(x = CI_HC20s[4], y = 0.2), shape=18, color = "orange", show.legend = FALSE, inherit.aes = FALSE) +
  geom_line(aes(y = 0.2), color = "red", alpha = 0.5) +
  geom_line(stat = "smooth", method = "nls",
            # x is mapped to EC10eq
            # y is mapped to y_rank
            formula = y ~ cum_norm_dist_function(x, mu, sig),
            se =  FALSE, # this is important
            method.args = list(start = c(mu = full_df_mu, sig = full_df_sig)),
            linewidth = 0.5) +
  geom_ribbon(aes(x=sp_mean, ymin=model_lower_fit, ymax=model_upper_fit), alpha=0.5 ) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  xlab("mean(log) toxicity - ec10eq") +
  theme(
    legend.position = "none",
    
    #axis.text.y = element_blank()
  )




#Examples
plot(wtloss$Days, wtloss$Weight)

expn1 <- deriv(y ~ b0 + b1 * 2^(-x/th), c("b0", "b1", "th"),
               function(b0, b1, th, x) {})

wtloss.gr <- nls(Weight ~ expn1(b0, b1, th, Days),
                 data = wtloss, start = c(b0=90, b1=95, th=120))
expn2 <- deriv(~b0 + b1*((w0 - b0)/b1)^(x/d0),
               c("b0","b1","d0"), function(b0, b1, d0, x, w0) {})

wtloss.init <- function(obj, w0) {
  p <- coef(obj)
  d0 <-  - log((w0 - p["b0"])/p["b1"])/log(2) * p["th"]
  c(p[c("b0", "b1")], d0 = as.vector(d0))
}

out <- NULL
w0s <- c(110, 100, 90)
for(w0 in w0s) {
  fm <- nls(Weight ~ expn2(b0, b1, d0, Days, w0),
            wtloss, start = wtloss.init(wtloss.gr, w0))
  out <- rbind(out, c(coef(fm)["d0"], confint(fm, "d0")))
}
dimnames(out) <- list(paste(w0s, "kg:"),  c("d0", "low", "high"))
out


# getInitial -Example
PurTrt <- Puromycin[ Puromycin$state == "treated", ]
print(getInitial( rate ~ SSmicmen( conc, Vm, K ), PurTrt ), digits = 3)

# testing out the SSlogis
Chick.1 <- ChickWeight[ChickWeight$Chick == 1, ]
SSlogis(Chick.1$Time, 368, 14, 6)  # response only
local({
  Asym <- 368; xmid <- 14; scal <- 6
  SSlogis(Chick.1$Time, Asym, xmid, scal) # response _and_ gradient
})
getInitial(weight ~ SSlogis(Time, Asym, xmid, scal), data = Chick.1)
## Initial values are in fact the converged one here, "Number of iter...: 0" :
fm1 <- nls(weight ~ SSlogis(Time, Asym, xmid, scal), data = Chick.1)
summary(fm1)
## but are slightly improved here:
fm2 <- update(fm1, control=nls.control(tol = 1e-9, warnOnly=TRUE), trace = TRUE)
all.equal(coef(fm1), coef(fm2)) # "Mean relative difference: 9.6e-6"
str(fm2$convInfo) # 3 iterations

erf(1)

nls_results <-  nls(y_rank ~ cum_norm_dist_function(sp_mean, mu, sig), 
                    data = EC10eq_df, 
                    start = list(mu = full_df_mu, sig = full_df_sig)

getInitial( y_rank ~ SSlogis(sp_mean, Asym, mu, sig), EC10eq_df)

nls_results_sslogis <-  nls(y_rank ~ SSlogis(sp_mean, Asym, mu, sig), 
                        data = EC10eq_df)
summary(nls_results_sslogis)                    

nls_results_sslogis_2 <- update(nls_results_sslogis, control=nls.control(tol = 1e-9, warnOnly=TRUE), trace = TRUE)
all.equal(coef(nls_results_sslogis), coef(nls_results_sslogis_2)) 
str(nls_results_sslogis_2$convInfo)



