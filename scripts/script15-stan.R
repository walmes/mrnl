# Mais código em: https://gitlab.c3sl.ufpr.br/walmes/stan

#-----------------------------------------------------------------------
# Code copied from Mages' Blog.
# http://www.magesblog.com/2015/10/non-linear-growth-curves-with-stan.html

# ATTENTION! TODO created in 08/01/2016, check this link.
# http://datascienceplus.com/bayesian-regression-with-stan-part-1-normal-regression/

#-----------------------------------------------------------------------
# Packages.

library(lattice)
library(latticeExtra)

library(rstan)      # Stan: http://mc-stan.org/.
library(tidyverse)  # Tidyverse: https://www.tidyverse.org/.

# rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())

# Please, visit the GIFs in the page below.
browseURL("http://leg.ufpr.br/~walmes/ensino/ce089-2014-02/ce089-2014-02-aula06.html")

#-----------------------------------------------------------------------
# Data.

dat <- list(
    "N" = 27,
    "n_new" = 50,
    "x" = c(1, 1.5, 1.5, 1.5, 2.5, 4, 5, 5, 7, 8, 8.5, 9, 9.5, 9.5, 10,
            12, 12, 13, 13, 14.5, 15.5, 15.5, 16.5, 17, 22.5, 29, 31.5),
    "Y" = c(1.8, 1.85, 1.87, 1.77, 2.02, 2.27, 2.15, 2.26, 2.47, 2.19,
            2.26, 2.4, 2.39, 2.41, 2.5, 2.32, 2.32, 2.43, 2.47, 2.56,
            2.65, 2.47, 2.64, 2.56, 2.7, 2.72, 2.57))
dat$x_new <- with(dat, seq(min(x), max(x), length.out = n_new))

xyplot(Y ~ x,
       data = dat,
       type = c("p", "smooth"))

#-----------------------------------------------------------------------
# Fitting a nonlinear regression model.

nlm <- nls(Y ~ alpha - beta * lambda^x,
           data = dat[c("x", "Y")],
           start = list(alpha = 1, beta = 1, lambda = 0.9))

summary(nlm)

# Estimate and CI.
cbind(summary(nlm)$coefficients[, 1:2],
      confint.default(nlm))

# Fitted curve over the points.
xyplot(Y ~ x,
       data = dat) +
    latticeExtra::layer({
        panel.curve(alpha - beta * lambda^x)
    },
    data = as.list(coef(nlm)))

#-----------------------------------------------------------------------
# Defining a Stan model.

stanmodel <- "
data {
  int<lower=0> N;     // Sample size.
  int<lower=0> n_new; // Prediction vector size.
  real x[N];          // Vector of independent variable.
  real Y[N];          // Vector of dependent variable.
  real x_new[n_new];  // Vector to predict the response.
}
parameters {
  real alpha;
  real beta;
  real<lower=.5,upper=1> lambda;
  real<lower=0> tau;
}
transformed parameters {
  real sigma;
  sigma = 1/sqrt(tau);
}
model {
  real m[N];
  for (i in 1:N) {
    m[i] = alpha - beta * pow(lambda, x[i]);
  }
  // Likelihood.
  Y ~ normal(m, sigma);
  // Priors.
  alpha ~ normal(0.0, 10); // Variance not so big.
  beta ~ normal(0.0, 10);
  lambda ~ uniform(.5, 1);
  tau ~ gamma(.0001, .0001);
}
generated quantities {
  real Y_mean[n_new]; // Mean.
  real Y_pred[n_new]; // Future observations.
  for (i in 1:n_new) {
    // Posterior parameter distribution of the mean.
    Y_mean[i] = alpha - beta * pow(lambda, x_new[i]);
    // Posterior predictive distribution.
    Y_pred[i] = normal_rng(Y_mean[i], sigma);
  }
}
"

# Fitting the bayesian model.
system.time({
    fit <- stan(model_code = stanmodel,
                model_name = "GrowthCurve",
                data = dat)
})
#   user  system elapsed
# 39.652   0.848  40.551
# 38.868   0.984  39.947
# 35.696   0.872  36.894

#-----------------------------------------------------------------------
# Exploring the results.

# Select the parameters to show the results.
p <- head(names(fit@sim$samples[[1]]), n = 4)
p

#  Overall results.
print(fit, pars = p)

# Traces for each parameter.
traceplot(fit, pars = p)

# Scatterplot matrix.
pairs(fit, pars = p)

#-----------------------------------------------------------------------
# Regression parameters.

# Function to resume the posterior.
stat_fun <- function(x) {
    stat <- c(mean(x),
              quantile(x, probs = c(0.025, 0.975)))
    hpd <- coda::HPDinterval(coda::as.mcmc(x),
                             prob = 0.95)
    stat <- append(stat, as.vector(hpd))
    names(stat) <- c("mean",
                     "q_lwr", "q_upr",
                     "hpd_lwr", "hpd_upr")
    return(stat)
}

# Extract the samples for each model parameter.
L <- sapply(p, FUN = function(x) rstan::extract(fit, x)[[1]])

# Sample summaries.
s <- apply(L, MARGIN = 2, FUN = stat_fun) %>% as.data.frame()
s

# Stack the values.
L <- L %>% as.data.frame() %>% stack()
str(L)

densityplot(~values | ind,
            data = L,
            as.table = TRUE,
            scales = "free") +
    latticeExtra::layer({
        panel.abline(v = s[, packet.number()],
                     lty = c(1, 2, 2, 3, 3))
    })

#-----------------------------------------------------------------------
# Regression confidence bands.

# List with values.
L <- list(mean = rstan::extract(fit, "Y_mean")[[1]],
          pred = rstan::extract(fit, "Y_pred")[[1]])
str(L)

#  Get the summaries.
L <- lapply(L,
       FUN = function(x) {
           apply(x, MARGIN = 2, FUN = stat_fun) %>%
               t() %>%
               as.data.frame() %>%
               mutate(x = dat$x_new)
       })
str(L)

# Legend.
leg <- list(corner = c(0.95, 0.05),
            lines = list(lty = c(0, 1, 2, 3),
                         pch = c(1, NA, NA, NA)),
            text = list(c("Observed",
                          "Posterior mean",
                          "Quantile CI (95%)",
                          "HDP CI (95%)")),
            divide = 1,
            type = "o")

# For the mean.
xyplot(Y ~ x,
       data = dat,
       key = leg) +
    latticeExtra::as.layer({
        xyplot(mean + q_lwr + q_upr + hpd_lwr + hpd_upr ~ x,
               data = L$mean,
               col = c(1, 2, 2, 3, 3),
               lty = c(1, 2, 2, 3, 3),
               type = "l")
    })

# For future observations.
xyplot(Y ~ x,
       data = dat,
       key = leg) +
    latticeExtra::as.layer({
        xyplot(mean + q_lwr + q_upr + hpd_lwr + hpd_upr ~ x,
               data = L$pred,
               col = c(1, 2, 2, 3, 3),
               lty = c(1, 2, 2, 3, 3),
               type = "l")
    })

#-----------------------------------------------------------------------
# Confidence bands using delta method (asymptotic approximation).

L$dm <- data.frame(x = dat$x_new)

formula(nlm)

f <- function(theta, x_new){
    with(as.list(theta), alpha - beta * lambda^x_new)
}

# Gradient of f for each x_new value at MLE.
X <- rootSolve::gradient(f, x = coef(nlm), x_new = L$dm$x)

# Variance-covariance matrix.
U <- chol(vcov(nlm))

# Fitted values.
L$dm$fit <- predict(nlm, newdata = L$dm)

# Standard error for the fitted value.
L$dm$se <- sqrt(apply(X %*% t(U), 1, function(x) sum(x^2)))

#  Get the IC (95%).
tval <- qt(p = c(lwr = 0.025, upr = 0.975),
           df = df.residual(nlm))
L$dm$me <- outer(L$dm$se, tval, "*")
L$dm <- cbind(L$dm, sweep(L$dm$me, 1, L$dm$fit, "+"))
str(L$dm)

xyplot(Y ~ x,
       data = dat) +
    latticeExtra::as.layer({
        xyplot(mean + q_lwr + q_upr + hpd_lwr + hpd_upr ~ x,
               data = L$mean,
               col = c(1, 2, 2, 3, 3),
               lty = c(1, 2, 2, 3, 3),
               type = "l")
    }) +
    latticeExtra::as.layer({
        xyplot(fit + lwr + upr ~ x,
               data = L$dm,
               col = c(4, 4, 4),
               lty = c(2, 4, 4),
               type = "l")
    })

#-----------------------------------------------------------------------
# Inspecting the confidence regions of the parameters.

library(nlstools)
ls("package:nlstools")

# `nlstools` methods needs a `data.frame` not a `list`.
da <- as.data.frame(dat[c("x", "Y")])

nlm <- nls(Y ~ alpha - beta * lambda^x,
           data = da,
           start = list(alpha = 1, beta = 1, lambda = 0.9))

# Residual sum of squares contours.
crss <- nlsContourRSS(nlm)
plot(crss)

# Confidence regions by Monte Carlo simulation.
creg <- nlsConfRegions(nlm)
plot(creg)

#-----------------------------------------------------------------------
