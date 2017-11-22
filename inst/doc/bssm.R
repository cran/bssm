## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----UKgas---------------------------------------------------------------
library("bssm")
set.seed(123)

init_sd <- 0.1 * sd(log10(UKgas))
prior <- halfnormal(init_sd, 1)
model <- bsm(log10(UKgas), sd_y = prior, sd_level = prior,
  sd_slope = prior, sd_seasonal =  prior)

## ----mcmc----------------------------------------------------------------
mcmc_out <- run_mcmc(model, n_iter = 1e5)
mcmc_out

## ----plot----------------------------------------------------------------
theta <- expand_sample(mcmc_out, "theta") ## until bayesplot is updated
library("bayesplot")
mcmc_areas(theta, bw = 0.001)
level <- expand_sample(mcmc_out, "alpha",  times = 101:108, states = 1)
mcmc_areas(level)

# posterior mode estimates
mcmc_out$theta[which.max(mcmc_out$posterior), ]


## ----trend, dev.args=list(pointsize = 10), fig.cap="Smoothed trend component with 95% intervals."----
level <- expand_sample(mcmc_out, "alpha", states = 1)
# or using summary method:
sumr <- summary(mcmc_out)
level <- sumr$states$Mean[, 1]
lwr <- level - 1.96 * sumr$states$SD[, 1]
upr <- level + 1.96 * sumr$states$SD[, 1]
ts.plot(model$y, cbind(level, lwr, upr), col = c(1, 2, 2, 2), lty = c(1, 1, 2, 2))

## ----predict, dev.args=list(pointsize = 10), fig.cap="Mean predictions and prediction intervals."----
future_model <- model
future_model$y <- ts(rep(NA, 24), start = end(model$y) + c(0, 1), frequency = 4)
pred <- predict(mcmc_out, future_model, probs = c(0.025, 0.1, 0.9, 0.975))
ts.plot(log10(UKgas), pred$mean, pred$intervals[, -3],
  col = c(1, 2, c(3, 4, 4, 3)), lty = c(1, 1, rep(2, 4)))

## ----predict2, dev.args=list(pointsize = 10), fig.cap="Prediction plots with ggplot2."----
require("ggplot2")
level_fit <- ts(colMeans(expand_sample(mcmc_out, "alpha")$level), start = start(model$y), 
  frequency = 4)
autoplot(pred, y = model$y, fit = level_fit, interval_color = "red", alpha_fill = 0.2)

## ----predict3, dev.args=list(pointsize = 10), fig.cap="State prediction."----
pred_state <- predict(mcmc_out, future_model, probs = c(0.025, 0.1, 0.9, 0.975), type = "state")
ts.plot(log10(UKgas), level_fit, pred_state$mean[,"level"], pred_state$intervals$level[, -3],
  col = c(1, 2, 2, c(3, 4, 4, 3)), lty = c(1, 1, 1, rep(2, 4)))

