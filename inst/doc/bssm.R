## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----nhtemp-------------------------------------------------------------------
library("bssm")
data("nhtemp", package = "datasets")
prior <- halfnormal(1, 10)
bsm_model <- bsm_lg(y = nhtemp, sd_y = prior, sd_level = prior,
  sd_slope = prior)

## -----------------------------------------------------------------------------
set.seed(1)
suppressMessages(library("sde"))
x <- sde.sim(t0 = 0, T = 100, X0 = 1, N = 100,
  drift = expression(0.5 * (2 - x)),
  sigma = expression(1),
  sigma.x = expression(0))
y <- rpois(100, exp(x[-1]))

## ---- eval = FALSE------------------------------------------------------------
#  Rcpp::sourceCpp("ssm_sde_template.cpp")
#  pntrs <- create_xptrs()
#  sde_model <- ssm_sde(y, pntrs$drift, pntrs$diffusion,
#    pntrs$ddiffusion, pntrs$obs_density, pntrs$prior, c(0.5, 2, 1), 1, FALSE)

## ----mcmc_bsm-----------------------------------------------------------------
prior <- halfnormal(0.1, 1)
UKgas_model <- bsm_lg(log10(UKgas), sd_y = prior, sd_level = prior,
  sd_slope = prior, sd_seasonal =  prior)

mcmc_bsm <- run_mcmc(UKgas_model, iter = 4e4, seed = 1)
mcmc_bsm

## ----plot---------------------------------------------------------------------
suppressMessages(library("ggplot2"))
d <- as.data.frame(mcmc_bsm, variable = "theta")
ggplot(d, aes(x = value)) + 
  geom_density(adjust = 3, fill = "#92f0a8") + 
  facet_wrap(~ variable, scales = "free") + 
  theme_bw()

## ----trend, dev.args=list(pointsize = 10), fig.cap="Smoothed trend component with 95% intervals."----
suppressMessages(library("dplyr"))
d <- as.data.frame(mcmc_bsm, variable = "states")
level_fit <- d %>% 
  filter(variable == "level") %>%
  group_by(time) %>%
  summarise(consumption = mean(value), 
    lwr = quantile(value, 0.025), 
    upr = quantile(value, 0.975))

ggplot(level_fit, aes(x = time, y = consumption)) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), 
    fill = "#92f0a8", alpha = 0.25) +
  geom_line(colour = "#92f0a8") +
  geom_line(data = data.frame(
    consumption = log10(UKgas), 
    time = time(UKgas)), 
    colour = "grey30", linetype = "dashed") +
  theme_bw()

