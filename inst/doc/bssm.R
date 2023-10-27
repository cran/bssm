## ----echo = FALSE-------------------------------------------------------------
Sys.setenv("OMP_NUM_THREADS" = 2) # For CRAN
if (!requireNamespace("rmarkdown") ||
    !rmarkdown::pandoc_available("1.12.3")) {
  warning(call. = FALSE, "These vignettes assume rmarkdown and pandoc version",
    "1.12.3. These were not found. Older versions will not work.")
  knitr::knit_exit()
}

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----nhtemp-------------------------------------------------------------------
library("bssm")
data("nhtemp", package = "datasets")
prior <- halfnormal(1, 10)
bsm_model <- bsm_lg(y = nhtemp, sd_y = prior, sd_level = prior,
  sd_slope = prior)

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
level_fit <- d |> 
  filter(variable == "level") |>
  group_by(time) |>
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

