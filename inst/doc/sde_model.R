## ---- echo = FALSE------------------------------------------------------------
if (!requireNamespace("rmarkdown") ||
    !rmarkdown::pandoc_available("1.12.3")) {
  warning(call. = FALSE, "These vignettes assume rmarkdown and pandoc version 1.12.3. These were not found. Older versions will not work.")
  knitr::knit_exit()
}

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
#options(tinytex.verbose = TRUE)

## -----------------------------------------------------------------------------
set.seed(1)
n <- 40
suppressMessages(library("sde"))
x <- sde.sim(t0 = 0, T = n, X0 = 1, N = n * 2^5,
  drift = expression(0.5 * (2 - x)),
  sigma = expression(0.2),
  sigma.x = expression(0),
  method = "milstein")
integer_x <- x[seq(frequency(x) + 1, length(x), frequency(x))]
y <- rpois(n, exp(integer_x))

## -----------------------------------------------------------------------------
library("bssm")
Rcpp::sourceCpp("ssm_sde_template.cpp")
pntrs <- create_xptrs()
sde_model <- ssm_sde(y, pntrs$drift, pntrs$diffusion, 
  pntrs$ddiffusion, pntrs$obs_density, pntrs$prior, 
  theta = c(log_rho = log(0.5), mu = 2, log_sigma = log(0.2)),
  x0 = 1, positive = FALSE)

## -----------------------------------------------------------------------------
out <- run_mcmc(sde_model, iter = 2e4, particles = 30, L_c = 2, L_f = 5, threads = 2)

## -----------------------------------------------------------------------------
suppressMessages(library("ggplot2"))
suppressMessages(library("dplyr"))
suppressMessages(library("diagis"))

d <- as.data.frame(out, variable = "states")

state_fit <- d %>% 
  group_by(time) %>%
  summarise(state = weighted_mean(value, weight), 
    lwr = weighted_quantile(value, weight, 0.025), 
    upr = weighted_quantile(value, weight, 0.975))

ggplot(state_fit, aes(x = time, y = state)) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr), 
    fill = "#7570b3", alpha = 0.25) +
  geom_line(data = data.frame(
    state = x, 
    time = time(x)), 
    colour = "#d95f02", linetype = "dashed") +
  geom_line(colour = "#7570b3") +
  geom_point(colour = "#7570b3") +
  geom_point(data=data.frame(state=integer_x,time=1:n), colour = "#d95f02") +
  theme_bw()

