## ----echo = FALSE-------------------------------------------------------------
Sys.setenv("OMP_NUM_THREADS" = 2) # For CRAN
if (!requireNamespace("rmarkdown") ||
    !rmarkdown::pandoc_available("1.12.3")) {
  warning(call. = FALSE, "These vignettes assume rmarkdown and pandoc version 1.12.3. These were not found. Older versions will not work.")
  knitr::knit_exit()
}

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----data---------------------------------------------------------------------
set.seed(1)

p1 <- 50 # population size at t = 1
K <- 500 # carrying capacity
H <- 1 # standard deviation of obs noise
R_1 <- 0.05 # standard deviation of the noise on logit-growth
R_2 <- 1 # standard deviation of the noise in population level
#sample time
dT <- .1

#observation times
t <- seq(0.1, 30, dT)
n <- length(t)
r <- plogis(cumsum(c(-1.5, rnorm(n - 1, sd = R_1))))
p <- numeric(n)
p[1] <- p1
for(i in 2:n)
  p[i] <- rnorm(1, K * p[i-1] * exp(r[i-1] * dT) / (K + p[i-1] * (exp(r[i-1] * dT) - 1)), R_2)
# observations
y <- p + rnorm(n, 0, H)

## ----pointers-----------------------------------------------------------------
Rcpp::sourceCpp("ssm_nlg_template.cpp")
pntrs <- create_xptrs()

## ----theta--------------------------------------------------------------------
initial_theta <- c(log_H = 0, log_R1 = log(0.05), log_R2 = 0)

# dT, K, a1 and the prior variances of 1st and 2nd state (logit r and and p)
known_params <- c(dT = dT, K = K, a11 = -1, a12 = 50, P11 = 1, P12 = 100)

## ----test---------------------------------------------------------------------
T_fn(0, c(100, 200), initial_theta, known_params, matrix(1))

## ----model--------------------------------------------------------------------
library("bssm")
model <- ssm_nlg(y = y, a1=pntrs$a1, P1 = pntrs$P1, 
  Z = pntrs$Z_fn, H = pntrs$H_fn, T = pntrs$T_fn, R = pntrs$R_fn, 
  Z_gn = pntrs$Z_gn, T_gn = pntrs$T_gn,
  theta = initial_theta, log_prior_pdf = pntrs$log_prior_pdf,
  known_params = known_params, known_tv_params = matrix(1),
  n_states = 2, n_etas = 2, state_names = c("logit_r", "p"))

## ----ekf----------------------------------------------------------------------
out_filter <- ekf(model)
out_smoother <- ekf_smoother(model)
ts.plot(cbind(y, out_filter$att[, 2], 
  out_smoother$alphahat[, 2]), col = 1:3)
ts.plot(plogis(cbind(out_filter$att[, 1], 
  out_smoother$alphahat[, 1])), col = 1:2)

## ----mcmc---------------------------------------------------------------------
# Cholesky of initial proposal matrix (taken from previous runs)
# used here to speed up convergence due to the small number of iterations
S <- matrix(c(0.13, 0.13, -0.11, 0, 0.82, -0.04, 0, 0, 0.16), 3, 3) 
mcmc_is <- run_mcmc(model, iter = 6000, burnin = 1000, particles = 10, 
  mcmc_type = "is2", sampling_method = "psi", S = S)
mcmc_ekf <- run_mcmc(model, iter = 6000, burnin = 1000, 
  mcmc_type = "ekf", S = S)
summary(mcmc_is, return_se = TRUE)
summary(mcmc_ekf, return_se = TRUE)

## ----summaries----------------------------------------------------------------
library("dplyr")
library("diagis")
d1 <- as.data.frame(mcmc_is, variable = "states")
d2 <- as.data.frame(mcmc_ekf, variable = "states")
d1$method <- "is2-psi"
d2$method <- "approx ekf"

r_summary <- rbind(d1, d2) |> 
  filter(variable == "logit_r") |>
  group_by(time, method) |>
  summarise(
    mean = weighted_mean(plogis(value), weight), 
    lwr = weighted_quantile(plogis(value), weight, 0.025), 
    upr = weighted_quantile(plogis(value), weight, 0.975))

p_summary <- rbind(d1, d2) |> 
  filter(variable == "p") |>
  group_by(time, method) |>
  summarise(  
    mean = weighted_mean(value, weight), 
    lwr = weighted_quantile(value, weight, 0.025), 
    upr = weighted_quantile(value, weight, 0.975))

## ----figures------------------------------------------------------------------
library("ggplot2")
ggplot(r_summary, aes(x = time, y = mean)) + 
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = method), 
    colour = NA, alpha = 0.25) +
  geom_line(aes(colour = method)) +
  geom_line(data = data.frame(mean = r, time = seq_along(r))) +
  theme_bw()

p_summary$cut <- cut(p_summary$time, c(0, 100, 200, 301))
ggplot(p_summary, aes(x = time, y = mean)) + 
  geom_point(data = data.frame(
    mean = y, time = seq_along(y),
    cut = cut(seq_along(y), c(0, 100, 200, 301))), alpha = 0.1) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = method), 
    colour = NA, alpha = 0.25) +
  geom_line(aes(colour = method)) +
  theme_bw() + facet_wrap(~ cut, scales = "free")

## -----------------------------------------------------------------------------
mcmc_is$time
mcmc_ekf$time

