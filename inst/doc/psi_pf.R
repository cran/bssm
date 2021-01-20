## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
options(dplyr.summarise.inform = FALSE)

## ---- echo = FALSE, message=FALSE, warning=FALSE------------------------------
library(dplyr)

## ---- echo = FALSE, cache = FALSE, eval = FALSE-------------------------------
#  library("bssm")
#  library("foreach")
#  library("doParallel")
#  
#  growth_model_experiment <- function(n_cores, nsim, particles) {
#  
#    set.seed(1)
#  
#    p1 <- 50 # population size at t = 1
#    K <- 500 # carrying capacity
#  
#    #sample time
#    dT <- .1
#  
#    #observation times
#    t <- seq(0.1, 30, dT)
#    n <- length(t)
#    r <- plogis(cumsum(c(-1.5, rnorm(n - 1, sd = 0.05))))
#    p <- numeric(n)
#    p[1] <- p1
#    for(i in 2:n)
#      p[i] <- rnorm(1, K * p[i-1] * exp(r[i-1] * dT) / (K + p[i-1] * (exp(r[i-1] * dT) - 1)), 1)
#    # observations
#    y <- p + rnorm(n, 0, 1)
#  
#    initial_theta <- c(H = 1, R1 = 0.05, R2 = 1)
#  
#    # dT, K, a1 and the prior variances
#    known_params <- c(dT = dT, K = K, a11 = -1.5, a12 = 50, P11 = 1, P12 = 100)
#  
#    cl<-makeCluster(n_cores)
#    registerDoParallel(cl)
#  
#    results <- foreach (j = 1:n_cores, .combine = "rbind", .packages = "bssm") %dopar% {
#  
#      Rcpp::sourceCpp("growth_model_functions.cpp")
#      pntrs <- create_xptrs()
#      model <- ssm_nlg(y = y, a1=pntrs$a1, P1 = pntrs$P1,
#        Z = pntrs$Z_fn, H = pntrs$H_fn, T = pntrs$T_fn, R = pntrs$R_fn,
#        Z_gn = pntrs$Z_gn, T_gn = pntrs$T_gn,
#        theta = initial_theta, log_prior_pdf = pntrs$log_prior_pdf,
#        known_params = known_params, known_tv_params = matrix(1),
#        n_states = 2, n_etas = 2)
#  
#      bsf <- ekpf <- psi <- matrix(NA, 10, nsim / n_cores)
#  
#      for(i in 1:ncol(bsf)) {
#  
#        time <- system.time(out <- particle_smoother(model, particles = particles, method = "bsf"))[3]
#        bsf[, i] <- c(out$logLik, out$alphahat[1, ], diag(out$Vt[, , 1]),
#          out$alphahat[n, ], diag(out$Vt[, , n]), time)
#  
#        time <- system.time(out <- particle_smoother(model, particles = particles, method = "psi"))[3]
#        psi[, i] <- c(out$logLik, out$alphahat[1, ], diag(out$Vt[, , 1]),
#          out$alphahat[n, ], diag(out$Vt[, , n]), time)
#  
#        time <- system.time(out <- particle_smoother(model, particles = particles, method = "ekf"))[3]
#        ekpf[, i] <- c(out$logLik, out$alphahat[1, ], diag(out$Vt[, , 1]),
#          out$alphahat[n, ], diag(out$Vt[, , n]), time)
#      }
#      x <- t(cbind(bsf, ekpf, psi))
#      colnames(x) <- c("logLik", "alpha_11", "alpha_21", "V_11", "V_21",
#        "alpha_1n", "alpha_2n", "V_1n", "V_2n", "time")
#  
#      data.frame(x,
#        method = rep(factor(c("BSF", "EKPF", "PSI")), each = ncol(bsf)),
#        N = particles)
#    }
#    stopCluster(cl)
#    results
#  }
#  
#  gm_result_10 <- growth_model_experiment(1, 10000, 10)
#  saveRDS(gm_result_10, file = "gm_result_10.rds")
#  
#  gm_result_100 <- growth_model_experiment(1, 10000, 100)
#  saveRDS(gm_result_100, file = "gm_result_100.rds")
#  
#  gm_result_1000 <- growth_model_experiment(1, 10000, 1000)
#  saveRDS(gm_result_1000, file = "gm_result_1000.rds")
#  
#  # ground truth
#  set.seed(1)
#  
#  p1 <- 50 # population size at t = 1
#  K <- 500 # carrying capacity
#  
#  #sample time
#  dT <- .1
#  
#  #observation times
#  t <- seq(0.1, 30, dT)
#  n <- length(t)
#  r <- plogis(cumsum(c(-1.5, rnorm(n - 1, sd = 0.05))))
#  p <- numeric(n)
#  p[1] <- p1
#  for(i in 2:n)
#    p[i] <- rnorm(1, K * p[i-1] * exp(r[i-1] * dT) / (K + p[i-1] * (exp(r[i-1] * dT) - 1)), 1)
#  # observations
#  y <- p + rnorm(n, 0, 1)
#  
#  initial_theta <- c(H = 1, R1 = 0.05, R2 = 1)
#  
#  # dT, K, a1 and the prior variances
#  known_params <- c(dT = dT, K = K, a11 = -1.5, a12 = 50, P11 = 1, P12 = 100)
#  
#  Rcpp::sourceCpp("growth_model_functions.cpp")
#  
#  pntrs <- create_xptrs()
#  model <- ssm_nlg(y = y, a1=pntrs$a1, P1 = pntrs$P1,
#    Z = pntrs$Z_fn, H = pntrs$H_fn, T = pntrs$T_fn, R = pntrs$R_fn,
#    Z_gn = pntrs$Z_gn, T_gn = pntrs$T_gn,
#    theta = initial_theta, log_prior_pdf = pntrs$log_prior_pdf,
#    known_params = known_params, known_tv_params = matrix(1),
#    n_states = 2, n_etas = 2)
#  
#  out <- particle_smoother(model, particles = 1e5, method = "bsf")
#  truth <- c(out$logLik, out$alphahat[1, ], diag(out$Vt[, , 1]),
#    out$alphahat[n, ], diag(out$Vt[, , n]))
#  names(truth) <- c("logLik", "alpha_11", "alpha_21", "V_11", "V_21", "alpha_1n",
#    "alpha_2n", "V_1n", "V_2n")
#  saveRDS(truth, file = "gm_truth.rds")

## ---- echo = FALSE, eval = FALSE----------------------------------------------
#  gm10 <- readRDS("psi_pf_experiments/gm_result_10.rds")
#  gm100 <- readRDS("psi_pf_experiments/gm_result_100.rds")
#  gm1000 <- readRDS("psi_pf_experiments/gm_result_1000.rds")
#  
#  results <- rbind(gm10, gm100, gm1000)

## ----loglik, echo = FALSE, eval = FALSE---------------------------------------
#  reference <- readRDS("psi_pf_experiments/gm_truth.rds")
#  
#  IRE <- function(x, time) {
#      mean((x - truth)^2) * mean(time)
#  }
#  truth <- reference["logLik"]
#  sumr <- results %>% group_by(method, N) %>%
#    summarise(mean = mean(logLik), SD = sd(logLik),
#      IRE = IRE(logLik, time), time = mean(time))
#  table1 <- sumr %>% arrange(N) %>% knitr::kable(digit = 4,
#               caption = "Results for the log-likelihood estimates of the growth model. ")
#  saveRDS(table1, file = "psi_pf_experiments/table1.rds")

## ----tabl21, echo = FALSE-----------------------------------------------------
readRDS("psi_pf_experiments/table1.rds")

## ----alpha, echo = FALSE, eval = FALSE----------------------------------------
#  truth <- reference["alpha_11"]
#  sumr <- results %>% group_by(method, N) %>%
#    summarise(mean = mean(alpha_11), SD = sd(alpha_11),
#    IRE = IRE(alpha_11, time), time = mean(time))
#  
#  table2 <- sumr %>% arrange(N) %>% knitr::kable(digit = 4,
#               caption = "Results for the p_1 estimates of the growth model. ")
#  saveRDS(table2, file = "psi_pf_experiments/table2.rds")

## ----table2, echo = FALSE-----------------------------------------------------
readRDS("psi_pf_experiments/table2.rds")

## ---- echo = FALSE, eval = FALSE----------------------------------------------
#  library("bssm")
#  library("foreach")
#  library("doParallel")
#  
#  ar_exp_model_experiment <- function(n_cores, nsim, particles, theta) {
#  
#    set.seed(1)
#    n <- 100
#    alpha <- arima.sim(n = n, list(ar = 0.95), sd = theta)
#    y <- rnorm(n, exp(alpha))
#  
#    cl<-makeCluster(n_cores)
#    registerDoParallel(cl)
#  
#    results <- foreach (j = 1:n_cores, .combine = "rbind", .packages = "bssm") %dopar% {
#  
#      Rcpp::sourceCpp("ar_exp_model_functions.cpp")
#      pntrs <- create_xptrs()
#      model <- ssm_nlg(y = y, a1=pntrs$a1, P1 = pntrs$P1,
#        Z = pntrs$Z_fn, H = pntrs$H_fn, T = pntrs$T_fn, R = pntrs$R_fn,
#        Z_gn = pntrs$Z_gn, T_gn = pntrs$T_gn,
#        theta = theta, log_prior_pdf = pntrs$log_prior_pdf,
#        known_params = 0, known_tv_params = matrix(1),
#        n_states = 1, n_etas = 1)
#  
#  
#      bsf <- ekpf <- psi <- matrix(NA, 6, nsim / n_cores)
#  
#  
#      for(i in 1:ncol(bsf)) {
#        time <- system.time(out <- particle_smoother(model, particles = particles, method = "bsf"))[3]
#        bsf[, i] <- c(out$logLik, out$alphahat[1, ], out$Vt[, , 1],
#          out$alphahat[n, ], out$Vt[, , n], time)
#  
#        time <- system.time(out <- particle_smoother(model, particles = particles, method = "psi"))[3]
#        psi[, i] <- c(out$logLik, out$alphahat[1, ], out$Vt[, , 1],
#          out$alphahat[n, ], out$Vt[, , n], time)
#  
#        time <- system.time(out <- particle_smoother(model, particles = particles, method = "ekf"))[3]
#        ekpf[, i] <- c(out$logLik, out$alphahat[1, ], out$Vt[, , 1],
#          out$alphahat[n, ], out$Vt[, , n], time)
#      }
#      x <- t(cbind(bsf, ekpf, psi))
#      colnames(x) <- c("logLik", "alpha_1", "V_1", "alpha_n", "V_n", "time")
#  
#      data.frame(x,
#        method = rep(factor(c("BSF", "EKPF", "PSI")), each = ncol(bsf)),
#        N = particles)
#    }
#    stopCluster(cl)
#    results
#  }
#  
#  
#  ### TRUTH
#  
#  set.seed(1)
#  n <- 100
#  alpha <- arima.sim(n = n, list(ar = 0.95), sd = 0.1)
#  y <- rnorm(n, exp(alpha), 1)
#  
#  Rcpp::sourceCpp("ar_exp_model_functions.cpp")
#  pntrs <- create_xptrs()
#  model <- ssm_nlg(y = y, a1=pntrs$a1, P1 = pntrs$P1,
#    Z = pntrs$Z_fn, H = pntrs$H_fn, T = pntrs$T_fn, R = pntrs$R_fn,
#    Z_gn = pntrs$Z_gn, T_gn = pntrs$T_gn,
#    theta = 0.1, log_prior_pdf = pntrs$log_prior_pdf,
#    known_params = 0, known_tv_params = matrix(1),
#    n_states = 1, n_etas = 1)
#  
#  
#  out <- particle_smoother(model, nsim = 1e6, method = "bsf")
#  reference <- c(logLik = out$logLik, alpha_1=out$alphahat[1], V_1 = out$Vt[1],
#    alpha_n = out$alphahat[n], V_n = out$Vt[n])
#  saveRDS(reference, file = "ar_truth.rds")
#  
#  print("Running with 10 particles")
#  ar_result_10 <- ar_exp_model_experiment(1, 10000, 10, 0.1)
#  saveRDS(ar_result_10, file = "ar_result_10.rds")
#  
#  print("Running with 100 particles")
#  ar_result_100 <- ar_exp_model_experiment(1, 10000, 100, 0.1)
#  saveRDS(ar_result_100, file = "ar_result_100.rds")
#  
#  print("Running with 1000 particles")
#  ar_result_1000 <- ar_exp_model_experiment(1, 10000, 1000, 0.1)
#  saveRDS(ar_result_1000, file = "ar_result_1000.rds")
#  

## ---- echo = FALSE, eval = FALSE----------------------------------------------
#  ar10 <- readRDS("psi_pf_experiments/ar_result_10.rds")
#  ar100 <- readRDS("psi_pf_experiments/ar_result_100.rds")
#  ar1000 <- readRDS("psi_pf_experiments/ar_result_1000.rds")
#  
#  results <- rbind(ar10, ar100, ar1000)

## ----loglik_ar, echo = FALSE, eval = FALSE------------------------------------
#  reference <- readRDS("psi_pf_experiments/ar_truth.rds")
#  truth <- reference["logLik"]
#  sumr <- results %>% group_by(method, N) %>%
#    summarise(mean = mean(logLik), SD = sd(logLik),
#      IRE = IRE(logLik, time), time = mean(time))
#  table3 <- sumr %>% arrange(N) %>% knitr::kable(digit = 4,
#               caption = "Results for the log-likelihood estimates of the AR model. ")
#  saveRDS(table3, file = "psi_pf_experiments/table3.rds")

## ----table3, echo = FALSE-----------------------------------------------------
readRDS("psi_pf_experiments/table3.rds")

## ----state1_ar, echo = FALSE, eval = FALSE------------------------------------
#  truth <- reference["alpha_1"]
#  sumr <- results %>% group_by(method, N) %>%
#    summarise(mean = mean(alpha_1), SD = sd(alpha_1),
#      IRE = IRE(alpha_1, time))
#  table4 <- sumr %>% arrange(N) %>% knitr::kable(digit = 4,
#               caption = "Results for the alpha_1 estimates of the AR model. ")
#  saveRDS(table4, file = "psi_pf_experiments/table4.rds")

## ----table5, echo = FALSE-----------------------------------------------------
readRDS("psi_pf_experiments/table4.rds")

