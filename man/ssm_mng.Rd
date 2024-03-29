% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/models.R
\name{ssm_mng}
\alias{ssm_mng}
\title{General Non-Gaussian State Space Model}
\usage{
ssm_mng(
  y,
  Z,
  T,
  R,
  a1 = NULL,
  P1 = NULL,
  distribution,
  phi = 1,
  u,
  init_theta = numeric(0),
  D = NULL,
  C = NULL,
  state_names,
  update_fn = default_update_fn,
  prior_fn = default_prior_fn
)
}
\arguments{
\item{y}{Observations as multivariate time series or matrix with dimensions
n x p.}

\item{Z}{System matrix Z of the observation equation as p x m matrix or
p x m x n array.}

\item{T}{System matrix T of the state equation. Either a m x m matrix or a
m x m x n array.}

\item{R}{Lower triangular matrix R the state equation. Either a m x k
matrix or a
m x k x n array.}

\item{a1}{Prior mean for the initial state as a vector of length m.}

\item{P1}{Prior covariance matrix for the initial state as m x m matrix.}

\item{distribution}{A vector of distributions of the observed series.
Possible choices are
\code{"poisson"}, \code{"binomial"}, \code{"negative binomial"},
\code{"gamma"}, and \code{"gaussian"}.}

\item{phi}{Additional parameters relating to the non-Gaussian distributions.
For negative binomial distribution this is the dispersion term, for
gamma distribution this is the shape parameter, for Gaussian this is
standard deviation, and for other distributions this is ignored.}

\item{u}{A matrix of positive constants for non-Gaussian models
(of same dimensions as y). For Poisson,  gamma, and negative binomial
distribution, this corresponds to the offset term. For binomial, this is the
number of trials (and as such should be integer(ish)).}

\item{init_theta}{Initial values for the unknown hyperparameters theta
(i.e. unknown variables excluding latent state variables).}

\item{D}{Intercept terms for observation equation, given as p x n matrix.}

\item{C}{Intercept terms for state equation, given as m x n matrix.}

\item{state_names}{A character vector defining the names of the states.}

\item{update_fn}{A function which returns list of updated model
components given input vector theta. This function should take only one
vector argument which is used to create list with elements named as
\code{Z}, \code{T}, \code{R}, \code{a1}, \code{P1}, \code{D}, \code{C}, and
\code{phi}, where each element matches the dimensions of the original model.
If any of these components is missing, it is assumed to be constant wrt.
theta. It's best to check the internal dimensions with
\code{str(model_object)} as the dimensions of input arguments can differ
from the final dimensions.}

\item{prior_fn}{A function which returns log of prior density
given input vector theta.}
}
\value{
An object of class \code{ssm_mng}.
}
\description{
Construct an object of class \code{ssm_mng} by directly defining the
corresponding terms of the model.
}
\details{
The general multivariate non-Gaussian model is defined using the following
observational and state equations:

\deqn{p^i(y^i_t | D_t + Z_t \alpha_t), (\textrm{observation equation})}
\deqn{\alpha_{t+1} = C_t + T_t \alpha_t + R_t \eta_t, 
(\textrm{transition equation})}

where \eqn{\eta_t \sim N(0, I_k)} and
\eqn{\alpha_1 \sim N(a_1, P_1)} independently of each other, and
\eqn{p^i(y_t | .)} is either Poisson, binomial, gamma, Gaussian, or
negative binomial distribution for each observation series \eqn{i=1,...,p}.
Here k is the number of disturbance terms (which can be less than m,
the number of states).
}
\examples{
 
set.seed(1)
n <- 20
x <- cumsum(rnorm(n, sd = 0.5))
phi <- 2
y <- cbind(
  rgamma(n, shape = phi, scale = exp(x) / phi),
  rbinom(n, 10, plogis(x)))

Z <- matrix(1, 2, 1)
T <- 1
R <- 0.5
a1 <- 0
P1 <- 1

update_fn <- function(theta) {
  list(R = array(theta[1], c(1, 1, 1)), phi = c(theta[2], 1))
}

prior_fn <- function(theta) {
  ifelse(all(theta > 0), sum(dnorm(theta, 0, 1, log = TRUE)), -Inf)
}

model <- ssm_mng(y, Z, T, R, a1, P1, phi = c(2, 1), 
  init_theta = c(0.5, 2), 
  distribution = c("gamma", "binomial"),
  u = cbind(1, rep(10, n)),
  update_fn = update_fn, prior_fn = prior_fn,
  state_names = "random_walk",
  # using default values, but being explicit for testing purposes
  D = matrix(0, 2, 1), C = matrix(0, 1, 1))

# smoothing based on approximating gaussian model
ts.plot(cbind(y, fast_smoother(model)), 
  col = 1:3, lty = c(1, 1, 2))

}
