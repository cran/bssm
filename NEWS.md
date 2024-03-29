bssm 2.0.2 (Release date: 2023-10-18)
=====================================
   * Switched to markdown NEWS with a plan to be more clear about the future
     changes in the package.
   * Added more details to the `?bssm` help page.
   * Added more details to the `?bssm_prior` help page.
   * Added option to extract only hyperparameters in `as_draws` method. Also
     fixed a bug in `as_draws` which caused the it to ignore `states` argument.
   * Added a default plot method for the `run_mcmc` output.
   * Fixed the aliases of the main help page to accomodate changes in roxygen2.
   * Removed explicit C++ version requirement as required by new CRAN policies.
   * Removed `magrittr` dependency and switched to native pipe, leading to 
     requirement for R 4.1.0+.
   * Added Sys.setenv("OMP_NUM_THREADS" = 2) to (partially) fix CRAN issues with 
     parallelisation on Debian.

bssm 2.0.1 (Release date: 2022-05-02)
==============
   * Fixed weights to one in case of non-linear model with mcmc_type="approx".
   * Adjusted tolerance of some testthat tests to comply with CRAN's MKL checks.

bssm 2.0.0 (Release date: 2021-11-26)
==============
   * Added a progress bar for run_mcmc.
   * Added a fitted method for extraction of summary statistics of posterior
     predictive distribution p(y_t | y_1, ..., y_n) for t = 1, ..., n.
   * Rewrote the summary method completely, which now returns data.frame. This
     also resulted in some changes in order of the function arguments.
   * The output of predict method is now a data frame with column weight
     corresponding to the IS-weights in case of IS-MCMC. Previously resampling
     was done internally, but now this is left for the user if needed
     (i.e. for drawing state trajectories).
   * The asymptotic_var and iact functions are now exported to users, and they
     also contain alternative methods based on the posterior package.
   * New function estimate_ess can be used to compute effective sample size
     from weighted MCMC.
   * Added compatibility with the posterior package by defining as_draws
     method for converting run_mcmc output to draws_df object.
   * New function check_diagnostics for quick glance of ESS and Rhat values.
   * Large number of new tests, and improved documentation with added examples.
   * Large number of internal tweaks so that the package complies with
     goodpractices package and Ropensci statistical software standards.

bssm 1.1.7-1 (Release date: 2021-09-21)
==============
   * Fixed an error in automatic tests due to lack of fixed RNG seed.

bssm 1.1.7 (Release date: 2021-09-20)
==============
   * Added a function cpp_example_model which can be used to extract and
     compile some non-linear and SDE models used in the examples and vignettes.
   * Added as_draws method for run_mcmc output so samples can be analysed using
     the posterior package.
   * Added more examples.
   * Fixed a tolerance of one MCMC test to pass the test on OSX as well.
   * Fixed a bug in iterated extended Kalman smoothing which resulted incorrect
     estimates.

bssm 1.1.6 (Release date: 2021-09-06)
==============
   * Cleaned some codes and added lots of tests in line with pkgcheck tests.
   * Fixed a bug in EKF-based particle filter which returned filtered estimates
     also in place of one-step ahead predictions.
   * Fixed a bug which caused an error in suggest_N for nlg_ssm.
   * Fixed a bug which caused incorrect sampling of smoothing distribution for
     ar1_lg model when predicting past or when using simulation smoother.
   * Fixed a bug which caused an error when predicting past values in
     multivariate time series case.
   * Fixed log-likelihood computation for gamma model with non-constant shape
     parameter when using (intermediate) Gaussian approximation.
   * Fixed sampling of negative binomial distribution in predict method, which
     used std::negative_binomial which converts non-integer phi to integer.
     Sampling now uses Gamma-Poisson mixture for simulation.


bssm 1.1.5 (Release date: 2021-06-14)
==============
   * Added explicit check for nsim > 0 in predict method as sample function
     works with missing argument causing crypting warnings later.
   * Updated drownings data until 2019 and changed the temperature variable
     to an average over three stations.
   * Improved checks for observations and distributions in model building.

bssm 1.1.4 (Release date: 2021-04-13)
==============
   * Better documentation for SV model, and changed ordering of arguments to
     emphasise the recommended parameterization.
   * Fixed predict method for SV model.
   * Removed parallelization in one example which failed on Solaris for some
     unknown reason.

bssm 1.1.3-2 (Release date: 2021-02-24)
==============
   * Fixed missing parenthesis causing compilation fail in case of no OpenMP
     support.
   * Added pandoc version >= 1.12.3 to system requirements.
   * Restructured C++ classes so no R structures are present in OpenMP regions.

bssm 1.1.3-1 (Release date: 2021-02-22)
==============
   * Fixed PM-MCMC and DA-MCMC for SDE models and added an example to `ssm_sde`.
   * Fixed the state covariance estimates of IS-MCMC, approx-MCMC, and
     Gaussian MCMC when output_type = "summary".
   * Fixed memory leaks due to uninitialized variables due to aborted particle
     filter.
   * Fixed numerical issues of multivariate normal density for nonlinear
     models.
   * Removed dependency on R::lchoose for safer parallel code.
   * Added vignette for SDE models.
   * Updated citation information and streamlined the main vignette.

bssm 1.1.2 (Release date: 2021-02-08)
==============
   * Changed the definition of D in ssm_ulg and ssm_ung, functions now accept
     D as scalar or vector as
     was originally intended.
   * Fixed a segfault issue with parallel state sampling in general
     ssm_ulg/mlg/ung/mng models caused by calls to R function inside parallel
     region.
   * Fixed a bug from version 1.0.0 in IS1 type sampling which actually lead
     to IS2 type sampling.
   * Fixed out-of-bounds error in IS3 sampling.
   * Fixed weight computations for multivariate nonlinear models in case of
     psi-APF in some border cases with non-standard H.
   * Removed Armadillo bound checks for efficiency gains.

bssm 1.1.1 (Release date: 2021-01-22)
==============

   * Added missing scaling for Gamma distribution in importance sampling
     weights for added numerical robustness.
   * Fixed sequential importance sampling for multivariate non-gaussian models.
   * Fixed simulation smoother for multivariate Gaussian models.

bssm 1.1.0 (Release date: 2021-01-19)
==============

   * Added function `suggest_N` which can be used to choose
     suitable number of particles for IS-MCMC.
   * Added function `post_correct` which can be used to update
     previous approximate MCMC with IS-weights.
   * Gamma priors are now supported in easy-to-use models such as `bsm_lg`.
   * The adaptation of the proposal distribution now continues also after the
     burn-in by default.
   * Changed default MCMC type to typically most efficient and robust IS2.
   * Renamed `nsim` argument to `particles` in most of the R functions (`nsim`
     also works with a warning).
   * Fixed a bug with bsm models with covariates, where all standard deviation
     parameters were fixed. This resulted error within MCMC algorithms.
   * Fixed a dimension drop bug in the predict method which caused error for
     univariate models.
   * Fixed some docs and added more examples.
   * Fixed few typos in vignette (thanks Kyle Hussman)
   * Reduced runtime of MCMC in growth model vignette as requested by CRAN.


bssm 1.0.1-1 (Release date: 2020-11-12)
==============

  * Added an argument `future` for predict method which allows
    predictions for current time points by supplying the original model
    (e.g., for posterior predictive checks).
    At the same time the argument name `future_model` was changed to `model`.
  * Fixed a bug in summary.mcmc_run which resulted error when
    trying to obtain summary for states only.
  * Added a check for Kalman filter for a degenerate case where all
    observational level and state level variances are zero.
  * Renamed argument `n_threads` to `threads` for consistency
    with `iter` and `burnin` arguments.
  * Improved documentation, added examples.
  * Added a vignette regarding psi-APF for non-linear models.

bssm 1.0.0 (Release date: 2020-06-09)
==============
Major update

  * Major changes for model definitions, now model updating and priors
    can be defined via R functions (non-linear and SDE models still rely on
    C++ snippets).
  * Added support for multivariate non-Gaussian models.
  * Added support for gamma distributions.
  * Added the function as.data.frame for mcmc output which converts the MCMC
    samples to data.frame format for easier post-processing.
  * Added truncated normal prior.
  * Many argument names and model building functions have been changed for
    clarity and consistency.
  * Major overhaul of C++ internals which can bring minor efficiency gains
    and smaller installation size.
  * Allow zero as initial value for positive-constrained parameters of bsm
    models.
  * Small changes to summary method which can now return also only summaries
    of the states.
  * Fixed a bug in initializing run_mcmc for negative binomial model.
  * Fixed a bug in phi-APF for non-linear models.
  * Reimplemented predict method which now always produces data frame of
    samples.

bssm 0.1.11 (Release date: 2020-02-25)
==============
  * Switched (back) to approximate posterior in RAM for PM-SPDK and PM-PSI,
    as it seems to work better with noisy likelihood estimates.
  * Print and summary methods for MCMC output are now coherent in their output.

bssm 0.1.10 (Release date: 2020-02-04)
==============
  * Fixed missing weight update for IS-SPDK without OPENMP flag.
  * Removed unused usage argument ... from expand_sample.

bssm 0.1.9 (Release date: 2020-01-27)
==============
  * Fixed state sampling for PM-MCMC with SPDK.
  * Added ts attribute for svm model.
  * Corrected asymptotic variance for summary methods.

bssm 0.1.8-1 (Release date: 2019-12-20)
==============
  * Tweaked tests in order to pass MKL case at CRAN.

bssm 0.1.8 (Release date: 2019-09-23)
==============
  * Fixed a bug in predict method which prevented the method working in case
    of ngssm models.
  * Fixed a bug in predict method which threw an error due to dimension drop of
    models with single state.
  * Fixed issues with the vignette.

bssm 0.1.7 (Release date: 2019-03-19)
==============
  * Fixed a bug in EKF smoother which resulted wrong smoothed state estimates
    in case of partially missing multivariate observations. Thanks for Santeri
    Karppinen for spotting the bug.
  * Added twisted SMC based simulation smoothing algorithm for Gaussian models,
    as an alternative to Kalman smoother based simulation.

bssm 0.1.6-1 (Release date: 2018-11-20)
==============
  * Fixed wrong dimension declarations in pseudo-marginal MCMC and logLik
    methods for SDE and ng_ar1 models.
  * Added a missing Jacobian for ng_bsm and bsm models using IS-correction.
  * Changed internal parameterization of ng_bsm and bsm models from
    log(1+theta) to log(theta).

bssm 0.1.5 (Release date: 2018-05-23)
==============
  * Fixed the Cholesky decomposition in filtering recursions of multivariate
    models.
  * as_gssm now works for multivariate Gaussian models of KFAS as well.
  * Fixed several issues regarding partially missing observations in
    multivariate models.
  * Added the MASS package to Suggests as it is used in some unit tests.
  * Added missing type argument to SDE MCMC call with delayed acceptance.

bssm 0.1.4-1 (Release date: 2018-02-04)
==============
  * Fixed the use of uninitialized values in psi-filter from version 0.1.3.

bssm 0.1.4 (Release date: 2018-02-04)
==============
  * MCMC output can now be defined with argument `type`. Instead of returning
    joint posterior samples, run_mcmc can now return only marginal samples of
    theta, or summary statistics of the states.
  * Due to the above change, argument `sim_states` was removed from the
    Gaussian MCMC methods.
  * MCMC functions are now less memory intensive, especially with
    `type="theta"`.


bssm 0.1.3 (Release date: 2018-01-07)
==============
  * Streamlined the output of the print method for MCMC results.
  * Fixed major bugs in predict method which caused wrong values for the
    prediction intervals.
  * Fixed some package dependencies.
  * Sampling for standard deviation parameters of BSM and their non-Gaussian
    counterparts is now done in logarithmic scale for slightly increased
    efficiency.
  * Added a new model class ar1 for univariate (possibly noisy) Gaussian AR(1)
    processes.
  * MCMC output now includes posterior predictive distribution of states for
    one step ahead to the future.

bssm 0.1.2 (Release date: 2017-11-21)
==============
  * API change for run_mcmc: All MCMC methods are now under the argument
    method, instead of having separate arguments for delayed acceptance and IS
    schemes.
  * summary method for MCMC output now omits the computation of SE and ESS in
    order to speed up the function.
  * Added new model class lgg_ssm, which is a linear-Gaussian model defined
    directly via C++ like non-linear ssm_nlg models. This allows more flexible
    prior definitions and complex system matrix constructions.
  * Added another new model class, ssm_sde, which is a model with continuous
    state dynamics defined as SDE. These too are defined via couple
    simple C++ functions.
  * Added non-gaussian AR(1) model class.
  * Added argument nsim for predict method, which allows multiple draws per
    MCMC iteration.
  * The noise multiplier matrices H and R in ssm_nlg models can now depend on
    states.

bssm 0.1.1-1 (Release date: 2017-06-27)
==============
  * Use byte compiler.
  * Skip tests relying in certain numerical precision on CRAN.

bssm 0.1.1 (Release date: 2017-06-27)
==============

  * Switched from C++11 PRNGs to sitmo.
  * Fixed some portability issues in C++ codes.

bssm 0.1.0 (Release date: 2017-06-24)
==============

  * Initial release.
