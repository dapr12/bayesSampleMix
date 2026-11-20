#' bayesSampleMix: Contamination-aware Bayesian integration of probability
#' and non-probability samples
#'
#' Fits Gaussian, logistic and Poisson regression models to a trusted
#' probability sample combined with a non-probability sample in which an
#' unknown subset of units may be incompatible with the target model. A
#' latent compatibility indicator per non-probability unit is sampled by data
#' augmentation, giving posterior inference for the regression coefficients,
#' the contamination fraction and unit-level filter scores.
#'
#' Main functions: [mixture_fit()], [baseline_fit()],
#' [overlap_diagnostics()], [simulate_mixture_data()].
#'
#' @keywords internal
#' @importFrom stats coef predict vcov nobs
"_PACKAGE"
