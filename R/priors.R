#' Prior hyperparameters for the contamination-mixture model
#'
#' Collects the hyperparameters of all prior distributions used by
#' [mixture_fit()] and [baseline_fit()]. The defaults are exactly those of the
#' original research code, so that results are reproducible.
#'
#' @param beta_scale Prior variance \eqn{b_0} of each regression coefficient:
#'   \eqn{\beta \sim N_p(0, b_0 I_p)}. Also used for the rival regression
#'   \eqn{\beta_{bad}} of model `"M1"`.
#' @param sigma_shape,sigma_rate Inverse-gamma prior for the Gaussian error
#'   variance, \eqn{\sigma^2 \sim IG(a_\sigma, b_\sigma)} (Gaussian family only;
#'   also used for the rival-regression variance in Gaussian `"M1"`).
#' @param pi_shape1,pi_shape2 Beta prior for the probability \eqn{\pi} that a
#'   non-probability observation is compatible ("good"):
#'   \eqn{\pi \sim Beta(a_\pi, b_\pi)}. The contamination fraction is
#'   \eqn{1 - \pi}.
#' @param tau_shape,tau_rate Inverse-gamma prior for the variance \eqn{\tau^2}
#'   of the covariate-free contaminating outcome density (Gaussian `"M0"` and
#'   `"Mjoint"`).
#' @param mu0_var Prior variance of the mean \eqn{\mu_0} of the contaminating
#'   outcome density (Gaussian `"M0"`/`"Mjoint"`). The prior mean is the
#'   sample mean of the non-probability outcomes (an empirical-Bayes choice
#'   inherited from the original code).
#' @param q_shape1,q_shape2 Beta prior for the base rate \eqn{q} of the
#'   contaminating Bernoulli density (binomial `"M0"`/`"Mjoint"`).
#' @param lambda_shape,lambda_rate Gamma prior (shape, rate) for the rate
#'   \eqn{\lambda_0} of the contaminating Poisson density (Poisson
#'   `"M0"`/`"Mjoint"`).
#'
#' @return An object of class `bsm_prior` (a named list).
#' @seealso [mixture_fit()], [mixture_control()]
#' @examples
#' mixture_prior()
#' mixture_prior(pi_shape1 = 2, pi_shape2 = 1)  # a priori favour compatibility
#' @export
mixture_prior <- function(beta_scale = 25,
                          sigma_shape = 2, sigma_rate = 1,
                          pi_shape1 = 1, pi_shape2 = 1,
                          tau_shape = 2, tau_rate = 2, mu0_var = 100,
                          q_shape1 = 1, q_shape2 = 1,
                          lambda_shape = 1, lambda_rate = 1) {
  args <- list(beta_scale = beta_scale, sigma_shape = sigma_shape,
               sigma_rate = sigma_rate, pi_shape1 = pi_shape1,
               pi_shape2 = pi_shape2, tau_shape = tau_shape,
               tau_rate = tau_rate, mu0_var = mu0_var, q_shape1 = q_shape1,
               q_shape2 = q_shape2, lambda_shape = lambda_shape,
               lambda_rate = lambda_rate)
  for (nm in names(args)) .check_positive(args[[nm]], nm)
  structure(args, class = "bsm_prior")
}

#' @export
print.bsm_prior <- function(x, ...) {
  cat("Prior hyperparameters (bayesSampleMix)\n")
  v <- unlist(unclass(x))
  print(v)
  invisible(x)
}

#' Algorithmic settings for the contamination-mixture samplers
#'
#' @param c_infl Variance-inflation factor \eqn{c > 1} of the covariate
#'   density of contaminated observations in model `"Mjoint"`:
#'   good observations have \eqn{x \sim N(\hat\mu_P, \hat\Sigma_P)}, contaminated
#'   ones \eqn{x \sim N(\hat\mu_P, c\,\hat\Sigma_P)}. `NULL` (default) uses the
#'   family-specific default of the original code: 4 for Gaussian, 2.5 for
#'   binomial and Poisson.
#' @param pg_terms Number of terms \eqn{K} in the truncated sum-of-gammas
#'   representation used to draw Polya-Gamma variables (binomial family).
#'   Default 128, as in the original code.
#' @param irls_steps Number of Newton (IRLS) steps used to locate the
#'   conditional mode for the Laplace proposal of the Poisson coefficient
#'   update. Default 3, as in the original code. Larger values make the
#'   independence Metropolis-Hastings step closer to exact (see Details in
#'   [mixture_fit()]).
#' @param m1_update How the rival-regression component of model `"M1"` is
#'   updated when few observations are allocated to it. `"exact"` (default)
#'   always uses the correct full conditional given the allocated
#'   observations (a proper prior makes this well defined for any number of
#'   observations). `"original"` reproduces the original code, which draws
#'   \eqn{\beta_{bad}} from its prior whenever \eqn{\le p + 1} observations are
#'   allocated to the contaminated component; this ignores those observations
#'   and is not a valid Gibbs step, so it is provided only for reproducing
#'   earlier results.
#'
#' @return An object of class `bsm_control`.
#' @examples
#' mixture_control(c_infl = 3)
#' @export
mixture_control <- function(c_infl = NULL, pg_terms = 128L, irls_steps = 3L,
                            m1_update = c("exact", "original")) {
  if (!is.null(c_infl)) {
    .check_positive(c_infl, "c_infl")
    if (c_infl <= 1) {
      stop("`c_infl` must be greater than 1 (contaminated covariates are more dispersed).",
           call. = FALSE)
    }
  }
  pg_terms <- .check_count(pg_terms, "pg_terms", 1L)
  irls_steps <- .check_count(irls_steps, "irls_steps", 1L)
  m1_update <- match.arg(m1_update)
  structure(list(c_infl = c_infl, pg_terms = pg_terms, irls_steps = irls_steps,
                 m1_update = m1_update), class = "bsm_control")
}

## family-specific default inflation factor (from the original scripts)
.default_c_infl <- function(family) {
  switch(family, gaussian = 4, binomial = 2.5, poisson = 2.5)
}
