#' Simulate probability and contaminated non-probability samples
#'
#' Generates example data with a known data-generating process for
#' demonstrations, tests and simulation studies.
#'
#' Covariates \eqn{x_1, \dots, x_{k}} (\eqn{k} = `length(beta) - 1`) are
#' independent standard normal in the probability sample and in the
#' compatible part of the non-probability sample. Exactly
#' `round(contamination * n_nonprob)` NP units are contaminated: their
#' covariates are drawn from \eqn{N(x_{shift}, 1)} and their outcomes from the
#' same family with linear predictor \eqn{x'\beta_{bad} + \delta}
#' (`bad_beta`, `outcome_shift`). NP rows are randomly permuted.
#'
#' @param family `"gaussian"`, `"binomial"` or `"poisson"`.
#' @param n_prob,n_nonprob Sample sizes.
#' @param contamination Fraction of contaminated NP units, in \[0, 1\].
#' @param beta True target coefficients (intercept first). Defaults:
#'   `c(2, 5)` (Gaussian), `c(0.2, 1, -0.8)` (binomial),
#'   `c(1.2, 0.5, -0.4)` (Poisson).
#' @param sigma Error standard deviation (Gaussian only).
#' @param bad_beta Coefficients of the contaminating mechanism. Default: the
#'   intercept of `beta` with all slopes sign-flipped.
#' @param x_shift Mean shift of contaminated units' covariates.
#' @param outcome_shift Constant added to contaminated units' linear
#'   predictor.
#' @param seed Optional seed (the global RNG state is restored afterwards).
#' @return A list with data frames `prob` (columns `y`, `x1`, ...) and
#'   `nonprob` (same columns plus `contaminated`, 1 for contaminated units),
#'   and the list `truth` of generating values.
#' @examples
#' d <- simulate_mixture_data("poisson", n_prob = 50, n_nonprob = 40,
#'                            contamination = 0.25, seed = 1)
#' str(d$nonprob)
#' table(d$nonprob$contaminated)
#' @export
simulate_mixture_data <- function(family = c("gaussian", "binomial", "poisson"),
                                  n_prob = 200, n_nonprob = 200, contamination = 0.25,
                                  beta = NULL, sigma = 1, bad_beta = NULL,
                                  x_shift = 0, outcome_shift = 0, seed = NULL) {
  if (is.character(family) && length(family) > 1L) family <- family[1L]
  family <- .match_family(family)
  n_prob <- .check_count(n_prob, "n_prob", 1L)
  n_nonprob <- .check_count(n_nonprob, "n_nonprob", 1L)
  if (!is.numeric(contamination) || length(contamination) != 1L ||
      contamination < 0 || contamination > 1) {
    stop("`contamination` must be a number in [0, 1].", call. = FALSE)
  }
  if (is.null(beta)) {
    beta <- switch(family, gaussian = c(2, 5), binomial = c(0.2, 1, -0.8),
                   poisson = c(1.2, 0.5, -0.4))
  }
  if (length(beta) < 2L) stop("`beta` needs an intercept and at least one slope.", call. = FALSE)
  if (is.null(bad_beta)) bad_beta <- c(beta[1], -beta[-1])
  if (length(bad_beta) != length(beta)) stop("`bad_beta` must have the same length as `beta`.", call. = FALSE)
  .check_positive(sigma, "sigma")
  k <- length(beta) - 1L
  n_bad <- round(contamination * n_nonprob)
  n_good <- n_nonprob - n_bad

  gen_y <- function(eta) {
    n <- length(eta)
    switch(family,
           gaussian = eta + stats::rnorm(n, 0, sigma),
           binomial = stats::rbinom(n, 1, .sigmoid(eta)),
           poisson = stats::rpois(n, exp(pmin(eta, 30))))
  }
  make_x <- function(n, shift) matrix(stats::rnorm(n * k, mean = shift), n, k)
  lin <- function(X, b) drop(cbind(rep(1, nrow(X)), X) %*% b)

  .with_seed(seed, {
    Xp <- make_x(n_prob, 0)
    yp <- gen_y(lin(Xp, beta))
    Xg <- make_x(n_good, 0)
    yg <- gen_y(lin(Xg, beta))
    Xb <- make_x(n_bad, x_shift)
    yb <- gen_y(lin(Xb, bad_beta) + outcome_shift)
    perm <- sample.int(n_nonprob)
    xn <- paste0("x", seq_len(k))
    prob <- data.frame(y = yp, Xp)
    names(prob) <- c("y", xn)
    nonprob <- data.frame(y = c(yg, yb), rbind(Xg, Xb),
                          contaminated = c(rep(0L, n_good), rep(1L, n_bad)))
    names(nonprob) <- c("y", xn, "contaminated")
    nonprob <- nonprob[perm, , drop = FALSE]
    rownames(nonprob) <- NULL
    list(prob = prob, nonprob = nonprob,
         truth = list(family = family, beta = beta, bad_beta = bad_beta, sigma = sigma,
                      contamination = n_bad / n_nonprob, x_shift = x_shift,
                      outcome_shift = outcome_shift))
  })
}
