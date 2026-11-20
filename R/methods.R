## S3 methods and accessors for bsm_fit and bsm_baseline objects.

.target_draws <- function(object) {
  object$draws[, object$coef_names, drop = FALSE]
}

.posterior_table <- function(d, probs) {
  ess <- effective_size(d)
  sds <- apply(d, 2, stats::sd)
  qs <- t(apply(d, 2, stats::quantile, probs = c(probs[1], 0.5, probs[2]), names = FALSE))
  tab <- cbind(mean = colMeans(d), sd = sds, qs, ess = ess, mcse = sds / sqrt(ess))
  colnames(tab)[3:5] <- c(sprintf("%.1f%%", 100 * probs[1]), "50%", sprintf("%.1f%%", 100 * probs[2]))
  tab
}

.check_probs <- function(level) {
  if (!is.numeric(level) || length(level) != 1L || level <= 0 || level >= 1) {
    stop("`level` must be a single number in (0, 1).", call. = FALSE)
  }
  c((1 - level) / 2, 1 - (1 - level) / 2)
}

.family_label <- function(f) {
  switch(f, gaussian = "Gaussian (identity link)", binomial = "Binomial (logit link)",
         poisson = "Poisson (log link)")
}

.model_label <- function(m) {
  switch(m, M0 = "M0 (covariate-free contaminating density)",
         Mjoint = "Mjoint (contaminating density + covariate leverage term)",
         M1 = "M1 (rival regression for contaminated units)")
}

#' Posterior summaries of a contamination-mixture fit
#'
#' @param object,x A `bsm_fit` object from [mixture_fit()] or
#'   [mixture_fit_matrix()] (or a `bsm_baseline` object for the methods that
#'   also apply to [baseline_fit()]).
#' @param level Probability mass of the equal-tailed credible intervals.
#' @param digits Number of significant digits to print.
#' @param ... Not used.
#'
#' @return `summary()` returns an object of class `summary.bsm_fit` with
#'   elements
#'   \item{coefficients}{posterior mean, standard deviation, credible
#'     interval, median, effective sample size and Monte Carlo standard error
#'     of each target-model coefficient;}
#'   \item{dispersion}{the same for \eqn{\sigma^2} (Gaussian family only);}
#'   \item{contamination}{the same for \eqn{1 - \pi} and the realised
#'     contamination fraction;}
#'   \item{filter}{number of NP units, posterior expected number of
#'     compatible units \eqn{\sum_j s_j}, number of units with
#'     \eqn{s_j < 0.5}, and mean indicator entropy;}
#'   plus model settings. All quantities are computed from the retained
#'   posterior draws; no asymptotic approximations are used.
#' @examples
#' d <- simulate_mixture_data("gaussian", n_prob = 100, n_nonprob = 100,
#'                            contamination = 0.25, seed = 4)
#' fit <- mixture_fit(y ~ x1, d$prob, d$nonprob, n_iter = 600, burnin = 200, seed = 1)
#' summary(fit)
#' coef(fit)
#' contamination(fit)
#' @export
summary.bsm_fit <- function(object, level = 0.95, ...) {
  probs <- .check_probs(level)
  d <- object$draws
  coef_tab <- .posterior_table(d[, object$coef_names, drop = FALSE], probs)
  disp_tab <- if (object$family == "gaussian") .posterior_table(d[, "sigma2", drop = FALSE], probs) else NULL
  cont_tab <- .posterior_table(d[, c("contamination", "contamination_realized"), drop = FALSE], probs)
  s <- object$filter_score
  filt <- c(n_nonprob = length(s), expected_compatible = sum(s),
            flagged_s_below_0.5 = sum(s < 0.5), mean_entropy = mean(object$indicator_entropy))
  structure(list(call = object$call, family = object$family, model = object$model,
                 coefficients = coef_tab, dispersion = disp_tab, contamination = cont_tab,
                 filter = filt, n = object$n, na_info = object$na_info,
                 mcmc = object$mcmc, n_draws = nrow(d), acceptance = object$acceptance,
                 leverage = object$leverage, c_infl = object$c_infl, level = level),
            class = "summary.bsm_fit")
}

#' @rdname summary.bsm_fit
#' @export
print.summary.bsm_fit <- function(x, digits = 3, ...) {
  cat("Bayesian contamination-mixture fit (bayesSampleMix)\n")
  if (!is.null(x$call)) {
    cat("Call: ")
    print(x$call)
  }
  cat("Family:  ", .family_label(x$family), "\n", sep = "")
  cat("Model:   ", .model_label(x$model), "\n", sep = "")
  if (!is.null(x$leverage)) {
    cat("Leverage covariates: ", paste(x$leverage, collapse = ", "),
        "  (inflation c = ", x$c_infl, ")\n", sep = "")
  }
  cat(sprintf("Samples: n_P = %d, n_NP = %d", x$n[["prob"]], x$n[["nonprob"]]))
  if (!is.null(x$na_info) && any(x$na_info > 0)) {
    cat(sprintf("  (dropped for missingness: %d P, %d NP)", x$na_info[1], x$na_info[2]))
  }
  cat("\n")
  cat(sprintf("MCMC:    %d iterations, burn-in %d, thin %d -> %d draws\n",
              x$mcmc$n_iter, x$mcmc$burnin, x$mcmc$thin, x$n_draws))
  if (!is.na(x$acceptance)) {
    cat(sprintf("Metropolis-Hastings acceptance rate (target block): %.3f\n", x$acceptance))
  }
  cat(sprintf("\nTarget-model coefficients (posterior; %g%% credible intervals):\n", 100 * x$level))
  print(signif(x$coefficients, digits))
  if (!is.null(x$dispersion)) {
    cat("\nError variance:\n")
    print(signif(x$dispersion, digits))
  }
  cat("\nContamination:\n")
  print(signif(x$contamination, digits))
  cat("\nFilter scores of NP units (s_j = posterior P(compatible)):\n")
  f <- x$filter
  cat(sprintf("  expected compatible units: %.1f of %d; units with s_j < 0.5: %d; mean entropy: %.3f\n",
              f[["expected_compatible"]], as.integer(f[["n_nonprob"]]),
              as.integer(f[["flagged_s_below_0.5"]]), f[["mean_entropy"]]))
  low_ess <- x$coefficients[, "ess"] < 100
  if (any(low_ess, na.rm = TRUE)) {
    cat("\nNote: effective sample size < 100 for some coefficients; consider more iterations.\n")
  }
  invisible(x)
}

#' @rdname summary.bsm_fit
#' @export
print.bsm_fit <- function(x, digits = 3, ...) {
  cat("Bayesian contamination-mixture fit (bayesSampleMix)\n")
  cat("Family: ", x$family, "   Model: ", x$model, "\n", sep = "")
  cat(sprintf("n_P = %d, n_NP = %d, %d posterior draws\n", x$n[["prob"]], x$n[["nonprob"]], nrow(x$draws)))
  cat("\nPosterior mean coefficients:\n")
  print(signif(coef(x), digits))
  cm <- mean(x$draws[, "contamination"])
  ci <- stats::quantile(x$draws[, "contamination"], c(0.025, 0.975), names = FALSE)
  cat(sprintf("\nContamination (1 - pi): %.3f  [95%% CrI %.3f, %.3f]\n", cm, ci[1], ci[2]))
  invisible(x)
}

#' @param type For `coef()`: posterior `"mean"` (default) or `"median"`.
#' @rdname summary.bsm_fit
#' @export
coef.bsm_fit <- function(object, type = c("mean", "median"), ...) {
  type <- match.arg(type)
  d <- .target_draws(object)
  if (type == "mean") colMeans(d) else apply(d, 2, stats::median)
}

#' @rdname summary.bsm_fit
#' @export
vcov.bsm_fit <- function(object, ...) {
  stats::cov(.target_draws(object))
}

#' @rdname summary.bsm_fit
#' @export
nobs.bsm_fit <- function(object, ...) {
  sum(object$n)
}

#' Posterior draws as a matrix
#'
#' @param x A `bsm_fit` or `bsm_baseline` object.
#' @param which For `bsm_fit`: `"target"` (target-model parameters and
#'   contamination quantities, default), `"bad"` (contaminating-component
#'   parameters) or `"all"`.
#' @param ... Not used.
#' @return A numeric matrix with one row per retained draw; suitable for
#'   `coda::mcmc()`.
#' @examples
#' d <- simulate_mixture_data("poisson", n_prob = 80, n_nonprob = 80,
#'                            contamination = 0.2, seed = 5)
#' fit <- mixture_fit(y ~ x1, d$prob, d$nonprob, family = "poisson",
#'                    n_iter = 400, burnin = 100, seed = 1)
#' head(as.matrix(fit))
#' @export
as.matrix.bsm_fit <- function(x, which = c("target", "bad", "all"), ...) {
  which <- match.arg(which)
  switch(which, target = x$draws, bad = x$bad_draws, all = cbind(x$draws, x$bad_draws))
}

#' Filter scores and contamination summaries
#'
#' `filter_scores()` returns the posterior probability \eqn{s_j} that each
#' non-probability unit is compatible with the target model (the
#' Rao-Blackwellised posterior mean of its latent indicator). Small values flag
#' likely contaminated units. `contamination()` summarises the posterior of the
#' contamination fraction \eqn{1 - \pi}.
#'
#' @param object A `bsm_fit` object.
#' @param entropy If `TRUE`, return a data frame that also contains the
#'   binary entropy (nats) of each score, a measure of classification
#'   uncertainty (maximal, \eqn{\log 2}, at \eqn{s_j = 0.5}).
#' @param level Credible-interval level.
#' @return `filter_scores()`: a numeric vector named by the rows of
#'   `nonprob_data` used in the fit, or a data frame if `entropy = TRUE`.
#'   `contamination()`: a named vector with posterior mean, median and
#'   credible limits of \eqn{1 - \pi}.
#' @examples
#' d <- simulate_mixture_data("gaussian", n_prob = 100, n_nonprob = 100,
#'                            contamination = 0.3, seed = 6)
#' fit <- mixture_fit(y ~ x1, d$prob, d$nonprob, n_iter = 600, burnin = 200, seed = 1)
#' s <- filter_scores(fit)
#' tapply(s, d$nonprob$contaminated, mean)   # low scores for contaminated units
#' contamination(fit)
#' @export
filter_scores <- function(object, entropy = FALSE) {
  if (!inherits(object, "bsm_fit")) stop("`object` must be a bsm_fit object.", call. = FALSE)
  s <- object$filter_score
  nm <- object$nonprob_rows
  if (!entropy) return(stats::setNames(s, nm))
  data.frame(row = nm, filter_score = s, entropy = object$indicator_entropy,
             stringsAsFactors = FALSE)
}

#' @rdname filter_scores
#' @export
contamination <- function(object, level = 0.95) {
  if (!inherits(object, "bsm_fit")) stop("`object` must be a bsm_fit object.", call. = FALSE)
  probs <- .check_probs(level)
  v <- object$draws[, "contamination"]
  q <- stats::quantile(v, c(probs[1], 0.5, probs[2]), names = FALSE)
  c(mean = mean(v), median = q[2], lower = q[1], upper = q[3])
}

## ---- bsm_baseline methods -------------------------------------------------

#' @rdname baseline_fit
#' @param object A `bsm_baseline` object.
#' @param level Credible-interval level for `summary()`.
#' @param digits Significant digits for printing.
#' @param ... Not used.
#' @export
summary.bsm_baseline <- function(object, level = 0.95, ...) {
  probs <- .check_probs(level)
  tab <- .posterior_table(object$draws, probs)
  structure(list(call = object$call, family = object$family, table = tab, n = object$n,
                 mcmc = object$mcmc, acceptance = object$acceptance, level = level),
            class = "summary.bsm_baseline")
}

#' @rdname baseline_fit
#' @export
print.summary.bsm_baseline <- function(x, digits = 3, ...) {
  cat("Single-sample Bayesian regression (bayesSampleMix baseline)\n")
  cat("Family: ", .family_label(x$family), ";  n = ", x$n, "\n", sep = "")
  cat(sprintf("MCMC: %d iterations, burn-in %d, thin %d\n", x$mcmc$n_iter, x$mcmc$burnin, x$mcmc$thin))
  if (!is.na(x$acceptance)) cat(sprintf("MH acceptance rate: %.3f\n", x$acceptance))
  cat(sprintf("\nPosterior summaries (%g%% credible intervals):\n", 100 * x$level))
  print(signif(x$table, digits))
  invisible(x)
}

#' @rdname baseline_fit
#' @export
print.bsm_baseline <- function(x, digits = 3, ...) {
  cat("Single-sample Bayesian regression (", x$family, ", n = ", x$n, ")\n", sep = "")
  print(signif(coef(x), digits))
  invisible(x)
}

#' @rdname baseline_fit
#' @export
coef.bsm_baseline <- function(object, ...) {
  colMeans(object$draws[, object$coef_names, drop = FALSE])
}

#' @rdname as.matrix.bsm_fit
#' @export
as.matrix.bsm_baseline <- function(x, ...) x$draws
