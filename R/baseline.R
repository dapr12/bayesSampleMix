#' Single-sample Bayesian regression (reference estimators)
#'
#' Fits the target regression model alone, with the same priors and the same
#' MCMC updates as the target component of [mixture_fit()], but without any
#' contamination component. It reproduces the original baseline functions
#' `fit_bayes_lm()`, `fit_bayes_logit()` and `fit_bayes_pois()` and is
#' useful to compare
#' \itemize{
#'   \item the *probability-only* estimator (`data = prob_data`), and
#'   \item the *naively pooled* estimator (`data = rbind(prob_data,
#'     nonprob_data)`), which treats every NP unit as compatible,
#' }
#' with the contamination-aware fit.
#'
#' @inheritParams mixture_fit
#' @param data Data frame containing the variables in `formula`.
#' @param x For `baseline_fit_matrix()`, the design matrix; for the `print`
#'   methods, a `bsm_baseline` or `summary.bsm_baseline` object.
#' @param y Response vector (matrix interface).
#' @param n_iter,burnin Total iterations (including burn-in) and burn-in; the
#'   defaults (2000, 1000) are those of the original baseline functions.
#'
#' @return An object of class `bsm_baseline` with elements `draws`,
#'   `coef_names`, `family`, `acceptance` (Poisson only), `mcmc` and the
#'   information needed by [predict()].
#' @examples
#' d <- simulate_mixture_data("gaussian", n_prob = 80, n_nonprob = 80,
#'                            contamination = 0.4, seed = 2)
#' prob_only <- baseline_fit(y ~ x1, d$prob, n_iter = 600, burnin = 200, seed = 1)
#' pooled <- baseline_fit(y ~ x1, rbind(d$prob, d$nonprob[names(d$prob)]),
#'                        n_iter = 600, burnin = 200, seed = 1)
#' rbind(prob_only = coef(prob_only), pooled = coef(pooled))
#' @export
baseline_fit <- function(formula, data, family = c("gaussian", "binomial", "poisson"),
                         n_iter = 2000, burnin = 1000, thin = 1, seed = NULL,
                         prior = mixture_prior(), control = mixture_control(),
                         na_action = c("omit", "fail")) {
  if (is.character(family) && length(family) > 1L) family <- family[1L]
  family <- .match_family(family)
  na_action <- match.arg(na_action)
  des <- .build_design(formula, data, NULL, na_action)
  fit <- baseline_fit_matrix(des$Xp, des$yp, family = family, n_iter = n_iter,
                             burnin = burnin, thin = thin, seed = seed,
                             prior = prior, control = control)
  fit$call <- match.call()
  fit$formula <- des$formula
  fit$terms <- des$terms
  fit$xlevels <- des$xlevels
  fit$contrasts <- des$contrasts
  fit
}

#' @rdname baseline_fit
#' @export
baseline_fit_matrix <- function(x, y, family = "gaussian", n_iter = 2000, burnin = 1000,
                                thin = 1, seed = NULL, prior = mixture_prior(),
                                control = mixture_control()) {
  family <- .match_family(family)
  mcmc <- .check_mcmc(n_iter, burnin, thin)
  seed <- .check_seed(seed)
  if (!inherits(prior, "bsm_prior")) stop("`prior` must be created by mixture_prior().", call. = FALSE)
  if (!inherits(control, "bsm_control")) stop("`control` must be created by mixture_control().", call. = FALSE)
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  .check_design(x, "x")
  x <- .fix_colnames(x)
  y <- .check_response(y, family, "`y`")
  if (length(y) != nrow(x)) stop("length(y) must equal nrow(x).", call. = FALSE)
  if (family == "gaussian" && (nrow(x) <= ncol(x) || qr(x)$rank < ncol(x))) {
    stop("family = \"gaussian\" needs a full-rank design with more rows than columns.",
         call. = FALSE)
  }
  empty_x <- x[0, , drop = FALSE]
  res <- .with_seed(seed, {
    engine <- .make_engine(family, Xp = x, yp = y, Xnp = empty_x, ynp = numeric(0),
                           model = "M0", prior = prior, control = control)
    .run_baseline(engine, mcmc)
  })
  structure(list(family = family, draws = res$draws, coef_names = colnames(x),
                 acceptance = res$acceptance, n = nrow(x), mcmc = mcmc, seed = seed,
                 prior = prior, control = control, call = match.call(),
                 terms = NULL, xlevels = NULL, contrasts = NULL),
            class = "bsm_baseline")
}
