#' Contamination-aware Bayesian integration of probability and
#' non-probability samples
#'
#' Fits a regression model to a trusted probability (P) sample and a
#' non-probability (NP) sample in which an unknown subset of observations may
#' be incompatible with the target model ("contaminated"). Each NP observation
#' carries a latent indicator of compatibility; compatible observations are
#' pooled with the P sample to estimate the target regression, while
#' contaminated observations are absorbed by a separate component.
#'
#' @section Model:
#' For P-sample units \eqn{i} the outcome follows the target model
#' \eqn{f(y_i \mid x_i, \theta)} (Gaussian, logistic or Poisson regression with
#' coefficients \eqn{\beta}). For NP units \eqn{j},
#' \deqn{z_j \sim \mathrm{Bernoulli}(\pi),}
#' \deqn{p(y_j, x_j \mid z_j = 1) = f(y_j \mid x_j, \theta)\, g_1(x_j),\qquad
#'       p(y_j, x_j \mid z_j = 0) = h(y_j \mid x_j, \phi)\, g_0(x_j).}
#' The contaminating outcome density \eqn{h} depends on `model`:
#' \describe{
#'   \item{`"M0"`}{a covariate-free density: \eqn{N(\mu_0, \tau^2)} (Gaussian),
#'     \eqn{Bernoulli(q)} (binomial), \eqn{Poisson(\lambda_0)} (Poisson); and
#'     \eqn{g_1 = g_0 = 1}.}
#'   \item{`"Mjoint"`}{as `"M0"`, plus a covariate ("leverage") term for the
#'     continuous predictors \eqn{x_c} selected by `leverage`:
#'     \eqn{g_1 = N(x_c; \hat\mu_P, \hat\Sigma_P)} and
#'     \eqn{g_0 = N(x_c; \hat\mu_P, c\,\hat\Sigma_P)}, where
#'     \eqn{\hat\mu_P, \hat\Sigma_P} are the P-sample mean and covariance
#'     (fixed plug-in values) and \eqn{c} is `control$c_infl`. NP points far from
#'     the P covariate distribution (in Mahalanobis distance) are thereby more
#'     likely to be flagged.}
#'   \item{`"M1"`}{a rival regression of the same family with its own
#'     coefficients \eqn{\beta_{bad}} (and variance, Gaussian case); no
#'     covariate term.}
#' }
#' Because the P sample contributes only to the target component, the labels
#' of the two mixture components are identified (no label switching).
#'
#' @section Algorithm:
#' A data-augmentation Gibbs sampler alternates between (1) drawing every
#' \eqn{z_j} from its Bernoulli full conditional, (2) drawing
#' \eqn{\pi \sim Beta(a_\pi + \sum z_j, b_\pi + n_{NP} - \sum z_j)}, (3) updating
#' the target-model parameters given the P sample plus the NP units with
#' \eqn{z_j = 1}, and (4) updating the contaminating component given the units
#' with \eqn{z_j = 0}. Step (3) is family specific and follows the original
#' implementation:
#' \itemize{
#'   \item Gaussian: conjugate normal / inverse-gamma Gibbs updates.
#'   \item Binomial: Polya-Gamma augmentation (Polson, Scott and Windle, 2013),
#'     with PG(1, z) variables drawn from a `pg_terms`-term truncation of their
#'     sum-of-gammas representation plus a tail-mean correction. This is a
#'     close approximation rather than an exact PG draw.
#'   \item Poisson: a Metropolis-Hastings step with a Laplace (Gaussian)
#'     proposal centred at the conditional mode, which is found by
#'     `irls_steps` Newton iterations from the current value. Once these
#'     iterations have converged the proposal is state independent and the
#'     step is an exact independence sampler; otherwise it is an approximation.
#'     The post-burn-in acceptance rate is reported by [summary.bsm_fit()].
#' }
#' The filter score of NP unit \eqn{j} is the Rao-Blackwellised posterior
#' probability of compatibility, \eqn{s_j = E[P(z_j = 1 \mid \cdot) \mid
#' data]}, averaged over all post-burn-in iterations (thinning affects only the
#' stored parameter draws).
#'
#' @section Contamination summaries:
#' Two posterior quantities are reported: `contamination` is \eqn{1 - \pi},
#' the population-level probability that an NP unit is contaminated (as in the
#' original code); `contamination_realized` is the fraction of the observed NP
#' units allocated to the contaminating component at each iteration,
#' \eqn{1 - \sum_j z_j / n_{NP}}.
#'
#' \eqn{1 - \pi} is the weight of the contaminating component, which is
#' identified only to the extent that the two components can be told apart.
#' With the covariate-free contaminating densities of `"M0"`/`"Mjoint"`,
#' compatible units that fit both components about equally well are partly
#' allocated to the contaminating one, and the `"Mjoint"` covariate term also
#' allocates compatible units in the tails of the P covariate distribution to
#' it. In simulations with binomial and Poisson outcomes this led to
#' posterior means of \eqn{1 - \pi} well above the true share of
#' contaminated units, while the target coefficients remained accurate (see
#' the package vignettes). Read \eqn{1 - \pi} as a conservative measure of
#' incompatibility rather than a count of contaminated units.
#'
#' @param formula A two-sided model formula, e.g. `y ~ x1 + x2`. Offsets are
#'   not supported.
#' @param prob_data Data frame with the probability (trusted) sample.
#' @param nonprob_data Data frame with the non-probability sample. Must contain
#'   all variables used in `formula`.
#' @param family One of `"gaussian"`, `"binomial"` (logistic regression for a
#'   0/1, logical or two-level factor response) or `"poisson"`. A `family`
#'   object such as `binomial()` is also accepted (only its name is used; the
#'   links are fixed to identity, logit and log respectively).
#' @param model Contamination model: `"Mjoint"` (default), `"M0"` or `"M1"`;
#'   see the Model section.
#' @param leverage Names of the design-matrix columns entering the covariate
#'   density of `"Mjoint"`. `NULL` (default) selects all non-constant columns
#'   taking more than two distinct values in the P sample (i.e. continuous
#'   predictors; intercept and dummy variables are excluded).
#' @param n_iter Total number of MCMC iterations, including burn-in.
#' @param burnin Number of initial iterations discarded.
#' @param thin Keep every `thin`-th post-burn-in draw.
#' @param seed Optional random seed. When supplied the fit is reproducible and
#'   the global random-number stream is restored afterwards; when `NULL` the
#'   current stream is used.
#' @param prior Prior hyperparameters, see [mixture_prior()].
#' @param control Algorithm settings, see [mixture_control()].
#' @param na_action `"omit"` (default) drops rows with missing values in any
#'   model variable, separately within each sample, and records how many were
#'   dropped; `"fail"` stops with an error if any are present.
#' @param verbose If `TRUE`, report progress every 10% of iterations.
#'
#' @return An object of class `bsm_fit`, a list including
#'   \item{draws}{matrix of retained posterior draws of the target-model
#'     coefficients (plus `sigma2` for Gaussian), `contamination` and
#'     `contamination_realized`;}
#'   \item{bad_draws}{draws of the contaminating-component parameters;}
#'   \item{filter_score}{posterior compatibility probabilities \eqn{s_j} of
#'     the NP observations used in the fit;}
#'   \item{indicator_entropy}{binary entropy (nats) of each \eqn{s_j};}
#'   \item{nonprob_rows}{row names of `nonprob_data` used in the fit (after
#'     removing missing values), aligned with `filter_score`;}
#' and settings needed by the methods [summary.bsm_fit()],
#' [predict.bsm_fit()], [plot.bsm_fit()], [coef.bsm_fit()] and
#' [as.matrix.bsm_fit()].
#'
#' @references
#' Polson, N. G., Scott, J. G. and Windle, J. (2013). Bayesian inference for
#' logistic models using Polya-Gamma latent variables. *Journal of the
#' American Statistical Association*, 108(504), 1339-1349.
#'
#' Diebolt, J. and Robert, C. P. (1994). Estimation of finite mixture
#' distributions through Bayesian sampling. *Journal of the Royal Statistical
#' Society: Series B*, 56(2), 363-375.
#'
#' @seealso [baseline_fit()] for probability-only and naively pooled
#'   reference fits, [overlap_diagnostics()] for covariate-overlap checks,
#'   [simulate_mixture_data()] for example data.
#' @examples
#' d <- simulate_mixture_data("gaussian", n_prob = 100, n_nonprob = 100,
#'                            contamination = 0.3, seed = 1)
#' fit <- mixture_fit(y ~ x1, d$prob, d$nonprob, family = "gaussian",
#'                    n_iter = 600, burnin = 200, seed = 42)
#' fit
#' summary(fit)
#' head(filter_scores(fit))
#' @export
mixture_fit <- function(formula, prob_data, nonprob_data,
                        family = c("gaussian", "binomial", "poisson"),
                        model = c("Mjoint", "M0", "M1"), leverage = NULL,
                        n_iter = 3000, burnin = 1500, thin = 1, seed = NULL,
                        prior = mixture_prior(), control = mixture_control(),
                        na_action = c("omit", "fail"), verbose = FALSE) {
  if (is.character(family) && length(family) > 1L) family <- family[1L]
  family <- .match_family(family)
  model <- match.arg(model)
  na_action <- match.arg(na_action)
  des <- .build_design(formula, prob_data, nonprob_data, na_action)

  lev <- NULL
  if (model == "Mjoint") {
    if (is.null(leverage)) {
      idx <- .auto_leverage_cols(des$Xp)
      if (length(idx) == 0L) {
        stop("No continuous predictors found for the covariate (leverage) density ",
             "of model = \"Mjoint\". Use model = \"M0\" or specify `leverage`.",
             call. = FALSE)
      }
      lev <- colnames(des$Xp)[idx]
    } else {
      if (!is.character(leverage)) stop("`leverage` must be a character vector of column names.", call. = FALSE)
      miss <- setdiff(leverage, colnames(des$Xp))
      if (length(miss)) {
        stop("`leverage` columns not in the design matrix: ", paste(miss, collapse = ", "),
             ". Available: ", paste(colnames(des$Xp), collapse = ", "), call. = FALSE)
      }
      lev <- leverage
    }
  }

  fit <- mixture_fit_matrix(des$Xp, des$yp, des$Xnp, des$ynp, family = family,
                            model = model, leverage = lev, n_iter = n_iter,
                            burnin = burnin, thin = thin, seed = seed,
                            prior = prior, control = control, verbose = verbose)
  fit$call <- match.call()
  fit$formula <- des$formula
  fit$terms <- des$terms
  fit$xlevels <- des$xlevels
  fit$contrasts <- des$contrasts
  fit$nonprob_rows <- des$np_rows
  fit$na_info <- des$na_info
  fit
}

## Build aligned design matrices for the P and NP samples from a formula.
## Factor levels are taken from the union of both samples so that both design
## matrices have identical columns.
.build_design <- function(formula, prob_data, nonprob_data = NULL, na_action = "omit") {
  single <- is.null(nonprob_data)
  if (single) nonprob_data <- prob_data[0, , drop = FALSE]
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("`formula` must be a two-sided formula such as y ~ x1 + x2.", call. = FALSE)
  }
  if (!is.data.frame(prob_data)) stop("`prob_data` must be a data frame.", call. = FALSE)
  if (!is.data.frame(nonprob_data)) stop("`nonprob_data` must be a data frame.", call. = FALSE)
  if (nrow(prob_data) == 0L) stop("`prob_data` has no rows.", call. = FALSE)
  if (!single && nrow(nonprob_data) == 0L) stop("`nonprob_data` has no rows.", call. = FALSE)
  tt0 <- stats::terms(formula, data = prob_data)
  if (!is.null(attr(tt0, "offset"))) stop("Offsets are not supported.", call. = FALSE)
  formula <- stats::formula(tt0)
  environment(formula) <- environment(tt0)
  vars <- all.vars(formula)
  for (nm in c("prob_data", "nonprob_data")) {
    d <- get(nm)
    miss <- setdiff(vars, names(d))
    if (length(miss)) {
      stop(sprintf("Variables not found in `%s`: %s", nm, paste(miss, collapse = ", ")),
           call. = FALSE)
    }
  }
  pdat <- as.data.frame(prob_data)[, vars, drop = FALSE]
  npdat <- as.data.frame(nonprob_data)[, vars, drop = FALSE]
  np_rows_all <- rownames(nonprob_data)
  cc_p <- stats::complete.cases(pdat)
  cc_np <- stats::complete.cases(npdat)
  na_info <- c(prob = sum(!cc_p), nonprob = sum(!cc_np))
  if (any(na_info > 0)) {
    if (na_action == "fail") {
      stop(sprintf("Missing values in model variables: %d row(s) in the probability/first data set, %d in nonprob_data.",
                   na_info[1], na_info[2]), call. = FALSE)
    }
    if (single) {
      message(sprintf("Dropping %d row(s) with missing values.", na_info[1]))
    } else {
      message(sprintf("Dropping rows with missing values: %d from prob_data, %d from nonprob_data.",
                      na_info[1], na_info[2]))
    }
  }
  pdat <- pdat[cc_p, , drop = FALSE]
  npdat <- npdat[cc_np, , drop = FALSE]
  if (nrow(pdat) == 0L) stop("No complete rows in `prob_data`.", call. = FALSE)
  if (!single && nrow(npdat) == 0L) stop("No complete rows in `nonprob_data`.", call. = FALSE)
  for (v in vars) {
    if (is.character(pdat[[v]]) || is.character(npdat[[v]])) {
      lv <- sort(unique(c(as.character(pdat[[v]]), as.character(npdat[[v]]))))
      pdat[[v]] <- factor(pdat[[v]], levels = lv)
      npdat[[v]] <- factor(npdat[[v]], levels = lv)
    }
  }
  combined <- rbind(pdat, npdat)
  rownames(combined) <- NULL
  mf <- stats::model.frame(formula, data = combined, na.action = stats::na.pass)
  tt <- attr(mf, "terms")
  X <- stats::model.matrix(tt, mf)
  y <- stats::model.response(mf)
  if (is.matrix(y)) stop("Matrix responses (e.g. cbind(successes, failures)) are not supported.", call. = FALSE)
  ip <- seq_len(nrow(pdat))
  inp <- nrow(pdat) + seq_len(nrow(npdat))
  list(Xp = X[ip, , drop = FALSE], yp = y[ip],
       Xnp = X[inp, , drop = FALSE], ynp = y[inp],
       formula = formula, terms = tt, xlevels = stats::.getXlevels(tt, mf),
       contrasts = attr(X, "contrasts"),
       np_rows = np_rows_all[cc_np], na_info = na_info)
}
