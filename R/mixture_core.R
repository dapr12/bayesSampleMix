## Shared machinery for all families: covariate (leverage) densities, the
## generic data-augmentation Gibbs driver, and assembly of the fitted object.

.make_engine <- function(family, ...) {
  switch(family,
         gaussian = .engine_gaussian(...),
         binomial = .engine_binomial(...),
         poisson = .engine_poisson(...))
}

## Plug-in covariate log densities for model "Mjoint".
## good:  x_c ~ N(mu_P, Sigma_P);  contaminated:  x_c ~ N(mu_P, c * Sigma_P),
## with (mu_P, Sigma_P) the probability-sample mean and covariance.
.leverage_logdens <- function(Xp, Xnp, idx, c_infl) {
  Xp_c <- Xp[, idx, drop = FALSE]
  Xnp_c <- Xnp[, idx, drop = FALSE]
  if (nrow(Xp_c) <= ncol(Xp_c)) {
    stop("model = \"Mjoint\" needs more probability-sample rows than leverage ",
         "covariates to estimate their covariance matrix.", call. = FALSE)
  }
  muP <- colMeans(Xp_c)
  SigP <- stats::cov(Xp_c)
  ok <- tryCatch({
    chol(SigP)
    TRUE
  }, error = function(e) FALSE)
  if (!ok || rcond(SigP) < 1e-10) {
    stop("The probability-sample covariance matrix of the leverage covariates (",
         paste(colnames(Xp)[idx], collapse = ", "), ") is singular or nearly ",
         "singular. Remove collinear/constant leverage covariates or use model = \"M0\".",
         call. = FALSE)
  }
  list(good = .dmvnorm_log(Xnp_c, muP, SigP),
       bad = .dmvnorm_log(Xnp_c, muP, c_infl * SigP),
       mu = muP, Sigma = SigP)
}

## Generic data-augmentation sampler. The order of random-number calls per
## iteration is: z (runif), pi (rbeta), target-model block, contaminating
## block; identical to the original scripts for all three families.
.run_mixture <- function(engine, nNP, logcov_good, logcov_bad, prior, mcmc,
                         verbose = FALSE) {
  n_iter <- mcmc$n_iter
  burnin <- mcmc$burnin
  thin <- mcmc$thin
  n_keep <- (n_iter - burnin) %/% thin
  good_names <- engine$good_names
  bad_names <- engine$bad_names
  n_good_par <- length(good_names)
  n_bad_par <- length(bad_names)

  draws <- matrix(NA_real_, n_keep, n_good_par + 2L,
                  dimnames = list(NULL, c(good_names, "contamination", "contamination_realized")))
  bad_draws <- matrix(NA_real_, n_keep, n_bad_par, dimnames = list(NULL, bad_names))
  pi_a <- prior$pi_shape1
  pi_b <- prior$pi_shape2

  pi_w <- 0.7
  resp_sum <- numeric(nNP)
  n_post <- 0L
  k <- 0L
  acc_at_burn <- NULL
  if (burnin == 0L) acc_at_burn <- engine$accept()
  report_every <- max(1L, n_iter %/% 10L)

  for (it in seq_len(n_iter)) {
    log_good <- engine$loglik_good() + logcov_good + log(pi_w)
    log_bad <- engine$loglik_bad() + logcov_bad + log(1 - pi_w)
    m <- pmax(log_good, log_bad)
    a <- exp(log_good - m)
    b <- exp(log_bad - m)
    r <- a / (a + b)
    if (anyNA(r)) {
      stop("Numerical failure while computing allocation probabilities at iteration ",
           it, ". Consider rescaling covariates or using stronger priors.", call. = FALSE)
    }
    z <- as.numeric(stats::runif(nNP) < r)
    n_good_np <- sum(z)
    pi_w <- stats::rbeta(1, pi_a + n_good_np, pi_b + nNP - n_good_np)

    good <- z > 0.5
    engine$update_good(good)
    engine$update_bad(!good)

    if (it == burnin) acc_at_burn <- engine$accept()
    if (it > burnin) {
      resp_sum <- resp_sum + r
      n_post <- n_post + 1L
      if ((it - burnin) %% thin == 0L) {
        k <- k + 1L
        draws[k, ] <- c(engine$good_params(), 1 - pi_w, 1 - n_good_np / nNP)
        bad_draws[k, ] <- engine$bad_params()
      }
    }
    if (verbose && it %% report_every == 0L) {
      message(sprintf("iteration %d / %d", it, n_iter))
    }
  }
  acc_end <- engine$accept()
  acc_rate <- if (is.null(acc_end)) NA_real_ else {
    d <- acc_end - acc_at_burn
    if (d[2] > 0) d[1] / d[2] else NA_real_
  }
  list(draws = draws, bad_draws = bad_draws, filter_score = resp_sum / n_post,
       acceptance = acc_rate)
}

## Single-sample Bayesian regression with the same target-model updates
## (used by baseline_fit()).
.run_baseline <- function(engine, mcmc) {
  n_iter <- mcmc$n_iter
  burnin <- mcmc$burnin
  thin <- mcmc$thin
  n_keep <- (n_iter - burnin) %/% thin
  draws <- matrix(NA_real_, n_keep, length(engine$good_names),
                  dimnames = list(NULL, engine$good_names))
  k <- 0L
  acc_at_burn <- if (burnin == 0L) engine$accept() else NULL
  for (it in seq_len(n_iter)) {
    engine$update_good(logical(0))
    if (it == burnin) acc_at_burn <- engine$accept()
    if (it > burnin && (it - burnin) %% thin == 0L) {
      k <- k + 1L
      draws[k, ] <- engine$good_params()
    }
  }
  acc_end <- engine$accept()
  acc_rate <- if (is.null(acc_end)) NA_real_ else {
    d <- acc_end - acc_at_burn
    if (d[2] > 0) d[1] / d[2] else NA_real_
  }
  list(draws = draws, acceptance = acc_rate)
}

.check_design <- function(X, label) {
  if (!is.matrix(X) || !is.numeric(X)) {
    stop(sprintf("`%s` must be a numeric matrix.", label), call. = FALSE)
  }
  if (any(!is.finite(X))) {
    stop(sprintf("`%s` contains missing or non-finite values.", label), call. = FALSE)
  }
  invisible(X)
}

#' Fit the contamination-mixture model from design matrices
#'
#' Low-level interface that takes design matrices directly, mirroring the
#' original research functions `fit_mixture()`, `fit_logit_mixture()` and
#' `fit_pois_mixture()`. Most users should call [mixture_fit()], which builds
#' the design matrices from a formula and calls this function.
#'
#' @param Xp,Xnp Numeric design matrices (including an intercept column, if
#'   wanted) for the probability (P) and non-probability (NP) samples. Both
#'   must have the same columns.
#' @param yp,ynp Outcome vectors for the P and NP samples.
#' @param leverage Integer column indices (or column names) of `Xp`/`Xnp`
#'   that enter the covariate density of model `"Mjoint"`; typically the
#'   continuous predictors. Equivalent to `cont_idx` in the original code.
#'   Ignored for models `"M0"` and `"M1"`.
#' @inheritParams mixture_fit
#'
#' @return An object of class `bsm_fit`; see [mixture_fit()].
#' @examples
#' set.seed(1)
#' Xp <- cbind(1, rnorm(100)); yp <- drop(Xp %*% c(2, 5)) + rnorm(100)
#' Xnp <- cbind(1, rnorm(100)); ynp <- drop(Xnp %*% c(2, 5)) + rnorm(100)
#' ynp[1:25] <- rnorm(25, -20)          # 25% contaminated outcomes
#' fit <- mixture_fit_matrix(Xp, yp, Xnp, ynp, leverage = 2,
#'                           n_iter = 600, burnin = 200, seed = 1)
#' coef(fit)
#' @export
mixture_fit_matrix <- function(Xp, yp, Xnp, ynp, family = "gaussian",
                               model = c("Mjoint", "M0", "M1"), leverage = NULL,
                               n_iter = 3000, burnin = 1500, thin = 1,
                               seed = NULL, prior = mixture_prior(),
                               control = mixture_control(), verbose = FALSE) {
  family <- .match_family(family)
  model <- match.arg(model)
  mcmc <- .check_mcmc(n_iter, burnin, thin)
  seed <- .check_seed(seed)
  if (!inherits(prior, "bsm_prior")) stop("`prior` must be created by mixture_prior().", call. = FALSE)
  if (!inherits(control, "bsm_control")) stop("`control` must be created by mixture_control().", call. = FALSE)

  Xp <- as.matrix(Xp)
  Xnp <- as.matrix(Xnp)
  storage.mode(Xp) <- "double"
  storage.mode(Xnp) <- "double"
  .check_design(Xp, "Xp")
  .check_design(Xnp, "Xnp")
  if (ncol(Xp) != ncol(Xnp)) stop("`Xp` and `Xnp` must have the same number of columns.", call. = FALSE)
  Xp <- .fix_colnames(Xp, fallback = colnames(Xnp))
  colnames(Xnp) <- colnames(Xp)
  yp <- .check_response(yp, family, "the probability sample")
  ynp <- .check_response(ynp, family, "the non-probability sample")
  if (length(yp) != nrow(Xp)) stop("length(yp) must equal nrow(Xp).", call. = FALSE)
  if (length(ynp) != nrow(Xnp)) stop("length(ynp) must equal nrow(Xnp).", call. = FALSE)
  p <- ncol(Xp)
  if (nrow(Xnp) < 1L) stop("The non-probability sample has no observations.", call. = FALSE)
  if (nrow(Xp) < 2L) stop("The probability sample needs at least 2 observations.", call. = FALSE)
  rk <- qr(Xp)$rank
  if (rk < p) {
    msg <- paste0("The probability-sample design matrix is rank deficient (rank ", rk,
                  " < ", p, " columns). The P sample anchors the target model; some ",
                  "coefficients are identified only through the prior and the NP sample.")
    if (family == "gaussian") {
      stop(msg, " The Gaussian sampler is initialised by least squares on the P sample ",
           "and requires full rank; remove redundant terms.", call. = FALSE)
    }
    warning(msg, call. = FALSE)
  }
  if (family == "gaussian" && nrow(Xp) <= p) {
    stop("family = \"gaussian\" needs more probability-sample observations than coefficients.",
         call. = FALSE)
  }

  c_infl <- if (is.null(control$c_infl)) .default_c_infl(family) else control$c_infl
  lev_idx <- NULL
  if (model == "Mjoint") {
    if (is.null(leverage) || length(leverage) == 0L) {
      stop("model = \"Mjoint\" requires at least one leverage (continuous) covariate; ",
           "supply `leverage` or use model = \"M0\".", call. = FALSE)
    }
    lev_idx <- if (is.character(leverage)) match(leverage, colnames(Xp)) else as.integer(leverage)
    if (anyNA(lev_idx) || any(lev_idx < 1L | lev_idx > p)) {
      stop("`leverage` refers to columns not present in the design matrix.", call. = FALSE)
    }
    lev <- .leverage_logdens(Xp, Xnp, lev_idx, c_infl)
    logcov_good <- lev$good
    logcov_bad <- lev$bad
  } else {
    logcov_good <- rep(0, nrow(Xnp))
    logcov_bad <- rep(0, nrow(Xnp))
    lev <- NULL
  }

  start <- proc.time()[["elapsed"]]
  res <- .with_seed(seed, {
    engine <- .make_engine(family, Xp = Xp, yp = yp, Xnp = Xnp, ynp = ynp,
                           model = model, prior = prior, control = control)
    .run_mixture(engine, nrow(Xnp), logcov_good, logcov_bad, prior, mcmc, verbose)
  })
  elapsed <- proc.time()[["elapsed"]] - start

  structure(list(
    family = family, model = model,
    draws = res$draws, bad_draws = res$bad_draws,
    coef_names = colnames(Xp),
    filter_score = res$filter_score,
    indicator_entropy = .binary_entropy(res$filter_score),
    acceptance = res$acceptance,
    leverage = if (is.null(lev_idx)) NULL else colnames(Xp)[lev_idx],
    leverage_moments = if (is.null(lev)) NULL else list(mean = lev$mu, cov = lev$Sigma),
    c_infl = if (model == "Mjoint") c_infl else NULL,
    n = c(prob = nrow(Xp), nonprob = nrow(Xnp)),
    mcmc = mcmc, seed = seed, prior = prior, control = control,
    elapsed = elapsed, call = match.call(),
    terms = NULL, xlevels = NULL, contrasts = NULL,
    nonprob_rows = seq_len(nrow(Xnp)), na_info = NULL
  ), class = "bsm_fit")
}
