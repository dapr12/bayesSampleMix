.newdata_matrix <- function(object, newdata) {
  if (missing(newdata) || is.null(newdata)) {
    stop("`newdata` is required: supply the covariate values at which to predict.", call. = FALSE)
  }
  if (!is.null(object$terms)) {
    if (!is.data.frame(newdata)) stop("`newdata` must be a data frame for formula-based fits.", call. = FALSE)
    tt <- stats::delete.response(object$terms)
    mf <- stats::model.frame(tt, newdata, xlev = object$xlevels, na.action = stats::na.pass)
    X <- stats::model.matrix(tt, mf, contrasts.arg = object$contrasts)
  } else {
    X <- as.matrix(newdata)
    if (ncol(X) != length(object$coef_names)) {
      stop("`newdata` must have ", length(object$coef_names), " columns (the design matrix).",
           call. = FALSE)
    }
  }
  if (!identical(colnames(X), object$coef_names) && !is.null(object$terms)) {
    stop("Design matrix of `newdata` does not match the fitted model.", call. = FALSE)
  }
  X
}

.predict_draws <- function(object, newdata, type) {
  X <- .newdata_matrix(object, newdata)
  B <- object$draws[, object$coef_names, drop = FALSE]
  eta <- B %*% t(X)
  if (type == "link") return(eta)
  switch(object$family, gaussian = eta, binomial = .sigmoid(eta), poisson = exp(eta))
}

#' Posterior predictions from the target regression model
#'
#' Computes, for each row of `newdata`, the posterior distribution of the
#' target-model linear predictor (`type = "link"`) or mean response
#' (`type = "response"`: \eqn{x'\beta}, \eqn{\mathrm{logit}^{-1}(x'\beta)} or
#' \eqn{\exp(x'\beta)}), using every retained posterior draw. This is the
#' regression function of the population represented by the probability
#' sample and the compatible NP units; it does not include observation-level
#' noise.
#'
#' @param object A `bsm_fit` or `bsm_baseline` object.
#' @param newdata Data frame (formula fits) or design matrix (matrix fits)
#'   at which to predict. Required.
#' @param type `"response"` (default) or `"link"`.
#' @param interval `"credible"` (default) to add equal-tailed credible
#'   limits, or `"none"`.
#' @param level Credible-interval level.
#' @param summary If `FALSE`, return the full matrix of posterior draws
#'   (draws in rows, observations in columns) instead of summaries.
#' @param ... Not used.
#' @return A data frame with columns `fit` (posterior mean) and, if
#'   requested, `lwr` and `upr`; or a matrix of draws when `summary = FALSE`.
#' @examples
#' d <- simulate_mixture_data("binomial", n_prob = 150, n_nonprob = 150,
#'                            contamination = 0.3, seed = 7)
#' fit <- mixture_fit(y ~ x1, d$prob, d$nonprob, family = "binomial",
#'                    n_iter = 500, burnin = 200, seed = 1)
#' predict(fit, newdata = data.frame(x1 = c(-1, 0, 1)))
#' @export
predict.bsm_fit <- function(object, newdata, type = c("response", "link"),
                            interval = c("credible", "none"), level = 0.95,
                            summary = TRUE, ...) {
  type <- match.arg(type)
  interval <- match.arg(interval)
  probs <- .check_probs(level)
  dr <- .predict_draws(object, newdata, type)
  if (!summary) return(dr)
  out <- data.frame(fit = colMeans(dr))
  if (interval == "credible") {
    q <- apply(dr, 2, stats::quantile, probs = probs, names = FALSE)
    out$lwr <- q[1, ]
    out$upr <- q[2, ]
  }
  rownames(out) <- rownames(newdata)
  out
}

#' @rdname predict.bsm_fit
#' @export
predict.bsm_baseline <- predict.bsm_fit
