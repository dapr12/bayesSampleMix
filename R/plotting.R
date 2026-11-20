#' Plot a contamination-mixture fit
#'
#' Base-graphics displays of the posterior draws and filter scores.
#'
#' @param x A `bsm_fit` object.
#' @param type `"trace"` (trace plots of the retained draws), `"density"`
#'   (kernel density estimates of the marginal posteriors), `"scores"`
#'   (histogram of NP filter scores) or `"contamination"` (trace and
#'   histogram of \eqn{1 - \pi}).
#' @param pars Character vector of parameter names (columns of
#'   `as.matrix(x)`) for `"trace"` and `"density"`. Defaults to the
#'   target-model coefficients plus `contamination`.
#' @param ... Further graphical arguments passed to the underlying plotting
#'   functions.
#' @return `x`, invisibly.
#' @examples
#' d <- simulate_mixture_data("gaussian", n_prob = 100, n_nonprob = 100,
#'                            contamination = 0.3, seed = 8)
#' fit <- mixture_fit(y ~ x1, d$prob, d$nonprob, n_iter = 600, burnin = 200, seed = 1)
#' plot(fit, type = "trace")
#' plot(fit, type = "scores")
#' @export
plot.bsm_fit <- function(x, type = c("trace", "density", "scores", "contamination"),
                         pars = NULL, ...) {
  type <- match.arg(type)
  d <- x$draws
  if (is.null(pars)) pars <- c(x$coef_names, "contamination")
  all_d <- cbind(x$draws, x$bad_draws)
  bad <- setdiff(pars, colnames(all_d))
  if (length(bad)) stop("Unknown parameter(s): ", paste(bad, collapse = ", "), call. = FALSE)
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  if (type %in% c("trace", "density")) {
    k <- length(pars)
    nc <- min(2L, k)
    graphics::par(mfrow = c(ceiling(k / nc), nc), mar = c(4, 4, 2, 1))
    for (nm in pars) {
      v <- all_d[, nm]
      if (type == "trace") {
        graphics::plot(v, type = "l", xlab = "Retained draw", ylab = nm, main = nm, ...)
        graphics::abline(h = mean(v), col = "grey50", lty = 2)
      } else {
        graphics::plot(stats::density(v), main = nm, xlab = nm, ...)
      }
    }
  } else if (type == "scores") {
    graphics::hist(x$filter_score, breaks = seq(0, 1, by = 0.05), xlab = "Filter score s_j",
                   main = "Posterior P(compatible) of NP units", col = "grey80", ...)
    graphics::abline(v = 0.5, lty = 2)
  } else {
    graphics::par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
    v <- d[, "contamination"]
    graphics::plot(v, type = "l", xlab = "Retained draw", ylab = "1 - pi",
                   main = "Trace: contamination", ...)
    graphics::hist(v, breaks = 30, xlab = "1 - pi", main = "Posterior: contamination",
                   col = "grey80")
  }
  invisible(x)
}
