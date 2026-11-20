## ---- MCMC diagnostics --------------------------------------------------------

## Autocorrelations via FFT (biased estimator, as is standard for ESS).
.acf_fft <- function(x) {
  n <- length(x)
  xc <- x - mean(x)
  m <- stats::nextn(2L * n)
  f <- stats::fft(c(xc, rep(0, m - n)))
  ac <- Re(stats::fft(Mod(f)^2, inverse = TRUE))[seq_len(n)] / m
  ac / ac[1L]
}

#' Effective sample size of MCMC draws
#'
#' Computes the effective sample size of each column of a matrix of draws
#' from a single chain, using Geyer's (1992) initial monotone sequence
#' estimator of the integrated autocorrelation time.
#'
#' @param x Numeric vector or matrix (one column per parameter).
#' @return Named numeric vector of effective sample sizes. Constant columns
#'   return `NA`.
#' @references Geyer, C. J. (1992). Practical Markov chain Monte Carlo.
#'   *Statistical Science*, 7(4), 473-483.
#' @examples
#' set.seed(1)
#' x <- as.numeric(arima.sim(list(ar = 0.5), n = 2000))
#' effective_size(x)   # theory: about 2000 * (1 - 0.5) / (1 + 0.5) = 667
#' @export
effective_size <- function(x) {
  x <- as.matrix(x)
  out <- apply(x, 2, function(v) {
    n <- length(v)
    if (n < 4L || !all(is.finite(v)) || stats::var(v) == 0) return(NA_real_)
    rho <- .acf_fft(v)
    n_pairs <- (n - 1L) %/% 2L
    gam <- rho[2L * seq_len(n_pairs) - 1L] + rho[2L * seq_len(n_pairs)]  # rho_{2k} + rho_{2k+1}, k = 0, 1, ...
    pos <- which(gam <= 0)
    K <- if (length(pos)) pos[1L] - 1L else length(gam)
    if (K < 1L) return(as.numeric(n))
    gam <- cummin(gam[seq_len(K)])
    tau <- -1 + 2 * sum(gam)
    n / max(tau, 1 / log10(max(n, 10)))
  })
  names(out) <- colnames(x)
  out
}

## ---- overlap diagnostics (original stage0) -----------------------------------

#' Covariate-overlap diagnostics between the two samples
#'
#' Pre-fitting ("Stage 0") diagnostics of how well the non-probability sample
#' overlaps the probability sample in the space of continuous covariates.
#' This is the `stage0()` function of the original code.
#'
#' Mahalanobis distances of all units are computed with respect to the
#' probability-sample mean and covariance. Reported are
#' \describe{
#'   \item{`OC`}{an overlap coefficient: the area shared by kernel density
#'     estimates of the P and NP distance distributions (on a 400-point
#'     grid), relative to the area of the P density. Values near 1 indicate
#'     good overlap.}
#'   \item{`EF`}{the exceedance fraction: share of NP units whose distance
#'     exceeds the 99th percentile of the P distances (about 0.01 under
#'     perfect overlap).}
#'   \item{`AUC_PS`}{the area under the ROC curve of a logistic propensity
#'     model for sample membership (NP vs P) on the covariates; 0.5 indicates
#'     indistinguishable samples, values near 1 strong separation.}
#' }
#'
#' @param prob,nonprob Numeric matrices or data frames of the continuous
#'   covariates for the P and NP samples (with matching column names).
#' @param vars Optional character vector of column names to use; by default
#'   all columns (which must be numeric).
#' @return A named list with elements `OC`, `EF` and `AUC_PS`, of class
#'   `bsm_overlap`.
#' @examples
#' set.seed(3)
#' P <- data.frame(x1 = rnorm(200), x2 = rnorm(200))
#' NP <- data.frame(x1 = rnorm(150, 0.5), x2 = rnorm(150))
#' overlap_diagnostics(P, NP)
#' @export
overlap_diagnostics <- function(prob, nonprob, vars = NULL) {
  prob <- as.data.frame(prob)
  nonprob <- as.data.frame(nonprob)
  if (is.null(vars)) vars <- names(prob)
  miss <- setdiff(vars, c(intersect(names(prob), names(nonprob))))
  if (length(miss)) stop("Columns missing from one of the samples: ", paste(miss, collapse = ", "), call. = FALSE)
  Xpc <- as.matrix(prob[, vars, drop = FALSE])
  Xnc <- as.matrix(nonprob[, vars, drop = FALSE])
  if (!is.numeric(Xpc) || !is.numeric(Xnc)) stop("All selected columns must be numeric.", call. = FALSE)
  if (anyNA(Xpc) || anyNA(Xnc)) stop("Missing values are not allowed.", call. = FALSE)
  if (nrow(Xpc) <= ncol(Xpc)) stop("Need more probability-sample rows than covariates.", call. = FALSE)
  mu <- colMeans(Xpc)
  S <- stats::cov(Xpc)
  Si <- tryCatch(solve(S), error = function(e) {
    stop("The probability-sample covariance matrix is singular.", call. = FALSE)
  })
  dM <- function(M) {
    Mc <- sweep(M, 2, mu, "-")
    sqrt(rowSums((Mc %*% Si) * Mc))
  }
  dP <- dM(Xpc)
  dN <- dM(Xnc)
  q99 <- stats::quantile(dP, 0.99)
  EF <- mean(dN > q99)
  lo <- min(c(dP, dN))
  hi <- max(c(dP, dN))
  fp <- stats::density(dP, from = lo, to = hi, n = 400)$y
  fn <- stats::density(dN, from = lo, to = hi, n = 400)$y
  OC <- sum(pmin(fp, fn)) / sum(fp)
  memb <- c(rep(0, nrow(Xpc)), rep(1, nrow(Xnc)))
  Xall <- rbind(Xpc, Xnc)
  gl <- suppressWarnings(stats::glm(memb ~ Xall, family = stats::binomial()))
  ph <- stats::predict(gl, type = "response")
  auc <- .auc_rank(ph[memb == 1], ph[memb == 0])
  structure(list(OC = as.numeric(OC), EF = as.numeric(EF), AUC_PS = as.numeric(auc)),
            class = "bsm_overlap")
}

#' @export
print.bsm_overlap <- function(x, digits = 3, ...) {
  cat("Covariate overlap diagnostics (P vs NP)\n")
  cat(sprintf("  Overlap coefficient (OC):       %s\n", format(round(x$OC, digits), nsmall = digits)))
  cat(sprintf("  Exceedance fraction (EF, >q99): %s\n", format(round(x$EF, digits), nsmall = digits)))
  cat(sprintf("  Propensity-score AUC:           %s\n", format(round(x$AUC_PS, digits), nsmall = digits)))
  invisible(x)
}

#' Area under the ROC curve of filter scores against known labels
#'
#' For simulation studies in which the true compatibility labels are known,
#' measures how well the filter scores rank compatible (label 1) above
#' contaminated (label 0) units. Ties count one half.
#'
#' @param labels Vector of true labels (1 = compatible, 0 = contaminated).
#' @param scores Numeric scores, e.g. [filter_scores()] of a fit.
#' @return A single number in \[0, 1\], or `NA` if either class is absent.
#' @examples
#' filter_auc(c(1, 1, 0, 0), c(0.9, 0.7, 0.8, 0.1))
#' @export
filter_auc <- function(labels, scores) {
  if (length(labels) != length(scores)) stop("`labels` and `scores` must have equal length.", call. = FALSE)
  .auc_rank(scores[labels == 1], scores[labels == 0])
}
