## Internal helpers: linear algebra, numerics, validation, RNG handling.
## None of these are exported.

## ---- multivariate normal --------------------------------------------------

## One draw from N(mu, Sigma). Kept bit-for-bit identical to the original
## scripts (upper Cholesky factor, one call to rnorm(p)) so that seeded runs
## reproduce the original samplers exactly.
.rmvnorm <- function(mu, Sigma) {
  L <- chol(Sigma)
  mu + as.vector(t(L) %*% stats::rnorm(length(mu)))
}

## Log density of each row of X under N(mu, Sigma).
.dmvnorm_log <- function(X, mu, Sigma) {
  p <- ncol(X)
  L <- chol(Sigma)
  Xc <- sweep(X, 2, mu, "-")
  z <- backsolve(L, t(Xc), transpose = TRUE)
  -0.5 * (p * log(2 * pi) + 2 * sum(log(diag(L))) + colSums(z^2))
}

## Log density of a single vector x under N(mu, Sigma).
.dmvnorm_log1 <- function(x, mu, Sigma) {
  L <- chol(Sigma)
  z <- backsolve(L, x - mu, transpose = TRUE)
  -0.5 * (length(x) * log(2 * pi) + 2 * sum(log(diag(L))) + sum(z^2))
}

.dnorm_log <- function(y, mean, var) {
  -0.5 * (log(2 * pi * var) + (y - mean)^2 / var)
}

.sigmoid <- function(x) 1 / (1 + exp(-x))

.clip <- function(x, lo, hi) pmin(pmax(x, lo), hi)

## Binary entropy (nats) of posterior inclusion probabilities.
.binary_entropy <- function(s) {
  -(ifelse(s > 0, s * log(s), 0) + ifelse(s < 1, (1 - s) * log(1 - s), 0))
}

## Rank-based (Mann-Whitney) AUC. Mathematically identical to
## mean(outer(pos, neg, ">")) + 0.5 * mean(outer(pos, neg, "==")) used in the
## original code, but O(n log n) in time and O(n) in memory.
.auc_rank <- function(pos, neg) {
  n1 <- length(pos)
  n0 <- length(neg)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  rk <- rank(c(pos, neg))
  (sum(rk[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

## ---- reproducible RNG without permanently altering the user's stream -------

## Evaluate `expr` after set.seed(seed); restore the caller's RNG state on exit.
## With seed = NULL the current stream is used and advanced as usual.
.with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  genv <- globalenv()
  had_seed <- exists(".Random.seed", envir = genv, inherits = FALSE)
  if (had_seed) old <- get(".Random.seed", envir = genv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old, envir = genv)
    } else if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
      rm(".Random.seed", envir = genv)
    }
  }, add = TRUE)
  set.seed(seed)
  expr
}

## ---- argument validation ---------------------------------------------------

.check_count <- function(x, name, min = 1L) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != round(x) || x < min) {
    stop(sprintf("`%s` must be a single integer >= %d.", name, min), call. = FALSE)
  }
  as.integer(x)
}

.check_positive <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0) {
    stop(sprintf("`%s` must be a single positive finite number.", name), call. = FALSE)
  }
  invisible(x)
}

.check_mcmc <- function(n_iter, burnin, thin) {
  n_iter <- .check_count(n_iter, "n_iter", 2L)
  burnin <- .check_count(burnin, "burnin", 0L)
  thin <- .check_count(thin, "thin", 1L)
  if (burnin >= n_iter) {
    stop("`burnin` must be smaller than `n_iter` (n_iter is the total number of ",
         "iterations, including burn-in).", call. = FALSE)
  }
  if ((n_iter - burnin) %/% thin < 2L) {
    stop("Fewer than 2 draws would be retained; increase `n_iter` or reduce ",
         "`burnin`/`thin`.", call. = FALSE)
  }
  list(n_iter = n_iter, burnin = burnin, thin = thin)
}

.check_seed <- function(seed) {
  if (is.null(seed)) return(NULL)
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || !is.finite(seed)) {
    stop("`seed` must be NULL or a single finite number.", call. = FALSE)
  }
  seed
}

## Validate outcome values for a family; returns numeric vector.
.check_response <- function(y, family, label) {
  if (is.factor(y)) {
    if (family != "binomial") {
      stop(sprintf("Response in %s is a factor; only allowed for family = \"binomial\".",
                   label), call. = FALSE)
    }
    if (nlevels(y) != 2L) {
      stop(sprintf("Factor response in %s must have exactly two levels.", label),
           call. = FALSE)
    }
    y <- as.numeric(y == levels(y)[2L])
  } else if (is.logical(y)) {
    y <- as.numeric(y)
  }
  if (!is.numeric(y) || is.matrix(y) && ncol(y) != 1L) {
    stop(sprintf("Response in %s must be a numeric vector.", label), call. = FALSE)
  }
  y <- as.numeric(y)
  if (any(!is.finite(y))) {
    stop(sprintf("Response in %s contains non-finite values.", label), call. = FALSE)
  }
  if (family == "binomial" && any(!(y %in% c(0, 1)))) {
    stop(sprintf("family = \"binomial\" requires a 0/1 response; %s contains other values.",
                 label), call. = FALSE)
  }
  if (family == "poisson" && (any(y < 0) || any(y != round(y)))) {
    stop(sprintf("family = \"poisson\" requires non-negative integer counts; %s contains other values.",
                 label), call. = FALSE)
  }
  y
}

.match_family <- function(family) {
  if (inherits(family, "family")) family <- family$family
  if (!is.character(family) || length(family) != 1L) {
    stop("`family` must be one of \"gaussian\", \"binomial\", \"poisson\".", call. = FALSE)
  }
  family <- tolower(family)
  if (family == "logistic") family <- "binomial"
  if (!family %in% c("gaussian", "binomial", "poisson")) {
    stop("`family` must be one of \"gaussian\", \"binomial\", \"poisson\".", call. = FALSE)
  }
  family
}

## Columns of a design matrix that are "continuous" (not constant, not 0/1).
.auto_leverage_cols <- function(X) {
  which(apply(X, 2, function(v) {
    u <- unique(v)
    length(u) > 2L && !all(u %in% c(0, 1))
  }))
}

## Ensure a design matrix has non-empty, unique column names.
.fix_colnames <- function(X, fallback = NULL) {
  cn <- colnames(X)
  if (is.null(cn)) cn <- fallback
  if (is.null(cn)) cn <- rep("", ncol(X))
  empty <- is.na(cn) | cn == ""
  cn[empty] <- paste0("x", which(empty))
  colnames(X) <- make.unique(cn)
  X
}
