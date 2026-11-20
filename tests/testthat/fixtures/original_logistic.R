## Verbatim function definitions from the original mixture_filter_glm.R (simulation driver and demo code removed).
## mixture_filter_glm.R
## ----------------------------------------------------------------------------
## Logistic-GLM extension of the fully Bayesian contamination-mixture filter.
##
##   Good outcome density:  y_j ~ Bernoulli(sigmoid(x_j' beta))
##   Bad  outcome density:  M0  y_j ~ Bernoulli(q)                (max-ent base rate)
##                          M1  y_j ~ Bernoulli(sigmoid(x_j' beta_bad))
##   Covariate (leverage) term for Mjoint is IDENTICAL to the Gaussian case
##   (Proposition 1 is outcome-agnostic).
##
## Conjugacy is restored with Polya-Gamma augmentation (Polson, Scott & Windle 2013):
##   omega_i ~ PG(1, x_i' beta) makes beta | omega Gaussian. PG(1,z) is drawn by its
##   truncated convolution representation (with a deterministic tail-mean correction).
##
## Line-for-line port of glm_mixture.py. Base R only.
## ----------------------------------------------------------------------------

.rmvnorm <- function(mu, Sigma) { L <- chol(Sigma); mu + as.vector(t(L) %*% rnorm(length(mu))) }
.dmvnorm_log <- function(X, mu, Sigma) {
  p <- ncol(X); L <- chol(Sigma); Xc <- sweep(X, 2, mu, "-")
  z <- backsolve(L, t(Xc), transpose = TRUE)
  -0.5 * (p * log(2 * pi) + 2 * sum(log(diag(L))) + colSums(z^2))
}
.sigmoid <- function(x) 1 / (1 + exp(-x))

## Sample PG(1, z) elementwise via truncated convolution (K terms) + tail correction.
rpg1 <- function(z, K = 128) {
  n <- length(z)
  ck <- (seq_len(K) - 0.5)^2                        # length K
  denom <- outer((z / (2 * pi))^2, ck, "+")         # n x K
  g <- matrix(rexp(n * K), n, K)
  rowSums(g / denom) / (2 * pi^2) + (1 / (2 * pi^2)) * (1 / K)
}

## Bayesian logistic regression via PG Gibbs (baseline estimator).
fit_bayes_logit <- function(X, y, iters = 2000, burn = 1000, B0_scale = 25, seed = 0, K = 128) {
  set.seed(seed); X <- as.matrix(X); y <- as.numeric(y)
  n <- nrow(X); p <- ncol(X); B0i <- diag(1 / B0_scale, p); kappa <- y - 0.5
  beta <- rep(0, p); draws <- matrix(NA_real_, iters - burn, p)
  for (it in seq_len(iters)) {
    omega <- rpg1(as.vector(X %*% beta), K)
    V <- solve(B0i + crossprod(X * omega, X))
    beta <- .rmvnorm(as.vector(V %*% crossprod(X, kappa)), V)
    if (it > burn) draws[it - burn, ] <- beta
  }
  list(beta_mean = colMeans(draws), beta_sd = apply(draws, 2, sd))
}

## Logistic contamination-mixture Gibbs sampler.
fit_logit_mixture <- function(Xp, yp, Xnp, ynp, cont_idx, model = "M0",
                              iters = 3000, burn = 1500, c_infl = 2.5, B0_scale = 25,
                              pi_a = 1, pi_b = 1, q_a = 1, q_b = 1, seed = 0, K = 128) {
  set.seed(seed)
  Xp <- as.matrix(Xp); yp <- as.numeric(yp); Xnp <- as.matrix(Xnp); ynp <- as.numeric(ynp)
  p <- ncol(Xp); nNP <- nrow(Xnp); B0i <- diag(1 / B0_scale, p)
  Xp_c <- Xp[, cont_idx, drop = FALSE]; Xnp_c <- Xnp[, cont_idx, drop = FALSE]
  muP <- colMeans(Xp_c); SigP <- cov(Xp_c)
  if (model == "Mjoint") {
    logcov_good <- .dmvnorm_log(Xnp_c, muP, SigP)
    logcov_bad  <- .dmvnorm_log(Xnp_c, muP, c_infl * SigP)
  } else { logcov_good <- rep(0, nNP); logcov_bad <- rep(0, nNP) }

  bern_ll <- function(y, pr) { pr <- pmin(pmax(pr, 1e-9), 1 - 1e-9); y * log(pr) + (1 - y) * log(1 - pr) }
  beta <- rep(0, p); pi <- 0.7; q <- min(max(mean(ynp), 0.05), 0.95); beta_bad <- rep(0, p)
  keep <- matrix(NA_real_, iters - burn, p + 2); resp <- numeric(nNP); nk <- 0

  for (it in seq_len(iters)) {
    pg <- .sigmoid(as.vector(Xnp %*% beta))
    log_good <- bern_ll(ynp, pg) + logcov_good + log(pi)
    if (model == "M1") {
      log_bad <- bern_ll(ynp, .sigmoid(as.vector(Xnp %*% beta_bad))) + logcov_bad + log(1 - pi)
    } else {
      log_bad <- bern_ll(ynp, q) + logcov_bad + log(1 - pi)
    }
    m <- pmax(log_good, log_bad); a <- exp(log_good - m); b <- exp(log_bad - m); r <- a / (a + b)
    z <- as.numeric(runif(nNP) < r)
    pi <- rbeta(1, pi_a + sum(z), pi_b + nNP - sum(z))

    good <- z > 0.5
    Xg <- rbind(Xp, Xnp[good, , drop = FALSE]); yg <- c(yp, ynp[good]); kappa <- yg - 0.5
    omega <- rpg1(as.vector(Xg %*% beta), K)
    V <- solve(B0i + crossprod(Xg * omega, Xg))
    beta <- .rmvnorm(as.vector(V %*% crossprod(Xg, kappa)), V)

    bad <- !good; nb <- sum(bad)
    if (model == "M1") {
      if (nb > p + 1) {
        Xb <- Xnp[bad, , drop = FALSE]; yb <- ynp[bad]; kb <- yb - 0.5
        ob <- rpg1(as.vector(Xb %*% beta_bad), K)
        Vb <- solve(B0i + crossprod(Xb * ob, Xb))
        beta_bad <- .rmvnorm(as.vector(Vb %*% crossprod(Xb, kb)), Vb)
      } else beta_bad <- .rmvnorm(rep(0, p), diag(B0_scale, p))
    } else {
      if (nb > 0) { yb <- ynp[bad]; q <- rbeta(1, q_a + sum(yb), q_b + nb - sum(yb)) }
      else q <- rbeta(1, q_a, q_b)
    }
    if (it > burn) { keep[it - burn, ] <- c(beta, pi, 1 - pi); resp <- resp + r; nk <- nk + 1 }
  }
  filt <- resp / nk
  ent <- -(ifelse(filt > 0, filt * log(filt), 0) + ifelse(filt < 1, (1 - filt) * log(1 - filt), 0))
  list(beta_mean = colMeans(keep[, 1:p, drop = FALSE]),
       beta_sd = apply(keep[, 1:p, drop = FALSE], 2, sd),
       pi_mean = mean(keep[, p + 1]), contam_mean = mean(keep[, p + 2]),
       contam_draws = keep[, p + 2], filter_score = filt, indicator_entropy = ent)
}

