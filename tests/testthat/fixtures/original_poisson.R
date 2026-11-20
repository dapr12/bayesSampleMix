## Verbatim function definitions from the original mixture_filter_pois.R (demo code removed).
## mixture_filter_pois.R
## ----------------------------------------------------------------------------
## Poisson-GLM extension of the fully Bayesian contamination-mixture filter.
##
##   Good outcome density:  y_j ~ Poisson(exp(x_j' beta))
##   Bad  outcome density:  M0  y_j ~ Poisson(lambda0)            (background rate)
##                          M1  y_j ~ Poisson(exp(x_j' beta_bad))
##   Covariate (leverage) term for Mjoint is IDENTICAL to the Gaussian/logistic
##   case (Proposition 1 is outcome-agnostic).
##
## The Poisson log-link has no Polya-Gamma-style conjugate augmentation, so the
## good-set coefficient block is updated by an independence Metropolis-Hastings step
## whose proposal is the Laplace (IRLS) Gaussian approximation to the conditional
## posterior. Near-i.i.d. draws, high acceptance. Line-for-line port of
## poisson_mixture.py. Base R only.
## ----------------------------------------------------------------------------

.rmvnorm <- function(mu, S) { L <- chol(S); mu + as.vector(t(L) %*% rnorm(length(mu))) }
.dmvnorm_log <- function(X, mu, S) {
  p <- ncol(X); L <- chol(S); Xc <- sweep(X, 2, mu, "-")
  z <- backsolve(L, t(Xc), transpose = TRUE)
  -0.5 * (p * log(2 * pi) + 2 * sum(log(diag(L))) + colSums(z^2))
}
.dmvnorm_log1 <- function(x, mu, S) {
  L <- chol(S); z <- backsolve(L, x - mu, transpose = TRUE)
  -0.5 * (length(x) * log(2 * pi) + 2 * sum(log(diag(L))) + sum(z^2))
}
.clip <- function(x, lo, hi) pmin(pmax(x, lo), hi)

## penalised IRLS to the conditional mode; returns list(mode, V = inverse neg-Hessian)
.irls_mode <- function(X, y, B0i, beta0, steps = 3) {
  beta <- beta0
  for (s in seq_len(steps)) {
    mu <- exp(.clip(X %*% beta, -30, 30))[, 1]
    H <- crossprod(X * mu, X) + B0i
    grad <- crossprod(X, (y - mu)) - B0i %*% beta
    beta <- beta + solve(H, grad)
  }
  mu <- exp(.clip(X %*% beta, -30, 30))[, 1]
  list(mode = as.vector(beta), V = solve(crossprod(X * mu, X) + B0i))
}
.log_post_pois <- function(X, y, beta, B0i) {
  eta <- as.vector(X %*% beta)
  sum(y * eta - exp(.clip(eta, -30, 30))) - 0.5 * as.numeric(beta %*% B0i %*% beta)
}
## independence MH step for the Poisson coefficient block
.update_beta_pois <- function(X, y, beta, B0i, steps = 3) {
  fit <- .irls_mode(X, y, B0i, beta, steps); mode <- fit$mode; V <- fit$V
  prop <- .rmvnorm(mode, V)
  la <- .log_post_pois(X, y, prop, B0i) - .log_post_pois(X, y, beta, B0i)
  lq <- .dmvnorm_log1(beta, mode, V) - .dmvnorm_log1(prop, mode, V)
  if (log(runif(1)) < la + lq) prop else beta
}

fit_bayes_pois <- function(X, y, iters = 2000, burn = 1000, B0_scale = 25, seed = 0) {
  set.seed(seed); X <- as.matrix(X); y <- as.numeric(y)
  p <- ncol(X); B0i <- diag(1 / B0_scale, p)
  beta <- .irls_mode(X, y, B0i, rep(0, p), 5)$mode
  draws <- matrix(NA_real_, iters - burn, p)
  for (it in seq_len(iters)) {
    beta <- .update_beta_pois(X, y, beta, B0i)
    if (it > burn) draws[it - burn, ] <- beta
  }
  list(beta_mean = colMeans(draws), beta_sd = apply(draws, 2, sd))
}

fit_pois_mixture <- function(Xp, yp, Xnp, ynp, cont_idx, model = "M0",
                             iters = 3000, burn = 1500, c_infl = 2.5, B0_scale = 25,
                             pi_a = 1, pi_b = 1, lam_a = 1, lam_b = 1, seed = 0) {
  set.seed(seed)
  Xp <- as.matrix(Xp); yp <- as.numeric(yp); Xnp <- as.matrix(Xnp); ynp <- as.numeric(ynp)
  p <- ncol(Xp); nNP <- nrow(Xnp); B0i <- diag(1 / B0_scale, p)
  Xp_c <- Xp[, cont_idx, drop = FALSE]; Xnp_c <- Xnp[, cont_idx, drop = FALSE]
  muP <- colMeans(Xp_c); SigP <- cov(Xp_c)
  if (model == "Mjoint") {
    logcov_good <- .dmvnorm_log(Xnp_c, muP, SigP); logcov_bad <- .dmvnorm_log(Xnp_c, muP, c_infl * SigP)
  } else { logcov_good <- rep(0, nNP); logcov_bad <- rep(0, nNP) }
  pois_ll <- function(y, mu) { mu <- pmax(mu, 1e-9); y * log(mu) - mu }

  beta <- .irls_mode(Xp, yp, B0i, rep(0, p), 5)$mode
  pi <- 0.7; lam0 <- max(mean(ynp), 0.1); beta_bad <- beta
  keep <- matrix(NA_real_, iters - burn, p + 2); resp <- numeric(nNP); nk <- 0

  for (it in seq_len(iters)) {
    mu_g <- exp(.clip(as.vector(Xnp %*% beta), -30, 30))
    log_good <- pois_ll(ynp, mu_g) + logcov_good + log(pi)
    if (model == "M1") {
      log_bad <- pois_ll(ynp, exp(.clip(as.vector(Xnp %*% beta_bad), -30, 30))) + logcov_bad + log(1 - pi)
    } else {
      log_bad <- pois_ll(ynp, lam0) + logcov_bad + log(1 - pi)
    }
    m <- pmax(log_good, log_bad); a <- exp(log_good - m); b <- exp(log_bad - m); r <- a / (a + b)
    z <- as.numeric(runif(nNP) < r)
    pi <- rbeta(1, pi_a + sum(z), pi_b + nNP - sum(z))

    good <- z > 0.5
    Xg <- rbind(Xp, Xnp[good, , drop = FALSE]); yg <- c(yp, ynp[good])
    beta <- .update_beta_pois(Xg, yg, beta, B0i)

    bad <- !good; nb <- sum(bad)
    if (model == "M1") {
      if (nb > p + 1) beta_bad <- .update_beta_pois(Xnp[bad, , drop = FALSE], ynp[bad], beta_bad, B0i)
      else beta_bad <- .rmvnorm(rep(0, p), diag(B0_scale, p))
    } else {
      if (nb > 0) { yb <- ynp[bad]; lam0 <- rgamma(1, lam_a + sum(yb), rate = lam_b + nb) }
      else lam0 <- rgamma(1, lam_a, rate = lam_b)
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

