## Binomial (logit-link) engine.
##
## Target outcome model:   y ~ Bernoulli(sigmoid(x' beta))
## Contaminating outcome model:
##   M0 / Mjoint:  y ~ Bernoulli(q),  q ~ Beta(q_a, q_b)
##   M1:           y ~ Bernoulli(sigmoid(x' beta_bad))
## Prior: beta, beta_bad ~ N(0, b0 I).
##
## The coefficient blocks are updated by Polya-Gamma data augmentation
## (Polson, Scott and Windle, 2013): given omega_i ~ PG(1, x_i' beta),
## beta | omega, y ~ N(V X' kappa, V), V = (X' Omega X + B0^{-1})^{-1},
## kappa = y - 1/2. As in the original code, PG(1, z) variables are drawn
## from a K-term truncation of their infinite sum-of-gammas representation
## plus a deterministic correction for the mean of the omitted tail. This is
## an (accurate) approximation, not an exact PG draw.

## Approximate PG(1, z) draws (verbatim from the original rpg1()).
.rpg1 <- function(z, K = 128) {
  n <- length(z)
  ck <- (seq_len(K) - 0.5)^2
  denom <- outer((z / (2 * pi))^2, ck, "+")
  g <- matrix(stats::rexp(n * K), n, K)
  rowSums(g / denom) / (2 * pi^2) + (1 / (2 * pi^2)) * (1 / K)
}

.bern_ll <- function(y, pr) {
  pr <- pmin(pmax(pr, 1e-9), 1 - 1e-9)
  y * log(pr) + (1 - y) * log(1 - pr)
}

.pg_update_beta <- function(X, y, beta, B0i, K) {
  kappa <- y - 0.5
  omega <- .rpg1(as.vector(X %*% beta), K)
  V <- solve(B0i + crossprod(X * omega, X))
  .rmvnorm(as.vector(V %*% crossprod(X, kappa)), V)
}

.engine_binomial <- function(Xp, yp, Xnp, ynp, model, prior, control) {
  p <- ncol(Xp)
  B0i <- diag(1 / prior$beta_scale, p)
  B0_scale <- prior$beta_scale
  q_a <- prior$q_shape1
  q_b <- prior$q_shape2
  K <- control$pg_terms
  exact_m1 <- identical(control$m1_update, "exact")

  ## initial values (original code)
  beta <- rep(0, p)
  q <- if (length(ynp)) min(max(mean(ynp), 0.05), 0.95) else 0.5
  beta_bad <- rep(0, p)

  loglik_good <- function() .bern_ll(ynp, .sigmoid(as.vector(Xnp %*% beta)))
  loglik_bad <- function() {
    if (model == "M1") .bern_ll(ynp, .sigmoid(as.vector(Xnp %*% beta_bad))) else .bern_ll(ynp, q)
  }
  update_good <- function(good) {
    Xg <- rbind(Xp, Xnp[good, , drop = FALSE])
    yg <- c(yp, ynp[good])
    beta <<- .pg_update_beta(Xg, yg, beta, B0i, K)
    invisible(NULL)
  }
  update_bad <- function(bad) {
    nb <- sum(bad)
    if (model == "M1") {
      use_data <- if (exact_m1) nb > 0 else nb > p + 1
      if (use_data) {
        beta_bad <<- .pg_update_beta(Xnp[bad, , drop = FALSE], ynp[bad], beta_bad, B0i, K)
      } else {
        beta_bad <<- .rmvnorm(rep(0, p), diag(B0_scale, p))
      }
    } else {
      if (nb > 0) {
        yb <- ynp[bad]
        q <<- stats::rbeta(1, q_a + sum(yb), q_b + nb - sum(yb))
      } else {
        q <<- stats::rbeta(1, q_a, q_b)
      }
    }
    invisible(NULL)
  }
  good_params <- function() beta
  bad_params <- function() if (model == "M1") beta_bad else c(q = q)
  bad_names <- if (model == "M1") paste0("bad_", colnames(Xp)) else "bad_q"

  list(family = "binomial", p = p,
       loglik_good = loglik_good, loglik_bad = loglik_bad,
       update_good = update_good, update_bad = update_bad,
       good_params = good_params, bad_params = bad_params,
       good_names = colnames(Xp), bad_names = bad_names,
       accept = function() NULL)
}
