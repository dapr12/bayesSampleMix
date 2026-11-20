## Poisson (log-link) engine.
##
## Target outcome model:   y ~ Poisson(exp(x' beta))
## Contaminating outcome model:
##   M0 / Mjoint:  y ~ Poisson(lambda0),  lambda0 ~ Gamma(lam_a, lam_b)
##   M1:           y ~ Poisson(exp(x' beta_bad))
## Prior: beta, beta_bad ~ N(0, b0 I).
##
## The log link admits no convenient conjugate augmentation, so (as in the
## original code) each coefficient block is updated by a Metropolis-Hastings
## step whose proposal is the Laplace approximation N(mode, V) to the
## conditional posterior; the mode is located by `irls_steps` Newton steps
## started at the current value. When the Newton iterations have converged the
## proposal does not depend on the current state and the step is an exact
## independence sampler; with very few steps from a distant start it is an
## approximation (see ?mixture_fit).

## Penalised IRLS / Newton iterations towards the conditional mode.
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

## One MH step; returns list(beta, accepted).
.update_beta_pois <- function(X, y, beta, B0i, steps = 3) {
  fit <- .irls_mode(X, y, B0i, beta, steps)
  mode <- fit$mode
  V <- fit$V
  prop <- .rmvnorm(mode, V)
  la <- .log_post_pois(X, y, prop, B0i) - .log_post_pois(X, y, beta, B0i)
  lq <- .dmvnorm_log1(beta, mode, V) - .dmvnorm_log1(prop, mode, V)
  if (log(stats::runif(1)) < la + lq) list(beta = prop, accepted = TRUE) else list(beta = beta, accepted = FALSE)
}

.pois_ll <- function(y, mu) {
  mu <- pmax(mu, 1e-9)
  y * log(mu) - mu
}

.engine_poisson <- function(Xp, yp, Xnp, ynp, model, prior, control) {
  p <- ncol(Xp)
  B0i <- diag(1 / prior$beta_scale, p)
  B0_scale <- prior$beta_scale
  lam_a <- prior$lambda_shape
  lam_b <- prior$lambda_rate
  steps <- control$irls_steps
  exact_m1 <- identical(control$m1_update, "exact")

  ## initial values (original code)
  beta <- .irls_mode(Xp, yp, B0i, rep(0, p), 5)$mode
  lam0 <- if (length(ynp)) max(mean(ynp), 0.1) else 1
  beta_bad <- beta
  n_acc <- 0L
  n_prop <- 0L

  loglik_good <- function() .pois_ll(ynp, exp(.clip(as.vector(Xnp %*% beta), -30, 30)))
  loglik_bad <- function() {
    if (model == "M1") {
      .pois_ll(ynp, exp(.clip(as.vector(Xnp %*% beta_bad), -30, 30)))
    } else {
      .pois_ll(ynp, lam0)
    }
  }
  update_good <- function(good) {
    Xg <- rbind(Xp, Xnp[good, , drop = FALSE])
    yg <- c(yp, ynp[good])
    u <- .update_beta_pois(Xg, yg, beta, B0i, steps)
    beta <<- u$beta
    n_acc <<- n_acc + u$accepted
    n_prop <<- n_prop + 1L
    invisible(NULL)
  }
  update_bad <- function(bad) {
    nb <- sum(bad)
    if (model == "M1") {
      use_data <- if (exact_m1) nb > 0 else nb > p + 1
      if (use_data) {
        beta_bad <<- .update_beta_pois(Xnp[bad, , drop = FALSE], ynp[bad], beta_bad, B0i, steps)$beta
      } else {
        beta_bad <<- .rmvnorm(rep(0, p), diag(B0_scale, p))
      }
    } else {
      if (nb > 0) {
        yb <- ynp[bad]
        lam0 <<- stats::rgamma(1, lam_a + sum(yb), rate = lam_b + nb)
      } else {
        lam0 <<- stats::rgamma(1, lam_a, rate = lam_b)
      }
    }
    invisible(NULL)
  }
  good_params <- function() beta
  bad_params <- function() if (model == "M1") beta_bad else c(lambda0 = lam0)
  bad_names <- if (model == "M1") paste0("bad_", colnames(Xp)) else "bad_lambda0"

  list(family = "poisson", p = p,
       loglik_good = loglik_good, loglik_bad = loglik_bad,
       update_good = update_good, update_bad = update_bad,
       good_params = good_params, bad_params = bad_params,
       good_names = colnames(Xp), bad_names = bad_names,
       accept = function() c(n_acc, n_prop))
}
