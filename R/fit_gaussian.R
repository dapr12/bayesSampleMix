## Gaussian (identity-link) engine.
##
## Target ("good") outcome model:   y ~ N(x' beta, sigma2)
## Contaminating outcome model:
##   M0 / Mjoint:  y ~ N(mu0, tau2)
##   M1:           y ~ N(x' beta_bad, s2_bad)
## Priors: beta ~ N(0, b0 I), sigma2 ~ IG(a_sig, b_sig),
##         mu0 ~ N(mean(y_NP), mu0_var), tau2 ~ IG(a_tau, b_tau),
##         beta_bad ~ N(0, b0 I), s2_bad ~ IG(a_sig, b_sig).
## All full conditionals are conjugate (Gibbs). The code below follows the
## original mixture_filter.R statement by statement so that the sequence of
## random numbers, and therefore every draw, is identical for a given seed.

.engine_gaussian <- function(Xp, yp, Xnp, ynp, model, prior, control) {
  p <- ncol(Xp)
  B0i <- diag(1 / prior$beta_scale, p)
  b0 <- rep(0, p)
  a_sig <- prior$sigma_shape
  b_sig <- prior$sigma_rate
  a_tau <- prior$tau_shape
  b_tau <- prior$tau_rate
  m0_var <- prior$mu0_var
  B0_scale <- prior$beta_scale
  m0_mean <- if (length(ynp)) mean(ynp) else NA_real_
  exact_m1 <- identical(control$m1_update, "exact")

  ## initial values (original code)
  beta <- as.vector(qr.solve(Xp, yp))
  sigma2 <- as.numeric(stats::var(yp - Xp %*% beta)) + 1e-6
  mu0 <- m0_mean
  tau2 <- 4 * sigma2
  beta_bad <- beta
  s2_bad <- tau2

  loglik_good <- function() {
    .dnorm_log(ynp, as.vector(Xnp %*% beta), sigma2)
  }
  loglik_bad <- function() {
    if (model == "M1") {
      .dnorm_log(ynp, as.vector(Xnp %*% beta_bad), s2_bad)
    } else {
      .dnorm_log(ynp, mu0, tau2)
    }
  }
  update_good <- function(good) {
    Xg <- rbind(Xp, Xnp[good, , drop = FALSE])
    yg <- c(yp, ynp[good])
    Vb <- solve(B0i + crossprod(Xg) / sigma2)
    mb <- Vb %*% (B0i %*% b0 + crossprod(Xg, yg) / sigma2)
    beta <<- .rmvnorm(as.vector(mb), Vb)
    resid <- yg - Xg %*% beta
    ng <- length(yg)
    sigma2 <<- 1 / stats::rgamma(1, a_sig + ng / 2, b_sig + 0.5 * sum(resid^2))
    invisible(NULL)
  }
  update_bad <- function(bad) {
    nb <- sum(bad)
    if (model == "M1") {
      use_data <- if (exact_m1) nb > 0 else nb > p + 1
      if (use_data) {
        Xb <- Xnp[bad, , drop = FALSE]
        yb <- ynp[bad]
        Vbb <- solve(B0i + crossprod(Xb) / s2_bad)
        mbb <- Vbb %*% (crossprod(Xb, yb) / s2_bad)
        beta_bad <<- .rmvnorm(as.vector(mbb), Vbb)
        rb <- yb - Xb %*% beta_bad
        s2_bad <<- 1 / stats::rgamma(1, a_sig + nb / 2, b_sig + 0.5 * sum(rb^2))
      } else {
        beta_bad <<- .rmvnorm(rep(0, p), diag(B0_scale, p))
        s2_bad <<- 1 / stats::rgamma(1, a_sig, b_sig)
      }
    } else {
      if (nb > 0) {
        yb <- ynp[bad]
        post_var <- 1 / (1 / m0_var + nb / tau2)
        post_mean <- post_var * (m0_mean / m0_var + sum(yb) / tau2)
        mu0 <<- stats::rnorm(1, post_mean, sqrt(post_var))
        tau2 <<- 1 / stats::rgamma(1, a_tau + nb / 2, b_tau + 0.5 * sum((yb - mu0)^2))
      } else {
        mu0 <<- stats::rnorm(1, m0_mean, sqrt(m0_var))
        tau2 <<- 1 / stats::rgamma(1, a_tau, b_tau)
      }
    }
    invisible(NULL)
  }
  good_params <- function() c(beta, sigma2 = sigma2)
  bad_params <- function() {
    if (model == "M1") c(beta_bad, s2_bad = s2_bad) else c(mu0 = mu0, tau2 = tau2)
  }
  bad_names <- if (model == "M1") c(paste0("bad_", colnames(Xp)), "bad_sigma2") else c("bad_mu0", "bad_tau2")

  list(family = "gaussian", p = p,
       loglik_good = loglik_good, loglik_bad = loglik_bad,
       update_good = update_good, update_bad = update_bad,
       good_params = good_params, bad_params = bad_params,
       good_names = c(colnames(Xp), "sigma2"), bad_names = bad_names,
       accept = function() NULL)
}
