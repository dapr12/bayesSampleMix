## Verbatim function definitions from the original mixture_filter.R (top-level demo code removed).
## Used only by test-equivalence.R to verify that the package reproduces the original samplers.
## mixture_filter.R
## ----------------------------------------------------------------------------
## Fully Bayesian contamination-mixture model for probability / non-probability
## (P / NP) data integration.
##
## Model
##   P sample (all trusted-good):  y_i ~ N(x_i' beta, sigma^2)
##   NP sample (latent mixture):   z_j ~ Bernoulli(pi)
##      good (z=1): y_j ~ N(x_j' beta, sigma^2)   [ * N(x_j; muP, SigP)   for Mjoint ]
##      bad  (z=0): y_j ~ N(mu0, tau^2)           [ * N(x_j; muP, c*SigP) for Mjoint ]
##      (variant M1: bad outcome density is a rival regression N(x_j' beta_bad, s2_bad))
##
##   The P sample anchors the "good" component, so there is no label switching.
##   The discrete indicator z is sampled directly (data augmentation), which makes
##   every full conditional conjugate -> a fast, dependency-free Gibbs sampler.
##
##   The posterior mean responsibility  s_j = E[z_j | data]  is the filter score
##   (analogue of the BayesianFilterScore); no proxy labels and no GAM are needed.
##
## This is a line-for-line port of the reference Python implementation used to
## produce the paper's numbers (mixmodel.py). Base R only.
## ----------------------------------------------------------------------------

## ---- small linear-algebra helpers (no external packages) -------------------
.rmvnorm <- function(mu, Sigma) {
  p <- length(mu)
  L <- chol(Sigma)                 # upper triangular, Sigma = t(L) %*% L
  mu + as.vector(t(L) %*% rnorm(p))
}

.dmvnorm_log <- function(X, mu, Sigma) {
  ## log density of rows of X under N(mu, Sigma)
  p <- ncol(X)
  L <- chol(Sigma)
  Xc <- sweep(X, 2, mu, "-")
  z  <- backsolve(L, t(Xc), transpose = TRUE)   # solves t(L) z = Xc'
  quad <- colSums(z^2)
  logdet <- 2 * sum(log(diag(L)))
  -0.5 * (p * log(2 * pi) + logdet + quad)
}

.dnorm_log <- function(y, mean, var) {
  -0.5 * (log(2 * pi * var) + (y - mean)^2 / var)
}

## ---- conjugate Bayesian linear regression (baseline estimators) ------------
fit_bayes_lm <- function(X, y, iters = 2000, burn = 1000, B0_scale = 25,
                         a_sig = 2, b_sig = 1, seed = 0) {
  set.seed(seed)
  X <- as.matrix(X); y <- as.numeric(y)
  n <- nrow(X); p <- ncol(X)
  B0i <- diag(1 / B0_scale, p)
  beta <- as.vector(qr.solve(X, y))
  sigma2 <- as.numeric(var(y - X %*% beta)) + 1e-6
  XtX <- crossprod(X); Xty <- crossprod(X, y)
  draws <- matrix(NA_real_, iters - burn, p + 1)
  for (it in seq_len(iters)) {
    Vb <- solve(B0i + XtX / sigma2)
    mb <- Vb %*% (Xty / sigma2)
    beta <- .rmvnorm(as.vector(mb), Vb)
    resid <- y - X %*% beta
    sigma2 <- 1 / rgamma(1, a_sig + n / 2, b_sig + 0.5 * sum(resid^2))
    if (it > burn) draws[it - burn, ] <- c(beta, sigma2)
  }
  list(beta_mean = colMeans(draws[, 1:p, drop = FALSE]),
       beta_sd   = apply(draws[, 1:p, drop = FALSE], 2, sd),
       sigma2_mean = mean(draws[, p + 1]))
}

## ---- contamination-mixture Gibbs sampler -----------------------------------
## cont_idx : column indices of the design (incl. intercept) that are the
##            CONTINUOUS predictors used for the covariate (leverage) density.
## model    : "M0" (covariate-free bad outcome density),
##            "Mjoint" (covariate-augmented; leverage-aware),
##            "M1" (rival-regression bad component).
fit_mixture <- function(Xp, yp, Xnp, ynp, cont_idx, model = "M0",
                        iters = 3000, burn = 1500, c_infl = 4,
                        B0_scale = 25, a_sig = 2, b_sig = 1,
                        pi_a = 1, pi_b = 1, a_tau = 2, b_tau = 2,
                        m0_var = 100, seed = 0) {
  set.seed(seed)
  Xp <- as.matrix(Xp); yp <- as.numeric(yp)
  Xnp <- as.matrix(Xnp); ynp <- as.numeric(ynp)
  nP <- nrow(Xp); p <- ncol(Xp); nNP <- nrow(Xnp)
  B0i <- diag(1 / B0_scale, p); b0 <- rep(0, p)
  m0_mean <- mean(ynp)

  ## reference covariate moments (empirical Bayes plug-in for leverage screen)
  Xp_c  <- Xp[, cont_idx, drop = FALSE]
  Xnp_c <- Xnp[, cont_idx, drop = FALSE]
  muP  <- colMeans(Xp_c)
  SigP <- cov(Xp_c)

  if (model == "Mjoint") {
    logcov_good <- .dmvnorm_log(Xnp_c, muP, SigP)
    logcov_bad  <- .dmvnorm_log(Xnp_c, muP, c_infl * SigP)
  } else {
    logcov_good <- rep(0, nNP); logcov_bad <- rep(0, nNP)
  }

  ## initialise
  beta   <- as.vector(qr.solve(Xp, yp))
  sigma2 <- as.numeric(var(yp - Xp %*% beta)) + 1e-6
  pi     <- 0.7
  mu0    <- m0_mean
  tau2   <- 4 * sigma2
  beta_bad <- beta
  s2_bad   <- tau2

  resp_sum <- numeric(nNP); n_keep <- 0
  keep <- matrix(NA_real_, iters - burn, p + 3)   # beta, sigma2, pi, contam

  for (it in seq_len(iters)) {
    ## ---- z_j responsibilities ----
    log_good <- .dnorm_log(ynp, as.vector(Xnp %*% beta), sigma2) + logcov_good + log(pi)
    if (model == "M1") {
      log_bad <- .dnorm_log(ynp, as.vector(Xnp %*% beta_bad), s2_bad) + logcov_bad + log(1 - pi)
    } else {
      log_bad <- .dnorm_log(ynp, mu0, tau2) + logcov_bad + log(1 - pi)
    }
    m  <- pmax(log_good, log_bad)
    pg <- exp(log_good - m); pb <- exp(log_bad - m)
    r  <- pg / (pg + pb)                    # P(z_j = 1)
    z  <- as.numeric(runif(nNP) < r)

    ## ---- pi ----
    n_good_np <- sum(z)
    pi <- rbeta(1, pi_a + n_good_np, pi_b + nNP - n_good_np)

    ## ---- beta, sigma2 on good set (P + good NP) ----
    good <- z > 0.5
    Xg <- rbind(Xp, Xnp[good, , drop = FALSE])
    yg <- c(yp, ynp[good])
    Vb <- solve(B0i + crossprod(Xg) / sigma2)
    mb <- Vb %*% (B0i %*% b0 + crossprod(Xg, yg) / sigma2)
    beta <- .rmvnorm(as.vector(mb), Vb)
    resid <- yg - Xg %*% beta
    ng <- length(yg)
    sigma2 <- 1 / rgamma(1, a_sig + ng / 2, b_sig + 0.5 * sum(resid^2))

    ## ---- bad component ----
    bad <- !good; nb <- sum(bad)
    if (model == "M1") {
      if (nb > p + 1) {
        Xb <- Xnp[bad, , drop = FALSE]; yb <- ynp[bad]
        Vbb <- solve(B0i + crossprod(Xb) / s2_bad)
        mbb <- Vbb %*% (crossprod(Xb, yb) / s2_bad)
        beta_bad <- .rmvnorm(as.vector(mbb), Vbb)
        rb <- yb - Xb %*% beta_bad
        s2_bad <- 1 / rgamma(1, a_sig + nb / 2, b_sig + 0.5 * sum(rb^2))
      } else {
        beta_bad <- .rmvnorm(rep(0, p), diag(B0_scale, p))
        s2_bad <- 1 / rgamma(1, a_sig, b_sig)
      }
    } else {
      if (nb > 0) {
        yb <- ynp[bad]
        post_var <- 1 / (1 / m0_var + nb / tau2)
        post_mean <- post_var * (m0_mean / m0_var + sum(yb) / tau2)
        mu0 <- rnorm(1, post_mean, sqrt(post_var))
        tau2 <- 1 / rgamma(1, a_tau + nb / 2, b_tau + 0.5 * sum((yb - mu0)^2))
      } else {
        mu0 <- rnorm(1, m0_mean, sqrt(m0_var))
        tau2 <- 1 / rgamma(1, a_tau, b_tau)
      }
    }

    if (it > burn) {
      keep[it - burn, ] <- c(beta, sigma2, pi, 1 - pi)
      resp_sum <- resp_sum + r; n_keep <- n_keep + 1
    }
  }

  filt <- resp_sum / n_keep
  ent  <- -(ifelse(filt > 0, filt * log(filt), 0) +
            ifelse(filt < 1, (1 - filt) * log(1 - filt), 0))
  list(beta_draws = keep[, 1:p, drop = FALSE],
       beta_mean  = colMeans(keep[, 1:p, drop = FALSE]),
       beta_sd    = apply(keep[, 1:p, drop = FALSE], 2, sd),
       sigma2_mean = mean(keep[, p + 1]),
       pi_mean     = mean(keep[, p + 2]),
       contam_mean = mean(keep[, p + 3]),
       contam_draws = keep[, p + 3],
       filter_score = filt,
       indicator_entropy = ent)
}

## ---- Stage 0 overlap diagnostics -------------------------------------------
stage0 <- function(Xp, Xnp, cont_idx) {
  Xpc <- Xp[, cont_idx, drop = FALSE]; Xnc <- Xnp[, cont_idx, drop = FALSE]
  mu <- colMeans(Xpc); S <- cov(Xpc); Si <- solve(S)
  dM <- function(M) sqrt(rowSums((sweep(M, 2, mu, "-") %*% Si) * sweep(M, 2, mu, "-")))
  dP <- dM(Xpc); dN <- dM(Xnc)
  q99 <- quantile(dP, 0.99); EF <- mean(dN > q99)
  gr <- seq(min(c(dP, dN)), max(c(dP, dN)), length.out = 400)
  fp <- density(dP, from = min(gr), to = max(gr), n = 400)$y
  fn <- density(dN, from = min(gr), to = max(gr), n = 400)$y
  OC <- sum(pmin(fp, fn)) / sum(fp)
  ## propensity-score AUC
  memb <- c(rep(0, nrow(Xpc)), rep(1, nrow(Xnc)))
  Xall <- rbind(Xpc, Xnc)
  gl <- glm(memb ~ Xall, family = binomial())
  ph <- predict(gl, type = "response")
  pos <- ph[memb == 1]; neg <- ph[memb == 0]
  auc <- mean(outer(pos, neg, ">")) + 0.5 * mean(outer(pos, neg, "=="))
  list(OC = as.numeric(OC), EF = as.numeric(EF), AUC_PS = as.numeric(auc))
}

## AUC helper (filter score vs known labels; simulation only)
auc_score <- function(labels, scores) {
  pos <- scores[labels == 1]; neg <- scores[labels == 0]
  if (length(pos) == 0 || length(neg) == 0) return(NA_real_)
  mean(outer(pos, neg, ">")) + 0.5 * mean(outer(pos, neg, "=="))
}

