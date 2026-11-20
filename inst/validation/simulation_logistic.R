## Factorial Monte Carlo study for the logistic contamination mixture.
## Port of run_glm_sim() from the original mixture_filter_glm.R to the package
## API, with the same data-generating process, seeds and settings (results
## are identical to the original, see tests/testthat/test-equivalence.R).
##
## Usage:  Rscript simulation_logistic.R [nreps] [outdir]   (default nreps = 30)

library(bayesSampleMix)

run_glm_sim <- function(nreps = 30, nPs = c(100, 300, 500), iters = 800, burn = 400,
                        K = 64, c_infl = 2.5, seed0 = 20260901, out = "results_glm") {
  dir.create(out, showWarnings = FALSE)
  BETA <- c(0.2, 1.0, -0.8, 0.9, -1.2)
  MU <- c(0, 0)
  SIG <- matrix(c(1, .5, .5, 1), 2)
  SEV <- c(High = 0.4, Moderate = 1.0, Low = 1.7)
  NG <- 150
  NB <- 100
  SHIFT <- c(2, -2)
  DB <- c(-0.5, -2.2, 1.8, 0, 2.6)
  sigmoid <- function(x) 1 / (1 + exp(-x))
  rmvn <- function(n, mu, S) matrix(rnorm(n * length(mu)), n) %*% chol(S) +
    matrix(mu, n, length(mu), byrow = TRUE)
  design <- function(a, b, c) cbind(1, a, b, as.numeric(c == 1), as.numeric(c == 2))
  gen <- function(nP, cond, seed) {
    set.seed(seed)
    s <- SEV[[cond]]
    good <- function(n) {
      XX <- rmvn(n, MU, SIG)
      ct <- sample(0:2, n, TRUE, c(.5, .3, .2))
      Xd <- design(XX[, 1], XX[, 2], ct)
      list(X = Xd, y = rbinom(n, 1, sigmoid(as.vector(Xd %*% BETA))))
    }
    mub <- MU + s * SHIFT
    w <- min(s / 1.7, 1) * 0.9
    catp <- (1 - w) * c(.5, .3, .2) + w * c(.05, .05, .90)
    P <- good(nP)
    G <- good(NG)
    XXb <- rmvn(NB, mub, SIG)
    cb <- sample(0:2, NB, TRUE, catp)
    Xb <- design(XXb[, 1], XXb[, 2], cb)
    yb <- rbinom(NB, 1, sigmoid(as.vector(Xb %*% (BETA + s * DB))))
    list(Xp = P$X, yp = P$y, Xnp = rbind(G$X, Xb), ynp = c(G$y, yb),
         z = c(rep(1, NG), rep(0, NB)))
  }
  ctrl <- mixture_control(pg_terms = K, c_infl = c_infl)
  bl <- function(X, y, s) coef(baseline_fit_matrix(X, y, "binomial", n_iter = iters,
                                                   burnin = burn, seed = s,
                                                   control = ctrl))
  ests <- c("Prob", "Naive", "Oracle", "M0", "Mjoint")
  rows <- list()
  ci <- 0
  for (cond in names(SEV)) for (nP in nPs) {
    ci <- ci + 1
    B <- setNames(lapply(ests, function(.) matrix(NA_real_, nreps, 5)), ests)
    C <- A <- numeric(nreps)
    for (r in seq_len(nreps)) {
      s <- seed0 + 7919 * ci + r
      d <- gen(nP, cond, s)
      B$Prob[r, ] <- bl(d$Xp, d$yp, s)
      B$Naive[r, ] <- bl(rbind(d$Xp, d$Xnp), c(d$yp, d$ynp), s)
      gm <- d$z > 0.5
      B$Oracle[r, ] <- bl(rbind(d$Xp, d$Xnp[gm, ]), c(d$yp, d$ynp[gm]), s)
      B$M0[r, ] <- coef(mixture_fit_matrix(d$Xp, d$yp, d$Xnp, d$ynp, "binomial", "M0",
                                           c(2, 3), iters, burn, seed = s,
                                           control = mixture_control(pg_terms = K)))
      mj <- mixture_fit_matrix(d$Xp, d$yp, d$Xnp, d$ynp, "binomial", "Mjoint", c(2, 3),
                               iters, burn, seed = s, control = ctrl)
      B$Mjoint[r, ] <- coef(mj)
      C[r] <- contamination(mj)[["mean"]]
      A[r] <- filter_auc(d$z, mj$filter_score)
    }
    for (k in ests) {
      arr <- B[[k]]
      mse <- (colMeans(arr) - BETA)^2 + apply(arr, 2, var)
      rows[[length(rows) + 1]] <- data.frame(cond, nP, est = k, totMSE = sum(mse),
                                             contam = mean(C), AUC = mean(A, na.rm = TRUE))
    }
    cat(sprintf("%-9s nP=%d done\n", cond, nP))
  }
  res <- do.call(rbind, rows)
  write.csv(res, file.path(out, "glm_results.csv"), row.names = FALSE)
  res
}

if (sys.nframe() == 0) {
  a <- commandArgs(trailingOnly = TRUE)
  nr <- if (length(a) >= 1) as.integer(a[1]) else 30
  out <- if (length(a) >= 2) a[2] else "results_glm"
  run_glm_sim(nreps = nr, out = out)
}
