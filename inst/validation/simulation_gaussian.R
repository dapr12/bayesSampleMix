## Factorial Monte Carlo study for the Gaussian contamination mixture.
## Port of the original simulation.R to the package API. The data-generating
## process, seeds and estimator settings are unchanged, and because the
## package reproduces the original samplers draw for draw (see
## tests/testthat/test-equivalence.R), the results are identical to those of
## the original script.
##
## Usage (from any directory):  Rscript simulation_gaussian.R [nreps] [outdir]
## Default nreps = 80 as in the original; expect several hours.

library(bayesSampleMix)

BETA <- c(2, 1.5, -1, 3, -2)          # intercept, X1, X2, catB, catC
MU <- c(1, 2)
SIG <- matrix(c(1, .5, .5, 1), 2)
SDE <- 1.5
SEV <- c(High = 0.4, Moderate = 1.0, Low = 1.7)   # contamination severity
NG <- 150
NB <- 100                             # good / biased NP counts (contam = 0.40)
PARAM_NAMES <- c("Intercept", "X1", "X2", "catB", "catC")
CONT_IDX <- c(2, 3)

rmvn <- function(n, mu, S) {
  L <- chol(S)
  matrix(rnorm(n * length(mu)), nrow = n, ncol = length(mu)) %*% L +
    matrix(mu, n, length(mu), byrow = TRUE)
}
design <- function(a, b, cat) cbind(1, a, b, as.numeric(cat == 1), as.numeric(cat == 2))

gen <- function(nP, cond, seed, NG = 150, NB = 100) {
  set.seed(seed)
  s <- SEV[[cond]]
  good <- function(n) {
    XX <- rmvn(n, MU, SIG)
    cat <- sample(0:2, n, TRUE, c(.5, .3, .2))
    Xd <- design(XX[, 1], XX[, 2], cat)
    list(X = Xd, y = as.vector(Xd %*% BETA) + rnorm(n, 0, SDE))
  }
  mub <- MU + s * c(2, -2)
  w <- min(s / 1.7, 1) * 0.9
  catp <- (1 - w) * c(.5, .3, .2) + w * c(.05, .05, .90)
  dbeta <- s * c(-0.6, 0, 0, 0, 2.5)
  P <- good(nP)
  G <- good(NG)
  XXb <- rmvn(NB, mub, SIG)
  catb <- sample(0:2, NB, TRUE, catp)
  Xb <- design(XXb[, 1], XXb[, 2], catb)
  yb <- as.vector(Xb %*% (BETA + dbeta)) + rnorm(NB, 0, SDE)
  list(Xp = P$X, yp = P$y, Xnp = rbind(G$X, Xb), ynp = c(G$y, yb),
       z = c(rep(1, NG), rep(0, NB)))
}

run_sim <- function(nreps = 80, nPs = c(100, 300, 500), iters = 1600, burn = 700,
                    seed0 = 20260828, out = "results") {
  dir.create(out, showWarnings = FALSE)
  ests <- c("Prob", "Naive", "Oracle", "M0", "Mjoint")
  mse_rows <- bias_rows <- auc_rows <- contam_rows <- list()
  cells <- expand.grid(cond = names(SEV), nP = nPs, stringsAsFactors = FALSE)
  bl <- function(X, y, s) coef(baseline_fit_matrix(X, y, "gaussian", n_iter = iters, burnin = burn, seed = s))
  for (ci in seq_len(nrow(cells))) {
    cond <- cells$cond[ci]
    nP <- cells$nP[ci]
    B <- setNames(lapply(ests, function(.) matrix(NA_real_, nreps, length(BETA))), ests)
    contam <- aucs <- oc <- ef <- numeric(nreps)
    for (r in seq_len(nreps)) {
      s <- seed0 + 7919 * ci + r
      d <- gen(nP, cond, s)
      s0 <- overlap_diagnostics(d$Xp[, CONT_IDX], d$Xnp[, CONT_IDX])
      oc[r] <- s0$OC
      ef[r] <- s0$EF
      B$Prob[r, ] <- bl(d$Xp, d$yp, s)
      B$Naive[r, ] <- bl(rbind(d$Xp, d$Xnp), c(d$yp, d$ynp), s)
      gm <- d$z > 0.5
      B$Oracle[r, ] <- bl(rbind(d$Xp, d$Xnp[gm, ]), c(d$yp, d$ynp[gm]), s)
      B$M0[r, ] <- coef(mixture_fit_matrix(d$Xp, d$yp, d$Xnp, d$ynp, "gaussian", "M0",
                                           CONT_IDX, iters, burn, seed = s))
      mj <- mixture_fit_matrix(d$Xp, d$yp, d$Xnp, d$ynp, "gaussian", "Mjoint",
                               CONT_IDX, iters, burn, seed = s)
      B$Mjoint[r, ] <- coef(mj)
      contam[r] <- contamination(mj)[["mean"]]
      aucs[r] <- filter_auc(d$z, mj$filter_score)
    }
    for (k in ests) {
      arr <- B[[k]]
      bias <- colMeans(arr) - BETA
      mse <- bias^2 + apply(arr, 2, var)
      mse_rows[[length(mse_rows) + 1]] <- data.frame(cond, nP, est = k, totMSE = sum(mse),
                                                     mse_catC = mse[5])
      for (j in seq_along(PARAM_NAMES)) {
        bias_rows[[length(bias_rows) + 1]] <- data.frame(cond, nP, est = k,
                                                         param = PARAM_NAMES[j], bias = bias[j])
      }
    }
    auc_rows[[length(auc_rows) + 1]] <- data.frame(cond, nP, meanAUC = mean(aucs, na.rm = TRUE),
                                                   sdAUC = sd(aucs, na.rm = TRUE))
    contam_rows[[length(contam_rows) + 1]] <- data.frame(cond, nP, contam_mean = mean(contam),
                                                         contam_sd = sd(contam), OC = mean(oc),
                                                         EF = mean(ef))
    cat(sprintf("%-9s nP=%d done\n", cond, nP))
  }
  write.csv(do.call(rbind, mse_rows), file.path(out, "mse.csv"), row.names = FALSE)
  write.csv(do.call(rbind, bias_rows), file.path(out, "bias.csv"), row.names = FALSE)
  write.csv(do.call(rbind, auc_rows), file.path(out, "auc.csv"), row.names = FALSE)
  write.csv(do.call(rbind, contam_rows), file.path(out, "contam.csv"), row.names = FALSE)
  cat("results written to", out, "\n")
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  nreps <- if (length(args) >= 1) as.integer(args[1]) else 80
  out <- if (length(args) >= 2) args[2] else "results"
  run_sim(nreps = nreps, out = out)
}
