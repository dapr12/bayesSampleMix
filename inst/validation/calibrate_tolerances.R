## Calibration of the thresholds used in tests/testthat/test-gaussian.R,
## test-logistic.R and test-poisson.R.
##
## For each test design the data are regenerated with many seeds and the
## statistic used in the test is recorded. A threshold is acceptable when it
## is satisfied in every replicate (and, for the coefficient z-scores, the
## empirical distribution is close to N(0, 1), confirming that |z| < 4 is a
## tolerance on statistical variation rather than a loose bound).
##
## Run from the package root with:  Rscript inst/validation/calibrate_tolerances.R [nrep]
## Runtime: roughly 10-20 minutes for nrep = 50 on a laptop.

library(bayesSampleMix)
args <- commandArgs(trailingOnly = TRUE)
nrep <- if (length(args)) as.integer(args[1]) else 50L

zsc <- function(fit, truth) {
  d <- as.matrix(fit)[, fit$coef_names, drop = FALSE]
  (colMeans(d) - truth) / apply(d, 2, sd)
}

designs <- list(
  gaussian_A = list(family = "gaussian", truth = c(2, 5), formula = y ~ x1, n = c(200, 200),
                    sim = list(contamination = 0), iter = c(2000, 500)),
  gaussian_B = list(family = "gaussian", truth = c(2, 5), formula = y ~ x1, n = c(200, 200),
                    sim = list(contamination = 0.25, x_shift = 1.5), iter = c(2000, 500)),
  gaussian_C = list(family = "gaussian", truth = c(2, 5), formula = y ~ x1, n = c(100, 200),
                    sim = list(contamination = 0.5, bad_beta = c(-50, 0), x_shift = 10),
                    iter = c(2000, 1000)),
  logistic_D = list(family = "binomial", truth = c(0.2, 1, -0.8), formula = y ~ x1 + x2,
                    n = c(400, 400), sim = list(contamination = 0.3, x_shift = 2),
                    iter = c(1500, 500)),
  logistic_D0 = list(family = "binomial", truth = c(0.2, 1, -0.8), formula = y ~ x1 + x2,
                     n = c(300, 300), sim = list(contamination = 0), iter = c(1500, 500)),
  poisson_E = list(family = "poisson", truth = c(1.2, 0.5, -0.4), formula = y ~ x1 + x2,
                   n = c(300, 300), sim = list(contamination = 0.25, x_shift = 2),
                   iter = c(1500, 500)),
  poisson_E0 = list(family = "poisson", truth = c(1.2, 0.5, -0.4), formula = y ~ x1 + x2,
                    n = c(300, 300), sim = list(contamination = 0), iter = c(1500, 500))
)

results <- list()
for (nm in names(designs)) {
  ds <- designs[[nm]]
  rows <- vector("list", nrep)
  for (r in seq_len(nrep)) {
    d <- do.call(simulate_mixture_data,
                 c(list(family = ds$family, n_prob = ds$n[1], n_nonprob = ds$n[2],
                        beta = ds$truth, seed = 10000 + r), ds$sim))
    fit <- mixture_fit(ds$formula, d$prob, d$nonprob, family = ds$family,
                       n_iter = ds$iter[1], burnin = ds$iter[2], seed = r)
    ct <- contamination(fit)
    true_c <- ds$sim$contamination
    rows[[r]] <- data.frame(
      design = nm, rep = r, max_abs_z = max(abs(zsc(fit, ds$truth))),
      z1 = zsc(fit, ds$truth)[2],
      contam_mean = ct[["mean"]],
      covered = ct[["lower"]] < true_c && true_c < ct[["upper"]],
      auc = if (true_c > 0) filter_auc(1 - d$nonprob$contaminated, filter_scores(fit)) else NA,
      mean_score = mean(filter_scores(fit)),
      acceptance = fit$acceptance)
  }
  results[[nm]] <- do.call(rbind, rows)
  res <- results[[nm]]
  cat(sprintf(paste0("%-12s max|z|: max %.2f, 99%% %.2f | sd(z_slope) %.2f | contam mean range [%.3f, %.3f] ",
                     "| coverage %.2f | min AUC %s | min mean score %.3f | min acc %s\n"),
              nm, max(res$max_abs_z), quantile(res$max_abs_z, 0.99), sd(res$z1),
              min(res$contam_mean), max(res$contam_mean), mean(res$covered),
              format(round(min(res$auc), 3)), min(res$mean_score), format(round(min(res$acceptance), 3))))
}
out <- do.call(rbind, results)
dir.create("validation_results", showWarnings = FALSE)
write.csv(out, file.path("validation_results", "calibration.csv"), row.names = FALSE)
