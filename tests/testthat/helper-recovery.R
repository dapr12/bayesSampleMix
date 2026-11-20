## Helpers shared by the statistical-recovery tests.
##
## Tolerance rationale. For a correctly specified model and a sample of
## moderate size, (posterior mean - true value) / posterior sd is
## approximately standard normal over repeated samples (Bernstein-von Mises).
## We therefore require |z| < 4 for every coefficient (two-sided frequentist
## failure probability about 6e-5 per coefficient if the software is
## correct), while Monte Carlo error of the posterior mean is negligible in
## comparison (MCSE << posterior sd for the chain lengths used). A software
## error that biases estimates by a few posterior standard deviations is
## detected, whereas sampling variability is not mistaken for an error.
## Fixed seeds make every test deterministic; the calibration script in
## inst/validation/calibrate_tolerances.R checks the thresholds across many
## seeds.

z_scores <- function(fit, truth) {
  d <- as.matrix(fit)[, fit$coef_names, drop = FALSE]
  (colMeans(d) - truth) / apply(d, 2, sd)
}

expect_recovers <- function(fit, truth, z_max = 4) {
  z <- z_scores(fit, truth)
  testthat::expect_true(all(abs(z) < z_max),
                        info = paste("z-scores:", paste(round(z, 2), collapse = ", ")))
}

skip_if_short <- function() {
  ## The longer recovery tests are skipped on CRAN to respect check-time
  ## limits; they run locally and on GitHub Actions (NOT_CRAN = "true").
  ## Set BSM_SKIP_LONG=true to skip them elsewhere as well.
  testthat::skip_on_cran()
  testthat::skip_if(identical(Sys.getenv("BSM_SKIP_LONG"), "true"), "long tests skipped")
}
