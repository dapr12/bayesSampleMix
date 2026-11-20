## Tests A-C: Gaussian regression recovery under 0%, 25% and 50% contamination.
## Data-generating process (simulate_mixture_data): x1 ~ N(0, 1) for the P
## sample and compatible NP units, y = 2 + 5 x1 + N(0, 1); contaminated NP
## units as described in each test.

beta_true <- c(2, 5)

test_that("Test A: no contamination - coefficients recovered, contamination small", {
  skip_if_short()
  d <- simulate_mixture_data("gaussian", n_prob = 200, n_nonprob = 200,
                             contamination = 0, beta = beta_true, seed = 101)
  fit <- mixture_fit(y ~ x1, d$prob, d$nonprob, family = "gaussian",
                     n_iter = 2000, burnin = 500, seed = 1)
  expect_recovers(fit, beta_true)
  ## Calibrated with inst/validation/calibrate_tolerances.R (50 data sets):
  ## posterior mean of 1 - pi ranged 0.005-0.012, mean filter score >= 0.993.
  ## With a Beta(1, 1) prior the posterior of 1 - pi cannot be exactly 0.
  expect_lt(contamination(fit)[["mean"]], 0.05)
  expect_gt(mean(filter_scores(fit)), 0.95)
})

test_that("Test B: 25% contamination - recovery, contamination estimate, classification", {
  skip_if_short()
  d <- simulate_mixture_data("gaussian", n_prob = 200, n_nonprob = 200,
                             contamination = 0.25, beta = beta_true, x_shift = 1.5,
                             seed = 102)
  fit <- mixture_fit(y ~ x1, d$prob, d$nonprob, family = "gaussian",
                     n_iter = 2000, burnin = 500, seed = 2)
  expect_recovers(fit, beta_true)
  ct <- contamination(fit)
  expect_true(ct[["lower"]] < 0.25 && 0.25 < ct[["upper"]])
  ## Calibration (50 data sets): interval covered 0.25 in 50/50, AUC >= 0.94.
  ## AUC < 1 is expected: the contaminating line y = 2 - 5 x crosses the
  ## target line at x = 0, where the two mechanisms are indistinguishable.
  expect_gt(filter_auc(1 - d$nonprob$contaminated, filter_scores(fit)), 0.9)

  ## The test is not vacuous: naive pooling is badly biased on the same data.
  pooled <- baseline_fit(y ~ x1, rbind(d$prob, d$nonprob[names(d$prob)]),
                         n_iter = 1500, burnin = 500, seed = 2)
  zp <- (coef(pooled) - beta_true) / apply(as.matrix(pooled)[, 1:2], 2, sd)
  expect_gt(max(abs(zp)), 10)
})

test_that("Test B (M0 and M1 variants) also recover the coefficients", {
  skip_if_short()
  d <- simulate_mixture_data("gaussian", n_prob = 200, n_nonprob = 200,
                             contamination = 0.25, beta = beta_true, x_shift = 1.5,
                             seed = 103)
  f0 <- mixture_fit(y ~ x1, d$prob, d$nonprob, model = "M0",
                    n_iter = 2000, burnin = 500, seed = 3)
  f1 <- mixture_fit(y ~ x1, d$prob, d$nonprob, model = "M1",
                    n_iter = 2000, burnin = 500, seed = 3)
  expect_recovers(f0, beta_true)
  expect_recovers(f1, beta_true)
})

test_that("Test C: 50% contamination (Constanza example) with beta = c(2, 5)", {
  skip_if_short()
  ## Contaminated units: covariates shifted to mean 10, outcomes ~ N(-50, 1).
  d <- simulate_mixture_data("gaussian", n_prob = 100, n_nonprob = 200,
                             contamination = 0.5, beta = beta_true,
                             bad_beta = c(-50, 0), x_shift = 10, seed = 123)
  fit <- mixture_fit(y ~ x1, d$prob, d$nonprob, family = "gaussian",
                     n_iter = 2000, burnin = 1000, seed = 4)
  expect_recovers(fit, beta_true)
  ct <- contamination(fit)
  expect_true(ct[["lower"]] < 0.5 && 0.5 < ct[["upper"]])
  ## Separation is extreme, so every unit should be classified correctly.
  s <- filter_scores(fit)
  expect_true(all(s[d$nonprob$contaminated == 1] < 0.01))
  expect_true(all(s[d$nonprob$contaminated == 0] > 0.99))
  ## realised contamination of the observed units is exactly 0.5
  expect_equal(mean(as.matrix(fit)[, "contamination_realized"]), 0.5, tolerance = 1e-8)
})
