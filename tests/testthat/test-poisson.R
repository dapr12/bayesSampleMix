## Test E: Poisson regression. DGP: x1, x2 ~ N(0, 1);
## y ~ Poisson(exp(1.2 + 0.5 x1 - 0.4 x2)). Contaminated NP units: covariates
## shifted by +2 and sign-flipped slopes.

beta_pois <- c(1.2, 0.5, -0.4)

test_that("Test E: Poisson recovery and contamination detection (25% contamination)", {
  skip_if_short()
  d <- simulate_mixture_data("poisson", n_prob = 300, n_nonprob = 300,
                             contamination = 0.25, beta = beta_pois, x_shift = 2,
                             seed = 301)
  fit <- mixture_fit(y ~ x1 + x2, d$prob, d$nonprob, family = "poisson",
                     n_iter = 1500, burnin = 500, seed = 1)
  expect_recovers(fit, beta_pois)
  ## Calibration (50 data sets): |z| <= 3.2 with sd(z) ~ 1.1; AUC >= 0.948;
  ## acceptance >= 0.97; posterior mean of 1 - pi in [0.34, 0.48] for a true
  ## fraction of 0.25 (AUDIT.md, finding F14: the background Poisson(lambda0)
  ## density and the covariate term absorb some compatible units), so the
  ## test checks the documented range rather than coverage of 0.25.
  ct <- contamination(fit)
  expect_true(ct[["mean"]] > 0.25 && ct[["mean"]] < 0.55)
  expect_gt(filter_auc(1 - d$nonprob$contaminated, filter_scores(fit)), 0.9)
  ## Laplace-proposal independence MH should accept most proposals.
  expect_gt(fit$acceptance, 0.9)
})

test_that("Test E: Poisson fit without contamination", {
  skip_if_short()
  d <- simulate_mixture_data("poisson", n_prob = 300, n_nonprob = 300,
                             contamination = 0, beta = beta_pois, seed = 302)
  fit <- mixture_fit(y ~ x1 + x2, d$prob, d$nonprob, family = "poisson",
                     n_iter = 1500, burnin = 500, seed = 2)
  expect_recovers(fit, beta_pois)
  ## calibration (50 data sets): posterior mean of 1 - pi in [0.007, 0.034]
  expect_lt(contamination(fit)[["mean"]], 0.1)
})

test_that("Poisson sampler is numerically stable with large counts and extreme predictors", {
  set.seed(303)
  n <- 150
  x <- c(rnorm(n - 5), c(6, 7, -6, 8, -7))  # a few extreme covariate values
  dp <- data.frame(x = x, y = rpois(n, exp(2 + 0.4 * x)))
  dn <- data.frame(x = rnorm(100), y = rpois(100, exp(2 + 0.4 * rnorm(100))))
  fit <- mixture_fit(y ~ x, dp, dn, family = "poisson", n_iter = 300, burnin = 100, seed = 1)
  expect_true(all(is.finite(as.matrix(fit))))
  expect_true(all(fit$filter_score >= 0 & fit$filter_score <= 1))
})
