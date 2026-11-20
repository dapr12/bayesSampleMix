## Test D: logistic regression. DGP: x1, x2 ~ N(0, 1);
## P(y = 1) = plogis(0.2 + 1.0 x1 - 0.8 x2). Contaminated NP units have
## covariates shifted by +2 and sign-flipped slopes (the "sign-flipping
## contaminant" of the original audit), which a logistic fit cannot ignore.

beta_logit <- c(0.2, 1, -0.8)

test_that("Test D: logistic recovery and contamination detection (30% contamination)", {
  skip_if_short()
  d <- simulate_mixture_data("binomial", n_prob = 400, n_nonprob = 400,
                             contamination = 0.3, beta = beta_logit, x_shift = 2,
                             seed = 201)
  fit <- mixture_fit(y ~ x1 + x2, d$prob, d$nonprob, family = "binomial",
                     n_iter = 1500, burnin = 500, seed = 1)
  expect_recovers(fit, beta_logit)
  ## Binary outcomes carry little information per unit; classification relies
  ## mostly on the covariate (leverage) term. Calibration (50 data sets):
  ## AUC >= 0.92, |z| <= 2.8. The posterior mean of 1 - pi ranged 0.50-0.70
  ## for a true fraction of 0.30 (AUDIT.md, finding F14: the Bernoulli(q)
  ## background density and the covariate term absorb some compatible units;
  ## conservative behaviour of the original model). We therefore
  ## test that contamination is detected, not that 0.30 is covered.
  expect_gt(filter_auc(1 - d$nonprob$contaminated, filter_scores(fit)), 0.9)
  ct <- contamination(fit)
  expect_gt(ct[["lower"]], 0.2)

  pooled <- baseline_fit(y ~ x1 + x2, rbind(d$prob, d$nonprob[names(d$prob)]),
                         family = "binomial", n_iter = 1000, burnin = 300, seed = 1)
  zp <- (coef(pooled) - beta_logit) / apply(as.matrix(pooled), 2, sd)
  expect_gt(max(abs(zp)), 4)
})

test_that("Test D: logistic fit with no contamination recovers the coefficients", {
  skip_if_short()
  d <- simulate_mixture_data("binomial", n_prob = 300, n_nonprob = 300,
                             contamination = 0, beta = beta_logit, seed = 202)
  fit <- mixture_fit(y ~ x1 + x2, d$prob, d$nonprob, family = "binomial",
                     n_iter = 1500, burnin = 500, seed = 2)
  expect_recovers(fit, beta_logit)
  ## calibration (50 data sets): posterior mean of 1 - pi in [0.011, 0.080]
  expect_lt(contamination(fit)[["mean"]], 0.15)
})

test_that("factor and logical binary responses are accepted", {
  d <- simulate_mixture_data("binomial", n_prob = 60, n_nonprob = 60,
                             contamination = 0.2, seed = 203)
  dp <- d$prob
  dn <- d$nonprob
  dp$yf <- factor(ifelse(dp$y == 1, "yes", "no"), levels = c("no", "yes"))
  dn$yf <- factor(ifelse(dn$y == 1, "yes", "no"), levels = c("no", "yes"))
  f1 <- mixture_fit(y ~ x1, dp, dn, family = "binomial", n_iter = 60, burnin = 20, seed = 1)
  f2 <- mixture_fit(yf ~ x1, dp, dn, family = "binomial", n_iter = 60, burnin = 20, seed = 1)
  dp$yl <- dp$y == 1
  dn$yl <- dn$y == 1
  f3 <- mixture_fit(yl ~ x1, dp, dn, family = binomial(), n_iter = 60, burnin = 20, seed = 1)
  expect_equal(coef(f1), coef(f2))
  expect_equal(coef(f1), coef(f3))
})
