## Test F: computational robustness - invalid inputs, missing values,
## singular designs, extreme predictors, small samples, MCMC settings and
## reproducibility.

small <- function(seed = 1, family = "gaussian") {
  simulate_mixture_data(family, n_prob = 40, n_nonprob = 30, contamination = 0.2, seed = seed)
}

test_that("invalid MCMC settings are rejected with informative errors", {
  d <- small()
  expect_error(mixture_fit(y ~ x1, d$prob, d$nonprob, n_iter = 100, burnin = 100), "burnin")
  expect_error(mixture_fit(y ~ x1, d$prob, d$nonprob, n_iter = 100, burnin = -1), "burnin")
  expect_error(mixture_fit(y ~ x1, d$prob, d$nonprob, n_iter = 10.5, burnin = 1), "n_iter")
  expect_error(mixture_fit(y ~ x1, d$prob, d$nonprob, n_iter = 100, burnin = 50, thin = 0), "thin")
  expect_error(mixture_fit(y ~ x1, d$prob, d$nonprob, n_iter = 100, burnin = 50, thin = 40),
               "Fewer than 2")
  expect_error(mixture_fit(y ~ x1, d$prob, d$nonprob, seed = "a"), "seed")
})

test_that("invalid data, family and formula inputs are rejected", {
  d <- small()
  expect_error(mixture_fit(y ~ x1, d$prob, d$nonprob, family = "gamma"), "family")
  expect_error(mixture_fit(~ x1, d$prob, d$nonprob), "two-sided")
  expect_error(mixture_fit(y ~ x9, d$prob, d$nonprob), "not found")
  expect_error(mixture_fit(y ~ x1, as.matrix(d$prob), d$nonprob), "data frame")
  expect_error(mixture_fit(y ~ x1 + offset(x1), d$prob, d$nonprob), "Offsets")
  expect_error(mixture_fit(y ~ x1, d$prob, d$nonprob[0, ]), "no rows")
  expect_error(mixture_fit(y ~ x1, d$prob, d$nonprob, family = "binomial"), "0/1")
  expect_error(mixture_fit(y ~ x1, d$prob, d$nonprob, family = "poisson"), "non-negative integer")
  expect_error(mixture_fit(y ~ x1, d$prob, d$nonprob, leverage = "zz"), "not in the design")
  expect_error(mixture_fit(y ~ x1, d$prob, d$nonprob, prior = list()), "mixture_prior")
  expect_error(mixture_prior(beta_scale = -1), "beta_scale")
  expect_error(mixture_control(c_infl = 0.5), "greater than 1")
})

test_that("missing values are dropped (na_action = 'omit') or rejected ('fail')", {
  d <- small()
  d$prob$x1[c(2, 5)] <- NA
  d$nonprob$y[3] <- NA
  expect_message(fit <- mixture_fit(y ~ x1, d$prob, d$nonprob, n_iter = 60, burnin = 20, seed = 1),
                 "2 from prob_data, 1 from nonprob_data")
  expect_equal(fit$n[["prob"]], 38)
  expect_equal(fit$n[["nonprob"]], 29)
  expect_length(filter_scores(fit), 29)
  expect_false("3" %in% names(filter_scores(fit)))
  expect_error(mixture_fit(y ~ x1, d$prob, d$nonprob, na_action = "fail"), "Missing values")
  expect_error(mixture_fit_matrix(cbind(1, c(NA, 1:39)), d$prob$y, cbind(1, 1:30),
                                  rnorm(30), leverage = 2), "non-finite")
})

test_that("singular and rank-deficient designs give clear errors", {
  d <- small()
  d$prob$x2 <- 2 * d$prob$x1
  d$nonprob$x2 <- 2 * d$nonprob$x1
  expect_error(mixture_fit(y ~ x1 + x2, d$prob, d$nonprob, model = "M0", n_iter = 50, burnin = 10),
               "rank deficient")
  ## leverage covariates collinear in the P sample
  expect_error(mixture_fit(y ~ x1 + x2, d$prob, d$nonprob, family = "gaussian",
                           n_iter = 50, burnin = 10),
               "rank deficient|singular")
  ## constant continuous covariate in P only -> singular leverage covariance
  d2 <- small(2)
  d2$prob$x3 <- 1.5
  d2$nonprob$x3 <- rnorm(30)
  expect_error(mixture_fit(y ~ x1 + x3, d2$prob, d2$nonprob, family = "gaussian",
                           leverage = c("x1", "x3"), n_iter = 50, burnin = 10),
               "rank deficient|singular")
  ## same for a GLM family, where rank deficiency is only a warning and the
  ## singular leverage covariance must be caught
  d2$prob$y <- rpois(40, 2)
  d2$nonprob$y <- rpois(30, 2)
  expect_error(suppressWarnings(mixture_fit(y ~ x1 + x3, d2$prob, d2$nonprob, family = "poisson",
                                            leverage = c("x1", "x3"), n_iter = 50, burnin = 10)),
               "singular")
  ## no continuous predictor available for Mjoint
  d3 <- small(3)
  d3$prob$g <- factor(rep(c("a", "b"), 20))
  d3$nonprob$g <- factor(rep(c("a", "b"), 15))
  expect_error(mixture_fit(y ~ g, d3$prob, d3$nonprob, n_iter = 50, burnin = 10), "M0")
  expect_s3_class(mixture_fit(y ~ g, d3$prob, d3$nonprob, model = "M0", n_iter = 50, burnin = 10),
                  "bsm_fit")
})

test_that("near-singular but valid designs and extreme predictors still run", {
  set.seed(4)
  x1 <- rnorm(60)
  x2 <- x1 + rnorm(60, sd = 1e-3)  # nearly collinear but full rank
  dp <- data.frame(y = 1 + x1 + rnorm(60), x1 = x1, x2 = x2)
  xn <- rnorm(40)
  dn <- data.frame(y = 1 + xn + rnorm(40), x1 = xn, x2 = xn + rnorm(40, sd = 1e-3))
  fit <- mixture_fit(y ~ x1 + x2, dp, dn, model = "M0", n_iter = 200, burnin = 50, seed = 1)
  expect_true(all(is.finite(as.matrix(fit))))
  ## extreme covariate values in the NP sample (Mahalanobis distances ~ 1e3)
  dn$x1[1:3] <- c(1e3, -1e3, 5e2)
  dn$x2[1:3] <- dn$x1[1:3]
  fit2 <- mixture_fit(y ~ x1, dp, dn, n_iter = 200, burnin = 50, seed = 1)
  expect_true(all(is.finite(fit2$filter_score)))
  expect_true(all(filter_scores(fit2)[1:3] < 1e-6))
})

test_that("small samples work", {
  d <- simulate_mixture_data("gaussian", n_prob = 5, n_nonprob = 3, contamination = 0, seed = 5)
  fit <- mixture_fit(y ~ x1, d$prob, d$nonprob, n_iter = 100, burnin = 20, seed = 1)
  expect_length(filter_scores(fit), 3)
  d1 <- simulate_mixture_data("binomial", n_prob = 8, n_nonprob = 1, contamination = 0, seed = 6)
  fit1 <- mixture_fit(y ~ x1, d1$prob, d1$nonprob, family = "binomial", model = "M1",
                      n_iter = 100, burnin = 20, seed = 1)
  expect_length(filter_scores(fit1), 1)
  ## Gaussian needs more P rows than coefficients
  expect_error(mixture_fit(y ~ x1, d$prob[1:2, ], d$nonprob, model = "M0",
                           n_iter = 50, burnin = 10),
               "more probability-sample")
})

test_that("fixed seeds give reproducible fits and restore the global RNG", {
  d <- small()
  set.seed(999)
  before <- .Random.seed
  f1 <- mixture_fit(y ~ x1, d$prob, d$nonprob, n_iter = 100, burnin = 20, seed = 77)
  expect_identical(.Random.seed, before)
  f2 <- mixture_fit(y ~ x1, d$prob, d$nonprob, n_iter = 100, burnin = 20, seed = 77)
  expect_identical(f1$draws, f2$draws)
  expect_identical(f1$filter_score, f2$filter_score)
  f3 <- mixture_fit(y ~ x1, d$prob, d$nonprob, n_iter = 100, burnin = 20, seed = 78)
  expect_false(identical(f1$draws, f3$draws))
  ## seed = NULL uses (and advances) the global stream
  set.seed(5)
  g1 <- mixture_fit(y ~ x1, d$prob, d$nonprob, n_iter = 50, burnin = 10)
  set.seed(5)
  g2 <- mixture_fit(y ~ x1, d$prob, d$nonprob, n_iter = 50, burnin = 10)
  expect_identical(g1$draws, g2$draws)
})

test_that("M1 'exact' and 'original' updates coincide when the bad set is large", {
  d <- simulate_mixture_data("gaussian", n_prob = 80, n_nonprob = 80, contamination = 0.5,
                             bad_beta = c(-20, 1), seed = 8)
  a <- mixture_fit(y ~ x1, d$prob, d$nonprob, model = "M1", n_iter = 200, burnin = 50, seed = 1)
  b <- mixture_fit(y ~ x1, d$prob, d$nonprob, model = "M1", n_iter = 200, burnin = 50, seed = 1,
                   control = mixture_control(m1_update = "original"))
  ## about 40 units are always allocated to the rival regression (> p + 1),
  ## so both rules perform the same conjugate update
  expect_identical(a$draws, b$draws)
})

test_that("matrix interface handles missing or empty column names", {
  set.seed(9)
  a <- rnorm(50)
  X <- cbind(1, a)                      # colnames c("", "a")
  y <- 1 + a + rnorm(50)
  b <- baseline_fit_matrix(X, y, n_iter = 50, burnin = 10, seed = 1)
  expect_equal(names(coef(b)), c("x1", "a"))
  f <- mixture_fit_matrix(unname(X), y, unname(X[1:20, ]), y[1:20], leverage = 2,
                          n_iter = 50, burnin = 10, seed = 1)
  expect_equal(names(coef(f)), c("x1", "x2"))
})
