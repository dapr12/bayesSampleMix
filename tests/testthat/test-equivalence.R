## The package must reproduce the original research code draw for draw.
## The original function definitions are stored verbatim in fixtures/ and
## are sourced into private environments. Identical seeds must give identical
## posterior summaries (up to floating-point noise), which shows that the
## refactoring did not change the sampler, the priors or the initialisation.

load_original <- function(file) {
  env <- new.env(parent = globalenv())
  sys.source(test_path("fixtures", file), envir = env)
  env
}

gauss_data <- function(seed = 11) {
  set.seed(seed)
  Xp <- cbind(1, rnorm(80), rnorm(80))
  yp <- drop(Xp %*% c(1, 2, -1)) + rnorm(80)
  Xnp <- cbind(1, rnorm(70, 0.5), rnorm(70))
  ynp <- drop(Xnp %*% c(1, 2, -1)) + rnorm(70)
  ynp[1:20] <- rnorm(20, 8, 2)
  list(Xp = Xp, yp = yp, Xnp = Xnp, ynp = ynp)
}

test_that("Gaussian sampler reproduces the original fit_mixture() exactly", {
  orig <- load_original("original_gaussian.R")
  d <- gauss_data()
  for (m in c("M0", "Mjoint", "M1")) {
    o <- orig$fit_mixture(d$Xp, d$yp, d$Xnp, d$ynp, cont_idx = 2:3, model = m,
                          iters = 300, burn = 100, seed = 5)
    f <- mixture_fit_matrix(d$Xp, d$yp, d$Xnp, d$ynp, family = "gaussian", model = m,
                            leverage = 2:3, n_iter = 300, burnin = 100, seed = 5,
                            control = mixture_control(m1_update = "original"))
    expect_equal(unname(f$draws[, 1:3]), unname(o$beta_draws), tolerance = 1e-10, info = m)
    expect_equal(mean(f$draws[, "sigma2"]), o$sigma2_mean, tolerance = 1e-10, info = m)
    expect_equal(mean(f$draws[, "contamination"]), o$contam_mean, tolerance = 1e-10, info = m)
    expect_equal(f$filter_score, o$filter_score, tolerance = 1e-10, info = m)
    expect_equal(f$indicator_entropy, o$indicator_entropy, tolerance = 1e-10, info = m)
  }
})

test_that("Gaussian baseline reproduces the original fit_bayes_lm() exactly", {
  orig <- load_original("original_gaussian.R")
  d <- gauss_data()
  o <- orig$fit_bayes_lm(d$Xp, d$yp, iters = 300, burn = 100, seed = 3)
  f <- baseline_fit_matrix(d$Xp, d$yp, "gaussian", n_iter = 300, burnin = 100, seed = 3)
  expect_equal(unname(coef(f)), o$beta_mean, tolerance = 1e-10)
  expect_equal(mean(f$draws[, "sigma2"]), o$sigma2_mean, tolerance = 1e-10)
})

test_that("overlap_diagnostics() reproduces the original stage0()", {
  orig <- load_original("original_gaussian.R")
  d <- gauss_data()
  o <- orig$stage0(d$Xp, d$Xnp, 2:3)
  f <- overlap_diagnostics(d$Xp[, 2:3], d$Xnp[, 2:3])
  expect_equal(f$OC, o$OC, tolerance = 1e-12)
  expect_equal(f$EF, o$EF, tolerance = 1e-12)
  expect_equal(f$AUC_PS, o$AUC_PS, tolerance = 1e-12)
  scores <- c(0.9, 0.8, 0.8, 0.3, 0.2, 0.8)
  labs <- c(1, 1, 0, 0, 0, 1)
  expect_equal(filter_auc(labs, scores), orig$auc_score(labs, scores))
})

test_that("Logistic sampler reproduces the original fit_logit_mixture() exactly", {
  orig <- load_original("original_logistic.R")
  set.seed(21)
  Xp <- cbind(1, matrix(rnorm(160), 80, 2))
  yp <- rbinom(80, 1, plogis(drop(Xp %*% c(0.2, 1, -0.8))))
  Xnp <- cbind(1, matrix(rnorm(140, 0.5), 70, 2))
  ynp <- rbinom(70, 1, plogis(drop(Xnp %*% c(0.2, -1, 0.8))))
  for (m in c("M0", "Mjoint", "M1")) {
    o <- orig$fit_logit_mixture(Xp, yp, Xnp, ynp, cont_idx = 2:3, model = m,
                                iters = 200, burn = 80, seed = 9, K = 32)
    f <- mixture_fit_matrix(Xp, yp, Xnp, ynp, family = "binomial", model = m,
                            leverage = 2:3, n_iter = 200, burnin = 80, seed = 9,
                            control = mixture_control(pg_terms = 32, m1_update = "original"))
    expect_equal(unname(coef(f)), o$beta_mean, tolerance = 1e-10, info = m)
    expect_equal(mean(f$draws[, "contamination"]), o$contam_mean, tolerance = 1e-10, info = m)
    expect_equal(f$filter_score, o$filter_score, tolerance = 1e-10, info = m)
  }
  ob <- orig$fit_bayes_logit(Xp, yp, iters = 200, burn = 80, seed = 2, K = 32)
  fb <- baseline_fit_matrix(Xp, yp, "binomial", n_iter = 200, burnin = 80, seed = 2,
                            control = mixture_control(pg_terms = 32))
  expect_equal(unname(coef(fb)), ob$beta_mean, tolerance = 1e-10)
})

test_that("Poisson sampler reproduces the original fit_pois_mixture() exactly", {
  orig <- load_original("original_poisson.R")
  set.seed(31)
  Xp <- cbind(1, matrix(rnorm(160), 80, 2))
  yp <- rpois(80, exp(drop(Xp %*% c(1.2, 0.5, -0.4))))
  Xnp <- cbind(1, matrix(rnorm(140, 0.5), 70, 2))
  ynp <- rpois(70, exp(drop(Xnp %*% c(1.2, -0.5, 0.4))))
  for (m in c("M0", "Mjoint", "M1")) {
    o <- orig$fit_pois_mixture(Xp, yp, Xnp, ynp, cont_idx = 2:3, model = m,
                               iters = 200, burn = 80, seed = 4)
    f <- mixture_fit_matrix(Xp, yp, Xnp, ynp, family = "poisson", model = m,
                            leverage = 2:3, n_iter = 200, burnin = 80, seed = 4,
                            control = mixture_control(m1_update = "original"))
    expect_equal(unname(coef(f)), o$beta_mean, tolerance = 1e-10, info = m)
    expect_equal(mean(f$draws[, "contamination"]), o$contam_mean, tolerance = 1e-10, info = m)
    expect_equal(f$filter_score, o$filter_score, tolerance = 1e-10, info = m)
  }
  ob <- orig$fit_bayes_pois(Xp, yp, iters = 200, burn = 80, seed = 2)
  fb <- baseline_fit_matrix(Xp, yp, "poisson", n_iter = 200, burnin = 80, seed = 2)
  expect_equal(unname(coef(fb)), ob$beta_mean, tolerance = 1e-10)
})

test_that("thinning keeps a subset of the unthinned chain", {
  d <- gauss_data()
  f1 <- mixture_fit_matrix(d$Xp, d$yp, d$Xnp, d$ynp, leverage = 2:3,
                           n_iter = 300, burnin = 100, thin = 1, seed = 8)
  f5 <- mixture_fit_matrix(d$Xp, d$yp, d$Xnp, d$ynp, leverage = 2:3,
                           n_iter = 300, burnin = 100, thin = 5, seed = 8)
  expect_equal(nrow(f5$draws), 40L)
  expect_equal(f5$draws, f1$draws[seq(5, 200, by = 5), ])
  ## filter scores average every post-burn-in iteration, independent of thin
  expect_equal(f5$filter_score, f1$filter_score)
})
