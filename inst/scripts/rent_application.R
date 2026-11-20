## Empirical illustration: Munich Rent Index 1999 (rent99, package gamlss.data).
## Port of the original rent_application.R to the package API.
##
## As in the original script, average-location flats are treated as the trusted
## probability sample (P) and good-location flats as the "non-probability"
## sample (NP). This is a covariate-shift illustration (good-location flats are
## newer), NOT a genuine self-selected NP source; interpret accordingly.
##
##   install.packages("gamlss.data")
##   Rscript rent_application.R [outdir]

library(bayesSampleMix)
if (!requireNamespace("gamlss.data", quietly = TRUE)) {
  stop("Please install.packages('gamlss.data') to run this script.")
}
args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args)) args[1] else tempdir()

data("rent99", package = "gamlss.data")
df <- rent99
loc <- as.integer(as.character(df$location))
df <- df[loc %in% c(1, 2), ]
loc <- loc[loc %in% c(1, 2)]
Pdat <- df[loc == 1, ]
NPdat <- df[loc == 2, ]
cat(sprintf("nP (average) = %d   nNP (good) = %d\n", nrow(Pdat), nrow(NPdat)))

## covariates standardised by P-sample moments, as in the original analysis
tob <- function(x) as.numeric(as.character(x))
mA <- mean(Pdat$area); sA <- sd(Pdat$area)
mY <- mean(Pdat$yearc); sY <- sd(Pdat$yearc)
prep <- function(d) data.frame(
  logrent = log(d$rent),
  area_z = (d$area - mA) / sA,
  yearc_z = (d$yearc - mY) / sY,
  bath = tob(d$bath), kitchen = tob(d$kitchen), cheating = tob(d$cheating)
)
P <- prep(Pdat)
NP <- prep(NPdat)
f <- logrent ~ area_z + yearc_z + bath + kitchen + cheating

## Stage 0: overlap
print(overlap_diagnostics(P[, c("area_z", "yearc_z")], NP[, c("area_z", "yearc_z")]))

## Mixture fit (leverage covariates: area_z, yearc_z - selected automatically)
fit <- mixture_fit(f, P, NP, family = "gaussian", model = "Mjoint",
                   n_iter = 6000, burnin = 3000, seed = 1)
print(summary(fit))

## Baselines
po <- baseline_fit(f, P, n_iter = 4000, burnin = 2000, seed = 1)
nv <- baseline_fit(f, rbind(P, NP), n_iter = 4000, burnin = 2000, seed = 1)
tab <- data.frame(ProbOnly = coef(po), NaivePool = coef(nv), Mixture = coef(fit),
                  Mixture_sd = apply(as.matrix(fit)[, fit$coef_names], 2, sd))
print(round(tab, 3))

out <- cbind(NPdat, filter_score = filter_scores(fit),
             indicator_entropy = fit$indicator_entropy)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
write.csv(out, file.path(outdir, "rent99_np_scored.csv"), row.names = FALSE)
cat("scored NP observations written to", file.path(outdir, "rent99_np_scored.csv"), "\n")
