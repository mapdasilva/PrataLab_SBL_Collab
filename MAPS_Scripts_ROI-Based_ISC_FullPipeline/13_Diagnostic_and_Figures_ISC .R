#!/usr/bin/env Rscript

library(brms)
library(bayesplot)

## 1) Load fitted model ----
load("ModelA.RData")   # must contain object 'fm'

## 2) Basic summaries + diagnostics ----
cat("\n===== brms summary =====\n")
print(fm)

cat("\n===== Fixed effects (GroupStr etc.) =====\n")
print(fixef(fm))

cat("\n===== Random effects (ROI, Subj1:Subj2, mm) =====\n")
print(VarCorr(fm))

# Extract Rhat and ESS for fixed effects
summ <- summary(fm)
cat("\n===== Rhat for fixed effects =====\n")
print(summ$fixed[,"Rhat"])

cat("\n===== Effective sample size (Eff.Sample) for fixed effects =====\n")
print(summ$fixed[,"Eff.Sample"])

## 3) Posterior predictive check ----
pdf("ModelA_pp_check.pdf", width = 6, height = 4)
pp_check(fm, ndraws = 50)
dev.off()

## 4) Posterior for GroupStr contrast (fixed effect) ----
# Adjust pattern if your GroupStr term is named differently
pars_group <- grep("^b_GroupStr", parnames(fm), value = TRUE)
if (length(pars_group) == 0) {
  stop("No parameter matching '^b_GroupStr' found. Check names(fixef(fm)).")
}
posterior_mat <- as.matrix(fm, pars = pars_group)

pdf("GroupStr_neutral_posterior.pdf", width = 6, height = 4)
mcmc_areas(
  posterior_mat,
  prob       = 0.80,   # inner 80% credible interval
  prob_outer = 0.95    # outer 95% credible interval
)
dev.off()

## 5) ROI-wise random slopes for GroupStr (caterpillar plot) ----
re_roi <- ranef(fm)$ROI  # 3D array: ROI x effect x {Estimate, Est.Error, Q2.5, ...}

# Inspect names once in the log (optional)
cat("\n===== ROI random-effect dimension names =====\n")
print(dimnames(re_roi))

# Pick the first random-effect term that contains "GroupStr"
eff_name <- grep("GroupStr", dimnames(re_roi)[[2]], value = TRUE)[1]
if (is.na(eff_name)) {
  stop("No ROI random slope containing 'GroupStr' found. Check dimnames(ranef(fm)$ROI).")
}

roi_labels <- rownames(re_roi)
roi_est    <- re_roi[, eff_name, "Estimate"]
roi_lower  <- re_roi[, eff_name, "Q2.5"]
roi_upper  <- re_roi[, eff_name, "Q97.5"]

roi_df <- data.frame(
  ROI   = roi_labels,
  est   = roi_est,
  lower = roi_lower,
  upper = roi_upper
)

# Order ROIs by effect size
roi_df <- roi_df[order(roi_df$est), ]

# Caterpillar / forest-style plot
pdf("ROI_GroupStr_caterpillar.pdf", width = 7, height = 8)
par(mar = c(5, 10, 4, 2))  # extra space for ROI labels on y-axis
plot(
  roi_df$est,
  seq_len(nrow(roi_df)),
  xlim = range(c(roi_df$lower, roi_df$upper)),
  yaxt = "n",
  xlab = "GroupStr effect (z-ISC)",
  ylab = "",
  pch  = 16
)
arrows(
  roi_df$lower, seq_len(nrow(roi_df)),
  roi_df$upper, seq_len(nrow(roi_df)),
  angle = 90, code = 3, length = 0.02
)
axis(2, at = seq_len(nrow(roi_df)), labels = roi_df$ROI, las = 2)
abline(v = 0, lty = 2)
title("ROI-wise GroupStr effect (neutral videos)")
dev.off()
