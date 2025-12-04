#!/usr/bin/env Rscript

# Step 14: Turn ROI effects into a Yeo7-based NIfTI brain map

library(brms)
library(RNifti)

# 1) Load fitted model (must contain object 'fm')
load("ModelA.RData")

# 2) Extract ROI-wise random slope for GroupStr -------------------------

re_roi <- ranef(fm)$ROI  # 3D array: ROI x effect x {Estimate, Est.Error, Q2.5, Q97.5, ...}

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

# 3) Optional: mask to ROIs with 95% CI not crossing 0 -----------------

roi_df$signif  <- with(roi_df, lower * upper > 0)   # TRUE if CI doesn't cross 0
roi_df$est_sig <- ifelse(roi_df$signif, roi_df$est, 0)

# Save table for inspection
write.csv(roi_df, "ROI_GroupStr_effects.csv", row.names = FALSE)

# 4) Load Yeo7 atlas and ROI label table -------------------------------

# Assumes:
#  - "Yeo7resampled.nii" is the atlas NIfTI (integer labels)
#  - "Yeo7.txt" has at least columns:
#        Name  = ROI name (matches roi_df$ROI)
#        Label = integer label in the atlas

atlas_nii <- readNifti("Yeo7resampled.nii")
atlas_vec <- as.vector(atlas_nii)

yeo_info <- read.table("Yeo7.txt", header = TRUE, sep = "\t")

merge_df <- merge(roi_df, yeo_info,
                  by.x = "ROI", by.y = "Name", all.x = TRUE)

if (any(is.na(merge_df$Label))) {
  warning("Some ROIs in roi_df did not find a matching Label in Yeo7.txt")
}

# 5) Create NIfTI map of GroupStr effect -------------------------------

beta_map <- rep(0, length(atlas_vec))

for (i in seq_len(nrow(merge_df))) {
  lab <- merge_df$Label[i]
  val <- merge_df$est_sig[i]   # use est for unmasked effects
  if (!is.na(lab)) {
    beta_map[atlas_vec == lab] <- val
  }
}

beta_array <- array(beta_map, dim = dim(atlas_nii))
beta_nii   <- updateNifti(beta_array, atlas_nii)

writeNifti(beta_nii, "GroupStr_effect_Yeo7.nii.gz")

cat("Saved GroupStr_effect_Yeo7.nii.gz\n")
