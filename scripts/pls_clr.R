# CLR transformation + partial least squares (PLS) regression:
# can the ASVs predict the environmental variables (pH, DO, Temp, Ammonia)?

library(phyloseq)
library(mixOmics)

setwd("C:/Users/Justine/Documents/PhD/Courses/Machine Learning for Life Sciences/project")
data_dir    <- "ml4ds_ampliseq_data"
results_dir <- "results"
fig_dir     <- "figures"
dir.create(results_dir, showWarnings = FALSE)
dir.create(fig_dir, showWarnings = FALSE)

ps <- readRDS(file.path(data_dir, "dada2_phyloseq.rds"))

ps_subset <- subset_samples(ps, trimws(District) != "Positive_control")
ps_subset <- prune_taxa(taxa_sums(ps_subset) > 0, ps_subset)  # drop ASVs only in the controls

# ---- Community matrix (samples in rows, ASVs in columns) ----
comm_subset <- as(otu_table(ps_subset), "matrix")
if (taxa_are_rows(ps_subset)) comm_subset <- t(comm_subset)

# ---- Remove rare ASVs ----
# We only have 17 samples, so keep ASVs present in at least 5 samples
# (rare ASVs are mostly zeros and add noise)
min_samples <- 5
keep_asvs   <- colSums(comm_subset > 0) >= min_samples
comm_filtered <- comm_subset[, keep_asvs]
cat("ASVs kept:", ncol(comm_filtered), "of", ncol(comm_subset), "\n")

# ---- CLR transformation ----
# clr(x) = log(x) - mean(log(x)) within each sample.
log_counts <- log(comm_filtered + 1)
clr_data   <- log_counts - rowMeans(log_counts)  # subtract each sample's mean log

# ---- Environmental variables ----
env_vars <- c("pH", "DO", "Temp", "Ammonia")
metadata <- read.delim(file.path(data_dir, "metadata.tsv"), row.names = "sampleID")
env <- metadata[rownames(clr_data), env_vars]  # align to the sample order of clr_data
stopifnot(!anyNA(env))

# District of each sample, with the same colours as the NMDS plots
district <- factor(trimws(metadata[rownames(clr_data), "District"]))
district_colors <- hcl.colors(4, palette = "Dark 3")[1:3]
names(district_colors) <- levels(district)

# ---- Leave-one-out cross-validation ----
max_comp  <- 5
n_samples <- nrow(clr_data)

cv_pred <- array(NA, dim = c(n_samples, length(env_vars), max_comp),
                 dimnames = list(rownames(clr_data), env_vars, paste0("comp", 1:max_comp)))

for (v in env_vars) {
  for (i in 1:n_samples) {
    model_cv <- pls(clr_data[-i, ], env[-i, v], ncomp = max_comp, mode = "regression")
    pred <- predict(model_cv, clr_data[i, , drop = FALSE])$predict  # one value per component
    cv_pred[i, v, ] <- pred[1, 1, ]
  }
}

# Q2 = cross-validated R2: 1 - (prediction error / total variance).
# Q2 near 1 = good predictions; Q2 <= 0 = no better than predicting the mean.
q2 <- matrix(NA, nrow = length(env_vars), ncol = max_comp,
             dimnames = list(env_vars, paste0("comp", 1:max_comp)))
for (v in env_vars) {
  for (k in 1:max_comp) {
    q2[v, k] <- 1 - sum((env[, v] - cv_pred[, v, k])^2) / sum((env[, v] - mean(env[, v]))^2)
  }
}
cat("\n==== Q2 (leave-one-out) by number of PLS components ====\n")
print(round(q2, 3))

# Best number of components for each variable = highest Q2
best_comp <- apply(q2, 1, which.max)
cat("\nBest number of components:\n")
print(best_comp)

# ---- Plot: observed vs cross-validated predictions ----
png(file.path(fig_dir, "pls_clr_observed_vs_predicted.png"),
    width = 9, height = 9, units = "in", res = 300)
par(mfrow = c(2, 2))
for (v in env_vars) {
  k <- best_comp[v]
  observed  <- env[, v]
  predicted <- cv_pred[, v, k]
  axis_range <- range(observed, predicted)
  axis_range <- axis_range + c(-0.05, 0.05) * diff(axis_range)  # a bit of space around the points
  plot(observed, predicted, pch = 19, cex = 1.3,
       col = district_colors[district],  # colour by district
       xlim = axis_range, ylim = axis_range,
       xlab = paste("Observed", v), ylab = paste("Predicted", v, "(leave-one-out)"),
       main = sprintf("%s  (%d comp., Q2 = %.2f)", v, k, q2[v, k]))
  abline(0, 1, lty = 2, col = "grey50")  # perfect prediction line
  # Sample number above each point ("sample7" -> "7"); xpd = TRUE so labels near the edge are not cut off
  text(observed, predicted, labels = sub("sample", "", rownames(env)),
       pos = 3, cex = 0.7, xpd = TRUE)
  legend("topleft", legend = levels(district), col = district_colors, pch = 19,
         title = "District", bty = "n", cex = 0.8)
}
dev.off()

# ---- Final models on all samples + most important ASVs (VIP) ----
n_top <- 10
for (v in env_vars) {
  k <- best_comp[v]
  model_final <- pls(clr_data, env[, v], ncomp = k, mode = "regression")

  vip_scores <- vip(model_final)[, k]
  # Correlation of each ASV (CLR values) with the variable:
  # positive = this ASV is relatively more abundant where v is high
  correlation <- cor(clr_data, env[, v])[, 1]

  vip_table <- data.frame(VIP = vip_scores, correlation = correlation)
  vip_table <- vip_table[order(vip_table$VIP, decreasing = TRUE), ]

  cat("\n==== Top", n_top, "ASVs by VIP:", v, "====\n")
  print(head(vip_table, n_top), digits = 3)

  write.csv(vip_table, file.path(results_dir, paste0("pls_clr_vip_", v, ".csv")))
}

# ---- Plot: top ASVs by VIP for DO and Temp ----
# (only these two models predict better than the mean, Q2 > 0)
n_plot <- 15

# Taxonomy label for each ASV: the deepest rank that has a name, with the rank in front
# (e.g. "species: Flavobacterium columnare", or "family: Methylophilaceae" if the genus is unknown)
tax <- as(tax_table(ps_subset), "matrix")
tax_label <- rep("Unassigned", nrow(tax))
names(tax_label) <- rownames(tax)
for (rank in c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")) {
  assigned <- !is.na(tax[, rank]) & tax[, rank] != ""
  tax_label[assigned] <- paste0(tolower(rank), ": ", tax[assigned, rank])
}
has_species <- !is.na(tax[, "Species"]) & tax[, "Species"] != ""
tax_label[has_species] <- paste0("species: ", tax[has_species, "Genus"], " ", tax[has_species, "Species"])

direction_colors <- c(positive = "#763CD9", negative = "#FD7F4A")  # purple, orange

png(file.path(fig_dir, "pls_clr_vip_DO_Temp.png"), width = 14, height = 6, units = "in", res = 300)
par(mfrow = c(1, 2), mar = c(8, 15, 4, 1))  # 2 panels side by side, wide left margin for names,
                                           # big bottom margin for the legend

for (v in c("DO", "Temp")) {
  vip_table <- read.csv(file.path(results_dir, paste0("pls_clr_vip_", v, ".csv")), row.names = 1)
  top <- head(vip_table, n_plot)

  # Colour by the sign of the correlation with the variable
  direction <- ifelse(top$correlation > 0, "positive", "negative")

  # rev() puts the most important ASV at the top
  barplot(rev(top$VIP), names.arg = rev(tax_label[rownames(top)]),
          horiz = TRUE, las = 1, border = NA,
          col = rev(direction_colors[direction]),
          xlim = c(0, 2.2),  # same scale on both panels
          xlab = "VIP score", main = paste("Top", n_plot, "ASVs for", v))
  abline(v = 1, lty = 2, col = "grey40")  # VIP = 1 threshold
  # Legend below the x-axis (xpd = TRUE allows drawing outside the plot area)
  legend("bottom", inset = c(0, -0.42), xpd = TRUE,
         legend = c(paste("more abundant at high", v), paste("more abundant at low", v)),
         fill = direction_colors, border = NA, bty = "n", cex = 0.9)
}
dev.off()


# ---- Batch effect check: PLS on run A samples only ----
run <- metadata[rownames(comm_subset), "run"]
comm_runA <- comm_subset[run == "A", ]
cat("\n\n######## Run A only:", nrow(comm_runA), "samples ########\n")

# Same steps as above: remove rare ASVs (here: present in at least 3 of the 11 samples), then CLR
keep_asvs_runA <- colSums(comm_runA > 0) >= 3
comm_runA <- comm_runA[, keep_asvs_runA]
cat("ASVs kept:", ncol(comm_runA), "\n")

log_counts_runA <- log(comm_runA + 1)
clr_runA <- log_counts_runA - rowMeans(log_counts_runA)

env_runA      <- metadata[rownames(clr_runA), env_vars]
district_runA <- factor(trimws(metadata[rownames(clr_runA), "District"]), levels = levels(district))

# Leave-one-out cross-validation, as above
n_runA <- nrow(clr_runA)
cv_pred_runA <- array(NA, dim = c(n_runA, length(env_vars), max_comp),
                      dimnames = list(rownames(clr_runA), env_vars, paste0("comp", 1:max_comp)))
for (v in env_vars) {
  for (i in 1:n_runA) {
    model_cv <- pls(clr_runA[-i, ], env_runA[-i, v], ncomp = max_comp, mode = "regression")
    pred <- predict(model_cv, clr_runA[i, , drop = FALSE])$predict
    cv_pred_runA[i, v, ] <- pred[1, 1, ]
  }
}

q2_runA <- matrix(NA, nrow = length(env_vars), ncol = max_comp,
                  dimnames = list(env_vars, paste0("comp", 1:max_comp)))
for (v in env_vars) {
  for (k in 1:max_comp) {
    q2_runA[v, k] <- 1 - sum((env_runA[, v] - cv_pred_runA[, v, k])^2) /
                         sum((env_runA[, v] - mean(env_runA[, v]))^2)
  }
}
cat("\n==== Q2 (leave-one-out), run A only ====\n")
print(round(q2_runA, 3))
best_comp_runA <- apply(q2_runA, 1, which.max)

# Observed vs predicted, run A only (same layout and colours as the plot above)
png(file.path(fig_dir, "pls_clr_observed_vs_predicted_runA.png"),
    width = 9, height = 9, units = "in", res = 300)
par(mfrow = c(2, 2))
for (v in env_vars) {
  k <- best_comp_runA[v]
  observed  <- env_runA[, v]
  predicted <- cv_pred_runA[, v, k]
  axis_range <- range(observed, predicted)
  axis_range <- axis_range + c(-0.05, 0.05) * diff(axis_range)
  plot(observed, predicted, pch = 19, cex = 1.3,
       col = district_colors[district_runA],
       xlim = axis_range, ylim = axis_range,
       xlab = paste("Observed", v), ylab = paste("Predicted", v, "(leave-one-out)"),
       main = sprintf("%s, run A only  (%d comp., Q2 = %.2f)", v, k, q2_runA[v, k]))
  abline(0, 1, lty = 2, col = "grey50")
  text(observed, predicted, labels = sub("sample", "", rownames(env_runA)),
       pos = 3, cex = 0.7, xpd = TRUE)
  legend("topleft", legend = levels(district), col = district_colors, pch = 19,
         title = "District", bty = "n", cex = 0.8)
}
dev.off()


n_top <- 10
for (v in env_vars) {
  k <- best_comp_runA[v]
  model_final_runA <- pls(clr_runA, env_runA[, v], ncomp = k, mode = "regression")

  vip_scores_runA <- vip(model_final_runA)[, k]
  # Correlation of each ASV (CLR values) with the variable:
  # positive = this ASV is relatively more abundant where v is high
  correlation <- cor(clr_runA, env_runA[, v])[, 1]

  vip_table_runA <- data.frame(VIP = vip_scores_runA, correlation = correlation)
  vip_table_runA <- vip_table_runA[order(vip_table_runA$VIP, decreasing = TRUE), ]

  cat("\n==== Top", n_top, "ASVs by VIP:", v, "====\n")
  print(head(vip_table_runA, n_top), digits = 3)

  write.csv(vip_table_runA, file.path(results_dir, paste0("pls_clr_runA_vip_", v, ".csv")))
}

png(file.path(fig_dir, "pls_clr_vip_runA_Temp.png"), width = 14, height = 6, units = "in", res = 300)
par(mar = c(8, 17, 4, 1))  # wide left margin for the taxonomy names,
                           # big bottom margin for the legend

for (v in c("Temp")) {
  vip_table_runA <- read.csv(file.path(results_dir, paste0("pls_clr_runA_vip_", v, ".csv")), row.names = 1)
  top <- head(vip_table_runA, n_plot)

  # Colour by the sign of the correlation with the variable
  direction <- ifelse(top$correlation > 0, "positive", "negative")

  # rev() puts the most important ASV at the top
  barplot(rev(top$VIP), names.arg = rev(tax_label[rownames(top)]),
          horiz = TRUE, las = 1, border = NA,
          col = rev(direction_colors[direction]),
          xlim = c(0, 2.2),  # same scale on both panels
          xlab = "VIP score", main = paste("Top", n_plot, "ASVs for", v, "run A only"))
  abline(v = 1, lty = 2, col = "grey40")  # VIP = 1 threshold
  # Legend below the x-axis (xpd = TRUE allows drawing outside the plot area)
  legend("bottom", inset = c(0, -0.42), xpd = TRUE,
         legend = c(paste("more abundant at high", v), paste("more abundant at low", v)),
         fill = direction_colors, border = NA, bty = "n", cex = 0.9)
}
dev.off()
