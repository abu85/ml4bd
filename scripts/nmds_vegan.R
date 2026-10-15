# NMDS on the ampliseq data with vegan (Bray-Curtis and Jaccard)

library(phyloseq)
library(vegan)

setwd("C:/Users/Justine/Documents/PhD/Courses/Machine Learning for Life Sciences/project")
data_dir <- "ml4ds_ampliseq_data"
fig_dir  <- "figures"
dir.create(fig_dir, showWarnings = FALSE)

ps <- readRDS(file.path(data_dir, "dada2_phyloseq.rds"))

# ---- Community matrix (samples in rows, ASVs in columns) ----
comm <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) comm <- t(comm)

meta <- as(sample_data(ps), "data.frame")
tax <- as(tax_table(ps), "matrix")

district <- factor(sample_data(ps)$District)
colors   <- hcl.colors(nlevels(district), palette = "Dark 3")

# ---- Distance matrices ----
# Bray-Curtis uses abundances; Jaccard is computed on presence/absence
distances <- list(
  "Bray-Curtis" = vegdist(comm, method = "bray"),
  "Jaccard"     = vegdist(comm, method = "jaccard", binary = TRUE)
)

# ---- Plotting ----
draw_nmds <- function(nmds, dist_name) {
  old_par <- par(mar = c(5, 4, 4, 10))  # extra right margin for the legend
  on.exit(par(old_par))

  ordiplot(nmds, display = "sites", type = "n",
           main = paste("NMDS of", dist_name, "Distances"))
  mtext(sprintf("Stress = %.3f", nmds$stress), side = 3, line = 0.3)
  points(nmds, display = "sites", pch = 19, cex = 1.3,
         col = colors[district])
  legend("topleft", inset = c(1.02, 0), legend = levels(district),
         col = colors, pch = 19, title = "District", bty = "n", xpd = TRUE)
}

# Close any open/stale graphics devices 
graphics.off()

nmds_results <- list()
for (dist_name in names(distances)) {
  set.seed(42)  # NMDS uses random starts; fix the seed for reproducibility
  nmds <- metaMDS(distances[[dist_name]], k = 2, trymax = 100, autotransform = TRUE, trace = FALSE)
  nmds_results[[dist_name]] <- nmds

  file_id <- tolower(gsub("[^A-Za-z]+", "_", dist_name))

  draw_nmds(nmds, dist_name)
  png(file.path(fig_dir, paste0("vegan_nmds_", file_id, ".png")),
      width = 8, height = 6, units = "in", res = 300)
  draw_nmds(nmds, dist_name)
  dev.off()

  stressplot(nmds, main = paste("Shepard diagram -", dist_name))
  png(file.path(fig_dir, paste0("vegan_shepard_", file_id, ".png")),
      width = 7, height = 6, units = "in", res = 300)
  stressplot(nmds, main = paste("Shepard diagram -", dist_name))
  dev.off()
}

# ---- Environmental fitting (envfit) ----
# Fit water-quality variables as vectors onto each NMDS; significance from permutations
env_vars <- c("pH", "DO", "Temp", "Ammonia")
metadata <- read.delim(file.path(data_dir, "metadata.tsv"), row.names = "sampleID")
env <- metadata[rownames(comm), env_vars]  # align to the sample order of the NMDS
stopifnot(!anyNA(env))

draw_envfit <- function(nmds, ef, dist_name) {
  draw_nmds(nmds, dist_name)
  plot(ef, p.max = 1, col = "grey20", cex = 0.9)  # show all vectors, significant or not
}

envfit_results <- list()
for (dist_name in names(nmds_results)) {
  nmds <- nmds_results[[dist_name]]
  set.seed(42)
  ef <- envfit(nmds, env, permutations = 999)
  envfit_results[[dist_name]] <- ef

  cat("\n==== envfit:", dist_name, "====\n")
  print(ef)

  file_id <- tolower(gsub("[^A-Za-z]+", "_", dist_name))
  draw_envfit(nmds, ef, dist_name)
  png(file.path(fig_dir, paste0("vegan_nmds_envfit_", file_id, ".png")),
      width = 8, height = 6, units = "in", res = 300)
  draw_envfit(nmds, ef, dist_name)
  dev.off()
}

# ---- Sub-NMDS without the positive controls (river water) ----
ps_subset <- subset_samples(ps, trimws(District) != "Positive_control")
ps_subset <- prune_taxa(taxa_sums(ps_subset) > 0, ps_subset)  # drop ASVs only in the controls

comm_subset <- as(otu_table(ps_subset), "matrix")
if (taxa_are_rows(ps_subset)) comm_subset <- t(comm_subset)

district_subset <- factor(sample_data(ps_subset)$District)
colors_subset   <- colors[1:nlevels(district_subset)]  # same colours as the full plot
env_subset      <- metadata[rownames(comm_subset), env_vars]

distances_subset <- list(
  "Bray-Curtis" = vegdist(comm_subset, method = "bray"),
  "Jaccard"     = vegdist(comm_subset, method = "jaccard", binary = TRUE)
)

nmds_results_subset   <- list()
envfit_results_subset <- list()
for (dist_name in names(distances_subset)) {
  set.seed(42)
  nmds <- metaMDS(distances_subset[[dist_name]], k = 2, trymax = 100, autotransform = TRUE, trace = FALSE)
  nmds_results_subset[[dist_name]] <- nmds

  ef <- envfit(nmds, env_subset, permutations = 999)
  envfit_results_subset[[dist_name]] <- ef

  print(ef)

  file_id <- tolower(gsub("[^A-Za-z]+", "_", dist_name))

  par(mar = c(5, 4, 4, 10)) 
  ordiplot(nmds, display = "sites", type = "n",
           main = paste("NMDS of", dist_name, "Distances (no controls)"))
  mtext(sprintf("Stress = %.3f", nmds$stress), side = 3, line = 0.3)
  points(nmds, display = "sites", pch = 19, cex = 1.3,
         col = colors_subset[district_subset])
  legend("topleft", inset = c(1.02, 0), legend = levels(district_subset),
         col = colors_subset, pch = 19, title = "District", bty = "n", xpd = TRUE)
  plot(ef, p.max = 1, col = "grey20", cex = 0.9)
  par(mar = c(5, 4, 4, 2) + 0.1)  

  png(file.path(fig_dir, paste0("vegan_nmds_envfit_", file_id, "_no_controls.png")),
      width = 8, height = 6, units = "in", res = 300)
  par(mar = c(5, 4, 4, 10))
  ordiplot(nmds, display = "sites", type = "n",
           main = paste("NMDS of", dist_name, "Distances (no controls)"))
  mtext(sprintf("Stress = %.3f", nmds$stress), side = 3, line = 0.3)
  points(nmds, display = "sites", pch = 19, cex = 1.3,
         col = colors_subset[district_subset])
  legend("topleft", inset = c(1.02, 0), legend = levels(district_subset),
         col = colors_subset, pch = 19, title = "District", bty = "n", xpd = TRUE)
  plot(ef, p.max = 1, col = "grey20", cex = 0.9)
  dev.off()

  stressplot(nmds, main = paste("Shepard diagram -", dist_name, "(no controls)"))
  png(file.path(fig_dir, paste0("vegan_shepard_", file_id, "_no_controls.png")),
      width = 7, height = 6, units = "in", res = 300)
  stressplot(nmds, main = paste("Shepard diagram -", dist_name, "(no controls)"))
  dev.off()
}


# ---- PERMANOVA ----
meta_subset <- as(sample_data(ps_subset), "data.frame")
meta_subset <- cbind(meta_subset, env_subset)  # add pH, DO, Temp, Ammonia

set.seed(42)
permanova_district_bray <- adonis2(distances_subset[["Bray-Curtis"]] ~ District,
                                   data = meta_subset, permutations = 999)
print(permanova_district_bray)

permanova_district_jaccard <- adonis2(distances_subset[["Jaccard"]] ~ District,
                                      data = meta_subset, permutations = 999)
print(permanova_district_jaccard)

dispersion_bray <- betadisper(distances_subset[["Bray-Curtis"]], meta_subset$District)
print(permutest(dispersion_bray, permutations = 999))

dispersion_jaccard <- betadisper(distances_subset[["Jaccard"]], meta_subset$District)
print(permutest(dispersion_jaccard, permutations = 999))

permanova_env_bray <- adonis2(distances_subset[["Bray-Curtis"]] ~ pH + DO + Temp + Ammonia,
                              data = meta_subset, permutations = 999, by = "margin")
print(permanova_env_bray)

permanova_env_jaccard <- adonis2(distances_subset[["Jaccard"]] ~ pH + DO + Temp + Ammonia,
                                 data = meta_subset, permutations = 999, by = "margin")
print(permanova_env_jaccard)


# ---- Batch effect: sequencing run ----
print(table(meta_subset$run, meta_subset$District))

for (dist_name in names(distances_subset)) {
  d <- distances_subset[[dist_name]]
  set.seed(42)

  print(adonis2(d ~ run, data = meta_subset, permutations = 999))

  print(adonis2(d ~ run + District, data = meta_subset, permutations = 999, by = "margin"))

  print(adonis2(d ~ run + pH + DO + Temp + Ammonia, data = meta_subset,
                permutations = 999, by = "margin"))

  within_run <- how(nperm = 999, blocks = meta_subset$run)
  print(adonis2(d ~ District, data = meta_subset, permutations = within_run))

  print(adonis2(d ~ pH + DO + Temp + Ammonia, data = meta_subset,
                permutations = within_run, by = "margin"))

  print(envfit(nmds_results_subset[[dist_name]], env_subset, permutations = within_run))
}
