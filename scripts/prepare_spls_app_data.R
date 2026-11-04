## =====================================================================
## prepare_spls_app_data.R
## Run ONCE, outside the Shiny app. Fits one sPLS1 model per
## environmental variable and saves only what the two app plots need:
##   (1) LOO-CV observed vs predicted
##   (2) top driver ASVs (Order-Genus)
## The app then just readRDS()es the result: no model fitting in Shiny.
##
## Runtime: 10-40 min with nperm = 100 (about a minute with nperm = 0).
## This only matters once: the app stays instant because it only reads
## the saved file.
## =====================================================================

## ---- Settings ---------------------------------------------------------
# load data - used for all plots
ps = readRDS("/Users/lovisalindquist/Documents/ML4LSProject/Data/dada2_phyloseq.rds")
meta = read.table("/Users/lovisalindquist/Documents/ML4LSProject/Data/metadata.tsv", header = T, row.names = "sampleID")

out_file  <- "spls_app_data.rds"          # the only file the app needs

batch_col  <- "run"
group_col  <- "District"
env_vars   <- c("pH", "DO", "Temp", "Ammonia")
bc_method  <- "limma"                     # "limma" (correct for run) or "none"
keepX_grid <- c(5, 10, 15, 20, 30, 50, 75, 100)
n_top      <- 20                          # driver ASVs kept for the bar plot
nperm      <- 100                         # permutations per variable (0 = skip test, p = NA)

suppressPackageStartupMessages({
  library(phyloseq)
  library(limma)
  library(mixOmics)
  library(permute)
})
set.seed(5249)


## ---- 1. Data ----------------------------------------------------------
ps <- readRDS(ps_path)
md <- as.data.frame(readr::read_delim(meta_path, delim = "\t",
                                      escape_double = FALSE, trim_ws = TRUE))
rownames(md) <- md$sampleID
sample_data(ps) <- sample_data(md)
ps <- subset_samples(ps, Group %in% c("1", "2", "3"))   # drop positive control
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
taxa_names(ps) <- paste0("ASV", seq_len(ntaxa(ps)))     # short IDs

X_raw <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps)) X_raw <- t(X_raw)                # samples as rows

meta <- data.frame(sample_data(ps), check.names = FALSE)[rownames(X_raw), , drop = FALSE]
meta[[batch_col]] <- droplevels(factor(meta[[batch_col]]))
meta[[group_col]] <- droplevels(factor(meta[[group_col]]))
for (v in env_vars) if (!is.numeric(meta[[v]]))         # "7,2" -> 7.2
  meta[[v]] <- as.numeric(gsub(",", ".", as.character(meta[[v]])))


## ---- 2. "Order-Genus" labels -----------------------------------------
tax <- as.data.frame(as(tax_table(ps), "matrix"), stringsAsFactors = FALSE)
names(tax) <- tools::toTitleCase(tolower(names(tax)))
clean <- function(v) {
  v <- trimws(sub("^[a-zA-Z]__", "", as.character(v)))
  v[is.na(v) | v == "" |
      grepl("^(uncultured|unclassified|unknown|metagenome|NA$)", v, ignore.case = TRUE)] <- NA
  v
}
getr <- function(r) if (r %in% names(tax)) clean(tax[[r]]) else rep(NA_character_, nrow(tax))
ord <- ifelse(!is.na(getr("Order")), getr("Order"),
              ifelse(!is.na(getr("Class")), paste0("c_", getr("Class")), "Unclassified"))
gen <- ifelse(!is.na(getr("Genus")), getr("Genus"),
              ifelse(!is.na(getr("Family")), paste0("f_", getr("Family")), "unclassified"))
tax$OrderGenus <- paste(ord, gen, sep = "-")            # e.g. "Bacillales-Bacillus"


## ---- 3. Filter rare ASVs + CLR ---------------------------------------
keep    <- colSums(X_raw) / sum(X_raw) > 1e-4 & colSums(X_raw > 0) >= 2
lx      <- log(X_raw[, keep] + 1)
clr_raw <- lx - rowMeans(lx)


## ---- 4. Helpers -------------------------------------------------------
## Learn the run effect on (training) samples, protecting y ...
fit_batch <- function(X, batch, y) {
  batch <- droplevels(factor(batch))
  if (bc_method == "none" || nlevels(batch) < 2) return(NULL)
  B <- outer(as.character(batch), levels(batch)[-1], "==") * 1
  colnames(B) <- paste0("run_", levels(batch)[-1])
  full <- cbind(model.matrix(~ y), B)
  if (qr(full)$rank < ncol(full)) stop("Run effect not estimable (no bridge samples).")
  coef <- limma::lmFit(t(X), full)$coefficients[, colnames(B), drop = FALSE]
  list(levels = levels(batch), coef = coef)
}

## ... and subtract it from any sample according to its run
apply_batch <- function(X, batch, bf) {
  if (is.null(bf)) return(X)
  b <- as.character(batch)
  stopifnot(all(b %in% bf$levels))
  X - (outer(b, bf$levels[-1], "==") * 1) %*% t(bf$coef)
}

## Leave-one-out folds, run correction learned inside each fold (no leakage)
loo_folds <- function(X, y, batch) lapply(seq_along(y), function(i) {
  bf  <- fit_batch(X[-i, ], batch[-i], y[-i])
  Xtr <- apply_batch(X[-i, ], batch[-i], bf)
  Xte <- apply_batch(X[i, , drop = FALSE], batch[i], bf)
  ok  <- apply(Xtr, 2, sd) > 1e-8                       # drop constant ASVs
  list(tr = -i, Xtr = Xtr[, ok], Xte = Xte[, ok, drop = FALSE])
})

## LOO predictions and Q2 for one ncomp/keepX
cv_spls1 <- function(folds, y, keepX, ncomp) {
  pred <- sapply(folds, function(f) {
    m <- mixOmics::spls(f$Xtr, y[f$tr], ncomp = ncomp,
                        keepX = rep(min(keepX, ncol(f$Xtr)), ncomp), mode = "regression")
    predict(m, newdata = f$Xte)$predict[, 1, ncomp]
  })
  list(Q2 = 1 - sum((y - pred)^2) / sum((y - mean(y))^2), pred = pred)
}

## Best ncomp/keepX by Q2
tune_spls1 <- function(folds, y) {
  grid <- expand.grid(ncomp = 1:2,
                      keepX = keepX_grid[keepX_grid <= ncol(folds[[1]]$Xtr)])
  grid$Q2 <- mapply(function(nc, kx) cv_spls1(folds, y, kx, nc)$Q2, grid$ncomp, grid$keepX)
  grid[which.max(grid$Q2), ]
}


## ---- 5. Fit each variable and keep only what the plots need ----------
spls_app <- lapply(env_vars, function(v) {
  message("sPLS1: ", v)
  cc <- !is.na(meta[[v]])
  y  <- meta[[v]][cc]
  X  <- clr_raw[cc, ]
  b  <- droplevels(meta[[batch_col]][cc])

  folds <- loo_folds(X, y, b)
  best  <- tune_spls1(folds, y)
  cv    <- cv_spls1(folds, y, best$keepX, best$ncomp)

  # Permutation test (as in the original analysis): shuffle y within run,
  # redo the in-fold correction AND the tuning, p = share of null Q2 >= observed
  p <- NA; null <- NULL                                 # only kept if nperm = 0
  if (nperm > 0) {
    null <- replicate(nperm, {
      yp <- y[permute::shuffle(length(y), control = how(blocks = b))]
      tune_spls1(loo_folds(X, yp, b), yp)$Q2
    })
    p <- (sum(null >= best$Q2) + 1) / (nperm + 1)
  }

  # Final model on all samples -> driver ASVs
  Xc  <- apply_batch(X, b, fit_batch(X, b, y))
  fit <- mixOmics::spls(Xc, y, ncomp = best$ncomp,
                        keepX = rep(best$keepX, best$ncomp), mode = "regression")
  sel <- mixOmics::selectVar(fit, comp = 1)$X$value
  asv <- rownames(sel)
  drv <- data.frame(ASV        = asv,
                    OrderGenus = tax[asv, "OrderGenus"],
                    loading    = sel$value.var,
                    spearman   = sapply(asv, function(a) cor(Xc[, a], y, method = "spearman")),
                    row.names  = NULL)
  drv <- head(drv[order(-abs(drv$loading)), ], n_top)
  drv$Label     <- paste0(drv$OrderGenus, " (", drv$ASV, ")")
  drv$direction <- ifelse(drv$spearman > 0, "positive", "negative")

  list(variable = v, Q2 = best$Q2, p = p, null = null,
       ncomp = best$ncomp, keepX = best$keepX,
       pred = data.frame(observed = y, predicted = cv$pred,
                         meta[cc, , drop = FALSE], check.names = FALSE),  # all metadata -> any colour-by
       drivers = drv)
})
names(spls_app) <- env_vars

saveRDS(spls_app, out_file)
message("Saved ", out_file)
