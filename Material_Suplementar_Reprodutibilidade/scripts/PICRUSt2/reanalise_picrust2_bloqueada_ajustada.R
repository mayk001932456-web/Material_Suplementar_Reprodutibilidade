# Independent corrected reanalysis of bacteria-only PICRUSt2 KO profiles.
# Model: cultivar + sequencing batch. Cultivar is tested conditional on batch.
# All distinct cultivar assignments preserving two samples per cultivar within
# each batch are enumerated (2,520 x 2,520 = 6,350,400 assignments).

.script_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(.script_args)) stop("Run this script with Rscript so its location can be resolved.")
.script_file <- normalizePath(sub("^--file=", "", .script_args[1]), winslash = "/", mustWork = TRUE)
source(file.path(dirname(dirname(.script_file)), "package_config.R"))
# This script does not use vegan and writes only to a new output directory.

valid <- INPUTS_ROOT
pic <- file.path(INPUTS_ROOT, "picrust2")
outdir <- file.path(RESULTS_ROOT, "PICRUSt2_KO")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

meta <- read.delim(file.path(valid, "metadados_16_raizes.tsv"),
                   stringsAsFactors = FALSE, check.names = FALSE)
required_meta <- c("sample_id", "sequencing_batch", "cultivar", "role")
stopifnot(all(required_meta %in% names(meta)))
meta <- meta[meta$role == "raiz", , drop = FALSE]
meta$cultivar <- factor(meta$cultivar,
                        levels = c("BL-65", "BL-69", "BL-74", "BL-183"))
meta$sequencing_batch <- factor(meta$sequencing_batch,
                                levels = c("Lote 1", "Lote 2"))
stopifnot(nrow(meta) == 16L, !anyNA(meta$cultivar), !anyNA(meta$sequencing_batch),
          !anyDuplicated(meta$sample_id), all(table(meta$sequencing_batch, meta$cultivar) == 2L))

# Metadata sample_id carries L1_/L2_; PICRUSt2 columns use Go IDs only.
ko_path <- file.path(pic, "KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz")
ko0 <- read.delim(gzfile(ko_path), check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(ncol(ko0) >= 17L)
pic_ids <- sub("^L[12]_", "", meta$sample_id)
ko_sample_cols <- names(ko0)[-(1:2)]
pic_sample_cols <- sub("^X", "", ko_sample_cols)
idx <- match(pic_ids, pic_sample_cols)
stopifnot(!anyNA(idx), !anyDuplicated(pic_ids), !anyDuplicated(pic_sample_cols),
          all(grepl("^Go[0-9]+$", pic_ids)))
ko <- as.matrix(ko0[, idx + 2L, drop = FALSE])
storage.mode(ko) <- "numeric"
stopifnot(all(is.finite(ko)), all(ko >= 0), all(colSums(ko) > 0))
rownames(ko) <- ko0[[1L]]
colnames(ko) <- pic_ids
ko <- sweep(ko, 2L, colSums(ko), "/")
X <- t(ko)                         # samples x predicted KOs
n <- nrow(X)
stopifnot(n == 16L, ncol(X) == nrow(ko), all(abs(rowSums(X) - 1) < 1e-10))

# Bray-Curtis on relative profiles; for rows summing to one this is Manhattan/2.
D <- as.matrix(stats::dist(X, method = "manhattan")) / 2
stopifnot(all(is.finite(D)), max(abs(D - t(D))) < 1e-12,
          max(abs(diag(D))) < 1e-12, all(D >= 0), all(D <= 1 + 1e-12))
J <- diag(n) - matrix(1 / n, n, n)
G <- -0.5 * J %*% (D^2) %*% J
G <- (G + t(G)) / 2
ss_total <- sum(diag(G))
stopifnot(is.finite(ss_total), ss_total > 0)

# Build all distinct assignments to four cultivar labels in one batch.
make_assignments <- function() {
  out <- matrix(0L, 2520L, 8L)
  r <- 0L
  for (a in combn(8L, 2L, simplify = FALSE)) {
    rem1 <- setdiff(seq_len(8L), a)
    for (b in combn(rem1, 2L, simplify = FALSE)) {
      rem2 <- setdiff(rem1, b)
      for (cc in combn(rem2, 2L, simplify = FALSE)) {
        rem3 <- setdiff(rem2, cc)
        r <- r + 1L
        out[r, a] <- 1L
        out[r, b] <- 2L
        out[r, cc] <- 3L
        out[r, rem3] <- 4L
      }
    }
  }
  stopifnot(r == 2520L,
            all(rowSums(out == 1L) == 2L), all(rowSums(out == 2L) == 2L),
            all(rowSums(out == 3L) == 2L), all(rowSums(out == 4L) == 2L),
            !anyDuplicated(as.data.frame(out)))
  out
}
A <- make_assignments()
i1 <- which(meta$sequencing_batch == "Lote 1")
i2 <- which(meta$sequencing_batch == "Lote 2")
stopifnot(length(i1) == 8L, length(i2) == 8L)

# Model-based observed partition: cultivar SS is the partial SS added after
# batch. QR projection computes trace(P G) and avoids relying on balanced-design
# shortcuts. Residual SS and df come from the full additive model.
X_batch <- stats::model.matrix(~ sequencing_batch, data = meta)
X_full <- stats::model.matrix(~ sequencing_batch + cultivar, data = meta)
stopifnot(qr(X_batch)$rank == 2L, qr(X_full)$rank == 5L)
Q_batch <- qr.Q(qr(X_batch))[, seq_len(qr(X_batch)$rank), drop = FALSE]
Q_full <- qr.Q(qr(X_full))[, seq_len(qr(X_full)$rank), drop = FALSE]
ss_batch <- sum(diag(crossprod(Q_batch, G %*% Q_batch)))
ss_full_obs <- sum(diag(crossprod(Q_full, G %*% Q_full)))
ss_cult_obs <- ss_full_obs - ss_batch
ss_resid_obs <- ss_total - ss_full_obs
df_cult <- qr(X_full)$rank - qr(X_batch)$rank
df_resid <- n - qr(X_full)$rank
stopifnot(df_cult == 3L, df_resid == 11L,
          is.finite(ss_cult_obs), ss_cult_obs > 0,
          is.finite(ss_resid_obs), ss_resid_obs > 0)
F_obs <- (ss_cult_obs / df_cult) / (ss_resid_obs / df_resid)
partial_R2 <- ss_cult_obs / (ss_cult_obs + ss_resid_obs)
R2_total <- ss_cult_obs / ss_total

# For each within-batch reassignment, SS(cultivar | batch) is computed from
# trace(P_full G)-trace(P_batch G). Since every cultivar has two samples in
# both batches, cultivar and batch contrasts are orthogonal. The following
# quadratic-form expression is exactly that partial SS for every assignment.
G1 <- G[i1, i1, drop = FALSE]
G2 <- G[i2, i2, drop = FALSE]
G12 <- G[i1, i2, drop = FALSE]
within_score <- function(Bmat) {
  v <- numeric(nrow(A))
  for (k in seq_len(4L)) {
    M <- 1.0 * (A == k)
    v <- v + rowSums((M %*% Bmat) * M) / 4
  }
  v
}
w1 <- within_score(G1)
w2 <- within_score(G2)
cross <- matrix(0, nrow(A), nrow(A))
for (k in seq_len(4L)) {
  M <- 1.0 * (A == k)
  cross <- cross + (M %*% G12) %*% t(M)
}
ss_cult_perm <- outer(w1, w2, "+") + 0.5 * cross
ss_resid_perm <- ss_total - ss_batch - ss_cult_perm
stopifnot(length(ss_cult_perm) == 6350400L,
          all(is.finite(ss_cult_perm)), all(is.finite(ss_resid_perm)),
          all(ss_cult_perm > 0), all(ss_resid_perm > 0))
F_perm <- (ss_cult_perm / df_cult) / (ss_resid_perm / df_resid)
stopifnot(all(is.finite(F_perm)))

# Identify the observed assignment independently in each batch and verify the
# quadratic-form SS equals the direct QR model SS before reporting exact p.
obs1 <- as.integer(meta$cultivar[i1])
obs2 <- as.integer(meta$cultivar[i2])
obs_row1 <- which(rowSums(sweep(A, 2L, obs1, `!=`)) == 0L)
obs_row2 <- which(rowSums(sweep(A, 2L, obs2, `!=`)) == 0L)
stopifnot(length(obs_row1) == 1L, length(obs_row2) == 1L)
ss_quad_obs <- ss_cult_perm[obs_row1, obs_row2]
stopifnot(abs(ss_quad_obs - ss_cult_obs) <= 1e-8 * max(1, abs(ss_cult_obs)),
          abs(F_perm[obs_row1, obs_row2] - F_obs) <= 1e-8 * max(1, abs(F_obs)))

tol <- sqrt(.Machine$double.eps)
exceed <- sum(F_perm >= F_obs - tol)
permutation_space <- length(F_perm)
p_exact <- exceed / permutation_space
result <- data.frame(
  analysis = "PICRUSt2 bacteria-only KO profile; Bray-Curtis PERMANOVA",
  model = "cultivar + sequencing_batch; cultivar partial effect conditional on batch",
  samples = n, features = ncol(X), df_cultivar = df_cult, df_residual = df_resid,
  SS_batch = ss_batch, SS_cultivar_partial = ss_cult_obs,
  SS_residual_full_model = ss_resid_obs, pseudo_F = F_obs,
  partial_R2 = partial_R2, R2_total = R2_total,
  exact_p_upper_tail = p_exact, exceedances = exceed,
  permutations = permutation_space,
  stringsAsFactors = FALSE, check.names = FALSE
)
write.csv(result, file.path(outdir, "PERMANOVA_KO_partial_cultivar_given_batch.csv"),
          row.names = FALSE, na = "")

write.csv(data.frame(
  analysis = "Exact within-batch permutation distribution",
  df_cultivar = df_cult, df_residual = df_resid,
  permutations = permutation_space, F_observed = F_obs,
  F_min = min(F_perm), F_median = median(F_perm), F_max = max(F_perm),
  exceedances = exceed, exact_p_upper_tail = p_exact,
  note = "F recomputed for every assignment using SS_residual = SS_total - SS_batch - SS_cultivar; no Monte Carlo sampling."
), file.path(outdir, "PERMANOVA_KO_permutation_summary.csv"), row.names = FALSE)

writeLines(c(
  "Corrected PICRUSt2 bacteria-only KO PERMANOVA",
  "Model: sequencing_batch + cultivar; cultivar is a partial effect conditional on batch.",
  sprintf("Samples: %d; predicted KO features: %d.", n, ncol(X)),
  sprintf("Degrees of freedom: cultivar=%d; full-model residual=%d.", df_cult, df_resid),
  sprintf("Exact assignments: %d x %d = %d.", nrow(A), nrow(A), permutation_space),
  sprintf("Pseudo-F=%.12g; partial R2=%.12g; total R2=%.12g; exact upper-tail p=%.12g (%d/%d).",
          F_obs, partial_R2, R2_total, p_exact, exceed, permutation_space),
  "SS definitions: SS_cultivar = trace((P_full - P_batch) G); SS_residual = trace((I - P_full) G).",
  "The permutation statistic recomputes F from the same partial SS and residual SS formula for every assignment.",
  "This script was not executed in the originating environment; verify execution and outputs on a machine with R.",
  ""
), file.path(outdir, "README_reanalise.txt"))

print(result, row.names = FALSE)
