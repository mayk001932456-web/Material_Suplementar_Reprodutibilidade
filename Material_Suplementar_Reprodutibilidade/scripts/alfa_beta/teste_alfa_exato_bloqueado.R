# Teste exato da associação entre cultivar e diversidade alfa, com lote como bloco.
# Unidade de análise: média das 100 rarefações por amostra (16 amostras independentes).
# Rótulos de cultivar são enumerados dentro de cada lote; quatro métricas, correção Holm.

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
source(file.path(dirname(dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE))), "package_config.R"))
out_dir <- file.path(RESULTS_ROOT, "alfa")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
root_dir <- INPUTS_ROOT
alpha_file <- file.path(RESULTS_ROOT, "alfa_bacteriana", "alfa_rarefeita_ASV_70554_100iter.csv")
meta_file <- file.path(root_dir, "metadados_16_raizes.tsv")

alpha <- read.csv(alpha_file, check.names = FALSE, stringsAsFactors = FALSE)
meta <- read.delim(meta_file, stringsAsFactors = FALSE)
meta <- meta[meta$role == "raiz", , drop = FALSE]
meta$label <- alpha$SampleID[match(meta$label, alpha$SampleID)]
alpha <- alpha[match(meta$label, alpha$SampleID), , drop = FALSE]
stopifnot(nrow(meta) == 16L, !anyNA(alpha$SampleID), all(meta$label == alpha$SampleID))

meta$cultivar <- factor(meta$cultivar, levels = c("BL-65", "BL-69", "BL-74", "BL-183"))
meta$sequencing_batch <- factor(meta$sequencing_batch, levels = c("Lote 1", "Lote 2"))
stopifnot(all(table(meta$cultivar, meta$sequencing_batch) == 2L))

make_assignments <- function() {
  A <- matrix(0L, nrow = 2520L, ncol = 8L)
  row <- 0L
  for (g1 in combn(8L, 2L, simplify = FALSE)) {
    rest1 <- setdiff(1:8, g1)
    for (g2 in combn(rest1, 2L, simplify = FALSE)) {
      rest2 <- setdiff(rest1, g2)
      for (g3 in combn(rest2, 2L, simplify = FALSE)) {
        g4 <- setdiff(rest2, g3)
        row <- row + 1L
        A[row, g1] <- 1L
        A[row, g2] <- 2L
        A[row, g3] <- 3L
        A[row, g4] <- 4L
      }
    }
  }
  stopifnot(row == 2520L, all(rowSums(A == 1L) == 2L), all(rowSums(A == 2L) == 2L),
            all(rowSums(A == 3L) == 2L), all(rowSums(A == 4L) == 2L))
  A
}
A <- make_assignments()
batch_idx <- split(seq_len(nrow(meta)), meta$sequencing_batch)
stopifnot(length(batch_idx) == 2L, all(lengths(batch_idx) == 8L))

exact_blocked_test <- function(y, metric) {
  total <- sum(y)
  grand_ss <- sum((y - mean(y))^2)
  batch_totals <- vapply(batch_idx, function(ix) sum(y[ix]), numeric(1))
  ss_batch <- sum(batch_totals^2 / 8) - total^2 / 16

  Q <- lapply(batch_idx, function(ix) {
    yy <- y[ix]
    sapply(1:4, function(g) as.numeric((A == g) %*% yy))
  })
  # Enumerate every pair of within-lot label assignments (2520 x 2520).
  w1 <- rowSums(Q[[1]]^2) / 4
  w2 <- rowSums(Q[[2]]^2) / 4
  ss_cult <- outer(w1, w2, "+") + 0.5 * (Q[[1]] %*% t(Q[[2]])) - total^2 / 16
  f_all <- (ss_cult / 3) / ((grand_ss - ss_batch - ss_cult) / 11)

  obs_idx <- lapply(batch_idx, function(ix) {
    codes <- as.integer(meta$cultivar[ix])
    which(apply(A, 1L, function(z) all(z == codes)))
  })
  stopifnot(all(lengths(obs_idx) == 1L))
  ss_obs <- ss_cult[obs_idx[[1]], obs_idx[[2]]]
  f_obs <- (ss_obs / 3) / ((grand_ss - ss_batch - ss_obs) / 11)
  exceed <- sum(f_all >= f_obs - sqrt(.Machine$double.eps))

  # Independently verify the observed partial F against nested linear models.
  d <- data.frame(y = y, batch = meta$sequencing_batch, cultivar = meta$cultivar)
  fit_reduced <- lm(y ~ batch, data = d)
  fit_full <- lm(y ~ batch + cultivar, data = d)
  nested <- anova(fit_reduced, fit_full)
  f_lm <- nested[2L, "F"]
  if (!isTRUE(all.equal(as.numeric(f_obs), as.numeric(f_lm), tolerance = 1e-9))) {
    stop(sprintf("Partial-F cross-check failed for %s: exact F %.12g; nested F %.12g; exact SS cultivar %.12g; nested SS cultivar %.12g; exact SSE %.12g; nested SSE %.12g; batch SS %.12g; labels=%s",
                 metric, f_obs, f_lm, ss_obs, nested$RSS[1] - nested$RSS[2],
                 grand_ss - ss_batch - ss_obs, nested$RSS[2], ss_batch,
                 paste(as.integer(meta$cultivar[batch_idx[[1]]]), collapse = ",")))
  }

  data.frame(metric = metric, F = f_obs, df_numerator = 3L, df_residual = 11L,
             partial_R2 = ss_obs / (ss_obs + grand_ss - ss_batch - ss_obs),
             p_exact = exceed / length(f_all), exceed = exceed,
             permutation_space = length(f_all), stringsAsFactors = FALSE)
}

metrics <- c(observed = "Riqueza observada", chao1 = "Chao1", shannon = "Shannon", pielou = "Pielou")
res <- do.call(rbind, lapply(names(metrics), function(m) exact_blocked_test(alpha[[m]], metrics[[m]])))
res$p_holm_4_metrics <- p.adjust(res$p_exact, method = "holm")
write.csv(res, file.path(out_dir, "teste_alfa_exato_bloqueado.csv"), row.names = FALSE)

means <- aggregate(alpha[names(metrics)], list(cultivar = meta$cultivar), mean)
write.csv(means, file.path(out_dir, "resumo_alfa_por_cultivar_para_teste_bloqueado.csv"), row.names = FALSE)
cat("R:", R.version.string, "\n")
cat("Lotes x cultivares:\n"); print(table(meta$sequencing_batch, meta$cultivar))
cat("Enumerações por lote:", nrow(A), "| espaço exato total:", nrow(A)^2, "\n")
print(res, row.names = FALSE)

