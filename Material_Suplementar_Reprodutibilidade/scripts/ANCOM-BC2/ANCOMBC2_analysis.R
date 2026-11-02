## ANCOM-BC2: triagem exploratoria de abundancia diferencial em familias e generos
## Requer ANCOMBC e phyloseq instalados na biblioteca do R.


.script_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(.script_args)) stop("Run this script with Rscript so its location can be resolved.")
.script_file <- normalizePath(sub("^--file=", "", .script_args[1]), winslash = "/", mustWork = TRUE)
source(file.path(dirname(dirname(.script_file)), "package_config.R"))

if (!requireNamespace("ANCOMBC", quietly = TRUE)) {
  stop("Pacote ANCOMBC ausente. Instale com BiocManager::install('ANCOMBC').")
}
if (!requireNamespace("phyloseq", quietly = TRUE)) {
  stop("Pacote phyloseq ausente. Ele e uma dependencia necessaria para este script.")
}

suppressPackageStartupMessages(library(ANCOMBC))
suppressPackageStartupMessages(library(phyloseq))

project_dir <- PACKAGE_ROOT
data_dir <- INPUTS_ROOT
counts_path <- file.path(data_dir, "matriz_bacterias_ASV_16_raizes.tsv")
metadata_path <- file.path(data_dir, "metadados_16_raizes.tsv")
out_dir <- file.path(RESULTS_ROOT, "ANCOMBC2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

counts <- read.delim(counts_path, check.names = FALSE, stringsAsFactors = FALSE)
metadata <- read.delim(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)

required_meta <- c("sample_id", "sequencing_batch", "cultivar", "role")
required_tax <- c("ASV", "Dominio", "Filo", "Classe", "Ordem", "Familia", "Genero", "Especie")
if (!all(required_meta %in% names(metadata))) {
  stop("Colunas de metadados ausentes: ", paste(setdiff(required_meta, names(metadata)), collapse = ", "))
}
if (!all(required_tax %in% names(counts))) {
  stop("Colunas taxonomicas ausentes: ", paste(setdiff(required_tax, names(counts)), collapse = ", "))
}

metadata <- metadata[metadata$role == "raiz", , drop = FALSE]
sample_ids <- metadata$sample_id
if (!all(sample_ids %in% names(counts))) {
  stop("Nem todas as amostras dos metadados existem como colunas na matriz de contagens.")
}

count_mat <- as.matrix(counts[, sample_ids, drop = FALSE])
storage.mode(count_mat) <- "numeric"
rownames(count_mat) <- counts$ASV
if (anyNA(count_mat) || any(count_mat < 0) || any(count_mat != floor(count_mat))) {
  stop("A matriz precisa conter contagens inteiras nao negativas e sem valores ausentes.")
}

tax_df <- counts[, c("Dominio", "Filo", "Classe", "Ordem", "Familia", "Genero", "Especie"), drop = FALSE]
names(tax_df) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax_mat <- as.matrix(tax_df)
rownames(tax_mat) <- counts$ASV
tax_mat[!is.na(tax_mat) & trimws(tax_mat) == ""] <- NA_character_

metadata$cultivar <- factor(metadata$cultivar, levels = c("BL-65", "BL-69", "BL-74", "BL-183"))
metadata$sequencing_batch <- factor(metadata$sequencing_batch)
meta_df <- as.data.frame(metadata)
rownames(meta_df) <- meta_df$sample_id

ps <- phyloseq(
  otu_table(count_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(meta_df)
)

if (nsamples(ps) != 16L || ntaxa(ps) != nrow(counts)) {
  stop("A verificacao de dimensoes da matriz/objeto phyloseq falhou.")
}

run_ancombc2 <- function(rank) {
  set.seed(20261102)
  ancombc2(
    data = ps,
    tax_level = rank,
    fix_formula = "sequencing_batch + cultivar",
    rand_formula = NULL,
    p_adj_method = "BH",
    prv_cut = 0.10,
    lib_cut = 1000,
    s0_perc = 0.05,
    group = "cultivar",
    struc_zero = TRUE,
    neg_lb = TRUE,
    alpha = 0.05,
    global = TRUE,
    pairwise = TRUE,
    pseudo_sens = TRUE,
    verbose = TRUE,
    n_cl = 1,
    mdfdr_control = list(fwer_ctrl_method = "holm", B = 100)
  )
}

save_result <- function(result, name) {
  tables <- c("res", "res_global", "res_pair", "zero_ind")
  for (tab in tables) {
    if (!is.null(result[[tab]])) {
      write.csv(
        as.data.frame(result[[tab]]),
        file.path(out_dir, paste0("ANCOMBC2_", name, "_", tab, ".csv")),
        row.names = FALSE,
        na = ""
      )
    }
  }
}

family_result <- run_ancombc2("Family")
save_result(family_result, "family")

genus_result <- run_ancombc2("Genus")
save_result(genus_result, "genus")

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
cat("ANCOM-BC2 concluido. Resultados salvos em: ", out_dir, "\n", sep = "")
