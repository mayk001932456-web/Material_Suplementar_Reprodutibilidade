# Exact pairwise PERMANOVA under the two sequencing-batch restrictions.
.script_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(.script_args)) stop("Run this script with Rscript so its location can be resolved.")
.script_file <- normalizePath(sub("^--file=", "", .script_args[1]), winslash = "/", mustWork = TRUE)
source(file.path(dirname(dirname(.script_file)), "package_config.R"))
valid <- INPUTS_ROOT

meta <- read.delim(file.path(valid, "metadados_16_raizes.tsv"), stringsAsFactors=FALSE)
meta$cultivar <- factor(meta$cultivar, levels=c("BL-65", "BL-69", "BL-74", "BL-183"))
meta$sequencing_batch <- factor(meta$sequencing_batch, levels=c("Lote 1", "Lote 2"))

choose2 <- combn(4, 2, simplify=FALSE)
pseudo_f <- function(D, group) {
  group <- as.character(group)
  n <- nrow(D); k <- length(unique(group))
  total <- sum(D^2)/n
  within <- 0
  for (g in unique(group)) {
    ix <- which(group == g)
    within <- within + sum(D[ix, ix, drop=FALSE]^2)/length(ix)
  }
  between <- total - within
  c(F=(between/(k-1))/(within/(n-k)), R2=between/total)
}
exact_pair <- function(D, m) {
  group <- droplevels(factor(as.character(m$cultivar)))
  batch <- m$sequencing_batch
  obs <- pseudo_f(D, group)
  F_perm <- numeric(36); h <- 0L
  for (a in choose2) for (b in choose2) {
    perm <- rep("B", nrow(D))
    for (bi in levels(batch)) {
      ix <- which(batch == bi)
      perm[ix[a]] <- "A"
    }
    h <- h + 1L
    F_perm[h] <- pseudo_f(D, factor(perm, levels=c("A", "B")))["F"]
  }
  c(F=unname(obs["F"]), R2=unname(obs["R2"]), p=sum(F_perm >= obs["F"]-sqrt(.Machine$double.eps))/length(F_perm),
    exceed=sum(F_perm >= obs["F"]-sqrt(.Machine$double.eps)), space=length(F_perm))
}

matrices <- list()
asvdir <- file.path(valid, "permanova_bloqueada")
for (metric in c("Bray_Curtis", "Jaccard", "Aitchison")) {
  matrices[[paste("ASV", metric, sep="_")]] <- as.matrix(read.delim(file.path(asvdir, paste0("distancia_", metric, ".tsv")), row.names=1, check.names=FALSE))
}
gdat <- read.delim(file.path(valid, "sensibilidade_genero", "matriz_genero_com_nao_classificados.tsv"), check.names=FALSE)
gcols <- grep("^L[12]_Go", names(gdat), value=TRUE)
gc <- t(as.matrix(gdat[, gcols, drop=FALSE])); storage.mode(gc) <- "numeric"
rownames(gc) <- gcols
rel <- gc/rowSums(gc)
bc <- as.matrix(dist(rel, method="manhattan"))/2
presence <- 1.0*(gc>0); inter <- tcrossprod(presence); nn <- rowSums(presence)
ja <- 1-inter/(outer(nn,nn,"+")-inter); diag(ja)<-0
logged <- log(gc+1); clr <- logged-rowMeans(logged); ai <- as.matrix(dist(clr))
gd <- list(Bray_Curtis=bc,Jaccard=ja,Aitchison=ai)
for (nm in names(gd)) matrices[[paste("Genero",nm,sep="_")]] <- gd[[nm]]

results <- list(); idx <- 0L
for (matrix_name in names(matrices)) {
  bits <- strsplit(matrix_name, "_")[[1]]
  kind <- bits[1]; metric <- paste(bits[-1], collapse="_")
  Dall <- matrices[[matrix_name]]
  ids <- rownames(Dall)
  if (kind=="ASV") md <- meta[match(ids, meta$sample_id), ] else {
    ids_go <- sub("^L[12]_", "", ids)
    md <- meta[match(ids_go, sub("^L[12]_", "", meta$sample_id)), ]
  }
  stopifnot(!anyNA(md$cultivar), !anyNA(md$sequencing_batch))
  for (pair in combn(levels(droplevels(md$cultivar)),2,simplify=FALSE)) {
    keep <- md$cultivar %in% pair
    m <- droplevels(md[keep,])
    D <- Dall[keep,keep,drop=FALSE]
    z <- exact_pair(D,m)
    idx <- idx+1L
    results[[idx]] <- data.frame(matrix=kind,metric=metric,cultivar_1=pair[1],cultivar_2=pair[2],
      F=unname(z["F"]),R2=unname(z["R2"]),p_exact=unname(z["p"]),
      exceed=unname(z["exceed"]),permutation_space=unname(z["space"]))
  }
}
out <- do.call(rbind,results)
out$p_holm_within_metric <- ave(out$p_exact, interaction(out$matrix,out$metric), FUN=function(p)p.adjust(p,"holm"))
out$p_holm_all_18 <- ave(out$p_exact, out$matrix, FUN=function(p)p.adjust(p,"holm"))
write.csv(out,file.path(RESULTS_ROOT,"permanova","PERMANOVA_pareada_exata_bloqueada.csv"),row.names=FALSE)
print(out,row.names=FALSE)
