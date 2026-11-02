
.script_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(.script_args)) stop("Run this script with Rscript so its location can be resolved.")
.script_file <- normalizePath(sub("^--file=", "", .script_args[1]), winslash = "/", mustWork = TRUE)
source(file.path(dirname(dirname(.script_file)), "package_config.R"))

# MetaCyc pathway tests for the bacteria-only PICRUSt2 predictions.
# Only the 16 root samples enter inference; labels are permuted within batch.

valid <- INPUTS_ROOT
pic <- file.path(INPUTS_ROOT, "picrust2")
outdir <- file.path(RESULTS_ROOT, "MetaCyc")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

meta <- read.delim(file.path(valid, "metadados_16_raizes.tsv"), stringsAsFactors=FALSE)
meta$go <- sub("^L[12]_", "", meta$sample_id)
meta$cultivar <- factor(meta$cultivar, levels=c("BL-65", "BL-69", "BL-74", "BL-183"))
meta$sequencing_batch <- factor(meta$sequencing_batch, levels=c("Lote 1", "Lote 2"))
stopifnot(nrow(meta)==16L, all(table(meta$sequencing_batch, meta$cultivar)==2L))
samples <- meta$go
read_gz <- function(rel) read.delim(gzfile(file.path(pic, rel)), check.names=FALSE, stringsAsFactors=FALSE)

# Build all 2,520 allowed assignments for one batch (two replicates per cultivar).
make_assignments <- function() {
  out <- matrix(0L, 2520L, 8L); r <- 0L
  for (a in combn(8,2,simplify=FALSE)) {
    rem1 <- setdiff(1:8,a)
    for (b in combn(rem1,2,simplify=FALSE)) {
      rem2 <- setdiff(rem1,b)
      for (cc in combn(rem2,2,simplify=FALSE)) {
        rem3 <- setdiff(rem2,cc); r <- r+1L
        out[r,a] <- 1L; out[r,b] <- 2L; out[r,cc] <- 3L; out[r,rem3] <- 4L
      }
    }
  }
  out
}
A <- make_assignments()
i1 <- which(meta$sequencing_batch==levels(meta$sequencing_batch)[1])
i2 <- which(meta$sequencing_batch==levels(meta$sequencing_batch)[2])
Mobs1 <- sapply(1:4,function(k)as.numeric(meta$cultivar[i1]==levels(meta$cultivar)[k]))
Mobs2 <- sapply(1:4,function(k)as.numeric(meta$cultivar[i2]==levels(meta$cultivar)[k]))

# One-way cultivar F after accounting for batch; permutation labels stay inside batch.
path0 <- read_gz("pathways_out/path_abun_unstrat_descrip.tsv.gz")
Y <- as.matrix(path0[,samples,drop=FALSE]); storage.mode(Y) <- "numeric"
rownames(Y) <- path0$pathway
Y <- sweep(Y, 2, colSums(Y), "/")
Y1 <- Y[,i1,drop=FALSE]; Y2 <- Y[,i2,drop=FALSE]
meanall <- rowMeans(Y)
ss_total <- rowSums((Y1-rowMeans(Y1))^2)+rowSums((Y2-rowMeans(Y2))^2)
Sobs <- Y1%*%Mobs1 + Y2%*%Mobs2
ss_obs <- rowSums((Sobs-4*meanall)^2)/4
F_obs <- (ss_obs/3)/((ss_total-ss_obs)/11)

screen_seed <- function(seed, B=99999L) {
  set.seed(seed); exceed <- integer(nrow(Y))
  for (bi in seq_len(B)) {
    lab1 <- sample(rep(1:4,each=2)); lab2 <- sample(rep(1:4,each=2))
    M1 <- sapply(1:4,function(k)as.numeric(lab1==k))
    M2 <- sapply(1:4,function(k)as.numeric(lab2==k))
    S <- Y1%*%M1 + Y2%*%M2
    ssb <- rowSums((S-4*meanall)^2)/4
    Fp <- (ssb/3)/((ss_total-ssb)/11)
    exceed <- exceed + (Fp >= F_obs-1e-14)
  }
  (exceed+1)/(B+1)
}

B <- 99999L
p1 <- screen_seed(20261021L,B)
p2 <- screen_seed(20261022L,B)
q1 <- p.adjust(p1,"BH"); q2 <- p.adjust(p2,"BH")
pathres <- data.frame(pathway=rownames(Y), description=path0$description,
  F=F_obs, mean_relative_abundance=rowMeans(Y), p_seed_20261021=p1,
  q_BH_seed_20261021=q1, p_seed_20261022=p2, q_BH_seed_20261022=q2,
  significant_both_seeds=q1<0.05 & q2<0.05,
  significant_either_seed=q1<0.05 | q2<0.05, check.names=FALSE)
cv_means <- sapply(levels(meta$cultivar), function(cv) rowMeans(Y[,meta$cultivar==cv,drop=FALSE]))
colnames(cv_means) <- paste0("mean_relative_", levels(meta$cultivar))
pathres <- cbind(pathres, cv_means)
pathres <- pathres[order(pathres$q_BH_seed_20261021,pathres$p_seed_20261021),]
write.csv(pathres,file.path(outdir,"MetaCyc_bacteria_only_blocked_2seeds.csv"),row.names=FALSE)
write.csv(pathres[pathres$significant_both_seeds,],file.path(outdir,"MetaCyc_significativas_duas_rodadas_q05.csv"),row.names=FALSE)
stability <- data.frame(
  bacteria_only=TRUE, root_samples=nrow(meta), pathways=nrow(Y),
  permutations_per_seed=B, seed1=20261021L, seed2=20261022L,
  significant_seed1=sum(q1<0.05), significant_seed2=sum(q2<0.05),
  significant_both_seeds=sum(q1<0.05&q2<0.05),
  significant_either_seed=sum(q1<0.05|q2<0.05),
  status_changes=sum((q1<0.05)!=(q2<0.05)),
  spearman_p=cor(p1,p2,method="spearman"),
  minimum_q_seed1=min(q1), minimum_q_seed2=min(q2))
write.csv(stability,file.path(outdir,"variacao_monte_carlo_MetaCyc.csv"),row.names=FALSE)

