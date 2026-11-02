# Figura 7: PCoA exploratoria das comunidades bacterianas em nivel de ASV.
args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
source(file.path(dirname(dirname(normalizePath(script_arg, winslash = "/", mustWork = TRUE))), "package_config.R"))
fig_dir <- file.path(RESULTS_ROOT, "figuras")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
coord_dir <- file.path(INPUTS_ROOT, "permanova_bloqueada", "PCoA")
metrics <- c("Bray_Curtis", "Jaccard", "Aitchison")
labels <- c(Bray_Curtis="Bray-Curtis", Jaccard="Jaccard", Aitchison="Aitchison")
cols <- c("BL-65"="#E74C3C", "BL-69"="#3498DB", "BL-74"="#2ECC71", "BL-183"="#F39C12")
out <- file.path(fig_dir, "Figura_07_PCoA_bacteriana_ASV_validacao.png")
png(out, width=2400, height=1150, res=220, type="cairo")
op <- par(mfrow=c(1,3), mar=c(4.5,4.6,4.2,1), oma=c(3,0,3,0), family="sans", las=1)
for (m in metrics) {
  f <- file.path(coord_dir, paste0("coordenadas_",m,".csv"))
  d <- read.csv(f, stringsAsFactors=FALSE)
  file.copy(f, fig_dir, overwrite=TRUE)
  plot(d$Axis1, d$Axis2, type="n", xlab=sprintf("PCoA1 (%.1f%%)",d$axis1_percent[1]),
       ylab=sprintf("PCoA2 (%.1f%%)",d$axis2_percent[1]), main=unname(labels[m]),
       cex.lab=0.9, cex.axis=0.8, cex.main=1.05)
  abline(h=0,v=0,col="#D9D9D9",lty=2)
  pch <- ifelse(d$sequencing_batch=="Lote 1",21,24)
  points(d$Axis1,d$Axis2,pch=pch,bg=cols[d$cultivar],col="black",cex=1.45)
  if (m == "Jaccard") legend("topright", bty="n", pch=21, pt.bg=cols, col="black",
                               legend=names(cols), cex=0.75, title="Cultivar")
}
mtext("PCoA da comunidade bacteriana em nivel de ASV", outer=TRUE, side=3, line=1,
      cex=1.15, font=2)
mtext("Cor = cultivar; circulo = lote 1; triangulo = lote 2", outer=TRUE, side=1,
      line=1.2, cex=0.9)
par(op)
dev.off()
cat("Figura salva em", fig_dir, "\n")


