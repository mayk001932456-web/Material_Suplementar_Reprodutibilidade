# Figura 6: diversidade alfa das amostras radiculares.
args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
source(file.path(dirname(dirname(normalizePath(script_arg, winslash = "/", mustWork = TRUE))), "package_config.R"))
fig_dir <- file.path(RESULTS_ROOT, "figuras")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
valid <- file.path(RESULTS_ROOT, "alfa_bacteriana")
csv <- file.path(valid, "alfa_rarefeita_ASV_70554_100iter.csv")
out <- file.path(fig_dir, "Figura_06_Diversidade_Alfa_ASV_validacao.png")
df <- read.csv(csv, check.names = FALSE)
file.copy(csv, file.path(fig_dir, basename(csv)), overwrite = TRUE)
file.copy(file.path(valid, "resumo_alfa_por_cultivar.csv"), fig_dir, overwrite = TRUE)

ordem <- c("BL-65", "BL-69", "BL-74", "BL-183")
cores <- c("BL-65"="#E74C3C", "BL-69"="#3498DB", "BL-74"="#2ECC71", "BL-183"="#F39C12")
metricas <- list(
  list(col="observed", sd="observed_sd", title="Riqueza observada", y="ASVs observados"),
  list(col="chao1", sd="chao1_sd", title="Riqueza estimada", y="Chao1"),
  list(col="shannon", sd="shannon_sd", title="Diversidade de Shannon", y="H'"),
  list(col="pielou", sd="pielou_sd", title="Equitabilidade de Pielou", y="J'")
)
png(out, width=2200, height=1750, res=220, type="cairo")
op <- par(mfrow=c(2,2), mar=c(4.2,4.6,3.2,1.2), oma=c(3.6,0,3.5,0), las=1,
          family="sans", cex.axis=0.9, cex.lab=1.0)
for (m in metricas) {
  y <- df[[m$col]]
  means <- tapply(y, factor(df$Cultivar, levels=ordem), mean)
  sds <- tapply(y, factor(df$Cultivar, levels=ordem), sd)
  plot(NA, xlim=c(0.5,4.5), ylim=range(c(0, y, means-sds, means+sds), finite=TRUE),
       xaxt="n", xlab="Cultivar", ylab=m$y, main=m$title, bty="l")
  axis(1, at=1:4, labels=ordem)
  abline(h=pretty(range(c(0,y), finite=TRUE)), col="#E6E6E6", lty=1)
  for (i in seq_along(ordem)) {
    ii <- which(df$Cultivar == ordem[i])
    points(rep(i,length(ii))+c(-0.08,0.08,0.13,-0.13)[seq_along(ii)], y[ii],
           pch=ifelse(df$Batch[ii]=="L1", 21, 24), bg=cores[ordem[i]], col="black", cex=1.2)
    segments(i, means[i]-sds[i], i, means[i]+sds[i], lwd=1.8)
    segments(i-0.09, means[i]-sds[i], i+0.09, means[i]-sds[i], lwd=1.8)
    segments(i-0.09, means[i]+sds[i], i+0.09, means[i]+sds[i], lwd=1.8)
    points(i, means[i], pch=23, bg="white", col="black", cex=1.25)
  }
}
mtext("Diversidade alfa bacteriana apos rarefacao a 70.554 leituras", outer=TRUE,
      side=3, line=1.2, cex=1.2, font=2)
mtext("Circulos = lote 1; triangulos = lote 2; losango = media; barras = desvio padrao entre as quatro amostras.",
      outer=TRUE, side=1, line=1.0, cex=0.85)
par(op)
dev.off()
cat("Imagem e dados salvos em", fig_dir, "\n")



