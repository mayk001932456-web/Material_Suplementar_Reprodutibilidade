# Figura 8: massa umida da amostra radicular (raizes + solo aderido) e diversidade alfa.
# Requer que recalcular_biomassa_alfa_ASV.R tenha sido executado antes.
args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
if (!length(script_arg)) stop("Execute com Rscript para localizar o pacote.")
script_file <- normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE)
source(file.path(dirname(dirname(script_file)), "package_config.R"))
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Instale ggplot2 antes de executar.")

input <- file.path(RESULTS_ROOT, "biomassa_alfa", "biomassa_e_alfa_por_amostra.csv")
outdir <- file.path(RESULTS_ROOT, "figuras")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input)) stop("Execute primeiro scripts/alfa_beta/recalcular_biomassa_alfa_ASV.R")
d <- read.csv(input, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(nrow(d) == 16L, all(c("biomass_g", "observed", "shannon", "cultivar", "sequencing_batch") %in% names(d)))
d$cultivar <- factor(d$cultivar, levels = c("BL-65", "BL-69", "BL-74", "BL-183"))
d$sequencing_batch <- factor(d$sequencing_batch, levels = c("Lote 1", "Lote 2"))
cols <- c("BL-65"="#E74C3C", "BL-69"="#3498DB", "BL-74"="#2ECC71", "BL-183"="#F39C12")

rho <- vapply(c("observed", "shannon"), function(v) cor.test(d$biomass_g, d[[v]], method="spearman", exact=FALSE)$estimate, numeric(1))
pv <- vapply(c("observed", "shannon"), function(v) cor.test(d$biomass_g, d[[v]], method="spearman", exact=FALSE)$p.value, numeric(1))
qv <- p.adjust(pv, method="holm")
subtxt <- sprintf("Spearman: riqueza ρ = %.3f, p-Holm = %.3f; Shannon ρ = %.3f, p-Holm = %.3f", rho[1], qv[1], rho[2], qv[2])
long <- rbind(data.frame(d, metric="Riqueza observada", value=d$observed),
              data.frame(d, metric="Shannon", value=d$shannon))
p <- ggplot2::ggplot(long, ggplot2::aes(x=biomass_g, y=value, color=cultivar, shape=sequencing_batch)) +
  ggplot2::geom_smooth(ggplot2::aes(group=1), method="lm", se=FALSE, color="#555555", linewidth=0.7) +
  ggplot2::geom_point(size=3) +
  ggplot2::facet_wrap(~metric, scales="free_y", nrow=1) +
  ggplot2::scale_color_manual(values=cols, name="Cultivar") +
  ggplot2::scale_shape_manual(values=c("Lote 1"=16, "Lote 2"=17), name="Lote") +
  ggplot2::labs(title="Massa umida das amostras radiculares e diversidade alfa", subtitle=subtxt,
       x="Massa umida (raizes com solo aderido; g)", y=NULL,
       caption="Linha: ajuste linear simples, apresentado como guia visual. A correlacao principal foi Spearman.") +
  ggplot2::theme_bw(base_size=12) + ggplot2::theme(plot.title=ggplot2::element_text(face="bold"),
       legend.position="right", strip.text=ggplot2::element_text(face="bold"))
ggplot2::ggsave(file.path(outdir, "Figura_08_Massa_Alfa.png"), p, width=10, height=5.2, dpi=600)
