# Curva de rarefação da riqueza observada para o material suplementar.
# Os pontos finais são identificados pelo código da amostra; a legenda mostra somente cultivar.
script_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_args)) stop("Execute com Rscript para resolver a pasta do pacote.")
script_file <- normalizePath(sub("^--file=", "", script_args[1]), winslash = "/", mustWork = TRUE)
source(file.path(dirname(dirname(script_file)), "package_config.R"))

required <- c("readxl", "ggplot2", "scales", "ggrepel", "dplyr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Instale os pacotes necessários: ", paste(missing, collapse = ", "))

library(readxl)
library(ggplot2)
library(ggrepel)
library(dplyr)

input_file <- file.path(INPUTS_ROOT, "dados_rarefacao.xlsx")
output_dir <- file.path(RESULTS_ROOT, "figuras")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_file <- file.path(output_dir, "Figura_02_Curva_Rarefacao.png")
rarefaction_depth <- 70554

colors_cultivar <- c(
  "Solo (Controle)" = "#777777",
  "BL-65" = "#E74C3C",
  "BL-69" = "#3498DB",
  "BL-74" = "#2ECC71",
  "BL-183" = "#F39C12"
)

df <- read_excel(input_file)
required_cols <- c("Amostra", "Cultivar", "Profundidade_Reads", "Riqueza_OTUs")
if (!all(required_cols %in% names(df))) {
  stop("Colunas ausentes na tabela de rarefação: ", paste(setdiff(required_cols, names(df)), collapse = ", "))
}
df$CultivarRotulo <- ifelse(df$Cultivar == "Solo", "Solo (Controle)", as.character(df$Cultivar))
df$CultivarRotulo <- factor(df$CultivarRotulo, levels = names(colors_cultivar))
endpoints <- df %>%
  group_by(Amostra, CultivarRotulo) %>%
  slice_max(order_by = Profundidade_Reads, n = 1, with_ties = FALSE) %>%
  ungroup()

p <- ggplot(df, aes(x = Profundidade_Reads, y = Riqueza_OTUs,
                    color = CultivarRotulo, group = Amostra)) +
  geom_line(linewidth = 0.95, alpha = 0.88) +
  geom_vline(xintercept = rarefaction_depth, linetype = "dotted",
             color = "#222222", linewidth = 0.85) +
  annotate("label", x = rarefaction_depth + 2200, y = 230,
           label = paste0("Rarefação: ", format(rarefaction_depth, big.mark = "."), " leituras"),
           hjust = 0, size = 3.6, fill = "white", color = "#222222", fontface = "bold") +
  geom_text_repel(data = endpoints,
                  aes(x = Profundidade_Reads, y = Riqueza_OTUs, label = Amostra,
                      color = CultivarRotulo),
                  inherit.aes = FALSE, hjust = 0, direction = "y", nudge_x = 3500,
                  min.segment.length = 0, segment.color = "#aaaaaa",
                  segment.size = 0.35, segment.alpha = 0.8,
                  box.padding = 0.22, point.padding = 0.3,
                  size = 3.5, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = colors_cultivar, name = "Cultivar / genótipo") +
  scale_x_continuous(labels = scales::label_comma(big.mark = ".", decimal.mark = ","),
                     limits = c(0, 165000), breaks = seq(0, 150000, by = 25000),
                     expand = c(0, 0)) +
  scale_y_continuous(labels = scales::label_comma(big.mark = ".", decimal.mark = ","),
                     limits = c(0, 2900), breaks = seq(0, 2800, by = 500),
                     expand = c(0, 0)) +
  labs(title = "Curvas de rarefação da riqueza de ASVs bacterianas",
       subtitle = "Amostras radiculares (n = 16) e solo-controle (n = 1)",
       x = "Número de leituras amostradas",
       y = "Riqueza observada de ASVs") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold", size = 11),
        legend.text = element_text(size = 10.5),
        legend.key.width = grid::unit(1.3, "cm"),
        plot.title = element_text(face = "bold", size = 13.5, hjust = 0),
        plot.subtitle = element_text(size = 10.5, color = "#444444", hjust = 0),
        panel.background = element_rect(fill = "#FAFAFA"),
        panel.grid.major = element_line(color = "#EEEEEE", linewidth = 0.6),
        panel.grid.minor = element_blank(),
        axis.title = element_text(face = "bold", size = 11),
        axis.text = element_text(size = 10.5),
        plot.margin = margin(15, 20, 15, 15))

ggsave(output_file, p, width = 11.5, height = 7, units = "in", dpi = 1200, bg = "white")
message("Figura salva em: ", output_file)
