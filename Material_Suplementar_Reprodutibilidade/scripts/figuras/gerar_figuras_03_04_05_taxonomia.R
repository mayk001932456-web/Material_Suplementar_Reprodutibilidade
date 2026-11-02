# =============================================================================
# Figuras 3, 4 e 5: abundancia relativa por filo, familia e genero.
# Geradas a partir dos dados integrados de outubro/2026
# Paleta Oficial do TCC
# =============================================================================

# === PACOTES ===
required_packages <- c("ggplot2", "dplyr", "tidyr", "scales", "stringr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Instale os pacotes necessários: ", paste(missing_packages, collapse = ", "))
}
invisible(lapply(required_packages, library, character.only = TRUE))


.script_args <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(.script_args)) stop("Run this script with Rscript so its location can be resolved.")
.script_file <- normalizePath(sub("^--file=", "", .script_args[1]), winslash = "/", mustWork = TRUE)
source(file.path(dirname(dirname(.script_file)), "package_config.R"))

# === CAMINHOS ===
DADOS_DIR <- INPUTS_ROOT
OUT_DIR <- file.path(RESULTS_ROOT, "figuras")

# === PALETA OFICIAL DO TCC ===
CORES_CULTIVAR <- c(
  "Solo (Controle)" = "#777777",
  "BL-65"           = "#E74C3C",
  "BL-69"           = "#3498DB",
  "BL-74"           = "#2ECC71",
  "BL-183"          = "#F39C12"
)

ORDEM_CULTIVARES <- c("Solo (Controle)", "BL-65", "BL-69", "BL-74", "BL-183")

# === TEMA PADRÃO ===
tema_tcc <- theme_bw(base_size = 12) +
  theme(
    plot.title        = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle     = element_text(size = 10, color = "#444444", hjust = 0),
    legend.position   = "right",
    legend.title      = element_text(face = "bold", size = 10),
    legend.text       = element_text(size = 9.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color = "#EEEEEE"),
    axis.text.x        = element_text(face = "bold", size = 10.5),
    axis.title         = element_text(face = "bold", size = 11),
    plot.margin        = margin(15, 15, 15, 15)
  )

# =============================================================================
# FUNÇÃO AUXILIAR: preparar dados para plot de barras empilhadas
# =============================================================================
preparar_dados <- function(csv_path, n_top, label_outros = "Demais táxons") {
  df <- read.csv(csv_path, stringsAsFactors = FALSE)

  # Limpar nomes de táxons
  df$taxon <- str_trim(df$taxon)
  # Substituir "Não classificado" pelo rótulo padrão
  df$taxon <- ifelse(grepl("class|classif", df$taxon, ignore.case = TRUE) &
                       !grepl("^[A-Z]", df$taxon),
                     "Sem classificação", df$taxon)

  # Selecionar top N pelos rizosfera (mean_all_roots)
  top_taxa <- df %>%
    arrange(desc(mean_all_roots)) %>%
    slice_head(n = n_top) %>%
    pull(taxon)

  # Montar tabela longa por cultivar + solo
  df_long <- df %>%
    select(taxon,
           `BL-65`   = `mean_BL.65`,
           `BL-69`   = `mean_BL.69`,
           `BL-74`   = `mean_BL.74`,
           `BL-183`  = `mean_BL.183`,
           `Solo (Controle)` = mean_soil) %>%
    pivot_longer(-taxon, names_to = "Cultivar", values_to = "media_pct") %>%
    mutate(
      taxon   = ifelse(taxon %in% top_taxa, taxon, label_outros),
      Cultivar = factor(Cultivar, levels = ORDEM_CULTIVARES)
    ) %>%
    group_by(taxon, Cultivar) %>%
    summarise(media_pct = sum(media_pct, na.rm = TRUE), .groups = "drop")

  # Ordenar fator de táxon: top_taxa na frente, "Demais" e "Sem classificação" no final
  niveis_taxon <- c(
    intersect(top_taxa, df_long$taxon),
    intersect(c("Sem classificação", label_outros),
              unique(df_long$taxon))
  )
  df_long$taxon <- factor(df_long$taxon, levels = niveis_taxon)

  return(df_long)
}

# =============================================================================
# PALETA DE PREENCHIMENTO (filos, famílias, gêneros)
# =============================================================================
gerar_paleta_fill <- function(niveis) {
  cores_base <- c(
    "#1565C0", "#E65100", "#827717", "#1B5E20", "#4A148C",
    "#3E2723", "#006064", "#004D40", "#880E4F", "#B71C1C",
    "#37474F", "#558B2F", "#6A1B9A", "#0277BD", "#E64A19",
    "#00695C", "#283593", "#AD1457", "#F57F17", "#4E342E"
  )
  especiais <- c(
    "Demais táxons"    = "#AAAAAA",
    "Sem classificação" = "#CCCCCC"
  )
  normais <- setdiff(niveis, names(especiais))
  pal <- setNames(cores_base[seq_along(normais)], normais)
  c(pal, especiais[intersect(names(especiais), niveis)])
}

# =============================================================================
# FIGURA 03 — FILOS
# =============================================================================
message("\n>>> Gerando Figura 03 — Filos...")

df_filo <- preparar_dados(
  file.path(DADOS_DIR, "abundancia_phylum_integrada.csv"),
  n_top = 10
)

pal_filo <- gerar_paleta_fill(levels(df_filo$taxon))

p2 <- ggplot(df_filo, aes(x = Cultivar, y = media_pct, fill = taxon)) +
  geom_bar(stat = "identity", position = "stack",
           color = "white", linewidth = 0.25) +
  scale_fill_manual(values = pal_filo, name = "Filo Bacteriano") +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(0, 100, by = 20),
    expand = c(0, 0),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    title    = "Abundância relativa dos principais filos bacterianos",
    subtitle = "Médias por cultivar e solo-controle (dados integrados, n = 16 raízes)",
    x = NULL,
    y = "Abundância relativa (%)"
  ) +
  tema_tcc

dir_filo <- file.path(OUT_DIR, "Figura_03 - Abundancia Filos")
dir.create(dir_filo, recursive = TRUE, showWarnings = FALSE)
out_filo <- file.path(dir_filo, "Figura_03_Abundancia_Filos.png")
png(out_filo, width = 10.5, height = 6.5, units = "in", res = 300)
print(p2)
dev.off()
message("Salvo: ", out_filo)

# =============================================================================
# FIGURA 04 — FAMÍLIAS
# =============================================================================
message("\n>>> Gerando Figura 04 — Famílias...")

df_fam <- preparar_dados(
  file.path(DADOS_DIR, "abundancia_family_integrada.csv"),
  n_top = 12
)

pal_fam <- gerar_paleta_fill(levels(df_fam$taxon))

p3 <- ggplot(df_fam, aes(x = Cultivar, y = media_pct, fill = taxon)) +
  geom_bar(stat = "identity", position = "stack",
           color = "white", linewidth = 0.25) +
  scale_fill_manual(values = pal_fam, name = "Família Bacteriana") +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(0, 100, by = 20),
    expand = c(0, 0),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    title    = "Abundância relativa das principais famílias bacterianas",
    subtitle = "Médias por cultivar e solo-controle (dados integrados, n = 16 raízes)",
    x = NULL,
    y = "Abundância relativa (%)"
  ) +
  tema_tcc

dir_fam <- file.path(OUT_DIR, "Figura_04 - Abundancia Familia")
dir.create(dir_fam, recursive = TRUE, showWarnings = FALSE)
out_fam <- file.path(dir_fam, "Figura_04_Abundancia_Familia_R.png")
png(out_fam, width = 10.5, height = 6.5, units = "in", res = 300)
print(p3)
dev.off()
message("Salvo: ", out_fam)

# =============================================================================
# FIGURA 05 — GÊNEROS
# =============================================================================
message("\n>>> Gerando Figura 05 — Gêneros...")

df_gen <- preparar_dados(
  file.path(DADOS_DIR, "abundancia_genus_integrada.csv"),
  n_top = 12
)

pal_gen <- gerar_paleta_fill(levels(df_gen$taxon))

p4 <- ggplot(df_gen, aes(x = Cultivar, y = media_pct, fill = taxon)) +
  geom_bar(stat = "identity", position = "stack",
           color = "white", linewidth = 0.25) +
  scale_fill_manual(values = pal_gen, name = "Gênero Bacteriano") +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(0, 100, by = 20),
    expand = c(0, 0),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    title    = "Abundância relativa dos principais gêneros bacterianos",
    subtitle = "Médias por cultivar e solo-controle (dados integrados, n = 16 raízes)",
    x = NULL,
    y = "Abundância relativa (%)"
  ) +
  tema_tcc

dir_gen <- file.path(OUT_DIR, "Figura_05 - Abundancia Genero")
dir.create(dir_gen, recursive = TRUE, showWarnings = FALSE)
out_gen <- file.path(dir_gen, "Figura_05_Abundancia_Genero_Definitiva.png")
png(out_gen, width = 10.5, height = 6.5, units = "in", res = 300)
print(p4)
dev.off()
message("Salvo: ", out_gen)

message("\n=== Todas as figuras geradas com sucesso! ===")



