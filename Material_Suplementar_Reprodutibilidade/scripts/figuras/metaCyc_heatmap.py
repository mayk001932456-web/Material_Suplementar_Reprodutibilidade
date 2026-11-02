from pathlib import Path
import html
import textwrap

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Raiz do pacote suplementar, independente do diretório de execução.
package_root = Path(__file__).resolve().parents[2]
csv_path = package_root / "resultados" / "MetaCyc" / "MetaCyc_bacteria_only_blocked_2seeds.csv"
output_dir = package_root / "resultados" / "figuras"
output_dir.mkdir(parents=True, exist_ok=True)

data = pd.read_csv(csv_path)

# Mantém as vias com q < 0,05 nas duas rodadas
q1 = pd.to_numeric(data["q_BH_seed_20261021"], errors="coerce")
q2 = pd.to_numeric(data["q_BH_seed_20261022"], errors="coerce")
f_stat = pd.to_numeric(data["F"], errors="coerce")

selected = data.loc[(q1 < 0.05) & (q2 < 0.05)].copy()
selected["F"] = f_stat.loc[selected.index]

# Seleciona as 20 vias com maiores estatísticas F
selected = selected.sort_values("F", ascending=False).head(20)

cultivars = ["BL-65", "BL-69", "BL-74", "BL-183"]
mean_columns = [f"mean_relative_{cultivar}" for cultivar in cultivars]

means = selected[mean_columns].apply(pd.to_numeric, errors="coerce")
means.columns = cultivars

# Usa log10 para tornar visíveis diferenças entre abundâncias de escalas distintas.
# Valores zero ficam como células cinzas.
log_means = np.log10(means.where(means > 0))

labels = []
for _, row in selected.iterrows():
    description = html.unescape(str(row["description"]))
    label = f'{row["pathway"]} — {description}'
    labels.append(textwrap.fill(label, width=52))

fig_height = max(7, 0.38 * len(selected))
fig, ax = plt.subplots(figsize=(11, fig_height))

sns.heatmap(
    log_means,
    ax=ax,
    cmap="viridis",
    mask=log_means.isna(),
    linewidths=0.4,
    linecolor="white",
    cbar_kws={"label": "log10 da abundância relativa média"},
    yticklabels=labels,
    xticklabels=cultivars,
)

ax.set_title(
    "20 vias MetaCyc com maiores estatísticas F\n"
    "entre as vias com q < 0,05 nas duas rodadas",
    pad=14,
)
ax.set_xlabel("Cultivar")
ax.set_ylabel("Via MetaCyc")
ax.tick_params(axis="y", labelsize=8)
ax.tick_params(axis="x", labelrotation=0)

fig.text(
    0.01,
    0.01,
    "Células cinzas indicam abundância relativa média igual a zero. "
    "As médias são descritivas; o teste F é global e não identifica pares de cultivares.",
    ha="left",
    fontsize=8,
)

fig.tight_layout(rect=[0, 0.04, 1, 1])

png_path = output_dir / "Figura_10_MetaCyc_top20_maior_F.png"
pdf_path = output_dir / "Figura_10_MetaCyc_top20_maior_F.pdf"
table_path = output_dir / "Figura_10_MetaCyc_vias_selecionadas.csv"

fig.savefig(png_path, dpi=300, bbox_inches="tight")
fig.savefig(pdf_path, bbox_inches="tight")
selected.to_csv(table_path, index=False, encoding="utf-8-sig")
plt.close(fig)

print(f"Figura PNG salva em: {png_path}")
print(f"Figura PDF salva em: {pdf_path}")
print(f"Tabela das vias selecionadas salva em: {table_path}")
