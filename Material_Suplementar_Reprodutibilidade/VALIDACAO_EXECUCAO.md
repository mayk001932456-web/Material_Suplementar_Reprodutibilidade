# Registro de validacao local dos scripts

Data da verificacao: 2 de novembro de 2026. Ambiente: Windows, R 4.6.1. Os scripts foram chamados a partir da pasta do pacote suplementar para conferir o uso de caminhos relativos e a producao dos resultados.

## Executado com sucesso

| Script | Resultado da execucao |
|---|---|
| `recalcular_alfa_bacteriana.R` | Gerou indices de 16 raizes a 70.554 leituras, com 100 subamostragens por amostra. |
| `teste_alfa_exato_bloqueado.R` | Gerou os testes exatos de quatro metricas. Os p ajustados de Holm foram 0,398 (riqueza), 0,398 (Chao1), 0,178 (Shannon) e 0,154 (Pielou), em concordancia com os valores informados no TCC, considerando arredondamento. |
| `permanova_exata_bloqueada.R` | Reproduziu os seis testes globais de ASVs e genero; pseudo-F, R2, p exato e p ajustado coincidiram com a tabela do TCC. |
| `permanova_pareada_exata_bloqueada.R` | Executou os 36 contrastes (18 ASV e 18 genero); os pseudo-F, R2 e p exatos coincidiram com as tabelas suplementares. |
| `permdisp_referencia_vegan.R` | Executou a reimplementacao em R base do algoritmo documentado de PERMDISP. Para ASVs, os p foram 0,7894 (Bray-Curtis), 0,9310 (Jaccard) e 0,5890 (Aitchison); coincidem com os valores arredondados do texto. |
| `screen_taxa_permutacao_lote.R` | Avaliou 535 familias (455 variaveis) e 1.166 generos (981 variaveis), sem q < 0,05, conforme relatado. |
| `recalcular_biomassa_alfa_ASV.R` | Reproduziu as correlacoes de Spearman (rho = 0,371 e 0,406; p-Holm = 0,238 para ambas) e os modelos de massa. |
| `reanalise_picrust2_bloqueada_ajustada.R` | Reproduziu F parcial = 3,9041, R2 parcial = 0,5157, p exato = 0,005155 e 32.736 de 6.350.400 atribuicoes extremas. |
| `MetaCyc_permutacao_duas_sementes.R` | Reproduziu exatamente os valores numericos do CSV de referencia para as 534 vias: estatistica F, p/q das duas rodadas e medias por cultivar. Obteve 106, 104 e 104 vias significativas na primeira, segunda e ambas as rodadas, respectivamente. |
| Figuras 6, 7 e 9 | Os scripts geraram PNGs a partir dos arquivos do pacote. As imagens foram abertas para inspecao visual. |

As tabelas de resultados geradas estao em `resultados/`.

## Bloqueios por dependencias ausentes

- `ANCOMBC2_analysis.R`: nao executou porque `phyloseq` nao esta instalado na biblioteca usada; a instalacao local do `ANCOMBC` tambem nao carrega por falta de `data.table`. Os CSVs dos resultados ANCOM-BC2 ja calculados continuam incluidos em `dados/`.
- Figuras 2, 3, 4, 5 e 8: nao executaram porque faltam pacotes R usados pelos scripts, como `readxl`, `ggplot2`, `scales`, `ggrepel`, `dplyr`, `tidyr` e `stringr`.
- Figura 10: nao executou porque nao ha uma instalacao de Python acessivel nesta sessao. A tabela completa de MetaCyc e o script de geracao do mapa de calor estao incluidos.

A predicao PICRUSt2 nao foi refeita: os scripts testam os perfis de KO e MetaCyc fornecidos como entradas. O processamento inicial pelo prestador tambem nao foi reexecutado, pois as leituras brutas nao fazem parte do pacote.

## Observacoes do ambiente

O R 4.6.1 emitiu avisos de configuracao de localidade (`LC_COLLATE`, `LC_CTYPE`, `LC_MONETARY` e `LC_TIME`) em cada inicializacao. Esses avisos nao impediram os scripts listados como executados com sucesso. A validacao aqui confirma a execucao local destes scripts e entradas; nao substitui a instalacao das dependencias e a reproducao integral em outro computador.