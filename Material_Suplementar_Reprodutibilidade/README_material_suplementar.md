# Material suplementar para reprodutibilidade

## Finalidade

Este pacote reúne cópias portáteis dos scripts diretamente ligados às análises e figuras apresentadas no TCC, os arquivos de entrada necessários e tabelas suplementares. A tabela completa das 534 vias MetaCyc está em `dados/MetaCyc_bacteria_only_blocked_2seeds.csv`. Os arquivos originais do projeto não foram alterados.

## Estrutura

- `scripts/`: scripts estatísticos e de geração das Figuras 2 a 10, além de `package_config.R` para resolver caminhos.
- `dados/entradas/`: matrizes e tabelas de entrada. A pasta `alfa_bacteriana/` também contém os resultados de referência usados na revisão.
- `dados/`: resultados completos de MetaCyc e tabelas de sinais robustos do ANCOM-BC2.
- `ambiente/`: registro `sessionInfo` da execução do ANCOM-BC2.
- `VALIDACAO_EXECUCAO.md`: scripts executados nesta sessão e dependências que bloquearam as demais figuras/análises.
- `resultados/`: pasta criada e preenchida pelos scripts quando executados.

## Scripts de análise

- `scripts/alfa_beta/recalcular_alfa_bacteriana.R`: rarefação sem reposição a 70.554 leituras, 100 subamostragens por raiz e cálculo de riqueza observada, Chao1, Shannon e Pielou. Gera os arquivos em `resultados/alfa_bacteriana/`.
- `scripts/alfa_beta/teste_alfa_exato_bloqueado.R`: teste exato das quatro métricas de diversidade alfa, com cultivar no modelo e permutações dentro de lote.
- `scripts/alfa_beta/permanova_exata_bloqueada.R`: PERMANOVA exata para ASVs e matriz agregada em gênero.
- `scripts/alfa_beta/permanova_pareada_exata_bloqueada.R`: comparações pareadas por distância e matriz.
- `scripts/alfa_beta/permdisp_referencia_vegan.R`: validação do PERMDISP.
- `scripts/alfa_beta/recalcular_biomassa_alfa_ASV.R`: correlações de Spearman (principal), Pearson (sensibilidade) e modelos ajustados para massa e diversidade alfa.
- `scripts/abundancia_taxonomica/screen_taxa_permutacao_lote.R`: triagem por estatística F das abundâncias relativas de famílias e gêneros, com permutações restritas por lote.
- `scripts/ANCOM-BC2/ANCOMBC2_analysis.R`: análise complementar com filtro de prevalência `prv_cut = 0.10`, lote e cultivar no modelo.
- `scripts/PICRUSt2/reanalise_picrust2_bloqueada_ajustada.R`: teste exato do efeito parcial de cultivar nos KOs, após considerar lote.
- `scripts/PICRUSt2/MetaCyc_permutacao_duas_sementes.R`: testes das 534 vias MetaCyc em duas rodadas de 99.999 permutações.

## Scripts de figuras

- `scripts/figuras/curva_rarefacao_final.R`: Figura 2, curva de rarefação em 1.200 dpi, com os códigos das amostras e sem legenda de repetições.
- `scripts/figuras/gerar_figuras_03_04_05_taxonomia.R`: Figuras 3, 4 e 5, composição taxonômica por filo, família e gênero.
- `scripts/figuras/gerar_figura_06_alfa_ASV_validacao.R`: Figura 6, diversidade alfa.
- `scripts/figuras/gerar_figura_07_pcoa_ASV.R`: Figura 7, PCoA das matrizes de ASVs.
- `scripts/figuras/gerar_figura_08_massa_alfa.R`: Figura 8, massa úmida das raízes com solo aderido e diversidade alfa.
- `scripts/figuras/gerar_figura_09_pcoa_KO.R`: Figura 9, PCoA dos perfis de KOs preditos; usa apenas as 16 raízes, embora a entrada do PICRUSt2 contenha 17 bibliotecas.
- `scripts/figuras/metaCyc_heatmap.py`: Figura 10, mapa de calor das 20 vias com maiores estatísticas F entre as vias com q < 0,05 nas duas rodadas.

## Dados incluídos

- `dados/MetaCyc_bacteria_only_blocked_2seeds.csv`: resultados das 534 vias, estatísticas F, valores de p e q nas duas rodadas e médias por cultivar.
- `dados/ANCOMBC2_apendice_family_sinais_robustos.csv` e `dados/ANCOMBC2_apendice_genus_sinais_robustos.csv`: sinais globais robustos e contrastes pareados robustos.
- `dados/entradas/metadados_16_raizes.tsv` e `dados/entradas/matriz_bacterias_ASV_16_raizes.tsv`: metadados e matriz de contagens para as 16 amostras radiculares.
- `dados/entradas/alfa_bacteriana/`: resultados de referência da rarefação alfa usados durante a revisão; a rarefação pode ser refeita a partir da matriz de contagens pelo script acima.
- `dados/entradas/permanova_bloqueada/`: matrizes de distância e coordenadas da PCoA para as análises de ASVs.
- `dados/entradas/sensibilidade_genero/`: matriz agregada em nível de gênero.
- `dados/entradas/abundancia_*_integrada.csv`: tabelas resumidas usadas nos gráficos taxonômicos.
- `dados/entradas/picrust2/`: abundâncias preditas de KOs e vias MetaCyc e NSTI ponderado.
- `dados/entradas/dados_rarefacao.xlsx`: dados da curva de rarefação.

As leituras brutas de sequenciamento não estão incluídas. As tabelas de abundância predita pelo PICRUSt2 estão incluídas como entradas: a execução completa da predição não é refeita por estes scripts. A predição fornecida contém 17 bibliotecas (solo-controle e 16 raízes); as análises inferenciais descritas no TCC utilizam somente as 16 raízes.

## Ordem de execução

1. Execute `scripts/alfa_beta/recalcular_alfa_bacteriana.R`; em seguida, `teste_alfa_exato_bloqueado.R` e a Figura 6.
2. Execute `scripts/alfa_beta/recalcular_biomassa_alfa_ASV.R` antes de gerar a Figura 8.
3. Execute os scripts PERMANOVA e PERMDISP para reproduzir os testes beta; as entradas são as matrizes de distância fornecidas.
4. Execute os scripts de análise PICRUSt2 para refazer os testes de KOs e MetaCyc a partir das abundâncias preditas incluídas; rode `metaCyc_heatmap.py` após a análise MetaCyc.
5. A curva de rarefação, gráficos taxonômicos e PCoA de KOs usam os arquivos de entrada fornecidos.

Execute scripts R com `Rscript caminho/do_script.R`. Os caminhos são resolvidos a partir da estrutura do pacote. O script Python da Figura 10 também resolve seus arquivos a partir da própria localização.

## Ambiente e estado atual

O arquivo `ambiente/sessionInfo_ANCOMBC2.txt` registra R 4.6.1, ANCOMBC 2.14.0, phyloseq 1.56.0, vegan 2.7-6 e permute 0.9-10 para a execução do ANCOM-BC2. A predição PICRUSt2 foi realizada com PICRUSt2 2.6.2 em ambiente conda/WSL com Python 3.9. Instale as dependências indicadas nos scripts antes de executá-los.

Parte dos scripts foi executada a partir desta pasta; os resultados coincidiram com as tabelas do TCC. As figuras 2–5 e 8 e o ANCOM-BC2 ficaram bloqueados por dependências ausentes; a Figura 10 ficou bloqueada pela ausência de Python. Consulte `VALIDACAO_EXECUCAO.md` para os detalhes. Portanto, a validacao ainda nao cobre todas as analises ou figuras. Antes de citá-lo como material suplementar público, conferir as saídas geradas, registrar as dependências completas e depositar o pacote em repositório estável (por exemplo, Zenodo ou GitHub com versão arquivada). Depois do depósito, substituir o caminho local por URL ou DOI permanente.

