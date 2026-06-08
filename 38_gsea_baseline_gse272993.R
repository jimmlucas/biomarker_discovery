# =========================
# 38_gsea_baseline_gse272993.R
# =========================
# GSEA sobre resultados DESeq2 Baseline CD8+ GSE272993
# Ranking: estadístico Wald de DESeq2
# Comparación: Responder vs Non_responder
# =========================

base_dir <- "/Users/jimmlucas/Downloads/R_MASTER_RNASEQ"
setwd(base_dir)

source(file.path(base_dir, "00_config.R"))

library(tidyverse)
library(data.table)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(ggplot2)

set.seed(123)

# -------------------------
# 1. Cargar resultados DESeq2
# -------------------------

deseq_file <- file.path(
  res_dir,
  "GSE272993_DESeq2_baseline_responder_vs_nonresponder.tsv"
)

if (!file.exists(deseq_file)) {
  stop("No existe el archivo DESeq2: ", deseq_file)
}

res_tbl <- fread(
  deseq_file,
  data.table = FALSE
) %>%
  as_tibble()

# -------------------------
# 2. Preparar ranking
# -------------------------

gene_ranking <- res_tbl %>%
  filter(!is.na(stat)) %>%
  filter(!duplicated(gene)) %>%
  arrange(desc(stat)) %>%
  select(gene, stat)

ranked_vector <- gene_ranking$stat
names(ranked_vector) <- gene_ranking$gene

ranked_vector <- sort(ranked_vector, decreasing = TRUE)

message("Genes en ranking GSEA:")
print(length(ranked_vector))

# -------------------------
# 3. Cargar MSigDB Hallmark
# -------------------------

hallmark_sets <- msigdbr(
  species = "Homo sapiens",
  category = "H"
) %>%
  select(gs_name, gene_symbol)

message("Colecciones Hallmark cargadas:")
print(length(unique(hallmark_sets$gs_name)))

# -------------------------
# 4. Ejecutar GSEA
# -------------------------

gsea_hallmark <- GSEA(
  geneList = ranked_vector,
  TERM2GENE = hallmark_sets,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  verbose = FALSE,
  seed = TRUE
)

gsea_tbl <- as.data.frame(gsea_hallmark) %>%
  as_tibble() %>%
  arrange(p.adjust) %>%
  mutate(
    direction = case_when(
      NES > 0 ~ "Enriched_in_Responder",
      NES < 0 ~ "Enriched_in_Non_responder",
      TRUE ~ "No_direction"
    )
  )

fwrite(
  gsea_tbl,
  file.path(res_dir, "GSE272993_GSEA_Hallmark_baseline.tsv"),
  sep = "\t"
)

# -------------------------
# 5. Resumen GSEA
# -------------------------

gsea_summary <- tibble(
  dataset = "GSE272993",
  analysis = "GSEA Hallmark CD8 pseudobulk Baseline Responder vs Non_responder",
  n_ranked_genes = length(ranked_vector),
  n_pathways_tested = nrow(gsea_tbl),
  n_pathways_FDR_0_05 = sum(gsea_tbl$p.adjust < 0.05, na.rm = TRUE),
  n_pathways_FDR_0_10 = sum(gsea_tbl$p.adjust < 0.10, na.rm = TRUE),
  n_enriched_responder_FDR_0_10 = sum(
    gsea_tbl$p.adjust < 0.10 & gsea_tbl$NES > 0,
    na.rm = TRUE
  ),
  n_enriched_non_responder_FDR_0_10 = sum(
    gsea_tbl$p.adjust < 0.10 & gsea_tbl$NES < 0,
    na.rm = TRUE
  )
)

fwrite(
  gsea_summary,
  file.path(res_dir, "GSE272993_GSEA_Hallmark_baseline_summary.tsv"),
  sep = "\t"
)

print(gsea_summary)

# -------------------------
# 6. Dotplot GSEA
# -------------------------

top_gsea <- gsea_tbl %>%
  filter(!is.na(p.adjust)) %>%
  arrange(p.adjust) %>%
  slice_head(n = 20)

p_dot <- top_gsea %>%
  mutate(
    Description = stringr::str_replace_all(ID, "HALLMARK_", ""),
    Description = stringr::str_replace_all(Description, "_", " "),
    Description = factor(Description, levels = rev(Description))
  ) %>%
  ggplot(aes(
    x = NES,
    y = Description,
    size = setSize,
    color = p.adjust
  )) +
  geom_point() +
  scale_color_viridis_c(direction = -1) +
  theme_classic() +
  labs(
    x = "NES",
    y = NULL,
    size = "Genes",
    color = "FDR",
    title = "GSEA Hallmark CD8+ Baseline - GSE272993"
  )

ggsave(
  file.path(fig_dir, "GSE272993_GSEA_Hallmark_baseline_dotplot.png"),
  p_dot,
  width = 8,
  height = 7,
  dpi = 300
)

# -------------------------
# 7. Ridgeplot / resumen enriquecimiento
# -------------------------

png(
  file.path(fig_dir, "GSE272993_GSEA_Hallmark_baseline_ridgeplot.png"),
  width = 2200,
  height = 1800,
  res = 300
)

print(
  ridgeplot(gsea_hallmark, showCategory = 20) +
    labs(title = "GSEA Hallmark CD8+ Baseline - GSE272993")
)

dev.off()

# -------------------------
# 8. Curvas GSEA individuales
# -------------------------

top_terms <- gsea_tbl %>%
  filter(!is.na(p.adjust)) %>%
  arrange(p.adjust) %>%
  slice_head(n = 6) %>%
  pull(ID)

for (term in top_terms) {
  
  safe_term <- gsub("[^A-Za-z0-9_]", "_", term)
  
  png(
    file.path(
      fig_dir,
      paste0("GSE272993_GSEA_curve_", safe_term, ".png")
    ),
    width = 2000,
    height = 1500,
    res = 300
  )
  
  print(
    gseaplot2(
      gsea_hallmark,
      geneSetID = term,
      title = term
    )
  )
  
  dev.off()
}

# -------------------------
# 9. Tabla compacta para resultados
# -------------------------

gsea_compact <- gsea_tbl %>%
  select(
    ID,
    Description,
    setSize,
    enrichmentScore,
    NES,
    pvalue,
    p.adjust,
    qvalue,
    direction,
    core_enrichment
  ) %>%
  arrange(p.adjust)

fwrite(
  gsea_compact,
  file.path(res_dir, "GSE272993_GSEA_Hallmark_baseline_compact_results.tsv"),
  sep = "\t"
)

message("GSEA Hallmark Baseline de GSE272993 completado correctamente.")