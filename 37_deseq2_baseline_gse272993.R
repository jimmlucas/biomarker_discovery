# =========================
# 37_deseq2_baseline_gse272993.R
# =========================
# Análisis diferencial pseudobulk CD8+ en GSE272993.
# Comparación principal: Responders vs Non-responders en Baseline.
# =========================

base_dir <- "/Users/jimmlucas/Downloads/R_MASTER_RNASEQ"
setwd(base_dir)

source(file.path(base_dir, "00_config.R"))

library(tidyverse)
library(data.table)
library(Matrix)
library(DESeq2)
library(ggrepel)

set.seed(123)

# -------------------------
# 1. Cargar datos pseudobulk
# -------------------------

counts_file <- file.path(proc_dir, "GSE272993_pseudobulk_counts_filtered.rds")
meta_file   <- file.path(res_dir, "GSE272993_pseudobulk_metadata.tsv")

if (!file.exists(counts_file)) {
  stop("No existe la matriz pseudobulk filtrada: ", counts_file)
}

if (!file.exists(meta_file)) {
  stop("No existe la metadata pseudobulk: ", meta_file)
}

pseudobulk_counts <- readRDS(counts_file)

pseudobulk_meta <- fread(
  meta_file,
  data.table = FALSE
) %>%
  as_tibble()

message("Matriz pseudobulk cargada:")
print(dim(pseudobulk_counts))

message("Metadata pseudobulk cargada:")
print(dim(pseudobulk_meta))

# -------------------------
# 2. Seleccionar muestras Baseline
# -------------------------

baseline_meta <- pseudobulk_meta %>%
  filter(timepoint == "Baseline") %>%
  mutate(
    response_binary = factor(
      response_binary,
      levels = c("Non_responder", "Responder")
    ),
    treatment = factor(treatment),
    patient_alias = factor(patient_alias)
  )

baseline_samples <- baseline_meta$pseudobulk_sample_id

baseline_counts <- pseudobulk_counts[, baseline_samples, drop = FALSE]

baseline_meta <- baseline_meta %>%
  arrange(match(pseudobulk_sample_id, colnames(baseline_counts)))

stopifnot(all(baseline_meta$pseudobulk_sample_id == colnames(baseline_counts)))

baseline_counts <- as.matrix(baseline_counts)
storage.mode(baseline_counts) <- "integer"

baseline_meta <- as.data.frame(baseline_meta)
rownames(baseline_meta) <- baseline_meta$pseudobulk_sample_id

message("Muestras Baseline:")
print(dplyr::count(as_tibble(baseline_meta), response_binary))

message("Tratamientos en Baseline:")
print(dplyr::count(as_tibble(baseline_meta), response_binary, treatment))

# -------------------------
# 3. Crear objeto DESeq2
# -------------------------

dds <- DESeqDataSetFromMatrix(
  countData = baseline_counts,
  colData = baseline_meta,
  design = ~ response_binary
)

# Filtrado adicional mínimo
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]

message("Genes retenidos para DESeq2:")
print(nrow(dds))

# -------------------------
# 4. Ejecutar DESeq2
# -------------------------

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("response_binary", "Responder", "Non_responder"),
  alpha = 0.05
)

res_tbl <- as.data.frame(res) %>%
  rownames_to_column("gene") %>%
  as_tibble() %>%
  arrange(padj) %>%
  mutate(
    significant_FDR_0_05 = ifelse(!is.na(padj) & padj < 0.05, TRUE, FALSE),
    significant_FDR_0_10 = ifelse(!is.na(padj) & padj < 0.10, TRUE, FALSE),
    direction = case_when(
      is.na(padj) ~ "NA",
      padj < 0.05 & log2FoldChange > 0 ~ "Higher_in_Responder",
      padj < 0.05 & log2FoldChange < 0 ~ "Higher_in_Non_responder",
      TRUE ~ "Not_significant"
    )
  )

fwrite(
  res_tbl,
  file.path(res_dir, "GSE272993_DESeq2_baseline_responder_vs_nonresponder.tsv"),
  sep = "\t"
)

# -------------------------
# 5. Resumen
# -------------------------

deseq_summary <- tibble(
  dataset = "GSE272993",
  analysis = "CD8 pseudobulk Baseline Responder vs Non_responder",
  n_baseline_samples = ncol(baseline_counts),
  n_baseline_patients = dplyr::n_distinct(baseline_meta$patient_alias),
  n_responders = sum(baseline_meta$response_binary == "Responder"),
  n_non_responders = sum(baseline_meta$response_binary == "Non_responder"),
  n_genes_input = nrow(pseudobulk_counts),
  n_genes_tested = nrow(dds),
  n_genes_FDR_0_05 = sum(res_tbl$significant_FDR_0_05, na.rm = TRUE),
  n_genes_FDR_0_10 = sum(res_tbl$significant_FDR_0_10, na.rm = TRUE),
  n_higher_in_responder_FDR_0_05 = sum(
    res_tbl$direction == "Higher_in_Responder",
    na.rm = TRUE
  ),
  n_higher_in_non_responder_FDR_0_05 = sum(
    res_tbl$direction == "Higher_in_Non_responder",
    na.rm = TRUE
  )
)

fwrite(
  deseq_summary,
  file.path(res_dir, "GSE272993_DESeq2_baseline_summary.tsv"),
  sep = "\t"
)

print(deseq_summary)

# -------------------------
# 6. PCA con VST
# -------------------------

vsd <- vst(dds, blind = TRUE)

pca_data <- plotPCA(
  vsd,
  intgroup = c("response_binary", "treatment"),
  returnData = TRUE
)

percent_var <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = response_binary,
    shape = treatment
  )
) +
  geom_point(size = 3, alpha = 0.85) +
  theme_classic() +
  labs(
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    color = "Respuesta clínica",
    shape = "Tratamiento",
    title = "PCA pseudobulk CD8+ Baseline - GSE272993"
  )

ggsave(
  file.path(fig_dir, "GSE272993_DESeq2_baseline_PCA.png"),
  p_pca,
  width = 7,
  height = 5,
  dpi = 300
)

# -------------------------
# 7. MA plot
# -------------------------

png(
  file.path(fig_dir, "GSE272993_DESeq2_baseline_MAplot.png"),
  width = 1800,
  height = 1400,
  res = 300
)

plotMA(
  res,
  ylim = c(-5, 5),
  main = "GSE272993 CD8+ Baseline: Responder vs Non-responder"
)

dev.off()

# -------------------------
# 8. Volcano plot
# -------------------------

volcano_tbl <- res_tbl %>%
  mutate(
    neg_log10_padj = -log10(padj),
    volcano_group = case_when(
      !is.na(padj) & padj < 0.05 & log2FoldChange > 0 ~ "Higher in Responder",
      !is.na(padj) & padj < 0.05 & log2FoldChange < 0 ~ "Higher in Non-responder",
      TRUE ~ "Not significant"
    )
  )

top_labels <- volcano_tbl %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  slice_head(n = 15)

p_volcano <- ggplot(
  volcano_tbl,
  aes(
    x = log2FoldChange,
    y = neg_log10_padj,
    color = volcano_group
  )
) +
  geom_point(alpha = 0.65, size = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_text_repel(
    data = top_labels,
    aes(label = gene),
    size = 3,
    max.overlaps = Inf
  ) +
  theme_classic() +
  labs(
    x = "log2 Fold Change: Responder vs Non-responder",
    y = "-log10(FDR)",
    color = "Grupo",
    title = "DESeq2 CD8+ Baseline - GSE272993"
  )

ggsave(
  file.path(fig_dir, "GSE272993_DESeq2_baseline_volcano.png"),
  p_volcano,
  width = 8,
  height = 6,
  dpi = 300
)

# -------------------------
# 9. Guardar objetos
# -------------------------

saveRDS(
  dds,
  file.path(proc_dir, "GSE272993_DESeq2_baseline_dds.rds")
)

saveRDS(
  vsd,
  file.path(proc_dir, "GSE272993_DESeq2_baseline_vsd.rds")
)

message("Análisis DESeq2 Baseline de GSE272993 completado correctamente.")