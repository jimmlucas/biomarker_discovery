# =========================
# 36_build_pseudobulk_gse272993.R
# =========================
# Construcción de matriz pseudobulk CD8+ por paciente/timepoint.
# Se agregan cuentas crudas de todas las células CD8+ de cada muestra clínica.

base_dir <- "/Users/jimmlucas/Downloads/R_MASTER_RNASEQ"
setwd(base_dir)

source(file.path(base_dir, "00_config.R"))

library(Seurat)
library(tidyverse)
library(data.table)
library(Matrix)

set.seed(123)

# -------------------------
# 1. Cargar objeto GSE272993
# -------------------------

obj_file_proc <- file.path(proc_dir, "GSE272993_cd8_nn_labeled_FINAL.rds")

obj_file_raw <- file.path(
  raw_dir,
  "GSE272993",
  "GSE272993_cd8_nn_labeled_FINAL.RDS"
)

if (file.exists(obj_file_proc)) {
  obj <- readRDS(obj_file_proc)
} else if (file.exists(obj_file_raw)) {
  obj <- readRDS(obj_file_raw)
} else {
  stop("No se encuentra el objeto Seurat de GSE272993.")
}

DefaultAssay(obj) <- "RNA"

# Para Seurat v5 con varias layers
obj <- tryCatch(
  JoinLayers(obj),
  error = function(e) {
    message("JoinLayers no aplicado: ", e$message)
    obj
  }
)

# -------------------------
# 2. Preparar metadata
# -------------------------

meta <- obj@meta.data %>%
  rownames_to_column("cell_id") %>%
  as_tibble()

if (!"response_binary" %in% colnames(meta)) {
  meta <- meta %>%
    mutate(
      response_binary = case_when(
        response %in% c("CR", "PR") ~ "Responder",
        response %in% c("PD", "SD") ~ "Non_responder",
        TRUE ~ NA_character_
      )
    )
}

meta <- meta %>%
  filter(!is.na(response_binary)) %>%
  mutate(
    pseudobulk_sample_id = paste(
      patient_alias,
      timepoint,
      treatment,
      response_binary,
      sep = "__"
    )
  )

obj$pseudobulk_sample_id <- meta$pseudobulk_sample_id[
  match(colnames(obj), meta$cell_id)
]

# -------------------------
# 3. Extraer matriz de cuentas crudas
# -------------------------

get_counts_matrix <- function(seurat_obj) {
  
  mat <- tryCatch(
    GetAssayData(
      seurat_obj,
      assay = DefaultAssay(seurat_obj),
      layer = "counts"
    ),
    error = function(e) {
      GetAssayData(
        seurat_obj,
        assay = DefaultAssay(seurat_obj),
        slot = "counts"
      )
    }
  )
  
  return(mat)
}

counts <- get_counts_matrix(obj)

valid_cells <- meta$cell_id[meta$cell_id %in% colnames(counts)]

counts <- counts[, valid_cells, drop = FALSE]

meta <- meta %>%
  filter(cell_id %in% valid_cells)

stopifnot(ncol(counts) == nrow(meta))

# -------------------------
# 4. Construir pseudobulk por paciente/timepoint
# -------------------------

sample_ids <- unique(meta$pseudobulk_sample_id)

pb_list <- vector("list", length(sample_ids))
names(pb_list) <- sample_ids

for (sid in sample_ids) {
  
  message("Agregando muestra: ", sid)
  
  cells_sid <- meta %>%
    filter(pseudobulk_sample_id == sid) %>%
    pull(cell_id)
  
  pb_list[[sid]] <- Matrix::rowSums(
    counts[, cells_sid, drop = FALSE]
  )
}

pseudobulk_counts <- do.call(cbind, pb_list)
pseudobulk_counts <- Matrix(pseudobulk_counts, sparse = TRUE)

colnames(pseudobulk_counts) <- sample_ids
rownames(pseudobulk_counts) <- rownames(counts)

# -------------------------
# 5. Metadata pseudobulk
# -------------------------

pseudobulk_meta <- meta %>%
  group_by(pseudobulk_sample_id) %>%
  summarise(
    patient_alias = first(patient_alias),
    response = first(response),
    response_binary = first(response_binary),
    treatment = first(treatment),
    timepoint = first(timepoint),
    n_cells = n(),
    median_nCount_RNA = median(nCount_RNA, na.rm = TRUE),
    median_nFeature_RNA = median(nFeature_RNA, na.rm = TRUE),
    median_percent_mt = median(percent.mt, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    total_counts_pseudobulk = Matrix::colSums(pseudobulk_counts)[pseudobulk_sample_id],
    n_genes_detected_pseudobulk = Matrix::colSums(pseudobulk_counts > 0)[pseudobulk_sample_id],
    y = ifelse(response_binary == "Responder", 1, 0)
  )

# Asegurar mismo orden
pseudobulk_meta <- pseudobulk_meta %>%
  slice(match(colnames(pseudobulk_counts), pseudobulk_sample_id))

stopifnot(all(pseudobulk_meta$pseudobulk_sample_id == colnames(pseudobulk_counts)))

# -------------------------
# 6. Filtrado básico de genes
# -------------------------

gene_detection <- Matrix::rowSums(pseudobulk_counts > 0)
gene_total_counts <- Matrix::rowSums(pseudobulk_counts)

gene_filter_summary <- tibble(
  gene = rownames(pseudobulk_counts),
  n_samples_detected = gene_detection,
  total_counts = gene_total_counts
)

keep_genes <- gene_filter_summary %>%
  filter(
    n_samples_detected >= 5,
    total_counts >= 20
  ) %>%
  pull(gene)

pseudobulk_counts_filtered <- pseudobulk_counts[keep_genes, , drop = FALSE]

gene_filter_report <- tibble(
  dataset = "GSE272993",
  n_genes_initial = nrow(pseudobulk_counts),
  n_genes_after_filter = nrow(pseudobulk_counts_filtered),
  n_pseudobulk_samples = ncol(pseudobulk_counts_filtered),
  n_patients = n_distinct(pseudobulk_meta$patient_alias),
  n_responders = sum(pseudobulk_meta$response_binary == "Responder"),
  n_non_responders = sum(pseudobulk_meta$response_binary == "Non_responder")
)

# -------------------------
# 7. Guardar resultados
# -------------------------

saveRDS(
  pseudobulk_counts,
  file.path(proc_dir, "GSE272993_pseudobulk_counts_raw.rds")
)

saveRDS(
  pseudobulk_counts_filtered,
  file.path(proc_dir, "GSE272993_pseudobulk_counts_filtered.rds")
)

fwrite(
  pseudobulk_meta,
  file.path(res_dir, "GSE272993_pseudobulk_metadata.tsv"),
  sep = "\t"
)

fwrite(
  gene_filter_summary,
  file.path(res_dir, "GSE272993_pseudobulk_gene_filter_summary.tsv"),
  sep = "\t"
)

fwrite(
  gene_filter_report,
  file.path(res_dir, "GSE272993_pseudobulk_gene_filter_report.tsv"),
  sep = "\t"
)

print(gene_filter_report)
print(pseudobulk_meta %>% count(response_binary, timepoint))

# -------------------------
# 8. Figuras QC pseudobulk
# -------------------------

p_cells <- pseudobulk_meta %>%
  ggplot(aes(x = response_binary, y = n_cells, fill = response_binary)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.6) +
  theme_classic() +
  labs(
    x = "Respuesta clínica",
    y = "Número de células agregadas",
    title = "Número de células CD8+ por pseudobulk"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(fig_dir, "GSE272993_pseudobulk_n_cells_by_response.png"),
  p_cells,
  width = 6,
  height = 5,
  dpi = 300
)

p_counts <- pseudobulk_meta %>%
  ggplot(aes(x = response_binary, y = total_counts_pseudobulk, fill = response_binary)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.6) +
  scale_y_log10() +
  theme_classic() +
  labs(
    x = "Respuesta clínica",
    y = "Cuentas totales pseudobulk",
    title = "Profundidad pseudobulk por respuesta"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(fig_dir, "GSE272993_pseudobulk_total_counts_by_response.png"),
  p_counts,
  width = 6,
  height = 5,
  dpi = 300
)

p_genes <- pseudobulk_meta %>%
  ggplot(aes(x = response_binary, y = n_genes_detected_pseudobulk, fill = response_binary)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.6) +
  theme_classic() +
  labs(
    x = "Respuesta clínica",
    y = "Genes detectados",
    title = "Genes detectados por pseudobulk"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(fig_dir, "GSE272993_pseudobulk_detected_genes_by_response.png"),
  p_genes,
  width = 6,
  height = 5,
  dpi = 300
)

message("Pseudobulk de GSE272993 construido correctamente.")