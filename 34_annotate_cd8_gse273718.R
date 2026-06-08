# =========================
# 34_annotate_cd8_gse273718.R
# =========================

base_dir <- "/Users/jimmlucas/Downloads/R_MASTER_RNASEQ"
setwd(base_dir)

source(file.path(base_dir, "00_config.R"))

library(Seurat)
library(tidyverse)
library(data.table)
library(Matrix)

set.seed(123)

# -------------------------
# 1. Cargar objeto
# -------------------------

obj <- readRDS(
  file.path(proc_dir, "GSE273718_seurat_processed.rds")
)

DefaultAssay(obj) <- "RNA"

# En Seurat v5 puede haber múltiples layers.
# JoinLayers reduce problemas de memoria y acceso.
obj <- JoinLayers(obj)

# -------------------------
# 2. Función ligera para score medio sparse
# -------------------------

add_sparse_score <- function(seurat_obj, genes, score_name) {
  
  genes <- intersect(genes, rownames(seurat_obj))
  
  if (length(genes) == 0) {
    warning("No genes found for ", score_name)
    seurat_obj[[score_name]] <- NA_real_
    return(seurat_obj)
  }
  
  mat <- GetAssayData(
    seurat_obj,
    assay = DefaultAssay(seurat_obj),
    layer = "data"
  )
  
  score <- Matrix::colMeans(mat[genes, , drop = FALSE])
  
  seurat_obj[[score_name]] <- as.numeric(score)
  
  return(seurat_obj)
}

# -------------------------
# 3. Marcadores generales
# -------------------------

marker_sets <- list(
  T_cell = c("CD3D", "CD3E", "CD3G", "TRAC"),
  CD8 = c("CD8A", "CD8B"),
  CD4 = c("CD4", "IL7R", "CCR7"),
  NK = c("NKG7", "GNLY", "KLRD1"),
  B_cell = c("MS4A1", "CD79A"),
  Myeloid = c("LYZ", "S100A8", "S100A9", "FCGR3A"),
  Tumor_or_epithelial = c("EPCAM", "MLANA", "PMEL", "TYR"),
  Exhaustion = c("PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX"),
  Effector = c("GZMB", "PRF1", "IFNG", "NKG7"),
  Memory = c("IL7R", "CCR7", "TCF7", "LEF1"),
  Cycling = c("MKI67", "TOP2A", "STMN1")
)

available_marker_sets <- lapply(marker_sets, function(x) intersect(x, rownames(obj)))
print(available_marker_sets)

for (nm in names(marker_sets)) {
  obj <- add_sparse_score(
    seurat_obj = obj,
    genes = marker_sets[[nm]],
    score_name = paste0("score_", nm)
  )
}

# -------------------------
# 4. Anotación celular amplia
# -------------------------

meta <- obj@meta.data %>%
  rownames_to_column("cell_id") %>%
  mutate(
    broad_celltype = case_when(
      score_T_cell > 0.20 & score_CD8 > 0.05 ~ "CD8_T",
      score_T_cell > 0.20 & score_CD4 > 0.05 ~ "CD4_T",
      score_NK > 0.25 ~ "NK",
      score_B_cell > 0.15 ~ "B_cell",
      score_Myeloid > 0.25 ~ "Myeloid",
      score_Tumor_or_epithelial > 0.15 ~ "Tumor_or_epithelial",
      score_T_cell > 0.20 ~ "Other_T",
      TRUE ~ "Other"
    )
  )

obj$broad_celltype <- meta$broad_celltype[match(colnames(obj), meta$cell_id)]

# -------------------------
# 5. Subset CD8
# -------------------------

cd8_obj <- subset(obj, subset = broad_celltype == "CD8_T")

message("Número de células CD8 detectadas: ", ncol(cd8_obj))

if (ncol(cd8_obj) < 50) {
  warning("Pocas células CD8 detectadas. Revisa umbrales.")
}

cd8_obj <- JoinLayers(cd8_obj)

cd8_obj <- NormalizeData(cd8_obj, verbose = FALSE)
cd8_obj <- FindVariableFeatures(cd8_obj, selection.method = "vst", nfeatures = 2000)
cd8_obj <- ScaleData(cd8_obj, verbose = FALSE)
cd8_obj <- RunPCA(cd8_obj, npcs = 20, verbose = FALSE)
cd8_obj <- FindNeighbors(cd8_obj, dims = 1:15)
cd8_obj <- FindClusters(cd8_obj, resolution = 0.4)
cd8_obj <- RunUMAP(cd8_obj, dims = 1:15)

# -------------------------
# 6. Programas CD8 comparables
# -------------------------

cd8_programs <- list(
  exhaustion = c("PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX"),
  effector = c("GZMB", "PRF1", "IFNG", "NKG7"),
  memory = c("IL7R", "CCR7", "TCF7", "LEF1"),
  activation = c("CD69", "HLA-DRA"),
  cycling = c("MKI67", "TOP2A", "STMN1")
)

available_cd8_programs <- lapply(cd8_programs, function(x) intersect(x, rownames(cd8_obj)))
print(available_cd8_programs)

for (nm in names(cd8_programs)) {
  cd8_obj <- add_sparse_score(
    seurat_obj = cd8_obj,
    genes = cd8_programs[[nm]],
    score_name = paste0("cd8_score_", nm)
  )
}

# -------------------------
# 7. Estado funcional CD8 dominante
# -------------------------

cd8_meta <- cd8_obj@meta.data %>%
  rownames_to_column("cell_id") %>%
  mutate(
    cd8_state = case_when(
      cd8_score_cycling > 0.20 ~ "Cycling",
      cd8_score_exhaustion > 0.12 ~ "Exhaustion",
      cd8_score_effector > 0.20 ~ "Effector",
      cd8_score_memory > 0.12 ~ "Memory",
      cd8_score_activation > 0.12 ~ "Activated",
      TRUE ~ "Other_CD8"
    )
  )

cd8_obj$cd8_state <- cd8_meta$cd8_state[match(colnames(cd8_obj), cd8_meta$cell_id)]

# -------------------------
# 8. Tablas resumen
# -------------------------

celltype_summary <- obj@meta.data %>%
  as_tibble() %>%
  count(tissue, broad_celltype, name = "n_cells") %>%
  group_by(tissue) %>%
  mutate(freq = n_cells / sum(n_cells)) %>%
  ungroup() %>%
  arrange(tissue, desc(freq))

cd8_state_summary <- cd8_obj@meta.data %>%
  as_tibble() %>%
  count(tissue, patient_alias, cd8_state, name = "n_cells") %>%
  group_by(tissue, patient_alias) %>%
  mutate(freq = n_cells / sum(n_cells)) %>%
  ungroup() %>%
  arrange(tissue, patient_alias, desc(freq))

cd8_program_patient <- cd8_obj@meta.data %>%
  as_tibble() %>%
  group_by(patient_alias, tissue, sample_id) %>%
  summarise(
    n_cd8_cells = n(),
    across(
      starts_with("cd8_score_"),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        median = ~ median(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

fwrite(
  celltype_summary,
  file.path(res_dir, "GSE273718_broad_celltype_summary.tsv"),
  sep = "\t"
)

fwrite(
  cd8_state_summary,
  file.path(res_dir, "GSE273718_cd8_state_summary.tsv"),
  sep = "\t"
)

fwrite(
  cd8_program_patient,
  file.path(res_dir, "GSE273718_cd8_program_patient_summary.tsv"),
  sep = "\t"
)

print(celltype_summary, n = Inf)
print(cd8_state_summary, n = Inf)
print(cd8_program_patient, width = Inf)

# -------------------------
# 9. Guardar objetos
# -------------------------

saveRDS(
  obj,
  file.path(proc_dir, "GSE273718_seurat_annotated.rds")
)

saveRDS(
  cd8_obj,
  file.path(proc_dir, "GSE273718_cd8_annotated.rds")
)

fwrite(
  obj@meta.data %>% rownames_to_column("cell_id"),
  file.path(res_dir, "GSE273718_metadata_annotated.tsv"),
  sep = "\t"
)

fwrite(
  cd8_obj@meta.data %>% rownames_to_column("cell_id"),
  file.path(res_dir, "GSE273718_cd8_metadata_annotated.tsv"),
  sep = "\t"
)

# -------------------------
# 10. Figuras generales
# -------------------------

p_umap_celltype <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "broad_celltype",
  pt.size = 0.2,
  label = TRUE
) +
  ggtitle("GSE273718 - anotación celular amplia")

ggsave(
  file.path(fig_dir, "GSE273718_umap_broad_celltype.png"),
  p_umap_celltype,
  width = 8,
  height = 6,
  dpi = 300
)

p_umap_tissue <- DimPlot(
  obj,
  reduction = "umap",
  group.by = "tissue",
  pt.size = 0.2
) +
  ggtitle("GSE273718 - PBMC vs TIL")

ggsave(
  file.path(fig_dir, "GSE273718_umap_tissue_annotated.png"),
  p_umap_tissue,
  width = 8,
  height = 6,
  dpi = 300
)

# -------------------------
# 11. Figuras CD8
# -------------------------

p_cd8_state <- DimPlot(
  cd8_obj,
  reduction = "umap",
  group.by = "cd8_state",
  label = TRUE,
  pt.size = 0.35
) +
  ggtitle("GSE273718 CD8+ - estados funcionales")

ggsave(
  file.path(fig_dir, "GSE273718_cd8_umap_states.png"),
  p_cd8_state,
  width = 8,
  height = 6,
  dpi = 300
)

p_cd8_tissue <- DimPlot(
  cd8_obj,
  reduction = "umap",
  group.by = "tissue",
  pt.size = 0.35
) +
  ggtitle("GSE273718 CD8+ - PBMC vs TIL")

ggsave(
  file.path(fig_dir, "GSE273718_cd8_umap_tissue.png"),
  p_cd8_tissue,
  width = 8,
  height = 6,
  dpi = 300
)

# -------------------------
# 12. Composición CD8 PBMC vs TIL
# -------------------------

p_cd8_comp <- cd8_state_summary %>%
  group_by(tissue, cd8_state) %>%
  summarise(
    mean_freq = mean(freq, na.rm = TRUE),
    sd_freq = sd(freq, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = cd8_state, y = mean_freq, fill = tissue)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  theme_classic() +
  labs(
    x = "Estado CD8+",
    y = "Frecuencia media",
    title = "Composición de estados CD8+ en PBMC y TIL"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(fig_dir, "GSE273718_cd8_state_composition_pbmc_til.png"),
  p_cd8_comp,
  width = 8,
  height = 5,
  dpi = 300
)

# -------------------------
# 13. DotPlot CD8
# -------------------------

marker_plot_genes <- unique(c(
  "CD3D", "CD3E", "CD8A", "CD8B",
  "IL7R", "CCR7", "TCF7",
  "GZMB", "PRF1", "NKG7",
  "PDCD1", "LAG3", "TIGIT", "HAVCR2",
  "MKI67", "TOP2A"
))

marker_plot_genes <- intersect(marker_plot_genes, rownames(cd8_obj))

p_dot <- DotPlot(
  cd8_obj,
  features = marker_plot_genes,
  group.by = "cd8_state"
) +
  RotatedAxis() +
  ggtitle("Marcadores funcionales de estados CD8+ en GSE273718")

ggsave(
  file.path(fig_dir, "GSE273718_cd8_marker_dotplot.png"),
  p_dot,
  width = 11,
  height = 5,
  dpi = 300
)

message("Anotación CD8 funcional de GSE273718 completada.")