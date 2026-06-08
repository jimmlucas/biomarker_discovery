# =========================
# 33_create_seurat_gse273718.R
# =========================
# Importa matrices 10x de GSE273718, crea objeto Seurat,
# añade metadata y genera UMAP básico para resultados.

base_dir <- "/Users/jimmlucas/Downloads/R_MASTER_RNASEQ"
setwd(base_dir)

source(file.path(base_dir, "00_config.R"))

library(Seurat)
library(tidyverse)
library(data.table)

set.seed(123)

gse_id <- "GSE273718"
extract_dir <- file.path(raw_dir, gse_id, "extracted")

stopifnot(dir.exists(extract_dir))

# -------------------------
# 1. Detectar matrices
# -------------------------

matrix_files <- list.files(
  extract_dir,
  pattern = "_matrix\\.mtx\\.gz$",
  full.names = TRUE
)

if (length(matrix_files) == 0) {
  stop("No se encontraron archivos matrix.mtx.gz en: ", extract_dir)
}

sample_info <- tibble(
  matrix_file = matrix_files,
  prefix = str_remove(basename(matrix_file), "_matrix\\.mtx\\.gz$"),
  barcode_file = file.path(extract_dir, paste0(prefix, "_barcodes.tsv.gz")),
  feature_file = file.path(extract_dir, paste0(prefix, "_features.tsv.gz"))
) %>%
  mutate(
    gsm = str_extract(prefix, "GSM[0-9]+"),
    patient_alias = paste0("P", str_extract(prefix, "(?<=_)[0-9]+(?=_pre)")),
    timepoint = "Baseline",
    tissue = case_when(
      str_detect(prefix, "_PBMC") ~ "PBMC",
      str_detect(prefix, "_TIL") ~ "TIL",
      TRUE ~ NA_character_
    ),
    sample_id = paste(patient_alias, timepoint, tissue, sep = "_")
  )

stopifnot(all(file.exists(sample_info$barcode_file)))
stopifnot(all(file.exists(sample_info$feature_file)))

fwrite(
  sample_info,
  file.path(res_dir, "GSE273718_sample_info.tsv"),
  sep = "\t"
)

print(sample_info, n = Inf)

# -------------------------
# 2. Leer cada muestra como Seurat
# -------------------------

seurat_list <- list()

for (i in seq_len(nrow(sample_info))) {
  
  message("Leyendo muestra: ", sample_info$sample_id[i])
  
  counts <- ReadMtx(
    mtx = sample_info$matrix_file[i],
    cells = sample_info$barcode_file[i],
    features = sample_info$feature_file[i],
    feature.column = 2,
    unique.features = TRUE
  )
  
  obj <- CreateSeuratObject(
    counts = counts,
    project = "GSE273718",
    min.cells = 3,
    min.features = 200
  )
  
  obj$dataset <- "GSE273718"
  obj$gsm <- sample_info$gsm[i]
  obj$patient_alias <- sample_info$patient_alias[i]
  obj$timepoint <- sample_info$timepoint[i]
  obj$tissue <- sample_info$tissue[i]
  obj$sample_id <- sample_info$sample_id[i]
  
  obj <- RenameCells(
    obj,
    add.cell.id = sample_info$sample_id[i]
  )
  
  seurat_list[[sample_info$sample_id[i]]] <- obj
  
  rm(counts, obj)
  gc()
}

# -------------------------
# 3. Unir muestras
# -------------------------

gse273718 <- merge(
  seurat_list[[1]],
  y = seurat_list[-1],
  project = "GSE273718"
)

rm(seurat_list)
gc()

# -------------------------
# 4. QC básico
# -------------------------

gse273718[["percent.mt"]] <- PercentageFeatureSet(
  gse273718,
  pattern = "^MT-"
)

qc_summary <- gse273718@meta.data %>%
  as_tibble() %>%
  summarise(
    dataset = "GSE273718",
    n_cells = n(),
    n_patients = n_distinct(patient_alias),
    n_samples = n_distinct(sample_id),
    median_nCount_RNA = median(nCount_RNA, na.rm = TRUE),
    median_nFeature_RNA = median(nFeature_RNA, na.rm = TRUE),
    median_percent_mt = median(percent.mt, na.rm = TRUE)
  )

fwrite(
  qc_summary,
  file.path(res_dir, "GSE273718_qc_summary.tsv"),
  sep = "\t"
)

print(qc_summary)

cells_by_sample <- gse273718@meta.data %>%
  as_tibble() %>%
  count(patient_alias, tissue, sample_id, name = "n_cells") %>%
  arrange(patient_alias, tissue)

fwrite(
  cells_by_sample,
  file.path(res_dir, "GSE273718_cells_by_sample.tsv"),
  sep = "\t"
)

print(cells_by_sample, n = Inf)

# -------------------------
# 5. Filtrado QC conservador
# -------------------------

gse273718 <- subset(
  gse273718,
  subset = nFeature_RNA >= 200 &
    nFeature_RNA <= 6000 &
    percent.mt <= 20
)

qc_summary_filtered <- gse273718@meta.data %>%
  as_tibble() %>%
  summarise(
    dataset = "GSE273718",
    n_cells_after_qc = n(),
    n_patients = n_distinct(patient_alias),
    n_samples = n_distinct(sample_id),
    median_nCount_RNA = median(nCount_RNA, na.rm = TRUE),
    median_nFeature_RNA = median(nFeature_RNA, na.rm = TRUE),
    median_percent_mt = median(percent.mt, na.rm = TRUE)
  )

fwrite(
  qc_summary_filtered,
  file.path(res_dir, "GSE273718_qc_summary_filtered.tsv"),
  sep = "\t"
)

print(qc_summary_filtered)

# -------------------------
# 6. Procesamiento Seurat básico
# -------------------------

gse273718 <- NormalizeData(gse273718)
gse273718 <- FindVariableFeatures(gse273718, selection.method = "vst", nfeatures = 3000)
gse273718 <- ScaleData(gse273718, verbose = FALSE)
gse273718 <- RunPCA(gse273718, npcs = 30, verbose = FALSE)
gse273718 <- FindNeighbors(gse273718, dims = 1:20)
gse273718 <- FindClusters(gse273718, resolution = 0.5)
gse273718 <- RunUMAP(gse273718, dims = 1:20)

# -------------------------
# 7. Guardar objeto
# -------------------------

saveRDS(
  gse273718,
  file.path(proc_dir, "GSE273718_seurat_processed.rds")
)

# -------------------------
# 8. Figuras QC y UMAP
# -------------------------

p_qc <- VlnPlot(
  gse273718,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "tissue",
  pt.size = 0.05,
  ncol = 3
)

ggsave(
  file.path(fig_dir, "GSE273718_qc_violin_by_tissue.png"),
  p_qc,
  width = 12,
  height = 5,
  dpi = 300
)

p_umap_tissue <- DimPlot(
  gse273718,
  reduction = "umap",
  group.by = "tissue",
  pt.size = 0.25
) +
  ggtitle("GSE273718 UMAP por tejido")

ggsave(
  file.path(fig_dir, "GSE273718_umap_by_tissue.png"),
  p_umap_tissue,
  width = 7,
  height = 6,
  dpi = 300
)

p_umap_patient <- DimPlot(
  gse273718,
  reduction = "umap",
  group.by = "patient_alias",
  pt.size = 0.25
) +
  ggtitle("GSE273718 UMAP por paciente")

ggsave(
  file.path(fig_dir, "GSE273718_umap_by_patient.png"),
  p_umap_patient,
  width = 7,
  height = 6,
  dpi = 300
)

p_umap_cluster <- DimPlot(
  gse273718,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  pt.size = 0.25
) +
  ggtitle("GSE273718 UMAP por clusters Seurat")

ggsave(
  file.path(fig_dir, "GSE273718_umap_by_cluster.png"),
  p_umap_cluster,
  width = 7,
  height = 6,
  dpi = 300
)

# -------------------------
# 9. Guardar metadata
# -------------------------

meta_273718 <- gse273718@meta.data %>%
  rownames_to_column("cell_id")

fwrite(
  meta_273718,
  file.path(res_dir, "GSE273718_metadata_processed.tsv"),
  sep = "\t"
)

saveRDS(
  meta_273718,
  file.path(proc_dir, "GSE273718_metadata_processed.rds")
)

message("Objeto Seurat de GSE273718 creado y guardado.")