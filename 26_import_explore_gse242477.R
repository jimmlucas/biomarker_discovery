##26_import_explore_gse_242477

base_dir <- "/Users/jimmlucas/Downloads/R_MASTER_RNASEQ"
setwd(base_dir)

source(file.path(base_dir, "00_config.R"))

library(tidyverse)
library(data.table)
library(Seurat)
library(stringr)

data_dir <- file.path(base_dir, "data/raw/GSE242477")
out_dir  <- file.path(base_dir, "data/processed/GSE242477")
res_dir  <- file.path(base_dir, "results/GSE242477")
fig_dir  <- file.path(base_dir, "figures/GSE242477")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Detectar y descomprimir TAR si existe
# ------------------------------------------------------------

files <- list.files(data_dir, recursive = TRUE, full.names = TRUE)

tar_files <- files[grepl("\\.tar$", files)]

if (length(tar_files) > 0) {
  for (tf in tar_files) {
    message("Untarring: ", tf)
    untar(tf, exdir = data_dir)
  }
}

files <- list.files(data_dir, recursive = TRUE, full.names = TRUE)

writeLines(
  files,
  file.path(res_dir, "files_detected_after_untar.txt")
)

# ------------------------------------------------------------
# 2. Localizar archivos 10x por tipo
# ------------------------------------------------------------

barcode_files <- files[grepl("barcodes\\.tsv\\.gz$", files)]
feature_files <- files[grepl("features\\.tsv\\.gz$", files)]
matrix_files  <- files[grepl("matrix\\.mtx\\.gz$", files)]

message("Barcode files: ", length(barcode_files))
message("Feature files: ", length(feature_files))
message("Matrix files: ", length(matrix_files))

if (length(matrix_files) == 0) {
  stop("No matrix.mtx.gz files found.")
}

# ------------------------------------------------------------
# 3. Extraer metadata de nombre de archivo
# Ejemplo:
# GSM7764406_1_MEL03_TIL_barcodes.tsv.gz
# GSM7764406_4_MEL03_TIL_matrix.mtx.gz
# ------------------------------------------------------------

parse_file_info <- function(path) {
  fname <- basename(path)
  
  tibble(
    path = path,
    file = fname,
    gsm = str_extract(fname, "GSM[0-9]+"),
    patient_id = str_extract(fname, "MEL[0-9]+"),
    compartment = case_when(
      str_detect(fname, "_TIL_") ~ "TIL",
      str_detect(fname, "_PBMC_") ~ "PBMC",
      TRUE ~ NA_character_
    ),
    file_type = case_when(
      str_detect(fname, "barcodes\\.tsv\\.gz$") ~ "barcodes",
      str_detect(fname, "features\\.tsv\\.gz$") ~ "features",
      str_detect(fname, "matrix\\.mtx\\.gz$") ~ "matrix",
      TRUE ~ NA_character_
    ),
    sample_id = paste(patient_id, compartment, sep = "_")
  )
}

file_info <- map_dfr(
  c(barcode_files, feature_files, matrix_files),
  parse_file_info
)

write.csv(
  file_info,
  file.path(res_dir, "gse242477_10x_file_info.csv"),
  row.names = FALSE
)

print(file_info)

# ------------------------------------------------------------
# 4. Construir carpetas temporales 10x estándar por muestra
# ------------------------------------------------------------

tmp_10x_dir <- file.path(out_dir, "tmp_10x_by_sample")
dir.create(tmp_10x_dir, recursive = TRUE, showWarnings = FALSE)

valid_samples <- file_info %>%
  filter(!is.na(sample_id)) %>%
  group_by(sample_id, patient_id, compartment) %>%
  summarise(
    n_files = n_distinct(file_type),
    has_barcodes = any(file_type == "barcodes"),
    has_features = any(file_type == "features"),
    has_matrix = any(file_type == "matrix"),
    .groups = "drop"
  ) %>%
  filter(has_barcodes, has_features, has_matrix)

write.csv(
  valid_samples,
  file.path(res_dir, "gse242477_valid_10x_samples.csv"),
  row.names = FALSE
)

print(valid_samples)

if (nrow(valid_samples) == 0) {
  stop("No valid 10x samples found with barcodes + features + matrix.")
}

# ------------------------------------------------------------
# 5. Leer cada muestra y crear objeto Seurat
# ------------------------------------------------------------

seurat_list <- list()

for (i in seq_len(nrow(valid_samples))) {
  
  sid <- valid_samples$sample_id[i]
  pid <- valid_samples$patient_id[i]
  comp <- valid_samples$compartment[i]
  
  message("Processing sample: ", sid)
  
  sample_files <- file_info %>%
    filter(sample_id == sid)
  
  barcode_path <- sample_files$path[sample_files$file_type == "barcodes"][1]
  feature_path <- sample_files$path[sample_files$file_type == "features"][1]
  matrix_path  <- sample_files$path[sample_files$file_type == "matrix"][1]
  
  sample_dir <- file.path(tmp_10x_dir, sid)
  dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
  
  file.copy(barcode_path, file.path(sample_dir, "barcodes.tsv.gz"), overwrite = TRUE)
  file.copy(feature_path, file.path(sample_dir, "features.tsv.gz"), overwrite = TRUE)
  file.copy(matrix_path,  file.path(sample_dir, "matrix.mtx.gz"), overwrite = TRUE)
  
  counts <- Read10X(data.dir = sample_dir)
  
  obj_i <- CreateSeuratObject(
    counts = counts,
    project = sid,
    min.cells = 3,
    min.features = 200
  )
  
  obj_i$sample_id <- sid
  obj_i$patient_id <- pid
  obj_i$compartment <- comp
  obj_i$dataset <- "GSE242477"
  
  obj_i <- RenameCells(
    obj_i,
    add.cell.id = sid
  )
  
  seurat_list[[sid]] <- obj_i
}

# ------------------------------------------------------------
# 6. Combinar objetos
# ------------------------------------------------------------

if (length(seurat_list) == 1) {
  obj <- seurat_list[[1]]
} else {
  obj <- merge(
    seurat_list[[1]],
    y = seurat_list[-1],
    project = "GSE242477"
  )
}

# ------------------------------------------------------------
# 7. Añadir métricas QC básicas
# ------------------------------------------------------------

obj[["percent.mt"]] <- PercentageFeatureSet(
  obj,
  pattern = "^MT-"
)

# ------------------------------------------------------------
# 8. Guardar objeto y metadata
# ------------------------------------------------------------

saveRDS(
  obj,
  file.path(out_dir, "gse242477_raw.rds")
)

write.csv(
  obj@meta.data,
  file.path(res_dir, "raw_metadata.csv")
)

metadata_summary <- data.frame(
  column = colnames(obj@meta.data),
  class = sapply(obj@meta.data, function(x) class(x)[1]),
  n_unique = sapply(obj@meta.data, function(x) length(unique(x))),
  example_values = sapply(
    obj@meta.data,
    function(x) paste(head(unique(x), 8), collapse = "; ")
  )
)

write.csv(
  metadata_summary,
  file.path(res_dir, "metadata_summary.csv"),
  row.names = FALSE
)

sample_summary <- obj@meta.data %>%
  group_by(patient_id, compartment, sample_id) %>%
  summarise(
    n_cells = n(),
    median_nFeature_RNA = median(nFeature_RNA),
    median_nCount_RNA = median(nCount_RNA),
    median_percent_mt = median(percent.mt),
    .groups = "drop"
  )

write.csv(
  sample_summary,
  file.path(res_dir, "sample_summary.csv"),
  row.names = FALSE
)

print(sample_summary)

message("Import OK")
message("Saved object: ", file.path(out_dir, "gse242477_raw.rds"))