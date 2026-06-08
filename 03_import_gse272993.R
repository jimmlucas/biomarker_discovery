# =========================
# 03_import_gse272993.R
# =========================
# Importación inicial de la cohorte GSE272993.
# Por limitaciones de memoria, se extrae y guarda principalmente la metadata,
# evitando guardar de nuevo el objeto Seurat completo.

source("scripts/00_config.R")

library(Seurat)
library(tidyverse)
library(data.table)

if (!requireNamespace("R.utils", quietly = TRUE)) {
  install.packages("R.utils", type = "binary")
}
library(R.utils)

gse272993_rds_gz <- file.path(
  raw_dir,
  "GSE272993",
  "GSE272993_cd8_nn_labeled_FINAL.RDS.gz"
)

gse272993_rds <- sub("\\.gz$", "", gse272993_rds_gz)

stopifnot(file.exists(gse272993_rds_gz))

# Descomprimir manteniendo el archivo .gz original.
gunzip(
  filename = gse272993_rds_gz,
  destname = gse272993_rds,
  overwrite = TRUE,
  remove = FALSE
)

# Comprobación básica del archivo descomprimido.
print(system2("file", gse272993_rds, stdout = TRUE))
print(file.info(gse272993_rds)$size)

# Cargar objeto Seurat.
obj_272993 <- readRDS(gse272993_rds)

# Inspección básica.
print(class(obj_272993))
print(dim(obj_272993))
print(colnames(obj_272993@meta.data))

# Extraer metadata celular.
meta_272993 <- obj_272993@meta.data %>%
  tibble::rownames_to_column("cell_id")

saveRDS(
  meta_272993,
  file.path(proc_dir, "GSE272993_metadata_only.rds")
)

fwrite(
  meta_272993,
  file.path(proc_dir, "GSE272993_metadata_only.tsv"),
  sep = "\t"
)

# Resumen del objeto Seurat.
object_summary <- tibble(
  dataset = "GSE272993",
  object_class = paste(class(obj_272993), collapse = ";"),
  n_features = nrow(obj_272993),
  n_cells = ncol(obj_272993),
  default_assay = DefaultAssay(obj_272993),
  assays = paste(Assays(obj_272993), collapse = ";"),
  reductions = paste(Reductions(obj_272993), collapse = ";")
)

fwrite(
  object_summary,
  file.path(proc_dir, "GSE272993_object_summary.tsv"),
  sep = "\t"
)

# Liberar memoria.
rm(obj_272993)
gc()