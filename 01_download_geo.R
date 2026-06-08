# -------------------------
# 1. Cargar objeto GSE272993
# -------------------------

obj_file_raw_gz <- file.path(
  raw_dir,
  "GSE272993",
  "GSE272993_cd8_nn_labeled_FINAL.RDS.gz"
)

obj_file_raw <- file.path(
  raw_dir,
  "GSE272993",
  "GSE272993_cd8_nn_labeled_FINAL.RDS"
)

if (!file.exists(obj_file_raw)) {
  
  if (!file.exists(obj_file_raw_gz)) {
    stop("No existe ni el archivo RDS ni el RDS.gz de GSE272993 en data/raw.")
  }
  
  if (!requireNamespace("R.utils", quietly = TRUE)) {
    install.packages("R.utils", type = "binary")
  }
  
  R.utils::gunzip(
    filename = obj_file_raw_gz,
    destname = obj_file_raw,
    overwrite = TRUE,
    remove = FALSE
  )
}

obj <- readRDS(obj_file_raw)

message("Objeto GSE272993 cargado correctamente desde data/raw.")
print(class(obj))
print(dim(obj))