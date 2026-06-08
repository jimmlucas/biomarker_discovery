# =========================
# 13_import_explore_gse120575.R
# =========================
# Exploración inicial de GSE120575 para validación externa OE4.
# Objetivo:
# 1) localizar archivos descargados
# 2) identificar si hay objetos RDS/RData, matrices o metadata
# 3) guardar resumen para decidir el siguiente paso

base_dir <- "/Users/jimmlucas/Downloads/R_MASTER_RNASEQ"
setwd(base_dir)

source(file.path(base_dir, "00_config.R"))

library(tidyverse)
library(data.table)

gse_id <- "GSE120575"

gse_dir <- file.path(raw_dir, gse_id)

if (!dir.exists(gse_dir)) {
  stop("No existe la carpeta: ", gse_dir)
}

# -------------------------
# 1. Listado de archivos
# -------------------------

files_gse120575 <- tibble(
  file_path = list.files(gse_dir, recursive = TRUE, full.names = TRUE),
  file_name = basename(file_path),
  extension = tools::file_ext(file_path),
  size_mb = file.info(file_path)$size / 1024^2
) %>%
  arrange(desc(size_mb))

print(files_gse120575, n = Inf)

fwrite(
  files_gse120575,
  file.path(res_dir, "GSE120575_available_files.tsv"),
  sep = "\t"
)

# -------------------------
# 2. Identificar posibles archivos útiles
# -------------------------

candidate_files <- files_gse120575 %>%
  filter(
    str_detect(
      file_name,
      regex("rds|RDS|rda|RData|csv|tsv|txt|mtx|metadata|meta|annotation|cell", ignore_case = TRUE)
    )
  )

print(candidate_files, n = Inf)

fwrite(
  candidate_files,
  file.path(res_dir, "GSE120575_candidate_files.tsv"),
  sep = "\t"
)

# -------------------------
# 3. Intentar leer archivos RDS
# -------------------------

rds_files <- files_gse120575 %>%
  filter(str_detect(file_name, regex("\\.rds$|\\.RDS$|\\.rds.gz$|\\.RDS.gz$", ignore_case = TRUE))) %>%
  pull(file_path)

rds_summary <- list()

if (length(rds_files) > 0) {
  
  for (i in seq_along(rds_files)) {
    
    f <- rds_files[i]
    message("Intentando leer RDS: ", f)
    
    obj <- tryCatch(
      readRDS(f),
      error = function(e) {
        message("No se pudo leer: ", f)
        message(e$message)
        return(NULL)
      }
    )
    
    if (!is.null(obj)) {
      
      obj_class <- paste(class(obj), collapse = ";")
      
      obj_dim <- tryCatch(
        paste(dim(obj), collapse = " x "),
        error = function(e) NA_character_
      )
      
      meta_cols <- NA_character_
      n_meta_rows <- NA_integer_
      
      if ("Seurat" %in% class(obj)) {
        meta_cols <- paste(colnames(obj@meta.data), collapse = ";")
        n_meta_rows <- nrow(obj@meta.data)
      }
      
      rds_summary[[i]] <- tibble(
        file_path = f,
        file_name = basename(f),
        object_class = obj_class,
        object_dim = obj_dim,
        n_meta_rows = n_meta_rows,
        meta_cols = meta_cols
      )
      
      rm(obj)
      gc()
    }
  }
}

rds_summary_df <- bind_rows(rds_summary)

if (nrow(rds_summary_df) > 0) {
  
  fwrite(
    rds_summary_df,
    file.path(res_dir, "GSE120575_rds_object_summary.tsv"),
    sep = "\t"
  )
  
  print(rds_summary_df)
  
} else {
  
  message("No se encontraron RDS legibles o no hay archivos RDS.")
}

# -------------------------
# 4. Intentar explorar archivos tabulares pequeños
# -------------------------

tabular_files <- files_gse120575 %>%
  filter(
    str_detect(file_name, regex("\\.csv$|\\.csv.gz$|\\.tsv$|\\.tsv.gz$|\\.txt$|\\.txt.gz$", ignore_case = TRUE)),
    size_mb < 200
  ) %>%
  pull(file_path)

tabular_summary <- list()

if (length(tabular_files) > 0) {
  
  for (i in seq_along(tabular_files)) {
    
    f <- tabular_files[i]
    message("Explorando archivo tabular: ", f)
    
    dt <- tryCatch(
      fread(f, nrows = 20),
      error = function(e) {
        message("No se pudo leer: ", f)
        message(e$message)
        return(NULL)
      }
    )
    
    if (!is.null(dt)) {
      
      tabular_summary[[i]] <- tibble(
        file_path = f,
        file_name = basename(f),
        n_preview_rows = nrow(dt),
        n_cols = ncol(dt),
        columns = paste(colnames(dt), collapse = ";")
      )
    }
  }
}

tabular_summary_df <- bind_rows(tabular_summary)

if (nrow(tabular_summary_df) > 0) {
  
  fwrite(
    tabular_summary_df,
    file.path(res_dir, "GSE120575_tabular_file_summary.tsv"),
    sep = "\t"
  )
  
  print(tabular_summary_df, n = Inf)
  
} else {
  
  message("No se encontraron archivos tabulares pequeños explorables.")
}

message("Exploración inicial de GSE120575 finalizada.")