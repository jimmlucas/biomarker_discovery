# =========================
# 31_import_explore_gse273718.R
# =========================
# Exploración inicial de GSE273718 para preparar la proyección/anotación.
# Objetivo:
# 1) listar archivos descargados
# 2) identificar objetos RDS/RData o matrices
# 3) explorar metadata/anotaciones disponibles
# 4) guardar resumen para decidir si se puede hacer Symphony o label transfer

base_dir <- "/Users/jimmlucas/Downloads/R_MASTER_RNASEQ"
setwd(base_dir)

source(file.path(base_dir, "00_config.R"))

library(tidyverse)
library(data.table)

gse_id <- "GSE273718"
gse_dir <- file.path(raw_dir, gse_id)

if (!dir.exists(gse_dir)) {
  stop("No existe la carpeta: ", gse_dir)
}

# -------------------------
# 1. Listado de archivos
# -------------------------

files_gse273718 <- tibble(
  file_path = list.files(gse_dir, recursive = TRUE, full.names = TRUE),
  file_name = basename(file_path),
  extension = tools::file_ext(file_path),
  size_mb = file.info(file_path)$size / 1024^2
) %>%
  arrange(desc(size_mb))

fwrite(
  files_gse273718,
  file.path(res_dir, "GSE273718_available_files.tsv"),
  sep = "\t"
)

print(files_gse273718, n = Inf)

# -------------------------
# 2. Archivos candidatos
# -------------------------

candidate_files <- files_gse273718 %>%
  filter(
    str_detect(
      file_name,
      regex(
        "rds|RDS|rda|RData|h5|h5ad|mtx|csv|tsv|txt|metadata|meta|annotation|cell|seurat",
        ignore_case = TRUE
      )
    )
  )

fwrite(
  candidate_files,
  file.path(res_dir, "GSE273718_candidate_files.tsv"),
  sep = "\t"
)

print(candidate_files, n = Inf)

# -------------------------
# 3. Intentar leer RDS/RDS.gz
# -------------------------

rds_files <- files_gse273718 %>%
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
      
      assays <- NA_character_
      reductions <- NA_character_
      meta_cols <- NA_character_
      n_meta_rows <- NA_integer_
      
      if ("Seurat" %in% class(obj)) {
        assays <- paste(Seurat::Assays(obj), collapse = ";")
        reductions <- paste(Seurat::Reductions(obj), collapse = ";")
        meta_cols <- paste(colnames(obj@meta.data), collapse = ";")
        n_meta_rows <- nrow(obj@meta.data)
      }
      
      rds_summary[[i]] <- tibble(
        file_path = f,
        file_name = basename(f),
        object_class = obj_class,
        object_dim = obj_dim,
        assays = assays,
        reductions = reductions,
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
    file.path(res_dir, "GSE273718_rds_object_summary.tsv"),
    sep = "\t"
  )
  
  print(rds_summary_df, width = Inf)
  
} else {
  message("No se encontraron RDS legibles o no hay archivos RDS.")
}

# -------------------------
# 4. Explorar archivos tabulares pequeños
# -------------------------

tabular_files <- files_gse273718 %>%
  filter(
    str_detect(file_name, regex("\\.csv$|\\.csv.gz$|\\.tsv$|\\.tsv.gz$|\\.txt$|\\.txt.gz$", ignore_case = TRUE)),
    size_mb < 300
  ) %>%
  pull(file_path)

tabular_summary <- list()

if (length(tabular_files) > 0) {
  
  for (i in seq_along(tabular_files)) {
    
    f <- tabular_files[i]
    message("Explorando archivo tabular: ", f)
    
    dt <- tryCatch(
      fread(f, nrows = 20, data.table = FALSE),
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
      
      preview_file <- file.path(
        res_dir,
        paste0("GSE273718_preview_", make.names(basename(f)), ".tsv")
      )
      
      fwrite(dt, preview_file, sep = "\t")
    }
  }
}

tabular_summary_df <- bind_rows(tabular_summary)

if (nrow(tabular_summary_df) > 0) {
  
  fwrite(
    tabular_summary_df,
    file.path(res_dir, "GSE273718_tabular_file_summary.tsv"),
    sep = "\t"
  )
  
  print(tabular_summary_df, n = Inf)
  
} else {
  message("No se encontraron archivos tabulares pequeños explorables.")
}

# -------------------------
# 5. Resumen final
# -------------------------

exploration_summary <- tibble(
  dataset = "GSE273718",
  n_files = nrow(files_gse273718),
  n_candidate_files = nrow(candidate_files),
  n_rds_files = length(rds_files),
  n_readable_rds = nrow(rds_summary_df),
  n_tabular_files_explored = nrow(tabular_summary_df)
)

fwrite(
  exploration_summary,
  file.path(res_dir, "GSE273718_exploration_summary.tsv"),
  sep = "\t"
)

print(exploration_summary)

message("Exploración inicial de GSE273718 finalizada.")