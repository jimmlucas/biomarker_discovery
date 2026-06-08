# =========================
# 32_extract_inspect_gse273718_tar.R
# =========================
# Extrae e inspecciona el .tar de GSE273718.

base_dir <- "/Users/jimmlucas/Downloads/R_MASTER_RNASEQ"
setwd(base_dir)

source(file.path(base_dir, "00_config.R"))

library(tidyverse)
library(data.table)

gse_id <- "GSE273718"
gse_dir <- file.path(raw_dir, gse_id)

tar_files <- list.files(
  gse_dir,
  pattern = "\\.tar$|\\.tar\\.gz$|\\.tgz$",
  full.names = TRUE
)

if (length(tar_files) == 0) {
  stop("No se encontró ningún archivo tar en: ", gse_dir)
}

tar_file <- tar_files[1]

message("Archivo TAR encontrado: ", tar_file)

# -------------------------
# 1. Ver contenido sin extraer
# -------------------------

tar_contents <- untar(tar_file, list = TRUE)

tar_contents_df <- tibble(
  file_inside_tar = tar_contents,
  file_name = basename(tar_contents),
  extension = tools::file_ext(tar_contents)
)

fwrite(
  tar_contents_df,
  file.path(res_dir, "GSE273718_tar_contents.tsv"),
  sep = "\t"
)

print(tar_contents_df, n = Inf)

# -------------------------
# 2. Extraer en subcarpeta
# -------------------------

extract_dir <- file.path(gse_dir, "extracted")

if (!dir.exists(extract_dir)) {
  dir.create(extract_dir, recursive = TRUE)
}

untar(
  tarfile = tar_file,
  exdir = extract_dir
)

message("Archivos extraídos en: ", extract_dir)

# -------------------------
# 3. Listar archivos extraídos
# -------------------------

extracted_files <- tibble(
  file_path = list.files(extract_dir, recursive = TRUE, full.names = TRUE),
  file_name = basename(file_path),
  extension = tools::file_ext(file_path),
  size_mb = file.info(file_path)$size / 1024^2
) %>%
  arrange(desc(size_mb))

fwrite(
  extracted_files,
  file.path(res_dir, "GSE273718_extracted_files.tsv"),
  sep = "\t"
)

print(extracted_files, n = Inf)

# -------------------------
# 4. Detectar candidatos útiles
# -------------------------

candidate_files <- extracted_files %>%
  filter(
    str_detect(
      file_name,
      regex(
        "rds|RDS|rda|RData|h5|h5ad|mtx|csv|tsv|txt|metadata|meta|annotation|cell|barcodes|features|genes|matrix",
        ignore_case = TRUE
      )
    )
  )

fwrite(
  candidate_files,
  file.path(res_dir, "GSE273718_extracted_candidate_files.tsv"),
  sep = "\t"
)

print(candidate_files, n = Inf)

# -------------------------
# 5. Explorar tabulares pequeños
# -------------------------

tabular_files <- extracted_files %>%
  filter(
    str_detect(
      file_name,
      regex("\\.csv$|\\.csv.gz$|\\.tsv$|\\.tsv.gz$|\\.txt$|\\.txt.gz$", ignore_case = TRUE)
    ),
    size_mb < 300
  ) %>%
  pull(file_path)

tabular_summary <- list()

if (length(tabular_files) > 0) {
  
  for (i in seq_along(tabular_files)) {
    
    f <- tabular_files[i]
    message("Explorando tabular: ", f)
    
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
        paste0("GSE273718_extracted_preview_", make.names(basename(f)), ".tsv")
      )
      
      fwrite(dt, preview_file, sep = "\t")
    }
  }
}

tabular_summary_df <- bind_rows(tabular_summary)

if (nrow(tabular_summary_df) > 0) {
  fwrite(
    tabular_summary_df,
    file.path(res_dir, "GSE273718_extracted_tabular_summary.tsv"),
    sep = "\t"
  )
  
  print(tabular_summary_df, n = Inf)
} else {
  message("No hay tabulares pequeños explorables.")
}

# -------------------------
# 6. Resumen final
# -------------------------

extract_summary <- tibble(
  dataset = "GSE273718",
  tar_file = basename(tar_file),
  n_files_inside_tar = length(tar_contents),
  n_extracted_files = nrow(extracted_files),
  n_candidate_files = nrow(candidate_files),
  n_tabular_files_explored = nrow(tabular_summary_df)
)

fwrite(
  extract_summary,
  file.path(res_dir, "GSE273718_extract_summary.tsv"),
  sep = "\t"
)

print(extract_summary)

message("Extracción e inspección de GSE273718 finalizada.")