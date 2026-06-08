# =========================
# 14_inspect_gse120575_structure.R
# =========================

base_dir <- "/Users/jimmlucas/Downloads/R_MASTER_RNASEQ"
setwd(base_dir)

source(file.path(base_dir, "00_config.R"))

library(tidyverse)
library(data.table)

gse_dir <- file.path(raw_dir, "GSE120575")

meta_file <- file.path(gse_dir, "GSE120575_patient_ID_single_cells.txt")
expr_file <- file.path(gse_dir, "GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz")

stopifnot(file.exists(meta_file))
stopifnot(file.exists(expr_file))

# -------------------------
# 1. Metadata externa
# -------------------------

meta_raw <- fread(
  meta_file,
  header = FALSE,
  fill = TRUE,
  data.table = FALSE
)

cat("Dim metadata:\n")
print(dim(meta_raw))

fwrite(
  meta_raw,
  file.path(res_dir, "GSE120575_patient_ID_single_cells_raw.tsv"),
  sep = "\t"
)

# -------------------------
# 2. Leer solo líneas iniciales del archivo TPM
# -------------------------

con <- gzfile(expr_file, open = "rt")
first_lines <- readLines(con, n = 5)
close(con)

# Línea 1: cell IDs
# Línea 2: paciente/timepoint
# Línea 3+: genes

cell_ids <- strsplit(first_lines[1], "\t")[[1]][-1]
sample_labels <- strsplit(first_lines[2], "\t")[[1]][-1]

expr_header_df <- tibble(
  cell_id = cell_ids,
  sample_label = sample_labels
) %>%
  mutate(
    timepoint = case_when(
      str_detect(sample_label, "^Pre") ~ "Pre",
      str_detect(sample_label, "^Post") ~ "Post",
      TRUE ~ NA_character_
    ),
    patient_id = str_extract(sample_label, "P[0-9]+")
  )

cat("Número de células en matriz TPM:\n")
print(nrow(expr_header_df))

cat("Resumen por timepoint:\n")
print(expr_header_df %>% count(timepoint))

cat("Resumen por paciente:\n")
print(expr_header_df %>% count(patient_id, timepoint) %>% arrange(patient_id, timepoint))

fwrite(
  expr_header_df,
  file.path(res_dir, "GSE120575_cell_metadata_from_expression_header.tsv"),
  sep = "\t"
)

# -------------------------
# 3. Genes preliminares sin cargar toda la matriz
# -------------------------

gene_preview <- first_lines[-c(1, 2)] %>%
  map_chr(~ strsplit(.x, "\t")[[1]][1])

gene_preview_df <- tibble(
  row_index = seq_along(gene_preview),
  gene = gene_preview
)

fwrite(
  gene_preview_df,
  file.path(res_dir, "GSE120575_gene_preview.tsv"),
  sep = "\t"
)

# -------------------------
# 4. Resumen final
# -------------------------

summary_check <- tibble(
  dataset = "GSE120575",
  n_cells_expression = nrow(expr_header_df),
  n_metadata_rows = nrow(meta_raw),
  n_patients_detected = n_distinct(expr_header_df$patient_id),
  n_pre_cells = sum(expr_header_df$timepoint == "Pre", na.rm = TRUE),
  n_post_cells = sum(expr_header_df$timepoint == "Post", na.rm = TRUE)
)

fwrite(
  summary_check,
  file.path(res_dir, "GSE120575_expression_metadata_dimension_check.tsv"),
  sep = "\t"
)

print(summary_check)

message("Inspección estructural de GSE120575 finalizada.")