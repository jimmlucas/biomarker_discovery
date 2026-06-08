
base_dir <- "/Users/jimmlucas/Downloads/R_MASTER_RNASEQ"
setwd(base_dir)

source(file.path(base_dir, "00_config.R"))

library(tidyverse)
library(data.table)
library(Seurat)
library(stringr)
library(tibble)

in_file <- file.path(
  base_dir,
  "data/processed/GSE242477/gse242477_raw.rds"
)

out_file <- file.path(
  base_dir,
  "data/processed/GSE242477/gse242477_clean.rds"
)

res_dir <- file.path(base_dir, "results/GSE242477")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(in_file)

meta <- obj@meta.data %>%
  rownames_to_column("cell_id")

clinical_map <- data.frame(
  patient_id = c("MEL03", "MEL09", "MEL10", "MEL11", "MEL12"),
  response_detail = c(
    "Non-responder_PD",
    "Mixed_responder_primary_CR_metastases_PR_PD",
    "Surgery_only_tumor_free",
    "Non-responder_PD",
    "Non-responder_PD"
  ),
  response_group = c(
    "Non-responder",
    "Mixed responder",
    "Excluded",
    "Non-responder",
    "Non-responder"
  ),
  stringsAsFactors = FALSE
)

meta <- meta %>%
  mutate(
    patient_id = toupper(as.character(patient_id)),
    sample_id = as.character(sample_id),
    compartment = as.character(compartment)
  ) %>%
  left_join(clinical_map, by = "patient_id") %>%
  column_to_rownames("cell_id")

# Mantener exactamente el orden de células del objeto
meta <- meta[colnames(obj), , drop = FALSE]

stopifnot(identical(rownames(meta), colnames(obj)))

obj@meta.data <- meta

saveRDS(obj, out_file)

write.csv(
  obj@meta.data,
  file.path(res_dir, "metadata_clean.csv")
)

patient_response_summary <- obj@meta.data %>%
  rownames_to_column("cell_id") %>%
  distinct(
    patient_id,
    compartment,
    sample_id,
    response_group,
    response_detail
  ) %>%
  arrange(patient_id, compartment)

write.csv(
  patient_response_summary,
  file.path(res_dir, "patient_response_summary.csv"),
  row.names = FALSE
)

print(patient_response_summary)

message("Metadata OK")
message("Saved object: ", out_file)