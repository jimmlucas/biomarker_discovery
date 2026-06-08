# =========================
# 15_clean_metadata_gse120575.R
# =========================

base_dir <- "/Users/jimmlucas/Downloads/R_MASTER_RNASEQ"
setwd(base_dir)

source(file.path(base_dir, "00_config.R"))

library(tidyverse)
library(data.table)

# -------------------------
# 1. Cargar metadata derivada del header de expresión
# -------------------------

cell_meta <- fread(
  file.path(res_dir, "GSE120575_cell_metadata_from_expression_header.tsv"),
  data.table = FALSE
) %>%
  as_tibble()

cell_meta <- cell_meta %>%
  mutate(
    patient_alias = patient_id,
    timepoint = case_when(
      timepoint == "Pre" ~ "Baseline",
      timepoint == "Post" ~ "Post",
      TRUE ~ timepoint
    )
  )

# -------------------------
# 2. Resumen por paciente/timepoint
# -------------------------

sample_summary <- cell_meta %>%
  count(patient_alias, timepoint, name = "n_cells") %>%
  arrange(patient_alias, timepoint)

fwrite(
  sample_summary,
  file.path(res_dir, "GSE120575_patient_timepoint_cell_counts.tsv"),
  sep = "\t"
)

print(sample_summary, n = Inf)

# -------------------------
# 3. Resumen global
# -------------------------

cohort_summary <- cell_meta %>%
  summarise(
    dataset = "GSE120575",
    n_cells = n(),
    n_patients = n_distinct(patient_alias),
    n_baseline_cells = sum(timepoint == "Baseline", na.rm = TRUE),
    n_post_cells = sum(timepoint == "Post", na.rm = TRUE)
  )

fwrite(
  cohort_summary,
  file.path(res_dir, "GSE120575_clean_metadata_summary.tsv"),
  sep = "\t"
)

print(cohort_summary)

# -------------------------
# 4. Guardar metadata limpia
# -------------------------

fwrite(
  cell_meta,
  file.path(res_dir, "GSE120575_cell_metadata_clean.tsv"),
  sep = "\t"
)

saveRDS(
  cell_meta,
  file.path(proc_dir, "GSE120575_cell_metadata_clean.rds")
)

message("Metadata limpia de GSE120575 guardada.")