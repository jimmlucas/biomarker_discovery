# =========================
# 17_resolve_response_gse120575.R
# =========================

base_dir <- "/Users/jimmlucas/Downloads/R_MASTER_RNASEQ"
setwd(base_dir)

source(file.path(base_dir, "00_config.R"))

library(tidyverse)
library(data.table)

cell_meta <- fread(
  file.path(res_dir, "GSE120575_cell_metadata_with_response.tsv"),
  data.table = FALSE
) %>%
  as_tibble()

# -------------------------
# 1. Detectar conflictos por paciente/timepoint
# -------------------------

response_resolution <- cell_meta %>%
  filter(!is.na(response_binary)) %>%
  count(patient_alias, timepoint, response_binary, name = "n_cells_response") %>%
  group_by(patient_alias, timepoint) %>%
  mutate(
    n_total = sum(n_cells_response),
    frac = n_cells_response / n_total,
    n_response_classes = n_distinct(response_binary)
  ) %>%
  ungroup()

fwrite(
  response_resolution,
  file.path(res_dir, "GSE120575_response_resolution_by_patient_timepoint.tsv"),
  sep = "\t"
)

print(response_resolution, n = Inf)

# -------------------------
# 2. Crear etiquetas estrictas
#    - si solo hay una clase: se conserva
#    - si hay más de una clase: se marca como Mixed y se excluye del modelado
# -------------------------

strict_labels <- response_resolution %>%
  group_by(patient_alias, timepoint) %>%
  summarise(
    n_response_classes = n_distinct(response_binary),
    dominant_response = response_binary[which.max(n_cells_response)],
    dominant_fraction = max(frac),
    n_cells = sum(n_cells_response),
    response_binary_strict = ifelse(
      n_response_classes == 1,
      dominant_response,
      "Mixed"
    ),
    .groups = "drop"
  )

fwrite(
  strict_labels,
  file.path(res_dir, "GSE120575_strict_patient_timepoint_labels.tsv"),
  sep = "\t"
)

print(strict_labels, n = Inf)

# -------------------------
# 3. Añadir etiqueta estricta a cada célula
# -------------------------

cell_meta_strict <- cell_meta %>%
  select(-any_of("response_binary_strict")) %>%
  left_join(
    strict_labels %>%
      select(patient_alias, timepoint, response_binary_strict),
    by = c("patient_alias", "timepoint")
  )

fwrite(
  cell_meta_strict,
  file.path(res_dir, "GSE120575_cell_metadata_with_response_strict.tsv"),
  sep = "\t"
)

saveRDS(
  cell_meta_strict,
  file.path(proc_dir, "GSE120575_cell_metadata_with_response_strict.rds")
)

# -------------------------
# 4. Resumen final para modelado
# -------------------------

strict_summary <- strict_labels %>%
  count(response_binary_strict, name = "n_patient_timepoints")

strict_patient_summary <- strict_labels %>%
  filter(response_binary_strict != "Mixed") %>%
  distinct(patient_alias, response_binary_strict) %>%
  count(response_binary_strict, name = "n_patients")

fwrite(
  strict_summary,
  file.path(res_dir, "GSE120575_strict_response_summary.tsv"),
  sep = "\t"
)

fwrite(
  strict_patient_summary,
  file.path(res_dir, "GSE120575_strict_patient_response_summary.tsv"),
  sep = "\t"
)

print(strict_summary)
print(strict_patient_summary)

message("Etiquetas clínicas estrictas de GSE120575 generadas.")