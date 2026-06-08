#30_gse242477_summary_for_results

base_dir <- "/Users/jimmlucas/Downloads/R_MASTER_RNASEQ"
setwd(base_dir)

source(file.path(base_dir, "00_config.R"))

library(tidyverse)
library(data.table)
library(dplyr)
library(tidyr)

res_dir <- file.path(base_dir, "results/GSE242477")

patient_scores <- read.csv(
  file.path(res_dir, "gse242477_patient_level_cd8_module_scores.csv")
)

score_cols <- grep(
  "_score1$",
  colnames(patient_scores),
  value = TRUE
)

# ------------------------------------------------------------
# 1. Resumen por grupo clínico
# ------------------------------------------------------------

summary_by_group <- patient_scores %>%
  group_by(response_group) %>%
  summarise(
    n_patients = n(),
    patients = paste(patient_id, collapse = "; "),
    median_cd8_like_cells = median(n_cd8_like_cells, na.rm = TRUE),
    across(
      all_of(score_cols),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        median = ~ median(.x, na.rm = TRUE),
        sd = ~ sd(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

write.csv(
  summary_by_group,
  file.path(res_dir, "gse242477_summary_by_response_group.csv"),
  row.names = FALSE
)

print(summary_by_group)

# ------------------------------------------------------------
# 2. Tabla direccional Mixed responder vs Non-responder
# No se calculan p-values por tamaño muestral insuficiente
# ------------------------------------------------------------

long_scores <- patient_scores %>%
  pivot_longer(
    cols = all_of(score_cols),
    names_to = "program",
    values_to = "score"
  ) %>%
  mutate(
    program = gsub("_score1", "", program)
  )

direction_table <- long_scores %>%
  group_by(program, response_group) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    median_score = median(score, na.rm = TRUE),
    n_patients = n(),
    patients = paste(patient_id, collapse = "; "),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = response_group,
    values_from = c(mean_score, median_score, n_patients, patients),
    names_sep = "__"
  )

direction_table <- direction_table %>%
  mutate(
    mean_difference_mixed_minus_NR =
      `mean_score__Mixed responder` - `mean_score__Non-responder`,
    median_difference_mixed_minus_NR =
      `median_score__Mixed responder` - `median_score__Non-responder`,
    direction = case_when(
      mean_difference_mixed_minus_NR > 0 ~ "Higher in mixed responder",
      mean_difference_mixed_minus_NR < 0 ~ "Higher in non-responders",
      TRUE ~ "No difference"
    )
  )

write.csv(
  direction_table,
  file.path(res_dir, "gse242477_direction_of_effect.csv"),
  row.names = FALSE
)

print(direction_table)

# ------------------------------------------------------------
# 3. Tabla compacta para la memoria
# ------------------------------------------------------------

compact_table <- direction_table %>%
  select(
    program,
    mean_score_NR = `mean_score__Non-responder`,
    mean_score_mixed = `mean_score__Mixed responder`,
    mean_difference_mixed_minus_NR,
    direction
  ) %>%
  arrange(program)

write.csv(
  compact_table,
  file.path(res_dir, "gse242477_compact_results_for_tfm.csv"),
  row.names = FALSE
)

print(compact_table)

message("GSE242477 summary OK")