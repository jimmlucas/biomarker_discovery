# =========================
# 19_score_cd8_programs_gse120575.R
# =========================

base_dir <- "/Users/jimmlucas/Downloads/R_MASTER_RNASEQ"
setwd(base_dir)

source(file.path(base_dir, "00_config.R"))

library(tidyverse)
library(data.table)

# -------------------------
# 1. Cargar expresión de marcadores y metadata clínica estricta
# -------------------------

marker_expr <- readRDS(
  file.path(proc_dir, "GSE120575_marker_expression_long.rds")
)

cell_meta <- readRDS(
  file.path(proc_dir, "GSE120575_cell_metadata_with_response_strict.rds")
)

# -------------------------
# 2. Pasar a formato ancho célula x gen
# -------------------------

marker_wide <- marker_expr %>%
  pivot_wider(
    names_from = gene,
    values_from = expression,
    values_fill = 0
  )

# -------------------------
# 3. Calcular scores funcionales simples
# -------------------------

score_df <- marker_wide %>%
  mutate(
    score_exhaustion = rowMeans(select(., any_of(c("PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX"))), na.rm = TRUE),
    score_effector   = rowMeans(select(., any_of(c("GZMB", "PRF1", "IFNG"))), na.rm = TRUE),
    score_memory     = rowMeans(select(., any_of(c("IL7R", "CCR7"))), na.rm = TRUE),
    score_activation = rowMeans(select(., any_of(c("CD69", "HLA-DRA"))), na.rm = TRUE),
    score_cycling    = rowMeans(select(., any_of(c("MKI67"))), na.rm = TRUE)
  )

# -------------------------
# 4. Unir con metadata
# -------------------------

score_meta <- score_df %>%
  left_join(
    cell_meta %>%
      select(cell_id, patient_alias, timepoint, response_binary_strict, therapy),
    by = "cell_id"
  )

# -------------------------
# 5. Agregar a paciente/timepoint
# -------------------------

patient_scores <- score_meta %>%
  filter(
    !is.na(response_binary_strict),
    response_binary_strict != "Mixed"
  ) %>%
  group_by(patient_alias, timepoint, response_binary_strict, therapy) %>%
  summarise(
    n_cells = n(),
    across(
      starts_with("score_"),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        median = ~ median(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

# -------------------------
# 6. Guardar resultados
# -------------------------

fwrite(
  score_meta,
  file.path(res_dir, "GSE120575_cell_level_cd8_program_scores.tsv"),
  sep = "\t"
)

fwrite(
  patient_scores,
  file.path(res_dir, "GSE120575_patient_timepoint_cd8_program_scores.tsv"),
  sep = "\t"
)

saveRDS(
  patient_scores,
  file.path(proc_dir, "GSE120575_patient_timepoint_cd8_program_scores.rds")
)

# -------------------------
# 7. Tests exploratorios
# -------------------------

score_long <- patient_scores %>%
  pivot_longer(
    cols = starts_with("score_"),
    names_to = "score",
    values_to = "value"
  )

wilcox_scores <- score_long %>%
  group_by(score) %>%
  summarise(
    p_value = wilcox.test(value ~ response_binary_strict)$p.value,
    median_responder = median(value[response_binary_strict == "Responder"], na.rm = TRUE),
    median_non_responder = median(value[response_binary_strict == "Non_responder"], na.rm = TRUE),
    delta_median = median_responder - median_non_responder,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    direction = case_when(
      delta_median > 0 ~ "Higher_in_responders",
      delta_median < 0 ~ "Lower_in_responders",
      TRUE ~ "No_difference"
    )
  ) %>%
  arrange(p_adj)

fwrite(
  wilcox_scores,
  file.path(res_dir, "GSE120575_wilcox_cd8_program_scores_response.tsv"),
  sep = "\t"
)

print(wilcox_scores)

# -------------------------
# 8. Figura exploratoria
# -------------------------

p_scores <- score_long %>%
  ggplot(aes(x = response_binary_strict, y = value, fill = response_binary_strict)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.6) +
  facet_wrap(~ score, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Respuesta clínica",
    y = "Score medio por paciente/timepoint",
    title = "Programas funcionales CD8+ en GSE120575"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(fig_dir, "GSE120575_cd8_program_scores_by_response.png"),
  p_scores,
  width = 12,
  height = 8,
  dpi = 300
)

message("Scores funcionales CD8+ de GSE120575 generados.")