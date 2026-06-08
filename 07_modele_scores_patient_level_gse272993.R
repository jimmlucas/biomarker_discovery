# =========================
# 07_module_scores_patient_level_gse272993.R
# =========================
# Agregación de module scores funcionales por paciente/timepoint.
# Estos scores se usan como biomarcadores candidatos interpretativos.

source("scripts/00_config.R")

library(tidyverse)
library(data.table)

meta <- readRDS(file.path(proc_dir, "GSE272993_metadata_with_response_binary.rds"))

score_cols <- grep("Jo_top250|module", colnames(meta), value = TRUE)

module_patient <- meta %>%
  filter(!is.na(response_binary)) %>%
  group_by(patient_alias, response, response_binary, treatment, timepoint) %>%
  summarise(
    across(
      all_of(score_cols),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        median = ~ median(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  mutate(
    across(
      where(is.numeric),
      ~ ifelse(is.finite(.x), .x, NA_real_)
    )
  )

fwrite(
  module_patient,
  file.path(res_dir, "GSE272993_patient_timepoint_module_scores.tsv"),
  sep = "\t"
)

comp_wide <- fread(
  file.path(res_dir, "GSE272993_patient_timepoint_metacluster_composition_wide.tsv")
)

model_df_modules <- comp_wide %>%
  left_join(
    module_patient,
    by = c("patient_alias", "response", "response_binary", "treatment", "timepoint")
  ) %>%
  filter(!is.na(response_binary)) %>%
  mutate(y = ifelse(response_binary == "Responder", 1, 0))

fwrite(
  model_df_modules,
  file.path(res_dir, "GSE272993_model_matrix_composition_plus_modules.tsv"),
  sep = "\t"
)

module_long <- module_patient %>%
  pivot_longer(
    cols = ends_with("_mean"),
    names_to = "score",
    values_to = "value"
  )

module_na_summary <- module_long %>%
  group_by(score) %>%
  summarise(
    n_total = n(),
    n_missing = sum(is.na(value) | !is.finite(value)),
    pct_missing = 100 * n_missing / n_total,
    .groups = "drop"
  ) %>%
  arrange(desc(n_missing))

fwrite(
  module_na_summary,
  file.path(res_dir, "GSE272993_module_score_missing_summary.tsv"),
  sep = "\t"
)

module_long_clean <- module_long %>%
  filter(!is.na(value), is.finite(value))

wilcox_modules <- module_long_clean %>%
  group_by(score) %>%
  summarise(
    p_value = wilcox.test(value ~ response_binary)$p.value,
    median_responder = median(value[response_binary == "Responder"], na.rm = TRUE),
    median_non_responder = median(value[response_binary == "Non_responder"], na.rm = TRUE),
    delta_median = median_responder - median_non_responder,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    direction = case_when(
      delta_median > 0 ~ "Higher_in_responders",
      delta_median < 0 ~ "Lower_in_responders",
      TRUE ~ "No_difference"
    ),
    evidence = case_when(
      p_adj < 0.05 ~ "FDR_significant",
      p_value < 0.05 & p_adj >= 0.05 ~ "Nominal_trend",
      TRUE ~ "No_evidence"
    )
  ) %>%
  arrange(p_adj)

fwrite(
  wilcox_modules,
  file.path(res_dir, "GSE272993_wilcox_module_scores_response_binary.tsv"),
  sep = "\t"
)

top_scores <- wilcox_modules %>%
  arrange(p_adj) %>%
  slice_head(n = 6) %>%
  pull(score)

p_modules <- module_long_clean %>%
  filter(score %in% top_scores) %>%
  ggplot(aes(x = response_binary, y = value, fill = response_binary)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
  facet_wrap(~ score, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Respuesta clínica",
    y = "Module score medio por paciente/timepoint",
    title = "Module scores CD8+ candidatos asociados a respuesta"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(fig_dir, "GSE272993_top_module_scores_by_response.png"),
  p_modules,
  width = 12,
  height = 8,
  dpi = 300
)

fwrite(
  module_long_clean,
  file.path(res_dir, "GSE272993_module_scores_long_clean.tsv"),
  sep = "\t"
)

print(wilcox_modules)