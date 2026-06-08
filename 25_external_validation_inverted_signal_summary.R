# =========================
# 25_external_validation_inverted_signal_summary.R
# =========================
# Consolidación final de OE4:
# validación externa, domain shift e inversión de señal.

base_dir <- "/Users/jimmlucas/Downloads/R_MASTER_RNASEQ"
setwd(base_dir)

source(file.path(base_dir, "00_config.R"))

library(tidyverse)
library(data.table)

# -------------------------
# 1. Cargar resultados principales
# -------------------------

external_perf <- fread(
  file.path(res_dir, "external_validation_marker_scores_performance_summary.tsv"),
  data.table = FALSE
) %>%
  as_tibble()

domain_norm_perf <- fread(
  file.path(res_dir, "domain_shift_corrected_external_validation_performance.tsv"),
  data.table = FALSE
) %>%
  as_tibble()

dataset_auc_norm <- fread(
  file.path(res_dir, "domain_shift_dataset_auc_by_normalization.tsv"),
  data.table = FALSE
) %>%
  as_tibble()

prediction_summary <- fread(
  file.path(res_dir, "domain_shift_corrected_prediction_summary.tsv"),
  data.table = FALSE
) %>%
  as_tibble()

therapy_pred <- fread(
  file.path(res_dir, "GSE120575_external_predictions_by_therapy.tsv"),
  data.table = FALSE
) %>%
  as_tibble()

wilcox_therapy <- fread(
  file.path(res_dir, "GSE120575_wilcox_cd8_scores_by_therapy.tsv"),
  data.table = FALSE
) %>%
  as_tibble()

auc_therapy <- fread(
  file.path(res_dir, "GSE120575_auc_cd8_scores_by_therapy.tsv"),
  data.table = FALSE
) %>%
  as_tibble()

# -------------------------
# 2. Tabla resumen de validación externa
# -------------------------

external_validation_summary <- external_perf %>%
  transmute(
    analysis = "Train GSE272993 -> External GSE120575",
    model,
    train_auc,
    train_auprc,
    external_auc,
    external_auprc,
    train_n_patients,
    external_n_patients,
    external_n_responders,
    external_n_non_responders,
    interpretation = "High external AUROC but inverted predicted probabilities"
  )

fwrite(
  external_validation_summary,
  file.path(res_dir, "OE4_external_validation_summary.tsv"),
  sep = "\t"
)

print(external_validation_summary, width = Inf)

# -------------------------
# 3. Tabla de domain shift y normalización
# -------------------------

domain_shift_summary <- dataset_auc_norm %>%
  mutate(
    interpretation = case_when(
      normalization == "raw" & dataset_auc > 0.9 ~
        "Strong dataset/domain shift",
      normalization != "raw" & dataset_auc < 0.6 ~
        "Dataset discrimination largely removed",
      TRUE ~
        "Intermediate domain effect"
    )
  )

fwrite(
  domain_shift_summary,
  file.path(res_dir, "OE4_domain_shift_summary.tsv"),
  sep = "\t"
)

print(domain_shift_summary)

# -------------------------
# 4. Tabla de inversión de señal
# -------------------------

signal_inversion_summary <- prediction_summary %>%
  select(
    normalization,
    response_binary_common,
    n,
    median_pred,
    mean_pred
  ) %>%
  pivot_wider(
    names_from = response_binary_common,
    values_from = c(n, median_pred, mean_pred)
  ) %>%
  mutate(
    median_difference_responder_minus_non_responder =
      median_pred_Responder - median_pred_Non_responder,
    mean_difference_responder_minus_non_responder =
      mean_pred_Responder - mean_pred_Non_responder,
    signal_direction = case_when(
      median_difference_responder_minus_non_responder > 0 ~
        "Expected_direction",
      median_difference_responder_minus_non_responder < 0 ~
        "Inverted_direction",
      TRUE ~
        "No_difference"
    )
  )

fwrite(
  signal_inversion_summary,
  file.path(res_dir, "OE4_signal_inversion_summary.tsv"),
  sep = "\t"
)

print(signal_inversion_summary, width = Inf)

# -------------------------
# 5. Tabla resumen por terapia
# -------------------------

therapy_signal_summary <- therapy_pred %>%
  pivot_wider(
    names_from = response_binary,
    values_from = c(n, median_pred, mean_pred)
  ) %>%
  mutate(
    median_difference_responder_minus_non_responder =
      median_pred_Responder - median_pred_Non_responder,
    signal_direction = case_when(
      median_difference_responder_minus_non_responder > 0 ~
        "Expected_direction",
      median_difference_responder_minus_non_responder < 0 ~
        "Inverted_direction",
      TRUE ~
        "No_difference"
    )
  )

fwrite(
  therapy_signal_summary,
  file.path(res_dir, "OE4_signal_inversion_by_therapy.tsv"),
  sep = "\t"
)

print(therapy_signal_summary, width = Inf)

# -------------------------
# 6. Top señales univariantes en GSE120575 anti-PD1
# -------------------------

top_antipd1_scores <- auc_therapy %>%
  filter(therapy_clean == "anti-PD1") %>%
  arrange(desc(auc_direction_corrected)) %>%
  select(
    therapy_clean,
    score,
    n,
    n_responders,
    n_non_responders,
    auc_raw,
    auc_direction_corrected,
    direction
  )

fwrite(
  top_antipd1_scores,
  file.path(res_dir, "OE4_top_GSE120575_antiPD1_univariate_scores.tsv"),
  sep = "\t"
)

print(top_antipd1_scores, n = Inf)

# -------------------------
# 7. Figura resumen de métricas
# -------------------------

plot_perf <- domain_norm_perf %>%
  select(normalization, train_auc, external_auc, train_auprc, external_auprc) %>%
  pivot_longer(
    cols = c(train_auc, external_auc, train_auprc, external_auprc),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      train_auc = "Train AUROC",
      external_auc = "External AUROC",
      train_auprc = "Train AUPRC",
      external_auprc = "External AUPRC"
    )
  )

p_perf <- plot_perf %>%
  ggplot(aes(x = normalization, y = value, fill = metric)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  theme_classic() +
  ylim(0, 1) +
  labs(
    x = "Normalización",
    y = "Métrica",
    title = "Rendimiento del modelo CD8+ en validación externa"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(fig_dir, "OE4_external_validation_performance_summary.png"),
  p_perf,
  width = 9,
  height = 5,
  dpi = 300
)

# -------------------------
# 8. Figura resumen de inversión de señal
# -------------------------

plot_inversion <- prediction_summary %>%
  mutate(
    response_binary_common = factor(
      response_binary_common,
      levels = c("Non_responder", "Responder")
    )
  )

p_inversion <- plot_inversion %>%
  ggplot(aes(x = response_binary_common, y = median_pred, fill = response_binary_common)) +
  geom_col(width = 0.65) +
  facet_wrap(~ normalization) +
  theme_classic() +
  labs(
    x = "Respuesta clínica externa",
    y = "Mediana de probabilidad predicha",
    title = "Inversión de la señal predictiva en GSE120575"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(fig_dir, "OE4_signal_inversion_by_normalization.png"),
  p_inversion,
  width = 9,
  height = 5,
  dpi = 300
)

# -------------------------
# 9. Figura de inversión por terapia
# -------------------------

plot_therapy <- therapy_pred %>%
  mutate(
    response_binary = factor(
      response_binary,
      levels = c("Non_responder", "Responder")
    )
  )

p_therapy <- plot_therapy %>%
  ggplot(aes(x = response_binary, y = median_pred, fill = response_binary)) +
  geom_col(width = 0.65) +
  facet_wrap(~ therapy_clean) +
  theme_classic() +
  labs(
    x = "Respuesta clínica",
    y = "Mediana de probabilidad predicha",
    title = "Inversión de señal predictiva por tipo de terapia"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(fig_dir, "OE4_signal_inversion_by_therapy.png"),
  p_therapy,
  width = 9,
  height = 5,
  dpi = 300
)

# -------------------------
# 10. Tabla final compacta para memoria
# -------------------------

OE4_final_table <- tibble(
  item = c(
    "External validation AUROC",
    "External validation AUPRC",
    "Raw dataset discrimination AUROC",
    "Dataset AUROC after z-score",
    "Dataset AUROC after rank normalization",
    "External signal direction",
    "Therapy-stratified interpretation"
  ),
  value = c(
    round(external_perf$external_auc[1], 3),
    round(external_perf$external_auprc[1], 3),
    round(dataset_auc_norm$dataset_auc[dataset_auc_norm$normalization == "raw"], 3),
    round(dataset_auc_norm$dataset_auc[dataset_auc_norm$normalization == "zscore_by_dataset"], 3),
    round(dataset_auc_norm$dataset_auc[dataset_auc_norm$normalization == "rank_by_dataset"], 3),
    "Inverted",
    "Inversion persists within anti-PD1"
  )
)

fwrite(
  OE4_final_table,
  file.path(res_dir, "OE4_final_compact_results_table.tsv"),
  sep = "\t"
)

print(OE4_final_table)

message("Consolidación final de OE4 completada.")