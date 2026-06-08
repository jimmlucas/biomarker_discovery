# =========================
# 24_therapy_stratified_gse120575.R
# =========================
# Análisis estratificado por terapia en GSE120575.
# Objetivo:
# comprobar si la señal de respuesta en scores CD8+ depende del tipo de tratamiento.

base_dir <- "/Users/jimmlucas/Downloads/R_MASTER_RNASEQ"
setwd(base_dir)

source(file.path(base_dir, "00_config.R"))

library(tidyverse)
library(data.table)
library(pROC)

set.seed(123)

# -------------------------
# 1. Cargar scores de GSE120575
# -------------------------

gse120575 <- fread(
  file.path(res_dir, "GSE120575_patient_timepoint_cd8_program_scores.tsv"),
  data.table = FALSE
) %>%
  as_tibble() %>%
  filter(
    !is.na(response_binary_strict),
    response_binary_strict != "Mixed"
  ) %>%
  mutate(
    response_binary = response_binary_strict,
    y = ifelse(response_binary == "Responder", 1, 0),
    therapy_clean = case_when(
      therapy %in% c("anti-PD1", "aPD1", "aPD-1") ~ "anti-PD1",
      therapy %in% c("anti-CTLA4", "aCTLA-4") ~ "anti-CTLA4",
      therapy %in% c("anti-CTLA4+PD1", "Combination") ~ "Combination",
      TRUE ~ therapy
    )
  )

feature_cols <- c(
  "score_exhaustion_mean",
  "score_effector_mean",
  "score_memory_mean",
  "score_activation_mean",
  "score_cycling_mean"
)

# -------------------------
# 2. Resumen por terapia y respuesta
# -------------------------

therapy_response_summary <- gse120575 %>%
  count(therapy_clean, response_binary, name = "n_patient_timepoints") %>%
  arrange(therapy_clean, response_binary)

fwrite(
  therapy_response_summary,
  file.path(res_dir, "GSE120575_therapy_response_summary.tsv"),
  sep = "\t"
)

print(therapy_response_summary)

# -------------------------
# 3. Tests Wilcoxon por terapia
# -------------------------

score_long <- gse120575 %>%
  pivot_longer(
    cols = all_of(feature_cols),
    names_to = "score",
    values_to = "value"
  )

wilcox_by_therapy <- score_long %>%
  group_by(therapy_clean, score) %>%
  filter(n_distinct(response_binary) == 2) %>%
  summarise(
    n = n(),
    n_responders = sum(response_binary == "Responder"),
    n_non_responders = sum(response_binary == "Non_responder"),
    p_value = wilcox.test(value ~ response_binary)$p.value,
    median_responder = median(value[response_binary == "Responder"], na.rm = TRUE),
    median_non_responder = median(value[response_binary == "Non_responder"], na.rm = TRUE),
    delta_median = median_responder - median_non_responder,
    .groups = "drop"
  ) %>%
  group_by(therapy_clean) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    direction = case_when(
      delta_median > 0 ~ "Higher_in_responders",
      delta_median < 0 ~ "Lower_in_responders",
      TRUE ~ "No_difference"
    )
  ) %>%
  ungroup() %>%
  arrange(therapy_clean, p_adj)

fwrite(
  wilcox_by_therapy,
  file.path(res_dir, "GSE120575_wilcox_cd8_scores_by_therapy.tsv"),
  sep = "\t"
)

print(wilcox_by_therapy, n = Inf)

# -------------------------
# 4. Figuras de scores por terapia y respuesta
# -------------------------

p_therapy_scores <- score_long %>%
  ggplot(aes(x = response_binary, y = value, fill = response_binary)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 0.9, alpha = 0.6) +
  facet_grid(therapy_clean ~ score, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Respuesta clínica",
    y = "Score medio por paciente/timepoint",
    title = "Programas CD8+ en GSE120575 estratificados por terapia"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(fig_dir, "GSE120575_cd8_scores_by_response_and_therapy.png"),
  p_therapy_scores,
  width = 16,
  height = 9,
  dpi = 300
)

# -------------------------
# 5. AUROC univariante por terapia y score
# -------------------------

auc_by_therapy <- score_long %>%
  group_by(therapy_clean, score) %>%
  filter(n_distinct(response_binary) == 2) %>%
  summarise(
    n = n(),
    n_responders = sum(response_binary == "Responder"),
    n_non_responders = sum(response_binary == "Non_responder"),
    auc_raw = as.numeric(auc(roc(y, value, quiet = TRUE))),
    auc_direction_corrected = max(auc_raw, 1 - auc_raw),
    direction = ifelse(auc_raw >= 0.5, "Higher_score_predicts_response", "Lower_score_predicts_response"),
    .groups = "drop"
  ) %>%
  arrange(therapy_clean, desc(auc_direction_corrected))

fwrite(
  auc_by_therapy,
  file.path(res_dir, "GSE120575_auc_cd8_scores_by_therapy.tsv"),
  sep = "\t"
)

print(auc_by_therapy, n = Inf)

# -------------------------
# 6. Modelo logístico simple solo anti-PD1
# -------------------------

anti_pd1_df <- gse120575 %>%
  filter(therapy_clean == "anti-PD1")

if (nrow(anti_pd1_df) >= 8 && n_distinct(anti_pd1_df$response_binary) == 2) {
  
  anti_pd1_model <- glm(
    y ~ score_memory_mean +
      score_exhaustion_mean +
      score_effector_mean +
      score_activation_mean +
      score_cycling_mean,
    data = anti_pd1_df,
    family = binomial()
  )
  
  anti_pd1_pred <- as.numeric(predict(anti_pd1_model, type = "response"))
  
  anti_pd1_roc <- roc(anti_pd1_df$y, anti_pd1_pred, quiet = TRUE)
  anti_pd1_auc <- as.numeric(auc(anti_pd1_roc))
  
  anti_pd1_coef <- broom::tidy(anti_pd1_model) %>%
    mutate(model = "logistic_all_scores_anti_PD1")
  
  anti_pd1_summary <- tibble(
    model = "logistic_all_scores_anti_PD1",
    dataset = "GSE120575",
    therapy = "anti-PD1",
    n_rows = nrow(anti_pd1_df),
    n_patients = n_distinct(anti_pd1_df$patient_alias),
    n_responders = sum(anti_pd1_df$y == 1),
    n_non_responders = sum(anti_pd1_df$y == 0),
    auc = anti_pd1_auc
  )
  
  fwrite(
    anti_pd1_coef,
    file.path(res_dir, "GSE120575_antiPD1_logistic_coefficients.tsv"),
    sep = "\t"
  )
  
  fwrite(
    anti_pd1_summary,
    file.path(res_dir, "GSE120575_antiPD1_logistic_summary.tsv"),
    sep = "\t"
  )
  
  print(anti_pd1_summary)
  print(anti_pd1_coef)
  
  png(
    file.path(fig_dir, "GSE120575_antiPD1_logistic_roc.png"),
    width = 1600,
    height = 1200,
    res = 200
  )
  
  plot(
    anti_pd1_roc,
    col = "#1B9E77",
    lwd = 3,
    main = paste0("GSE120575 anti-PD1 - AUROC = ", round(anti_pd1_auc, 3))
  )
  
  abline(a = 0, b = 1, lty = 2, col = "gray50")
  
  dev.off()
  
} else {
  
  message("No hay suficientes muestras anti-PD1 con ambas clases para modelo logístico.")
}

# -------------------------
# 7. Comparar predicciones del modelo externo por terapia
# -------------------------

external_predictions_file <- file.path(
  res_dir,
  "GSE120575_marker_scores_external_predictions.tsv"
)

if (file.exists(external_predictions_file)) {
  
  external_predictions <- fread(
    external_predictions_file,
    data.table = FALSE
  ) %>%
    as_tibble() %>%
    mutate(
      therapy_clean = case_when(
        treatment %in% c("anti-PD1", "aPD1", "aPD-1") ~ "anti-PD1",
        treatment %in% c("anti-CTLA4", "aCTLA-4") ~ "anti-CTLA4",
        treatment %in% c("anti-CTLA4+PD1", "Combination") ~ "Combination",
        TRUE ~ treatment
      )
    )
  
  pred_by_therapy <- external_predictions %>%
    group_by(therapy_clean, response_binary) %>%
    summarise(
      n = n(),
      median_pred = median(pred_prob, na.rm = TRUE),
      mean_pred = mean(pred_prob, na.rm = TRUE),
      .groups = "drop"
    )
  
  fwrite(
    pred_by_therapy,
    file.path(res_dir, "GSE120575_external_predictions_by_therapy.tsv"),
    sep = "\t"
  )
  
  print(pred_by_therapy, n = Inf)
  
  p_pred_therapy <- external_predictions %>%
    ggplot(aes(x = response_binary, y = pred_prob, fill = response_binary)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 1.1, alpha = 0.6) +
    facet_wrap(~ therapy_clean, scales = "free_y") +
    theme_classic() +
    labs(
      x = "Respuesta clínica",
      y = "Probabilidad predicha de respuesta",
      title = "Predicciones externas estratificadas por terapia"
    ) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  ggsave(
    file.path(fig_dir, "GSE120575_external_predictions_by_therapy.png"),
    p_pred_therapy,
    width = 10,
    height = 5,
    dpi = 300
  )
}

message("Análisis estratificado por terapia en GSE120575 completado.")