# =========================
# 22_domain_shift_analysis.R
# =========================

base_dir <- "/Users/jimmlucas/Downloads/R_MASTER_RNASEQ"
setwd(base_dir)

source(file.path(base_dir, "00_config.R"))

library(tidyverse)
library(data.table)
library(pROC)

set.seed(123)

# -------------------------
# 1. Cargar scores de ambas cohortes
# -------------------------

gse272993 <- fread(
  file.path(res_dir, "GSE272993_patient_timepoint_marker_program_scores.tsv"),
  data.table = FALSE
) %>%
  as_tibble() %>%
  mutate(
    dataset                = "GSE272993",
    response_binary_common = response_binary,
    therapy_common         = treatment
  )

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
    dataset                = "GSE120575",
    response_binary_common = response_binary_strict,
    therapy_common         = therapy
  )

feature_cols <- c(
  "score_exhaustion_mean",
  "score_effector_mean",
  "score_memory_mean",
  "score_activation_mean",
  "score_cycling_mean"
)

# -------------------------
# 2. Unir datasets
# -------------------------

combined_scores <- bind_rows(
  gse272993 %>%
    select(dataset, patient_alias, timepoint, response_binary_common,
           therapy_common, n_cells, all_of(feature_cols)),
  gse120575 %>%
    select(dataset, patient_alias, timepoint, response_binary_common,
           therapy_common, n_cells, all_of(feature_cols))
) %>%
  filter(!is.na(response_binary_common))

fwrite(
  combined_scores,
  file.path(res_dir, "combined_GSE272993_GSE120575_marker_scores.tsv"),
  sep = "\t"
)

# -------------------------
# 3. Comparación de distribución por dataset
# -------------------------

score_long <- combined_scores %>%
  pivot_longer(
    cols      = all_of(feature_cols),
    names_to  = "score",
    values_to = "value"
  )

dataset_tests <- score_long %>%
  group_by(score) %>%
  summarise(
    p_value_wilcox_dataset = wilcox.test(value ~ dataset)$p.value,
    median_GSE272993       = median(value[dataset == "GSE272993"], na.rm = TRUE),
    median_GSE120575       = median(value[dataset == "GSE120575"], na.rm = TRUE),
    delta_median           = median_GSE120575 - median_GSE272993,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj_dataset = p.adjust(p_value_wilcox_dataset, method = "BH"),
    direction = case_when(
      delta_median > 0 ~ "Higher_in_GSE120575",
      delta_median < 0 ~ "Higher_in_GSE272993",
      TRUE             ~ "No_difference"
    )
  ) %>%
  arrange(p_adj_dataset)

fwrite(
  dataset_tests,
  file.path(res_dir, "domain_shift_marker_scores_dataset_tests.tsv"),
  sep = "\t"
)

print(dataset_tests)

# -------------------------
# 4. Figura: scores por dataset
# -------------------------

p_dataset_scores <- score_long %>%
  ggplot(aes(x = dataset, y = value, fill = dataset)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
  facet_wrap(~ score, scales = "free_y") +
  theme_classic() +
  labs(
    x     = "Dataset",
    y     = "Score medio por paciente/timepoint",
    title = "Diferencias de distribución de programas CD8+ entre cohortes"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(fig_dir, "domain_shift_marker_scores_by_dataset.png"),
  p_dataset_scores,
  width = 12, height = 8, dpi = 300
)

# -------------------------
# 5. PCA conjunto
# -------------------------

pca_df <- combined_scores %>%
  select(all_of(feature_cols)) %>%
  mutate(across(everything(), ~ ifelse(is.na(.x) | !is.finite(.x), median(.x, na.rm = TRUE), .x)))

pca_scaled <- scale(pca_df)

pca <- prcomp(pca_scaled, center = FALSE, scale. = FALSE)

pca_scores <- as_tibble(pca$x[, 1:3]) %>%
  bind_cols(
    combined_scores %>%
      select(dataset, patient_alias, timepoint, response_binary_common, therapy_common, n_cells)
  )

pca_var <- tibble(
  PC                 = paste0("PC", seq_along(pca$sdev)),
  variance_explained = (pca$sdev^2) / sum(pca$sdev^2)
)

fwrite(pca_scores, file.path(res_dir, "domain_shift_pca_scores.tsv"),             sep = "\t")
fwrite(pca_var,    file.path(res_dir, "domain_shift_pca_variance_explained.tsv"), sep = "\t")

p_pca_dataset <- ggplot(
  pca_scores,
  aes(x = PC1, y = PC2, color = dataset, shape = response_binary_common)
) +
  geom_point(size = 3, alpha = 0.85) +
  theme_classic() +
  labs(
    title    = "PCA de programas funcionales CD8+",
    subtitle = "Color por dataset; forma por respuesta clínica",
    x        = paste0("PC1 (", round(100 * pca_var$variance_explained[1], 1), "%)"),
    y        = paste0("PC2 (", round(100 * pca_var$variance_explained[2], 1), "%)")
  )

ggsave(
  file.path(fig_dir, "domain_shift_pca_dataset_response.png"),
  p_pca_dataset,
  width = 8, height = 6, dpi = 300
)

# -------------------------
# 6. PCA coloreado por terapia
# -------------------------

p_pca_therapy <- ggplot(
  pca_scores,
  aes(x = PC1, y = PC2, color = therapy_common, shape = dataset)
) +
  geom_point(size = 3, alpha = 0.85) +
  theme_classic() +
  labs(
    title = "PCA de programas funcionales CD8+ por terapia",
    x     = paste0("PC1 (", round(100 * pca_var$variance_explained[1], 1), "%)"),
    y     = paste0("PC2 (", round(100 * pca_var$variance_explained[2], 1), "%)")
  )

ggsave(
  file.path(fig_dir, "domain_shift_pca_therapy.png"),
  p_pca_therapy,
  width = 8, height = 6, dpi = 300
)

# -------------------------
# 7. Asociación de PCs con dataset y respuesta
# -------------------------

pc_tests <- pca_scores %>%
  pivot_longer(
    cols      = c(PC1, PC2, PC3),
    names_to  = "PC",
    values_to = "value"
  )

pc_dataset_tests <- pc_tests %>%
  group_by(PC) %>%
  summarise(
    p_dataset        = wilcox.test(value ~ dataset)$p.value,
    median_GSE272993 = median(value[dataset == "GSE272993"], na.rm = TRUE),
    median_GSE120575 = median(value[dataset == "GSE120575"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(p_adj_dataset = p.adjust(p_dataset, method = "BH"))

pc_response_tests <- pc_tests %>%
  group_by(PC) %>%
  summarise(
    p_response           = wilcox.test(value ~ response_binary_common)$p.value,
    median_responder     = median(value[response_binary_common == "Responder"],     na.rm = TRUE),
    median_non_responder = median(value[response_binary_common == "Non_responder"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(p_adj_response = p.adjust(p_response, method = "BH"))

fwrite(pc_dataset_tests,  file.path(res_dir, "domain_shift_pca_dataset_tests.tsv"),  sep = "\t")
fwrite(pc_response_tests, file.path(res_dir, "domain_shift_pca_response_tests.tsv"), sep = "\t")

print(pc_dataset_tests)
print(pc_response_tests)

# -------------------------
# 8. Modelo simple: predecir dataset a partir de scores
# -------------------------

combined_model <- combined_scores %>%
  mutate(dataset_binary = ifelse(dataset == "GSE120575", 1, 0))

dataset_glm <- glm(
  dataset_binary ~ score_exhaustion_mean +
    score_effector_mean +
    score_memory_mean +
    score_activation_mean +
    score_cycling_mean,
  data   = combined_model,
  family = binomial()
)

dataset_pred <- as.numeric(predict(dataset_glm, type = "response"))
y_dataset    <- combined_model$dataset_binary

roc_dataset <- roc(y_dataset, dataset_pred, quiet = TRUE)
auc_dataset <- as.numeric(auc(roc_dataset))

domain_classifier_summary <- tibble(
  model                      = "Logistic_regression_predict_dataset",
  outcome                    = "GSE120575_vs_GSE272993",
  auc_dataset_discrimination = auc_dataset,
  n_rows                     = nrow(combined_model),
  n_GSE272993                = sum(combined_model$dataset == "GSE272993"),
  n_GSE120575                = sum(combined_model$dataset == "GSE120575")
)

fwrite(
  domain_classifier_summary,
  file.path(res_dir, "domain_shift_dataset_classifier_summary.tsv"),
  sep = "\t"
)

print(domain_classifier_summary)

# -------------------------
# 9. Curva ROC clasificador de dominio (ggplot2, ejes 0-1 estrictos)
# -------------------------

roc_dataset_df <- data.frame(
  specificity = roc_dataset$specificities,
  sensitivity = roc_dataset$sensitivities
)

p_roc_dataset <- ggplot(roc_dataset_df, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(color = "#7B3294", linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0), name = "1 - Especificidad (FPR)") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), name = "Sensibilidad (TPR)") +
  theme_classic() +
  labs(title = paste0("Discriminación de dataset por scores CD8+ — AUROC = ", round(auc_dataset, 3)))

ggsave(
  file.path(fig_dir, "domain_shift_dataset_classifier_roc.png"),
  p_roc_dataset,
  width = 6, height = 6, dpi = 300
)

message("Análisis de domain shift completado.")