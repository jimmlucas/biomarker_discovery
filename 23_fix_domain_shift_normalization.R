# =========================
# 23_fix_domain_shift_normalization.R
# =========================
# Corrección exploratoria del domain shift entre GSE272993 y GSE120575.
# Se comparan dos estrategias:
# 1) Z-score dentro de cada dataset
# 2) Rank/percentile normalization dentro de cada dataset

base_dir <- "/Users/jimmlucas/Downloads/R_MASTER_RNASEQ"
setwd(base_dir)

source(file.path(base_dir, "00_config.R"))

library(tidyverse)
library(data.table)
library(glmnet)
library(pROC)
library(PRROC)

set.seed(123)

# -------------------------
# 1. Cargar datos
# -------------------------

train_df <- fread(
  file.path(res_dir, "GSE272993_patient_timepoint_marker_program_scores.tsv"),
  data.table = FALSE
) %>%
  as_tibble() %>%
  mutate(
    dataset = "GSE272993",
    response_binary_common = response_binary,
    therapy_common = treatment
  )

external_df <- fread(
  file.path(res_dir, "GSE120575_patient_timepoint_cd8_program_scores.tsv"),
  data.table = FALSE
) %>%
  as_tibble() %>%
  filter(
    !is.na(response_binary_strict),
    response_binary_strict != "Mixed"
  ) %>%
  mutate(
    dataset = "GSE120575",
    response_binary_common = response_binary_strict,
    therapy_common = therapy
  )

feature_cols <- c(
  "score_exhaustion_mean",
  "score_effector_mean",
  "score_memory_mean",
  "score_activation_mean",
  "score_cycling_mean"
)

# -------------------------
# 2. Armonizar columnas
# -------------------------

train_base <- train_df %>%
  filter(!is.na(response_binary_common)) %>%
  mutate(y = ifelse(response_binary_common == "Responder", 1, 0)) %>%
  select(
    dataset,
    patient_alias,
    timepoint,
    response_binary_common,
    therapy_common,
    n_cells,
    y,
    all_of(feature_cols)
  )

external_base <- external_df %>%
  filter(!is.na(response_binary_common)) %>%
  mutate(y = ifelse(response_binary_common == "Responder", 1, 0)) %>%
  select(
    dataset,
    patient_alias,
    timepoint,
    response_binary_common,
    therapy_common,
    n_cells,
    y,
    all_of(feature_cols)
  )

combined_base <- bind_rows(train_base, external_base)

# -------------------------
# 3. Funciones de normalización
# -------------------------

zscore_by_dataset <- function(df, feature_cols) {
  
  df %>%
    group_by(dataset) %>%
    mutate(
      across(
        all_of(feature_cols),
        ~ {
          s <- sd(.x, na.rm = TRUE)
          m <- mean(.x, na.rm = TRUE)
          
          if (is.na(s) || s == 0) {
            rep(0, length(.x))
          } else {
            (.x - m) / s
          }
        }
      )
    ) %>%
    ungroup()
}

rank_by_dataset <- function(df, feature_cols) {
  
  df %>%
    group_by(dataset) %>%
    mutate(
      across(
        all_of(feature_cols),
        ~ percent_rank(.x)
      )
    ) %>%
    ungroup()
}

# -------------------------
# 4. Crear datasets transformados
# -------------------------

combined_raw <- combined_base %>%
  mutate(normalization = "raw")

combined_zscore <- zscore_by_dataset(combined_base, feature_cols) %>%
  mutate(normalization = "zscore_by_dataset")

combined_rank <- rank_by_dataset(combined_base, feature_cols) %>%
  mutate(normalization = "rank_by_dataset")

combined_all <- bind_rows(
  combined_raw,
  combined_zscore,
  combined_rank
)

fwrite(
  combined_all,
  file.path(res_dir, "domain_shift_normalized_combined_scores.tsv"),
  sep = "\t"
)

# -------------------------
# 5. Función de entrenamiento y validación externa
# -------------------------

fit_external_validation <- function(df, norm_name) {
  
  train_model_df <- df %>%
    filter(dataset == "GSE272993")
  
  external_model_df <- df %>%
    filter(dataset == "GSE120575")
  
  for (f in feature_cols) {
    
    med_train <- median(train_model_df[[f]], na.rm = TRUE)
    
    train_model_df[[f]][is.na(train_model_df[[f]]) | !is.finite(train_model_df[[f]])] <- med_train
    external_model_df[[f]][is.na(external_model_df[[f]]) | !is.finite(external_model_df[[f]])] <- med_train
  }
  
  X_train <- as.matrix(train_model_df %>% select(all_of(feature_cols)))
  y_train <- train_model_df$y
  
  X_external <- as.matrix(external_model_df %>% select(all_of(feature_cols)))
  y_external <- external_model_df$y
  
  cvfit <- cv.glmnet(
    x = X_train,
    y = y_train,
    family = "binomial",
    alpha = 0.5,
    type.measure = "auc"
  )
  
  pred_train <- as.numeric(
    predict(cvfit, newx = X_train, s = "lambda.min", type = "response")
  )
  
  pred_external <- as.numeric(
    predict(cvfit, newx = X_external, s = "lambda.min", type = "response")
  )
  
  roc_train <- roc(y_train, pred_train, quiet = TRUE)
  roc_external <- roc(y_external, pred_external, quiet = TRUE)
  
  auc_train <- as.numeric(auc(roc_train))
  auc_external <- as.numeric(auc(roc_external))
  
  pr_train <- pr.curve(
    scores.class0 = pred_train[y_train == 1],
    scores.class1 = pred_train[y_train == 0],
    curve = TRUE
  )
  
  pr_external <- pr.curve(
    scores.class0 = pred_external[y_external == 1],
    scores.class1 = pred_external[y_external == 0],
    curve = TRUE
  )
  
  coef_df <- as.matrix(coef(cvfit, s = "lambda.min")) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("feature") %>%
    rename(coef = 2) %>%
    mutate(normalization = norm_name) %>%
    arrange(desc(abs(coef)))
  
  train_predictions <- train_model_df %>%
    mutate(
      pred_prob = pred_train,
      normalization = norm_name,
      validation_type = "training"
    )
  
  external_predictions <- external_model_df %>%
    mutate(
      pred_prob = pred_external,
      normalization = norm_name,
      validation_type = "external"
    )
  
  prediction_summary <- external_predictions %>%
    group_by(response_binary_common) %>%
    summarise(
      n = n(),
      median_pred = median(pred_prob, na.rm = TRUE),
      mean_pred = mean(pred_prob, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(normalization = norm_name)
  
  performance <- tibble(
    normalization = norm_name,
    model = "ElasticNet_marker_program_scores",
    train_auc = auc_train,
    train_auprc = pr_train$auc.integral,
    external_auc = auc_external,
    external_auprc = pr_external$auc.integral,
    train_n_rows = nrow(train_model_df),
    external_n_rows = nrow(external_model_df),
    external_n_responders = sum(y_external == 1),
    external_n_non_responders = sum(y_external == 0),
    lambda_min = cvfit$lambda.min
  )
  
  list(
    performance = performance,
    coef = coef_df,
    train_predictions = train_predictions,
    external_predictions = external_predictions,
    prediction_summary = prediction_summary,
    roc_train = roc_train,
    roc_external = roc_external
  )
}

# -------------------------
# 6. Ejecutar validaciones
# -------------------------

results_raw <- fit_external_validation(
  combined_raw,
  "raw"
)

results_zscore <- fit_external_validation(
  combined_zscore,
  "zscore_by_dataset"
)

results_rank <- fit_external_validation(
  combined_rank,
  "rank_by_dataset"
)

performance_all <- bind_rows(
  results_raw$performance,
  results_zscore$performance,
  results_rank$performance
)

coef_all <- bind_rows(
  results_raw$coef,
  results_zscore$coef,
  results_rank$coef
)

predictions_all <- bind_rows(
  results_raw$train_predictions,
  results_raw$external_predictions,
  results_zscore$train_predictions,
  results_zscore$external_predictions,
  results_rank$train_predictions,
  results_rank$external_predictions
)

prediction_summary_all <- bind_rows(
  results_raw$prediction_summary,
  results_zscore$prediction_summary,
  results_rank$prediction_summary
)

# -------------------------
# 7. Guardar resultados
# -------------------------

fwrite(
  performance_all,
  file.path(res_dir, "domain_shift_corrected_external_validation_performance.tsv"),
  sep = "\t"
)

fwrite(
  coef_all,
  file.path(res_dir, "domain_shift_corrected_external_validation_coefficients.tsv"),
  sep = "\t"
)

fwrite(
  predictions_all,
  file.path(res_dir, "domain_shift_corrected_external_validation_predictions.tsv"),
  sep = "\t"
)

fwrite(
  prediction_summary_all,
  file.path(res_dir, "domain_shift_corrected_prediction_summary.tsv"),
  sep = "\t"
)

print(performance_all, width = Inf)
print(prediction_summary_all, width = Inf)
print(coef_all)

# -------------------------
# 8. Comparar distribución tras normalización
# -------------------------

score_long_all <- combined_all %>%
  pivot_longer(
    cols = all_of(feature_cols),
    names_to = "score",
    values_to = "value"
  )

p_norm_scores <- score_long_all %>%
  ggplot(aes(x = dataset, y = value, fill = dataset)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 0.6, alpha = 0.4) +
  facet_grid(normalization ~ score, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Dataset",
    y = "Valor transformado",
    title = "Efecto de las normalizaciones sobre los scores CD8+"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(fig_dir, "domain_shift_normalization_score_distributions.png"),
  p_norm_scores,
  width = 16,
  height = 9,
  dpi = 300
)

# -------------------------
# 9. PCA tras normalización
# -------------------------

run_pca_plot <- function(df, norm_name) {
  
  pca_input <- df %>%
    select(all_of(feature_cols)) %>%
    mutate(across(everything(), ~ ifelse(is.na(.x) | !is.finite(.x), median(.x, na.rm = TRUE), .x)))
  
  pca <- prcomp(
    pca_input,
    center = TRUE,
    scale. = TRUE
  )
  
  pca_var <- (pca$sdev^2) / sum(pca$sdev^2)
  
  pca_scores <- as_tibble(pca$x[, 1:2]) %>%
    bind_cols(
      df %>%
        select(dataset, patient_alias, response_binary_common, therapy_common)
    )
  
  p <- ggplot(
    pca_scores,
    aes(x = PC1, y = PC2, color = dataset, shape = response_binary_common)
  ) +
    geom_point(size = 3, alpha = 0.85) +
    theme_classic() +
    labs(
      title = paste0("PCA tras normalización: ", norm_name),
      x = paste0("PC1 (", round(100 * pca_var[1], 1), "%)"),
      y = paste0("PC2 (", round(100 * pca_var[2], 1), "%)")
    )
  
  ggsave(
    file.path(fig_dir, paste0("domain_shift_pca_", norm_name, ".png")),
    p,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  invisible(p)
}

run_pca_plot(combined_raw, "raw")
run_pca_plot(combined_zscore, "zscore_by_dataset")
run_pca_plot(combined_rank, "rank_by_dataset")

# -------------------------
# 10. Capacidad de predecir dataset tras normalización
# -------------------------

dataset_auc_by_normalization <- map_dfr(
  unique(combined_all$normalization),
  function(norm_name) {
    
    df <- combined_all %>%
      filter(normalization == norm_name) %>%
      mutate(dataset_binary = ifelse(dataset == "GSE120575", 1, 0))
    
    fit <- glm(
      dataset_binary ~ score_exhaustion_mean +
        score_effector_mean +
        score_memory_mean +
        score_activation_mean +
        score_cycling_mean,
      data = df,
      family = binomial()
    )
    
    pred <- as.numeric(predict(fit, type = "response"))
    
    roc_obj <- roc(df$dataset_binary, pred, quiet = TRUE)
    
    tibble(
      normalization = norm_name,
      dataset_auc = as.numeric(auc(roc_obj)),
      n_rows = nrow(df)
    )
  }
)

fwrite(
  dataset_auc_by_normalization,
  file.path(res_dir, "domain_shift_dataset_auc_by_normalization.tsv"),
  sep = "\t"
)

print(dataset_auc_by_normalization)

message("Corrección exploratoria de domain shift completada.")