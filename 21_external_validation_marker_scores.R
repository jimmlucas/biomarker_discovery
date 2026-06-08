# =========================
# 21_external_validation_marker_scores.R
# =========================
# Validación externa cross-platform:
# entrenamiento en GSE272993 y prueba externa en GSE120575
# usando scores funcionales CD8+ comparables.

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
  as_tibble()

external_df <- fread(
  file.path(res_dir, "GSE120575_patient_timepoint_cd8_program_scores.tsv"),
  data.table = FALSE
) %>%
  as_tibble()

# -------------------------
# 2. Definir variables comparables
# -------------------------

feature_cols <- c(
  "score_exhaustion_mean",
  "score_effector_mean",
  "score_memory_mean",
  "score_activation_mean",
  "score_cycling_mean"
)

stopifnot(all(feature_cols %in% colnames(train_df)))
stopifnot(all(feature_cols %in% colnames(external_df)))

# -------------------------
# 3. Preparar entrenamiento GSE272993
# -------------------------

train_model_df <- train_df %>%
  filter(!is.na(response_binary)) %>%
  mutate(
    y = ifelse(response_binary == "Responder", 1, 0),
    dataset = "GSE272993"
  ) %>%
  select(
    dataset,
    patient_alias,
    response,
    response_binary,
    treatment,
    timepoint,
    n_cells,
    y,
    all_of(feature_cols)
  )

# -------------------------
# 4. Preparar validación externa GSE120575
# -------------------------

external_model_df <- external_df %>%
  filter(
    !is.na(response_binary_strict),
    response_binary_strict != "Mixed"
  ) %>%
  mutate(
    y = ifelse(response_binary_strict == "Responder", 1, 0),
    dataset = "GSE120575",
    response = response_binary_strict,
    treatment = therapy
  ) %>%
  select(
    dataset,
    patient_alias,
    response,
    response_binary = response_binary_strict,
    treatment,
    timepoint,
    n_cells,
    y,
    all_of(feature_cols)
  )

# -------------------------
# 5. Imputación sencilla
# -------------------------

for (f in feature_cols) {
  med_train <- median(train_model_df[[f]], na.rm = TRUE)
  train_model_df[[f]][is.na(train_model_df[[f]]) | !is.finite(train_model_df[[f]])] <- med_train
  external_model_df[[f]][is.na(external_model_df[[f]]) | !is.finite(external_model_df[[f]])] <- med_train
}

# -------------------------
# 6. Estandarización usando SOLO entrenamiento
# -------------------------

train_means <- sapply(train_model_df[feature_cols], mean, na.rm = TRUE)
train_sds   <- sapply(train_model_df[feature_cols], sd, na.rm = TRUE)

train_sds[train_sds == 0 | is.na(train_sds)] <- 1

scale_with_train <- function(df, feature_cols, means, sds) {
  df_scaled <- df
  for (f in feature_cols) {
    df_scaled[[f]] <- (df_scaled[[f]] - means[[f]]) / sds[[f]]
  }
  df_scaled
}

train_scaled <- scale_with_train(train_model_df, feature_cols, train_means, train_sds)
external_scaled <- scale_with_train(external_model_df, feature_cols, train_means, train_sds)

X_train    <- as.matrix(train_scaled %>% select(all_of(feature_cols)))
y_train    <- train_scaled$y

X_external <- as.matrix(external_scaled %>% select(all_of(feature_cols)))
y_external <- external_scaled$y

# -------------------------
# 7. Entrenar ElasticNet en GSE272993
# -------------------------

cvfit <- cv.glmnet(
  x = X_train,
  y = y_train,
  family = "binomial",
  alpha = 0.5,
  type.measure = "auc"
)

best_lambda <- cvfit$lambda.min

pred_train    <- as.numeric(predict(cvfit, newx = X_train,    s = "lambda.min", type = "response"))
pred_external <- as.numeric(predict(cvfit, newx = X_external, s = "lambda.min", type = "response"))

# -------------------------
# 8. Métricas AUROC y AUPRC
# -------------------------

roc_train    <- roc(y_train,    pred_train,    quiet = TRUE)
roc_external <- roc(y_external, pred_external, quiet = TRUE)

auc_train    <- as.numeric(auc(roc_train))
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

auprc_train    <- pr_train$auc.integral
auprc_external <- pr_external$auc.integral

# -------------------------
# 9. Coeficientes
# -------------------------

coef_df <- as.matrix(coef(cvfit, s = "lambda.min")) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("feature") %>%
  rename(coef = 2) %>%
  arrange(desc(abs(coef)))

fwrite(
  coef_df,
  file.path(res_dir, "external_validation_marker_scores_coefficients.tsv"),
  sep = "\t"
)

# -------------------------
# 10. Guardar predicciones
# -------------------------

train_predictions <- train_model_df %>%
  mutate(pred_prob = pred_train, validation_type = "training_internal")

external_predictions <- external_model_df %>%
  mutate(pred_prob = pred_external, validation_type = "external_validation")

all_predictions <- bind_rows(train_predictions, external_predictions)

fwrite(train_predictions,  file.path(res_dir, "GSE272993_marker_scores_training_predictions.tsv"),  sep = "\t")
fwrite(external_predictions, file.path(res_dir, "GSE120575_marker_scores_external_predictions.tsv"), sep = "\t")
fwrite(all_predictions,    file.path(res_dir, "marker_scores_training_external_predictions.tsv"),    sep = "\t")

# -------------------------
# 11. Resumen de rendimiento
# -------------------------

performance_summary <- tibble(
  model                  = "ElasticNet_marker_program_scores",
  train_dataset          = "GSE272993",
  external_dataset       = "GSE120575",
  alpha                  = 0.5,
  lambda_min             = best_lambda,
  n_features             = length(feature_cols),
  train_n_rows           = nrow(train_model_df),
  train_n_patients       = n_distinct(train_model_df$patient_alias),
  train_n_responders     = sum(y_train == 1),
  train_n_non_responders = sum(y_train == 0),
  train_auc              = auc_train,
  train_auprc            = auprc_train,
  external_n_rows           = nrow(external_model_df),
  external_n_patients       = n_distinct(external_model_df$patient_alias),
  external_n_responders     = sum(y_external == 1),
  external_n_non_responders = sum(y_external == 0),
  external_auc              = auc_external,
  external_auprc            = auprc_external
)

fwrite(
  performance_summary,
  file.path(res_dir, "external_validation_marker_scores_performance_summary.tsv"),
  sep = "\t"
)

print(performance_summary)

# -------------------------
# 12. Curva ROC comparativa (ggplot2, ejes 0-1 estrictos)
# -------------------------

roc_train_df <- data.frame(
  specificity = roc_train$specificities,
  sensitivity = roc_train$sensitivities,
  cohort      = paste0("GSE272993 entrenamiento (AUROC = ", round(auc_train, 3), ")")
)

roc_external_df <- data.frame(
  specificity = roc_external$specificities,
  sensitivity = roc_external$sensitivities,
  cohort      = paste0("GSE120575 externo (AUROC = ", round(auc_external, 3), ")")
)

roc_combined_df <- bind_rows(roc_train_df, roc_external_df)

p_roc_combined <- ggplot(roc_combined_df, aes(x = 1 - specificity, y = sensitivity, color = cohort)) +
  geom_line(linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0), name = "1 - Especificidad (FPR)") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), name = "Sensibilidad (TPR)") +
  scale_color_manual(values = c("#2C7BB6", "#D7191C")) +
  theme_classic() +
  theme(
    legend.position  = "bottom",
    legend.title     = element_blank(),
    legend.text      = element_text(size = 10)
  ) +
  labs(title = "Validación externa con programas funcionales CD8+")

ggsave(
  file.path(fig_dir, "external_validation_marker_scores_roc.png"),
  p_roc_combined,
  width = 7, height = 6, dpi = 300
)

# -------------------------
# 13. Boxplot de predicciones externas
# -------------------------

p_external_pred <- external_predictions %>%
  ggplot(aes(x = response_binary, y = pred_prob, fill = response_binary)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.6) +
  theme_classic() +
  labs(
    x     = "Respuesta clínica GSE120575",
    y     = "Probabilidad predicha de respuesta",
    title = "Validación externa del modelo entrenado en GSE272993"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(fig_dir, "GSE120575_external_validation_predicted_probabilities.png"),
  p_external_pred,
  width = 7, height = 5, dpi = 300
)

# -------------------------
# 14. Guardar parámetros de escalado
# -------------------------

scaling_params <- tibble(
  feature    = feature_cols,
  train_mean = as.numeric(train_means),
  train_sd   = as.numeric(train_sds)
)

fwrite(
  scaling_params,
  file.path(res_dir, "external_validation_marker_scores_scaling_params.tsv"),
  sep = "\t"
)

message("Validación externa con marker scores completada.")