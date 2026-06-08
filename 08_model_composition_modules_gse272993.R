# =========================
# 08_model_composition_modules_gse272993.R
# =========================
# Comparación de modelos:
# 1) composición + todos los module scores
# 2) composición + módulos seleccionados
#
# El objetivo no es maximizar el AUC aparente, sino evaluar generalización
# mediante predicciones out-of-fold agrupadas por paciente.

source("scripts/00_config.R")

library(tidyverse)
library(data.table)
library(glmnet)
library(pROC)

set.seed(123)

model_df_modules <- fread(
  file.path(res_dir, "GSE272993_model_matrix_composition_plus_modules.tsv")
) %>%
  select(-any_of(c("fold", "fold.x", "fold.y"))) %>%
  filter(!is.na(response_binary)) %>%
  mutate(y = ifelse(response_binary == "Responder", 1, 0))

patient_fold_df <- fread(
  file.path(res_dir, "GSE272993_patient_folds_5fold.tsv")
)

model_df_modules <- model_df_modules %>%
  left_join(patient_fold_df, by = "patient_alias")

foldid <- model_df_modules$fold
y <- model_df_modules$y

fit_oof_elasticnet <- function(data, feature_cols, model_name, output_prefix) {
  
  data <- data %>%
    mutate(
      across(
        all_of(feature_cols),
        ~ ifelse(is.na(.x) | !is.finite(.x), median(.x, na.rm = TRUE), .x)
      )
    )
  
  X <- as.matrix(data %>% select(all_of(feature_cols)))
  y <- data$y
  foldid <- data$fold
  
  cvfit <- cv.glmnet(
    x = X,
    y = y,
    family = "binomial",
    alpha = 0.5,
    foldid = foldid,
    type.measure = "auc"
  )
  
  pred_prob <- as.numeric(
    predict(cvfit, newx = X, s = "lambda.min", type = "response")
  )
  
  auc_apparent <- as.numeric(auc(roc(y, pred_prob, quiet = TRUE)))
  
  pred_oof <- rep(NA_real_, nrow(data))
  
  for (k in sort(unique(foldid))) {
    train_idx <- which(foldid != k)
    test_idx  <- which(foldid == k)
    
    fit_k <- cv.glmnet(
      x = X[train_idx, ],
      y = y[train_idx],
      family = "binomial",
      alpha = 0.5,
      type.measure = "auc"
    )
    
    pred_oof[test_idx] <- as.numeric(
      predict(fit_k, newx = X[test_idx, ], s = "lambda.min", type = "response")
    )
  }
  
  auc_oof <- as.numeric(auc(roc(y, pred_oof, quiet = TRUE)))
  
  coef_df <- as.matrix(coef(cvfit, s = "lambda.min")) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("feature") %>%
    rename(coef = 2) %>%
    filter(coef != 0) %>%
    arrange(desc(abs(coef)))
  
  pred_df <- data %>%
    select(patient_alias, response, response_binary, treatment, timepoint, y, fold) %>%
    mutate(
      pred_prob = pred_prob,
      pred_prob_oof = pred_oof
    )
  
  summary_df <- tibble(
    dataset = "GSE272993",
    model = model_name,
    alpha = 0.5,
    lambda_min = cvfit$lambda.min,
    apparent_auc = auc_apparent,
    n_rows = nrow(data),
    n_patients = n_distinct(data$patient_alias),
    n_features = length(feature_cols),
    n_responders = sum(y == 1),
    n_non_responders = sum(y == 0),
    oof_auc = auc_oof
  )
  
  fwrite(
    coef_df,
    file.path(res_dir, paste0(output_prefix, "_coefficients.tsv")),
    sep = "\t"
  )
  
  fwrite(
    pred_df,
    file.path(res_dir, paste0(output_prefix, "_predictions.tsv")),
    sep = "\t"
  )
  
  fwrite(
    summary_df,
    file.path(res_dir, paste0(output_prefix, "_model_summary.tsv")),
    sep = "\t"
  )
  
  return(summary_df)
}

# Modelo con todos los module scores disponibles.
all_feature_cols <- grep(
  "^freq_|^Jo_top250_|^module\\d+_mean",
  colnames(model_df_modules),
  value = TRUE
)

all_feature_cols <- all_feature_cols[
  sapply(model_df_modules %>% select(all_of(all_feature_cols)), is.numeric)
]

summary_modules <- fit_oof_elasticnet(
  data = model_df_modules,
  feature_cols = all_feature_cols,
  model_name = "ElasticNet_composition_plus_modules",
  output_prefix = "GSE272993_elasticnet_composition_modules"
)

# Modelo restringido guiado por resultados univariantes.
selected_features <- c(
  grep("^freq_", colnames(model_df_modules), value = TRUE),
  "Jo_top250_cm1_mean",
  "module21_mean"
)

selected_features <- selected_features[
  selected_features %in% colnames(model_df_modules)
]

summary_restricted <- fit_oof_elasticnet(
  data = model_df_modules,
  feature_cols = selected_features,
  model_name = "ElasticNet_composition_plus_selected_modules",
  output_prefix = "GSE272993_elasticnet_restricted"
)

# Comparación final de modelos.
model_comparison <- bind_rows(
  fread(file.path(res_dir, "GSE272993_elasticnet_composition_model_summary_with_oof.tsv")),
  summary_modules,
  summary_restricted
)

fwrite(
  model_comparison,
  file.path(res_dir, "GSE272993_model_comparison_summary.tsv"),
  sep = "\t"
)

print(model_comparison)