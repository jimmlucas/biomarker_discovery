# =========================
# 06_model_composition_gse272993.R
# =========================
# Modelo basal: ElasticNet usando únicamente composición de subestados CD8+.
# La validación se hace agrupando por paciente para evitar leakage entre timepoints.

source("00_config.R")

library(tidyverse)
library(data.table)
library(glmnet)
library(pROC)

set.seed(123)

comp_wide <- fread(
  file.path(res_dir, "GSE272993_patient_timepoint_metacluster_composition_wide.tsv")
)

model_df <- comp_wide %>%
  filter(!is.na(response_binary)) %>%
  mutate(y = ifelse(response_binary == "Responder", 1, 0)) %>%
  select(-any_of(c("fold", "fold.x", "fold.y")))

feature_cols <- grep("^freq_", colnames(model_df), value = TRUE)

X <- as.matrix(model_df %>% select(all_of(feature_cols)))
y <- model_df$y

# Crear folds agrupados por paciente
patients <- unique(model_df$patient_alias)

patient_fold_df <- tibble(
  patient_alias = patients,
  fold          = sample(rep(1:5, length.out = length(patients)))
)

fwrite(
  patient_fold_df,
  file.path(res_dir, "GSE272993_patient_folds_5fold.tsv"),
  sep = "\t"
)

model_df <- model_df %>%
  left_join(patient_fold_df, by = "patient_alias")

foldid <- model_df$fold

stopifnot(length(foldid) == nrow(model_df))
stopifnot(nrow(X) == nrow(model_df))
stopifnot(length(y) == nrow(model_df))

print(table(foldid, model_df$response_binary))

# ElasticNet con alpha = 0.5
cvfit <- cv.glmnet(
  x            = X,
  y            = y,
  family       = "binomial",
  alpha        = 0.5,
  foldid       = foldid,
  type.measure = "auc"
)

best_lambda <- cvfit$lambda.min

png(
  file.path(fig_dir, "GSE272993_elasticnet_composition_cv_auc.png"),
  width = 1600, height = 1200, res = 200
)
plot(cvfit)
dev.off()

coef_df <- as.matrix(coef(cvfit, s = "lambda.min")) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("feature") %>%
  rename(coef = 2) %>%
  filter(coef != 0) %>%
  arrange(desc(abs(coef)))

fwrite(
  coef_df,
  file.path(res_dir, "GSE272993_elasticnet_composition_coefficients.tsv"),
  sep = "\t"
)

# AUROC aparente
pred_prob    <- as.numeric(predict(cvfit, newx = X, s = "lambda.min", type = "response"))
auc_apparent <- as.numeric(auc(roc(y, pred_prob, quiet = TRUE)))

# Predicciones out-of-fold
pred_oof <- rep(NA_real_, nrow(model_df))

for (k in sort(unique(foldid))) {
  train_idx <- which(foldid != k)
  test_idx  <- which(foldid == k)
  
  fit_k <- cv.glmnet(
    x            = X[train_idx, ],
    y            = y[train_idx],
    family       = "binomial",
    alpha        = 0.5,
    type.measure = "auc"
  )
  
  pred_oof[test_idx] <- as.numeric(
    predict(fit_k, newx = X[test_idx, ], s = "lambda.min", type = "response")
  )
}

roc_oof <- roc(y, pred_oof, quiet = TRUE)
auc_oof <- as.numeric(auc(roc_oof))

# Curva ROC OOF (ggplot2, ejes 0-1 estrictos)
roc_oof_df <- data.frame(
  specificity = roc_oof$specificities,
  sensitivity = roc_oof$sensitivities
)

p_roc_oof <- ggplot(roc_oof_df, aes(x = 1 - specificity, y = sensitivity)) +
  geom_line(color = "#2C7BB6", linewidth = 1.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0), name = "1 - Especificidad (FPR)") +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0), name = "Sensibilidad (TPR)") +
  theme_classic() +
  labs(title = paste0("ElasticNet composición CD8+ — OOF AUROC = ", round(auc_oof, 3)))

ggsave(
  file.path(fig_dir, "GSE272993_elasticnet_composition_oof_roc.png"),
  p_roc_oof,
  width = 6, height = 6, dpi = 300
)

# Guardar predicciones
pred_df <- model_df %>%
  select(patient_alias, response, response_binary, treatment, timepoint, y, fold) %>%
  mutate(
    pred_prob     = pred_prob,
    pred_prob_oof = pred_oof
  )

fwrite(
  pred_df,
  file.path(res_dir, "GSE272993_elasticnet_composition_predictions.tsv"),
  sep = "\t"
)

model_summary <- tibble(
  dataset          = "GSE272993",
  model            = "ElasticNet_composition_only",
  alpha            = 0.5,
  lambda_min       = best_lambda,
  apparent_auc     = auc_apparent,
  n_rows           = nrow(model_df),
  n_patients       = n_distinct(model_df$patient_alias),
  n_features       = length(feature_cols),
  n_responders     = sum(y == 1),
  n_non_responders = sum(y == 0),
  oof_auc          = auc_oof
)

fwrite(
  model_summary,
  file.path(res_dir, "GSE272993_elasticnet_composition_model_summary_with_oof.tsv"),
  sep = "\t"
)

print(model_summary, width = Inf)