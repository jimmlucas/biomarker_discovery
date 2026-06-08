# =========================
# 10_robustness_permutation_gse272993.R
# =========================

base_dir <- "/Users/jimmlucas/Downloads/R_MASTER_RNASEQ"
setwd(base_dir)

source(file.path(base_dir, "00_config.R"))

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

# -------------------------
# Recuperar o generar folds agrupados por paciente
# -------------------------

fold_file <- file.path(res_dir, "GSE272993_patient_folds_5fold.tsv")

if (file.exists(fold_file)) {
  
  message("Archivo de folds encontrado.")
  patient_fold_df <- fread(fold_file)
  
} else {
  
  message("No existe archivo de folds. Intentando reconstruir desde predicciones OOF...")
  
  oof_file <- file.path(res_dir, "GSE272993_elasticnet_composition_oof_predictions.tsv")
  
  if (file.exists(oof_file)) {
    
    oof_pred <- fread(oof_file)
    
    patient_fold_df <- oof_pred %>%
      distinct(patient_alias, fold)
    
    fwrite(patient_fold_df, fold_file, sep = "\t")
    
    message("Folds reconstruidos desde predicciones OOF.")
    
  } else {
    
    message("No existe archivo OOF. Generando folds nuevos.")
    
    patients <- unique(model_df$patient_alias)
    
    patient_fold_df <- tibble(
      patient_alias = patients,
      fold = sample(rep(1:5, length.out = length(patients)))
    )
    
    fwrite(patient_fold_df, fold_file, sep = "\t")
  }
}

model_df <- model_df %>%
  left_join(patient_fold_df, by = "patient_alias")

foldid <- model_df$fold

stopifnot(length(foldid) == nrow(model_df))
stopifnot(!any(is.na(foldid)))
stopifnot(nrow(X) == nrow(model_df))
stopifnot(length(y) == nrow(model_df))

print(table(foldid, model_df$response_binary))

# -------------------------
# Permutación de etiquetas
# -------------------------

n_perm <- 500
perm_results <- vector("list", n_perm)

for (i in seq_len(n_perm)) {
  
  if (i %% 50 == 0) {
    message("Permutación ", i, " / ", n_perm)
  }
  
  y_perm <- sample(y)
  pred_oof <- rep(NA_real_, length(y_perm))
  
  for (k in sort(unique(foldid))) {
    
    train_idx <- which(foldid != k)
    test_idx  <- which(foldid == k)
    
    fit_k <- cv.glmnet(
      x = X[train_idx, ],
      y = y_perm[train_idx],
      family = "binomial",
      alpha = 0.5,
      type.measure = "auc"
    )
    
    pred_oof[test_idx] <- as.numeric(
      predict(
        fit_k,
        newx = X[test_idx, ],
        s = "lambda.min",
        type = "response"
      )
    )
  }
  
  auc_perm <- as.numeric(
    auc(
      roc(y_perm, pred_oof, quiet = TRUE)
    )
  )
  
  perm_results[[i]] <- tibble(
    permutation = i,
    auc_oof_permuted = auc_perm
  )
}

perm_results <- bind_rows(perm_results)

fwrite(
  perm_results,
  file.path(res_dir, "GSE272993_permutation_auc_results.tsv"),
  sep = "\t"
)

# -------------------------
# Comparación con AUROC real
# -------------------------

real_summary <- fread(
  file.path(res_dir, "GSE272993_elasticnet_composition_model_summary_with_oof.tsv")
)

real_auc <- real_summary$oof_auc[1]

p_perm <- ggplot(perm_results, aes(x = auc_oof_permuted)) +
  geom_histogram(bins = 30, fill = "gray75", color = "white") +
  geom_vline(
    xintercept = real_auc,
    color = "#D7191C",
    linewidth = 1.2
  ) +
  theme_classic() +
  labs(
    x = "AUROC OOF con etiquetas permutadas",
    y = "Número de permutaciones",
    title = "Control negativo por permutación de etiquetas"
  )

ggsave(
  file.path(fig_dir, "GSE272993_permutation_auc_distribution.png"),
  p_perm,
  width = 7,
  height = 5,
  dpi = 300
)

perm_summary <- tibble(
  dataset = "GSE272993",
  model = "ElasticNet_composition_only",
  real_oof_auc = real_auc,
  mean_permuted_auc = mean(perm_results$auc_oof_permuted),
  sd_permuted_auc = sd(perm_results$auc_oof_permuted),
  empirical_p_value = mean(perm_results$auc_oof_permuted >= real_auc),
  n_permutations = n_perm
)

fwrite(
  perm_summary,
  file.path(res_dir, "GSE272993_permutation_summary.tsv"),
  sep = "\t"
)

print(perm_summary)