# =========================
# 11_bootstrap_stability_gse272993.R
# =========================

base_dir <- "/Users/jimmlucas/Downloads/R_MASTER_RNASEQ"
setwd(base_dir)

source(file.path(base_dir, "00_config.R"))

library(tidyverse)
library(data.table)
library(glmnet)

set.seed(123)

comp_wide <- fread(
  file.path(res_dir, "GSE272993_patient_timepoint_metacluster_composition_wide.tsv")
)

model_df <- comp_wide %>%
  filter(!is.na(response_binary)) %>%
  mutate(y = ifelse(response_binary == "Responder", 1, 0)) %>%
  select(-any_of(c("fold", "fold.x", "fold.y")))

feature_cols <- grep("^freq_", colnames(model_df), value = TRUE)

patients <- unique(model_df$patient_alias)

n_boot <- 500
boot_coef_list <- vector("list", n_boot)

for (b in seq_len(n_boot)) {
  
  if (b %% 50 == 0) {
    message("Bootstrap ", b, " / ", n_boot)
  }
  
  sampled_patients <- sample(
    patients,
    size = length(patients),
    replace = TRUE
  )
  
  boot_df <- purrr::map_dfr(sampled_patients, function(p) {
    model_df %>% filter(patient_alias == p)
  })
  
  X_boot <- as.matrix(boot_df %>% select(all_of(feature_cols)))
  y_boot <- boot_df$y
  
  if (length(unique(y_boot)) < 2) {
    next
  }
  
  fit <- cv.glmnet(
    x = X_boot,
    y = y_boot,
    family = "binomial",
    alpha = 0.5,
    type.measure = "auc"
  )
  
  coef_b <- as.matrix(coef(fit, s = "lambda.min")) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("feature") %>%
    rename(coef = 2) %>%
    filter(feature %in% feature_cols) %>%
    mutate(
      bootstrap = b,
      selected = coef != 0
    )
  
  boot_coef_list[[b]] <- coef_b
}

boot_coef <- bind_rows(boot_coef_list)

fwrite(
  boot_coef,
  file.path(res_dir, "GSE272993_bootstrap_coefficients_long.tsv"),
  sep = "\t"
)

coef_stability <- boot_coef %>%
  group_by(feature) %>%
  summarise(
    mean_coef = mean(coef, na.rm = TRUE),
    sd_coef = sd(coef, na.rm = TRUE),
    median_coef = median(coef, na.rm = TRUE),
    selection_frequency = mean(selected, na.rm = TRUE),
    n_boot_valid = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(selection_frequency), desc(abs(mean_coef)))

fwrite(
  coef_stability,
  file.path(res_dir, "GSE272993_bootstrap_feature_stability.tsv"),
  sep = "\t"
)

p_stability <- coef_stability %>%
  mutate(feature = reorder(feature, selection_frequency)) %>%
  ggplot(aes(
    x = feature,
    y = selection_frequency,
    fill = mean_coef > 0
  )) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(
    x = "Subestado CD8+",
    y = "Frecuencia de selección",
    fill = "Coeficiente positivo",
    title = "Estabilidad de selección de variables en bootstrap"
  )

ggsave(
  file.path(fig_dir, "GSE272993_bootstrap_feature_stability.png"),
  p_stability,
  width = 8,
  height = 6,
  dpi = 300
)

print(coef_stability)