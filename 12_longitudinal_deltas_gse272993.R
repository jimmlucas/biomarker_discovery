# =========================
# 12_longitudinal_deltas_gse272993.R
# =========================

base_dir <- "/Users/jimmlucas/Downloads/R_MASTER_RNASEQ"
setwd(base_dir)

source(file.path(base_dir, "00_config.R"))

suppressPackageStartupMessages({
  library(magrittr)   # %>%
  library(dplyr)
  library(tidyr)
  library(data.table)
  library(ggplot2)
  library(glmnet)
  library(pROC)
})

set.seed(123)

comp_wide <- fread(
  file.path(res_dir, "GSE272993_patient_timepoint_metacluster_composition_wide.tsv")
)

feature_cols <- grep("^freq_", colnames(comp_wide), value = TRUE)

baseline_names <- c("Baseline")
early_names    <- c("Follow Up 1")

# =============================================================
# ANÁLISIS PAIRED: Wilcoxon pareado baseline vs Follow Up 1
# por paciente — responde al comentario de Laura (M3)
# =============================================================

paired_data <- comp_wide %>%
  filter(timepoint %in% c(baseline_names, early_names)) %>%
  mutate(
    timepoint_simple = case_when(
      timepoint %in% baseline_names ~ "baseline",
      timepoint %in% early_names   ~ "early",
      TRUE ~ NA_character_
    )
  ) %>%
  select(patient_alias, response_binary, timepoint_simple, all_of(feature_cols)) %>%
  distinct()

# Solo pacientes con AMBOS puntos temporales
paired_patients <- paired_data %>%
  group_by(patient_alias) %>%
  filter(all(c("baseline", "early") %in% timepoint_simple)) %>%
  ungroup()

wilcox_paired <- lapply(feature_cols, function(f) {
  
  wide_f <- paired_patients %>%
    select(patient_alias, response_binary, timepoint_simple, all_of(f)) %>%
    pivot_wider(names_from = timepoint_simple, values_from = all_of(f)) %>%
    filter(!is.na(baseline), !is.na(early))
  
  if (nrow(wide_f) < 4) return(NULL)
  
  test_result <- wilcox.test(
    wide_f$early,
    wide_f$baseline,
    paired = TRUE,       # <-- PAIRED
    exact  = FALSE
  )
  
  tibble(
    feature               = f,
    n_paired_patients     = nrow(wide_f),
    median_baseline       = median(wide_f$baseline, na.rm = TRUE),
    median_early          = median(wide_f$early,    na.rm = TRUE),
    median_delta          = median(wide_f$early - wide_f$baseline, na.rm = TRUE),
    p_value_paired        = test_result$p.value
  )
}) %>%
  bind_rows() %>%
  mutate(
    p_adj_paired = p.adjust(p_value_paired, method = "BH"),
    direction = case_when(
      median_delta > 0 ~ "Aumenta_FU1_vs_Baseline",
      median_delta < 0 ~ "Disminuye_FU1_vs_Baseline",
      TRUE             ~ "Sin_cambio"
    )
  ) %>%
  arrange(p_adj_paired)

fwrite(
  wilcox_paired,
  file.path(res_dir, "GSE272993_wilcox_paired_longitudinal.tsv"),
  sep = "\t"
)

print(wilcox_paired)