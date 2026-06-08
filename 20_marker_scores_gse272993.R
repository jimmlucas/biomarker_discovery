# =========================
# 20_marker_scores_gse272993.R
# =========================

base_dir <- "/Users/jimmlucas/Downloads/R_MASTER_RNASEQ"
setwd(base_dir)

source(file.path(base_dir, "00_config.R"))

library(Seurat)
library(tidyverse)
library(data.table)

# -------------------------
# 1. Cargar objeto Seurat GSE272993
# -------------------------

rds_file <- file.path(
  raw_dir,
  "GSE272993",
  "GSE272993_cd8_nn_labeled_FINAL.RDS"
)

if (!file.exists(rds_file)) {
  rds_gz <- paste0(rds_file, ".gz")
  
  if (file.exists(rds_gz)) {
    R.utils::gunzip(
      filename = rds_gz,
      destname = rds_file,
      overwrite = TRUE,
      remove = FALSE
    )
  } else {
    stop("No se encuentra el objeto RDS ni RDS.gz de GSE272993.")
  }
}

obj <- readRDS(rds_file)

DefaultAssay(obj) <- "RNA"

# -------------------------
# 2. Definir genes marcadores comparables
# -------------------------

programs <- list(
  exhaustion = c("PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX"),
  effector   = c("GZMB", "PRF1", "IFNG"),
  memory     = c("IL7R", "CCR7"),
  activation = c("CD69", "HLA-DRA"),
  cycling    = c("MKI67")
)

genes_available <- rownames(obj)

programs_available <- lapply(programs, function(g) {
  intersect(g, genes_available)
})

print(programs_available)

# -------------------------
# 3. Calcular scores por célula
# -------------------------

for (program_name in names(programs_available)) {
  
  genes <- programs_available[[program_name]]
  
  if (length(genes) == 0) {
    warning("No hay genes disponibles para: ", program_name)
    next
  }
  
  obj <- AddModuleScore(
    object = obj,
    features = list(genes),
    name = paste0("score_", program_name, "_")
  )
}

# AddModuleScore crea columnas tipo score_exhaustion_1
meta <- obj@meta.data %>%
  tibble::rownames_to_column("cell_id")

score_cols_raw <- grep("^score_.*_1$", colnames(meta), value = TRUE)

meta_scores <- meta %>%
  rename(
    score_exhaustion = any_of("score_exhaustion_1"),
    score_effector   = any_of("score_effector_1"),
    score_memory     = any_of("score_memory_1"),
    score_activation = any_of("score_activation_1"),
    score_cycling    = any_of("score_cycling_1")
  )

# -------------------------
# 4. Añadir respuesta binaria si no existe
# -------------------------

if (!"response_binary" %in% colnames(meta_scores)) {
  meta_scores <- meta_scores %>%
    mutate(
      response_binary = case_when(
        response %in% c("CR", "PR") ~ "Responder",
        response %in% c("PD", "SD") ~ "Non_responder",
        TRUE ~ NA_character_
      )
    )
}

# -------------------------
# 5. Agregar a paciente/timepoint
# -------------------------

score_cols <- c(
  "score_exhaustion",
  "score_effector",
  "score_memory",
  "score_activation",
  "score_cycling"
)

patient_scores <- meta_scores %>%
  filter(!is.na(response_binary)) %>%
  group_by(patient_alias, response, response_binary, treatment, timepoint) %>%
  summarise(
    n_cells = n(),
    across(
      all_of(score_cols),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        median = ~ median(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

# -------------------------
# 6. Guardar
# -------------------------

fwrite(
  meta_scores %>% select(cell_id, patient_alias, response, response_binary, treatment, timepoint, all_of(score_cols)),
  file.path(res_dir, "GSE272993_cell_level_marker_program_scores.tsv"),
  sep = "\t"
)

fwrite(
  patient_scores,
  file.path(res_dir, "GSE272993_patient_timepoint_marker_program_scores.tsv"),
  sep = "\t"
)

saveRDS(
  patient_scores,
  file.path(proc_dir, "GSE272993_patient_timepoint_marker_program_scores.rds")
)

# -------------------------
# 7. Tests exploratorios
# -------------------------

score_long <- patient_scores %>%
  pivot_longer(
    cols = starts_with("score_"),
    names_to = "score",
    values_to = "value"
  )

wilcox_scores <- score_long %>%
  group_by(score) %>%
  summarise(
    p_value = wilcox.test(value ~ response_binary)$p.value,
    median_responder = median(value[response_binary == "Responder"], na.rm = TRUE),
    median_non_responder = median(value[response_binary == "Non_responder"], na.rm = TRUE),
    delta_median = median_responder - median_non_responder,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    direction = case_when(
      delta_median > 0 ~ "Higher_in_responders",
      delta_median < 0 ~ "Lower_in_responders",
      TRUE ~ "No_difference"
    )
  ) %>%
  arrange(p_adj)

fwrite(
  wilcox_scores,
  file.path(res_dir, "GSE272993_wilcox_marker_program_scores_response.tsv"),
  sep = "\t"
)

print(wilcox_scores)

# -------------------------
# 8. Figura
# -------------------------

p_scores <- score_long %>%
  ggplot(aes(x = response_binary, y = value, fill = response_binary)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
  facet_wrap(~ score, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Respuesta clínica",
    y = "Marker program score medio",
    title = "Programas funcionales CD8+ comparables en GSE272993"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(fig_dir, "GSE272993_marker_program_scores_by_response.png"),
  p_scores,
  width = 12,
  height = 8,
  dpi = 300
)

rm(obj)
gc()

message("Scores funcionales comparables de GSE272993 generados.")