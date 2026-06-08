# =========================
# 35_gse273718_summary_for_results.R
# =========================

base_dir <- "/Users/jimmlucas/Downloads/R_MASTER_RNASEQ"
setwd(base_dir)

source(file.path(base_dir, "00_config.R"))

library(tidyverse)
library(data.table)

# -------------------------
# 1. Cargar resultados
# -------------------------

celltype_summary <- fread(
  file.path(res_dir, "GSE273718_broad_celltype_summary.tsv"),
  data.table = FALSE
) %>% as_tibble()

cd8_state_summary <- fread(
  file.path(res_dir, "GSE273718_cd8_state_summary.tsv"),
  data.table = FALSE
) %>% as_tibble()

cd8_program_patient <- fread(
  file.path(res_dir, "GSE273718_cd8_program_patient_summary.tsv"),
  data.table = FALSE
) %>% as_tibble()

qc_summary <- fread(
  file.path(res_dir, "GSE273718_qc_summary_filtered.tsv"),
  data.table = FALSE
) %>% as_tibble()

# -------------------------
# 2. Resumen global
# -------------------------

gse273718_results_summary <- tibble(
  dataset = "GSE273718",
  n_cells_after_qc = qc_summary$n_cells_after_qc[1],
  n_patients = qc_summary$n_patients[1],
  n_samples = qc_summary$n_samples[1],
  n_cd8_cells = sum(celltype_summary$n_cells[celltype_summary$broad_celltype == "CD8_T"]),
  n_cd8_pbmc = sum(celltype_summary$n_cells[celltype_summary$broad_celltype == "CD8_T" & celltype_summary$tissue == "PBMC"]),
  n_cd8_til = sum(celltype_summary$n_cells[celltype_summary$broad_celltype == "CD8_T" & celltype_summary$tissue == "TIL"])
)

fwrite(
  gse273718_results_summary,
  file.path(res_dir, "GSE273718_results_summary.tsv"),
  sep = "\t"
)

print(gse273718_results_summary)

# -------------------------
# 3. Composición celular amplia
# -------------------------

broad_celltype_compact <- celltype_summary %>%
  select(tissue, broad_celltype, n_cells, freq) %>%
  arrange(tissue, desc(freq))

fwrite(
  broad_celltype_compact,
  file.path(res_dir, "GSE273718_broad_celltype_compact_results.tsv"),
  sep = "\t"
)

print(broad_celltype_compact, n = Inf)

# -------------------------
# 4. Composición CD8 media por tejido
# -------------------------

cd8_state_by_tissue <- cd8_state_summary %>%
  group_by(tissue, cd8_state) %>%
  summarise(
    mean_freq = mean(freq, na.rm = TRUE),
    sd_freq = sd(freq, na.rm = TRUE),
    median_freq = median(freq, na.rm = TRUE),
    total_cells = sum(n_cells),
    .groups = "drop"
  ) %>%
  arrange(tissue, desc(mean_freq))

fwrite(
  cd8_state_by_tissue,
  file.path(res_dir, "GSE273718_cd8_state_by_tissue_summary.tsv"),
  sep = "\t"
)

print(cd8_state_by_tissue, n = Inf)

# -------------------------
# 5. Diferencias PBMC vs TIL en programas CD8
# -------------------------

program_long <- cd8_program_patient %>%
  pivot_longer(
    cols = starts_with("cd8_score_"),
    names_to = "program",
    values_to = "value"
  )

program_tissue_summary <- program_long %>%
  group_by(program, tissue) %>%
  summarise(
    n_samples = n(),
    median_value = median(value, na.rm = TRUE),
    mean_value = mean(value, na.rm = TRUE),
    .groups = "drop"
  )

program_tissue_wide <- program_tissue_summary %>%
  select(program, tissue, median_value) %>%
  pivot_wider(
    names_from = tissue,
    values_from = median_value,
    names_prefix = "median_"
  ) %>%
  mutate(
    delta_TIL_minus_PBMC = median_TIL - median_PBMC,
    direction = case_when(
      delta_TIL_minus_PBMC > 0 ~ "Higher_in_TIL",
      delta_TIL_minus_PBMC < 0 ~ "Higher_in_PBMC",
      TRUE ~ "No_difference"
    )
  ) %>%
  arrange(desc(abs(delta_TIL_minus_PBMC)))

fwrite(
  program_tissue_summary,
  file.path(res_dir, "GSE273718_cd8_program_tissue_summary.tsv"),
  sep = "\t"
)

fwrite(
  program_tissue_wide,
  file.path(res_dir, "GSE273718_cd8_program_pbmc_til_comparison.tsv"),
  sep = "\t"
)

print(program_tissue_wide, n = Inf)

# -------------------------
# 6. Figura composición celular amplia
# -------------------------

p_broad <- celltype_summary %>%
  ggplot(aes(x = tissue, y = freq, fill = broad_celltype)) +
  geom_col(position = "fill") +
  theme_classic() +
  labs(
    x = "Compartimento",
    y = "Fracción celular",
    fill = "Tipo celular",
    title = "Composición celular amplia en GSE273718"
  )

ggsave(
  file.path(fig_dir, "GSE273718_broad_celltype_composition.png"),
  p_broad,
  width = 8,
  height = 5,
  dpi = 300
)

# -------------------------
# 7. Figura programas CD8 PBMC vs TIL
# -------------------------

p_programs <- program_long %>%
  ggplot(aes(x = tissue, y = value, fill = tissue)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.7) +
  facet_wrap(~ program, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Compartimento",
    y = "Score CD8",
    title = "Programas funcionales CD8+ en PBMC y TIL"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(fig_dir, "GSE273718_cd8_programs_pbmc_vs_til.png"),
  p_programs,
  width = 12,
  height = 7,
  dpi = 300
)

# -------------------------
# 8. Tabla final compacta
# -------------------------

GSE273718_final_table <- tibble(
  item = c(
    "Cells after QC",
    "Patients",
    "Samples",
    "CD8 cells",
    "CD8 PBMC cells",
    "CD8 TIL cells",
    "Main CD8 pattern in TIL",
    "Main biological interpretation"
  ),
  value = c(
    as.character(gse273718_results_summary$n_cells_after_qc),
    as.character(gse273718_results_summary$n_patients),
    as.character(gse273718_results_summary$n_samples),
    as.character(gse273718_results_summary$n_cd8_cells),
    as.character(gse273718_results_summary$n_cd8_pbmc),
    as.character(gse273718_results_summary$n_cd8_til),
    "Enrichment of exhaustion/activation-associated states",
    "GSE273718 supports compartment-dependent CD8 functional heterogeneity"
  )
)

fwrite(
  GSE273718_final_table,
  file.path(res_dir, "GSE273718_final_compact_results_table.tsv"),
  sep = "\t"
)

print(GSE273718_final_table)

message("Resumen final de GSE273718 para Resultados completado.")