# =========================
# 04_qc_explore_gse272993.R
# =========================
# Exploración inicial y control de calidad usando metadata.
# Este paso evita cargar el objeto Seurat completo.

source("scripts/00_config.R")

library(tidyverse)
library(data.table)

meta <- readRDS(file.path(proc_dir, "GSE272993_metadata_only.rds"))

# Resumen global de la cohorte.
qc_summary <- meta %>%
  summarise(
    n_cells = n(),
    n_patients = n_distinct(patient_alias),
    n_samples = n_distinct(sample),
    n_timepoints = n_distinct(timepoint),
    n_treatments = n_distinct(treatment),
    n_responses = n_distinct(response),
    median_nCount_RNA = median(nCount_RNA, na.rm = TRUE),
    median_nFeature_RNA = median(nFeature_RNA, na.rm = TRUE),
    median_percent_mt = median(percent.mt, na.rm = TRUE)
  )

print(qc_summary)

fwrite(
  qc_summary,
  file.path(res_dir, "GSE272993_qc_summary.tsv"),
  sep = "\t"
)

# Conteos descriptivos.
tab_patient <- meta %>%
  count(patient_alias, response, treatment, sort = TRUE)

tab_sample <- meta %>%
  count(sample, patient_alias, timepoint, response, treatment, sort = TRUE)

tab_cluster <- meta %>%
  count(metaclusters, response, sort = TRUE)

tab_time <- meta %>%
  count(timepoint, response, sort = TRUE)

fwrite(tab_patient, file.path(res_dir, "GSE272993_cells_by_patient.tsv"), sep = "\t")
fwrite(tab_sample,  file.path(res_dir, "GSE272993_cells_by_sample.tsv"), sep = "\t")
fwrite(tab_cluster, file.path(res_dir, "GSE272993_cells_by_metacluster_response.tsv"), sep = "\t")
fwrite(tab_time,    file.path(res_dir, "GSE272993_cells_by_timepoint_response.tsv"), sep = "\t")

# Figuras QC.
p_nfeature <- ggplot(meta, aes(x = response, y = nFeature_RNA, fill = response)) +
  geom_violin(trim = TRUE, scale = "width") +
  theme_classic() +
  labs(x = "Respuesta clínica", y = "nFeature_RNA")

ggsave(
  file.path(fig_dir, "GSE272993_nFeature_RNA_by_response.png"),
  p_nfeature,
  width = 6,
  height = 5,
  dpi = 300
)

p_mt <- ggplot(meta, aes(x = response, y = percent.mt, fill = response)) +
  geom_violin(trim = TRUE, scale = "width") +
  theme_classic() +
  labs(x = "Respuesta clínica", y = "Porcentaje mitocondrial")

ggsave(
  file.path(fig_dir, "GSE272993_percent_mt_by_response.png"),
  p_mt,
  width = 6,
  height = 5,
  dpi = 300
)

p_ncount <- ggplot(meta, aes(x = response, y = nCount_RNA, fill = response)) +
  geom_violin(trim = TRUE, scale = "width") +
  theme_classic() +
  labs(x = "Respuesta clínica", y = "nCount_RNA")

ggsave(
  file.path(fig_dir, "GSE272993_nCount_RNA_by_response.png"),
  p_ncount,
  width = 6,
  height = 5,
  dpi = 300
)