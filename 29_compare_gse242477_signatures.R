#29_compare_gse242477_signatures
#modules scores por pacientes/respuesta

base_dir <- "/Users/jimmlucas/Downloads/R_MASTER_RNASEQ"
setwd(base_dir)

source(file.path(base_dir, "00_config.R"))

library(tidyverse)
library(data.table)
library(Seurat)
library(ggplot2)
library(tidyr)

in_file <- file.path(
  base_dir,
  "data/processed/GSE242477/gse242477_til_cd8_scored.rds"
)

res_dir <- file.path(base_dir, "results/GSE242477")
fig_dir <- file.path(base_dir, "figures/GSE242477")

dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(in_file)

meta <- obj@meta.data

score_cols <- grep(
  "_score1$",
  colnames(meta),
  value = TRUE
)

score_cols <- score_cols[
  grepl(
    "Cytotoxicity|IFN|Activation|Exhaustion|Proliferation",
    score_cols
  )
]

print(score_cols)

if (length(score_cols) == 0) {
  stop("No module score columns found.")
}

patient_scores <- meta %>%
  filter(!is.na(response_group)) %>%
  filter(response_group != "Excluded") %>%
  group_by(
    patient_id,
    response_group
  ) %>%
  summarise(
    n_cd8_like_cells = n(),
    across(
      all_of(score_cols),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  arrange(response_group, patient_id)

write.csv(
  patient_scores,
  file.path(res_dir, "gse242477_patient_level_cd8_module_scores.csv"),
  row.names = FALSE
)

print(patient_scores)

long_scores <- patient_scores %>%
  pivot_longer(
    cols = all_of(score_cols),
    names_to = "program",
    values_to = "score"
  ) %>%
  mutate(
    program = gsub("_score1", "", program),
    response_group = factor(
      response_group,
      levels = c("Non-responder", "Mixed responder")
    )
  )

write.csv(
  long_scores,
  file.path(res_dir, "gse242477_long_cd8_module_scores.csv"),
  row.names = FALSE
)

p1 <- ggplot(
  long_scores,
  aes(
    x = response_group,
    y = score,
    color = response_group
  )
) +
  geom_point(
    size = 4,
    position = position_jitter(width = 0.08, height = 0)
  ) +
  facet_wrap(
    ~ program,
    scales = "free_y"
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = "GSE242477: CD8-like TIL functional module scores",
    subtitle = "Descriptive comparison using author-provided clinical response labels",
    x = "Clinical response group",
    y = "Mean module score per patient"
  ) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey90")
  )

ggsave(
  file.path(fig_dir, "gse242477_cd8_module_scores_by_response.png"),
  p1,
  width = 9,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(fig_dir, "gse242477_cd8_module_scores_by_response.pdf"),
  p1,
  width = 9,
  height = 5
)

p2 <- ggplot(
  long_scores,
  aes(
    x = patient_id,
    y = score,
    fill = response_group
  )
) +
  geom_col(width = 0.7) +
  facet_wrap(
    ~ program,
    scales = "free_y"
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = "GSE242477: patient-level CD8-like TIL module scores",
    x = "Patient",
    y = "Mean module score"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey90")
  )

ggsave(
  file.path(fig_dir, "gse242477_cd8_module_scores_by_patient.png"),
  p2,
  width = 9,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(fig_dir, "gse242477_cd8_module_scores_by_patient.pdf"),
  p2,
  width = 9,
  height = 5
)

message("Comparison OK")