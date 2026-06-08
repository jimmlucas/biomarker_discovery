# =========================
# 05_patient_level_features_gse272993.R
# =========================
# Construcción de representaciones a nivel paciente/timepoint.
# Aquí se obtiene la composición de subestados CD8+ por muestra clínica.

source("scripts/00_config.R")

library(tidyverse)
library(data.table)

meta <- readRDS(file.path(proc_dir, "GSE272993_metadata_only.rds"))

# Revisión de categorías clínicas y celulares.
print(meta %>% count(response, sort = TRUE))
print(meta %>% count(treatment, sort = TRUE))
print(meta %>% count(timepoint, sort = TRUE))
print(meta %>% count(metaclusters, sort = TRUE))

# Definición de respuesta binaria para modelado.
# CR/PR se consideran respondedores; PD/SD no respondedores.
meta <- meta %>%
  mutate(
    response_binary = case_when(
      response %in% c("CR", "PR") ~ "Responder",
      response %in% c("PD", "SD") ~ "Non_responder",
      TRUE ~ NA_character_
    )
  )

print(meta %>% count(response, response_binary))

saveRDS(
  meta,
  file.path(proc_dir, "GSE272993_metadata_with_response_binary.rds")
)

fwrite(
  meta,
  file.path(proc_dir, "GSE272993_metadata_with_response_binary.tsv"),
  sep = "\t"
)

# Composición de metaclusters por paciente, tratamiento y timepoint.
comp_patient <- meta %>%
  filter(!is.na(response_binary)) %>%
  count(
    patient_alias,
    response,
    response_binary,
    treatment,
    timepoint,
    metaclusters
  ) %>%
  group_by(patient_alias, response, response_binary, treatment, timepoint) %>%
  mutate(freq = n / sum(n)) %>%
  ungroup()

fwrite(
  comp_patient,
  file.path(res_dir, "GSE272993_patient_timepoint_metacluster_composition_long.tsv"),
  sep = "\t"
)

# Formato ancho: una columna por subestado.
comp_patient_wide <- comp_patient %>%
  select(patient_alias, response, response_binary, treatment, timepoint, metaclusters, freq) %>%
  pivot_wider(
    names_from = metaclusters,
    values_from = freq,
    values_fill = 0,
    names_prefix = "freq_"
  )

fwrite(
  comp_patient_wide,
  file.path(res_dir, "GSE272993_patient_timepoint_metacluster_composition_wide.tsv"),
  sep = "\t"
)

# Figura global de composición.
p_comp <- comp_patient %>%
  ggplot(aes(x = response_binary, y = freq, fill = response_binary)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
  facet_wrap(~ metaclusters, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Respuesta clínica",
    y = "Frecuencia relativa",
    title = "Composición de subestados CD8+ por respuesta clínica"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(fig_dir, "GSE272993_metacluster_composition_by_response_binary.png"),
  p_comp,
  width = 12,
  height = 8,
  dpi = 300
)

# Comparación exploratoria por Wilcoxon.
wilcox_metaclusters <- comp_patient %>%
  group_by(metaclusters) %>%
  summarise(
    p_value = wilcox.test(freq ~ response_binary)$p.value,
    median_responder = median(freq[response_binary == "Responder"], na.rm = TRUE),
    median_non_responder = median(freq[response_binary == "Non_responder"], na.rm = TRUE),
    delta_median = median_responder - median_non_responder,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    direction = case_when(
      delta_median > 0 ~ "Higher_in_responders",
      delta_median < 0 ~ "Lower_in_responders",
      TRUE ~ "No_difference"
    ),
    evidence = case_when(
      p_adj < 0.05 ~ "FDR_significant",
      p_value < 0.05 & p_adj >= 0.05 ~ "Nominal_trend",
      TRUE ~ "No_evidence"
    )
  ) %>%
  arrange(p_adj)

fwrite(
  wilcox_metaclusters,
  file.path(res_dir, "GSE272993_metacluster_response_interpretation.tsv"),
  sep = "\t"
)

print(wilcox_metaclusters)

# Figura con los tres subestados candidatos principales.
top_clusters <- wilcox_metaclusters %>%
  arrange(p_adj) %>%
  slice_head(n = 3) %>%
  pull(metaclusters)

p_top <- comp_patient %>%
  filter(metaclusters %in% top_clusters) %>%
  ggplot(aes(x = response_binary, y = freq, fill = response_binary)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 1, alpha = 0.6) +
  facet_wrap(~ metaclusters, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Respuesta clínica",
    y = "Frecuencia relativa",
    title = "Subestados CD8+ candidatos asociados a respuesta"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(fig_dir, "GSE272993_top_metacluster_candidates_by_response.png"),
  p_top,
  width = 9,
  height = 4,
  dpi = 300
)