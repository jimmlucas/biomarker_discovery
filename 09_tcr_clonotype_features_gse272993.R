# =========================
# 09_tcr_clonotype_features_gse272993.R
# =========================
# Análisis exploratorio de clonotipos TCR.
# Se calculan métricas de diversidad y expansión clonal por paciente/timepoint.
# Este bloque se interpreta con cautela porque algunas métricas pueden estar
# influidas por profundidad de captura o número de células con TCR recuperado.

source("scripts/00_config.R")

library(tidyverse)
library(data.table)

meta <- readRDS(file.path(proc_dir, "GSE272993_metadata_with_response_binary.rds"))

tcr_meta <- meta %>%
  filter(!is.na(response_binary)) %>%
  mutate(
    clonotype_trb = TRB_chain_cdr3,
    clonotype_tra_trb = paste(TRA_chain_cdr3, TRB_chain_cdr3, sep = "_")
  ) %>%
  filter(!is.na(clonotype_trb), clonotype_trb != "")

clono_patient <- tcr_meta %>%
  count(patient_alias, response, response_binary, treatment, timepoint, clonotype_trb) %>%
  group_by(patient_alias, response, response_binary, treatment, timepoint) %>%
  summarise(
    n_cells_tcr = sum(n),
    n_clonotypes = n_distinct(clonotype_trb),
    max_clone_size = max(n),
    top_clone_fraction = max(n) / sum(n),
    shannon_entropy = -sum((n / sum(n)) * log(n / sum(n))),
    clonality = 1 - shannon_entropy / log(n_clonotypes),
    .groups = "drop"
  ) %>%
  mutate(
    clonality = ifelse(is.nan(clonality) | is.infinite(clonality), NA_real_, clonality)
  )

fwrite(
  clono_patient,
  file.path(res_dir, "GSE272993_patient_timepoint_tcr_clonality.tsv"),
  sep = "\t"
)

clono_long <- clono_patient %>%
  pivot_longer(
    cols = c(
      n_cells_tcr,
      n_clonotypes,
      max_clone_size,
      top_clone_fraction,
      shannon_entropy,
      clonality
    ),
    names_to = "metric",
    values_to = "value"
  ) %>%
  filter(!is.na(value), is.finite(value))

wilcox_tcr <- clono_long %>%
  group_by(metric) %>%
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
    ),
    evidence = case_when(
      p_adj < 0.05 ~ "FDR_significant",
      p_value < 0.05 & p_adj >= 0.05 ~ "Nominal_trend",
      TRUE ~ "No_evidence"
    )
  ) %>%
  arrange(p_adj)

fwrite(
  wilcox_tcr,
  file.path(res_dir, "GSE272993_wilcox_tcr_clonality_response.tsv"),
  sep = "\t"
)

p_tcr <- clono_long %>%
  ggplot(aes(x = response_binary, y = value, fill = response_binary)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
  facet_wrap(~ metric, scales = "free_y") +
  theme_classic() +
  labs(
    x = "Respuesta clínica",
    y = "Valor",
    title = "Métricas TCR/clonalidad CD8+ por respuesta"
  ) +
  theme(legend.position = "none")

ggsave(
  file.path(fig_dir, "GSE272993_tcr_clonality_by_response.png"),
  p_tcr,
  width = 12,
  height = 8,
  dpi = 300
)

print(wilcox_tcr)