# =========================
# 16_add_clinical_response_gse120575.R
# =========================

base_dir <- "/Users/jimmlucas/Downloads/R_MASTER_RNASEQ"
setwd(base_dir)

source(file.path(base_dir, "00_config.R"))

library(tidyverse)
library(data.table)

# -------------------------
# 1. Cargar metadata original GEO
# -------------------------

meta_raw <- fread(
  file.path(res_dir, "GSE120575_patient_ID_single_cells_raw.tsv"),
  data.table = FALSE
) %>%
  as_tibble() %>%
  mutate(across(everything(), as.character))

# -------------------------
# 2. Extraer tabla clínica desde fila 20
# -------------------------

header_row <- 20

header <- meta_raw[header_row, ] %>%
  unlist(use.names = FALSE)

header <- ifelse(
  is.na(header) | header == "",
  paste0("extra_", seq_along(header)),
  header
)

header <- make.names(header, unique = TRUE)

clinical_raw <- meta_raw[(header_row + 1):nrow(meta_raw), ]
colnames(clinical_raw) <- header

clinical_raw <- clinical_raw %>%
  as_tibble() %>%
  filter(!is.na(Sample.name), Sample.name != "")

# -------------------------
# 3. Limpiar variables clínicas
# -------------------------

clinical_clean <- clinical_raw %>%
  transmute(
    cell_id = title,
    sample_label = characteristics..patinet.ID..Pre.baseline..Post..on.treatment.,
    response = characteristics..response,
    therapy = characteristics..therapy,
    patient_alias = str_extract(sample_label, "P[0-9]+"),
    timepoint = case_when(
      str_detect(sample_label, "^Pre") ~ "Baseline",
      str_detect(sample_label, "^Post") ~ "Post",
      TRUE ~ NA_character_
    ),
    response_binary = case_when(
      response == "Responder" ~ "Responder",
      response == "Non-responder" ~ "Non_responder",
      response == "Non_responder" ~ "Non_responder",
      TRUE ~ NA_character_
    )
  )

# -------------------------
# 4. Cargar metadata limpia previa y unir
# -------------------------

cell_meta <- fread(
  file.path(res_dir, "GSE120575_cell_metadata_clean.tsv"),
  data.table = FALSE
) %>%
  as_tibble()

cell_meta_response <- cell_meta %>%
  left_join(
    clinical_clean %>%
      select(cell_id, response, response_binary, therapy),
    by = "cell_id"
  )

# -------------------------
# 5. Comprobaciones
# -------------------------

response_summary <- cell_meta_response %>%
  count(patient_alias, timepoint, response_binary, therapy, name = "n_cells") %>%
  arrange(patient_alias, timepoint)

cohort_response_summary <- cell_meta_response %>%
  summarise(
    dataset = "GSE120575",
    n_cells = n(),
    n_patients = n_distinct(patient_alias),
    n_cells_with_response = sum(!is.na(response_binary)),
    n_responders_cells = sum(response_binary == "Responder", na.rm = TRUE),
    n_non_responders_cells = sum(response_binary == "Non_responder", na.rm = TRUE),
    n_patients_with_response = n_distinct(patient_alias[!is.na(response_binary)])
  )

patient_response_summary <- cell_meta_response %>%
  filter(!is.na(response_binary)) %>%
  distinct(patient_alias, response_binary, therapy) %>%
  count(response_binary, therapy, name = "n_patients")

# -------------------------
# 6. Guardar resultados
# -------------------------

fwrite(
  clinical_clean,
  file.path(res_dir, "GSE120575_clinical_metadata_clean.tsv"),
  sep = "\t"
)

fwrite(
  cell_meta_response,
  file.path(res_dir, "GSE120575_cell_metadata_with_response.tsv"),
  sep = "\t"
)

saveRDS(
  cell_meta_response,
  file.path(proc_dir, "GSE120575_cell_metadata_with_response.rds")
)

fwrite(
  response_summary,
  file.path(res_dir, "GSE120575_patient_timepoint_response_summary.tsv"),
  sep = "\t"
)

fwrite(
  cohort_response_summary,
  file.path(res_dir, "GSE120575_cohort_response_summary.tsv"),
  sep = "\t"
)

fwrite(
  patient_response_summary,
  file.path(res_dir, "GSE120575_patient_response_summary.tsv"),
  sep = "\t"
)

print(cohort_response_summary)
print(patient_response_summary)
print(response_summary, n = Inf)

message("Respuesta clínica añadida a GSE120575.")