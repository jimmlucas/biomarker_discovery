#28_score_cd8_programs_gse242477

base_dir <- "/Users/jimmlucas/Downloads/R_MASTER_RNASEQ"
setwd(base_dir)

source(file.path(base_dir, "00_config.R"))

library(tidyverse)
library(data.table)
library(Seurat)
library(ggplot2)

in_file <- file.path(
  base_dir,
  "data/processed/GSE242477/gse242477_clean.rds"
)

out_all_file <- file.path(
  base_dir,
  "data/processed/GSE242477/gse242477_scored_all_cells.rds"
)

out_cd8_file <- file.path(
  base_dir,
  "data/processed/GSE242477/gse242477_til_cd8_scored.rds"
)

res_dir <- file.path(base_dir, "results/GSE242477")

obj <- readRDS(in_file)

DefaultAssay(obj) <- "RNA"

# ------------------------------------------------------------
# 1. Mantener TIL para análisis principal
# ------------------------------------------------------------

obj <- subset(
  obj,
  subset = compartment == "TIL"
)

# ------------------------------------------------------------
# 2. Excluir MEL10 porque no tiene inmunoterapia evaluable
# ------------------------------------------------------------

obj <- subset(
  obj,
  subset = response_group != "Excluded"
)

# ------------------------------------------------------------
# 3. Normalizar si no está normalizado
# ------------------------------------------------------------

has_data <- FALSE

try({
  data_mat <- GetAssayData(
    obj,
    assay = "RNA",
    slot = "data"
  )
  
  has_data <- nrow(data_mat) > 0 && ncol(data_mat) > 0
}, silent = TRUE)

if (!has_data) {
  message("Normalizing RNA data...")
  obj <- NormalizeData(obj)
} else {
  message("RNA data slot already present.")
}

# ------------------------------------------------------------
# 4. Definir programas funcionales CD8+
# Mantener genes simples y comparables entre cohortes
# ------------------------------------------------------------

programs <- list(
  Cytotoxicity = c(
    "GZMB", "GZMA", "GZMH", "PRF1", "NKG7", "GNLY", "CTSW"
  ),
  IFN = c(
    "IFIT1", "IFIT2", "IFIT3", "ISG15", "MX1", "MX2", "STAT1", "IRF7"
  ),
  Activation = c(
    "CD69", "IL2RA", "TNFRSF9", "ICOS", "CD40LG", "HLA-DRA", "HLA-DRB1"
  ),
  Exhaustion = c(
    "PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "CTLA4", "ENTPD1"
  ),
  Proliferation = c(
    "MKI67", "TOP2A", "STMN1", "TYMS", "PCNA", "UBE2C"
  )
)

programs_filtered <- lapply(
  programs,
  function(g) intersect(g, rownames(obj))
)

program_genes_used <- data.frame(
  program = names(programs_filtered),
  n_genes_available = sapply(programs_filtered, length),
  genes_used = sapply(programs_filtered, paste, collapse = ";")
)

write.csv(
  program_genes_used,
  file.path(res_dir, "gse242477_program_genes_used.csv"),
  row.names = FALSE
)

print(program_genes_used)

# ------------------------------------------------------------
# 5. Calcular module scores
# ------------------------------------------------------------

for (nm in names(programs_filtered)) {
  
  genes <- programs_filtered[[nm]]
  
  if (length(genes) >= 3) {
    
    message("Scoring program: ", nm)
    
    obj <- AddModuleScore(
      object = obj,
      features = list(genes),
      name = paste0(nm, "_score")
    )
    
  } else {
    
    warning("Skipping ", nm, ": fewer than 3 genes available.")
  }
}

# ------------------------------------------------------------
# 6. Definir CD8-like por marcador
# Dataset TIL enriquecido, pero se restringe a señal CD3/CD8
# ------------------------------------------------------------

cd8_marker_genes <- intersect(
  c("CD3D", "CD3E", "CD8A", "CD8B"),
  rownames(obj)
)

if (length(cd8_marker_genes) < 3) {
  stop("Not enough CD8 marker genes available.")
}

obj <- AddModuleScore(
  object = obj,
  features = list(cd8_marker_genes),
  name = "CD8_marker_score"
)

# Umbral conservador: top 50% dentro de TIL tratados/evaluables
# Evita perder demasiadas células en un dataset ya enriquecido en T cells.
threshold <- quantile(
  obj$CD8_marker_score1,
  0.50,
  na.rm = TRUE
)

obj$cd8_status <- ifelse(
  obj$CD8_marker_score1 >= threshold,
  "CD8_like",
  "Other_TIL"
)

obj_cd8 <- subset(
  obj,
  subset = cd8_status == "CD8_like"
)

# ------------------------------------------------------------
# 7. Guardar objetos y metadata
# ------------------------------------------------------------

saveRDS(
  obj,
  out_all_file
)

saveRDS(
  obj_cd8,
  out_cd8_file
)

write.csv(
  obj@meta.data,
  file.path(res_dir, "gse242477_til_all_cells_metadata_with_scores.csv")
)

write.csv(
  obj_cd8@meta.data,
  file.path(res_dir, "gse242477_til_cd8_metadata_with_scores.csv")
)

cell_summary <- obj@meta.data %>%
  group_by(patient_id, response_group, cd8_status) %>%
  summarise(
    n_cells = n(),
    .groups = "drop"
  )

write.csv(
  cell_summary,
  file.path(res_dir, "gse242477_cd8_like_cell_summary.csv"),
  row.names = FALSE
)

print(cell_summary)

message("Scoring OK")
message("Saved CD8-like object: ", out_cd8_file)