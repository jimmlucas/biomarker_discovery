# =========================
# 18_extract_marker_expression_gse120575.R
# =========================

base_dir <- "/Users/jimmlucas/Downloads/R_MASTER_RNASEQ"
setwd(base_dir)

source(file.path(base_dir, "00_config.R"))

library(tidyverse)
library(data.table)

gse_dir <- file.path(raw_dir, "GSE120575")

expr_file <- file.path(
  gse_dir,
  "GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz"
)

stopifnot(file.exists(expr_file))

# -------------------------
# 1. Definir genes de interés (CD8 states)
# -------------------------

marker_genes <- c(
  # Exhaustion
  "PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX",
  
  # Effector
  "GZMB", "PRF1", "IFNG",
  
  # Memory
  "IL7R", "CCR7",
  
  # Activation
  "CD69", "HLA-DRA",
  
  # Proliferation
  "MKI67"
)

marker_genes <- unique(marker_genes)

# -------------------------
# 2. Leer cabecera (células)
# -------------------------

con <- gzfile(expr_file, open = "rt")

header_line <- readLines(con, n = 1)
sample_line <- readLines(con, n = 1)

cell_ids <- strsplit(header_line, "\t")[[1]][-1]

message("Número de células detectadas: ", length(cell_ids))

# -------------------------
# 3. Iterar por genes sin cargar todo
# -------------------------

chunk_size <- 500
marker_data_list <- list()

i <- 0

repeat {
  
  lines <- readLines(con, n = chunk_size)
  
  if (length(lines) == 0) break
  
  for (line in lines) {
    
    parts <- strsplit(line, "\t")[[1]]
    gene <- parts[1]
    
    if (gene %in% marker_genes) {
      
      expr_values <- as.numeric(parts[-1])
      
      df <- tibble(
        gene = gene,
        cell_id = cell_ids,
        expression = expr_values
      )
      
      marker_data_list[[length(marker_data_list) + 1]] <- df
      
      message("Encontrado gen: ", gene)
    }
  }
  
  i <- i + 1
  if (i %% 20 == 0) message("Chunks procesados: ", i)
}

close(con)

# -------------------------
# 4. Combinar resultados
# -------------------------

marker_expr <- bind_rows(marker_data_list)

# -------------------------
# 5. Guardar
# -------------------------

fwrite(
  marker_expr,
  file.path(res_dir, "GSE120575_marker_expression_long.tsv"),
  sep = "\t"
)

saveRDS(
  marker_expr,
  file.path(proc_dir, "GSE120575_marker_expression_long.rds")
)

# -------------------------
# 6. Chequeo rápido
# -------------------------

summary_expr <- marker_expr %>%
  group_by(gene) %>%
  summarise(
    mean_expr = mean(expression, na.rm = TRUE),
    pct_nonzero = mean(expression > 0, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_expr)

fwrite(
  summary_expr,
  file.path(res_dir, "GSE120575_marker_expression_summary.tsv"),
  sep = "\t"
)

message("Extracción de genes marcadores completada.")