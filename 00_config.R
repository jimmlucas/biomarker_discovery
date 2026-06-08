# =========================
# 00_config.R
# =========================
# Configuración general del proyecto.
# Este script define rutas y crea carpetas necesarias.
# Se carga al inicio de los demás scripts.

base_dir <- "/Users/jimmlucas/Downloads/R_MASTER_RNASEQ"

raw_dir  <- file.path(base_dir, "data", "raw")
proc_dir <- file.path(base_dir, "data", "processed")
res_dir  <- file.path(base_dir, "results")
fig_dir  <- file.path(base_dir, "figures")

dir.create(raw_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(proc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir,  recursive = TRUE, showWarnings = FALSE)

setwd(base_dir)