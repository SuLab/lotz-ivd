#!/usr/bin/env Rscript
# ============================================================================
# Helper: Convert selected workflow's RDS files to h5ad for downstream Python
#
# Converts Seurat RDS objects from CCA or STACAS workflows into AnnData h5ad
# format, placing them at the canonical data/integrated/{object}.h5ad path
# expected by Modules 06-12.
#
# Usage:
#   Rscript scripts/05_select_workflow.R --workflow cca
#   Rscript scripts/05_select_workflow.R --workflow stacas --object NP
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(SeuratDisk)
  library(argparse)
})

# Detect script directory robustly (works with Rscript, nohup, source)
.get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(dirname(sub("^--file=", "", file_arg[1]))))
  }
  if (!is.null(sys.frame(1)$ofile)) {
    return(normalizePath(dirname(sys.frame(1)$ofile)))
  }
  return(normalizePath("."))
}
BASE <- normalizePath(file.path(.get_script_dir(), ".."), mustWork = FALSE)

parser <- ArgumentParser(description = "Convert selected workflow RDS to h5ad")
parser$add_argument("--workflow", required = TRUE,
                    choices = c("cca", "stacas"),
                    help = "Which workflow's RDS files to convert")
parser$add_argument("--object", default = NULL,
                    help = "Single object to convert (default: all)")
args <- parser$parse_args()

WF_DIR <- file.path(BASE, "data", "integrated", args$workflow)
OUT_DIR <- file.path(BASE, "data", "integrated")
OBJECTS <- c("NP", "AF", "CEP", "all_cells")

if (!is.null(args$object)) {
  OBJECTS <- args$object
}

for (obj_name in OBJECTS) {
  rds_path <- file.path(WF_DIR, paste0(obj_name, ".rds"))
  h5ad_path <- file.path(OUT_DIR, paste0(obj_name, ".h5ad"))
  h5seurat_path <- file.path(OUT_DIR, paste0(obj_name, ".h5seurat"))

  if (!file.exists(rds_path)) {
    message("  SKIP: ", rds_path, " not found")
    next
  }

  message("\n  Converting ", args$workflow, "/", obj_name, ".rds -> ", obj_name, ".h5ad ...")
  obj <- readRDS(rds_path)

  # Ensure default assay has data for conversion
  message("    Cells: ", ncol(obj), ", Assays: ", paste(names(obj@assays), collapse = ", "))
  message("    Reductions: ", paste(names(obj@reductions), collapse = ", "))

  # Save as h5seurat (intermediate), then convert to h5ad
  if (file.exists(h5seurat_path)) file.remove(h5seurat_path)
  SaveH5Seurat(obj, filename = h5seurat_path, overwrite = TRUE)
  Convert(h5seurat_path, dest = "h5ad", overwrite = TRUE)

  # Clean up intermediate h5seurat file
  if (file.exists(h5seurat_path)) file.remove(h5seurat_path)

  if (file.exists(h5ad_path)) {
    message("    OK: ", h5ad_path)
  } else {
    message("    ERROR: conversion failed for ", obj_name)
  }

  rm(obj)
  gc()
}

message("\nDone. Converted ", args$workflow, " outputs to h5ad at ", OUT_DIR)
