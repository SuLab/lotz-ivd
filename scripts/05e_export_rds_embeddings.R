#!/usr/bin/env Rscript
# ============================================================================
# Export CCA/STACAS embeddings + metadata from RDS for Python metric computation
#
# Extracts integrated PCA embedding and metadata (study, condition_harmonized)
# from each RDS file and saves as CSV for the Python comparison script.
#
# Usage:
#   Rscript scripts/05e_export_rds_embeddings.R
#   Rscript scripts/05e_export_rds_embeddings.R --workflow cca
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
})

args <- commandArgs(trailingOnly = TRUE)

BASE_DIR <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
INT_DIR <- file.path(BASE_DIR, "data", "integrated")

workflows <- c("cca", "stacas")
if ("--workflow" %in% args) {
  wf_idx <- which(args == "--workflow") + 1
  if (wf_idx <= length(args)) {
    workflows <- args[wf_idx]
  }
}

objects <- c("NP", "AF", "CEP", "all_cells")

for (wf in workflows) {
  wf_dir <- file.path(INT_DIR, wf)
  export_dir <- file.path(wf_dir, "embeddings_export")
  dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)

  for (obj_name in objects) {
    rds_path <- file.path(wf_dir, paste0(obj_name, ".rds"))
    if (!file.exists(rds_path)) {
      message("  SKIP: ", rds_path, " not found")
      next
    }

    message("  Loading ", wf, "/", obj_name, "...")
    obj <- readRDS(rds_path)
    message("    Cells: ", ncol(obj))

    # Find the integrated reduction
    reductions <- names(obj@reductions)
    message("    Available reductions: ", paste(reductions, collapse = ", "))

    # For CCA: integrated.cca or integrated.dr; for STACAS: pca
    emb_key <- NULL
    for (key in c("integrated.cca", "integrated.dr", "pca")) {
      if (key %in% reductions) {
        emb_key <- key
        break
      }
    }

    if (is.null(emb_key)) {
      message("    WARNING: No suitable embedding found, skipping")
      next
    }

    message("    Using embedding: ", emb_key)
    emb <- Embeddings(obj, reduction = emb_key)

    # Export embedding (first 50 dims or all available)
    n_dims <- min(50, ncol(emb))
    emb_df <- as.data.frame(emb[, 1:n_dims])

    # Add metadata
    meta <- obj@meta.data
    emb_df$study <- if ("study" %in% colnames(meta)) meta$study else NA
    emb_df$condition_harmonized <- if ("condition_harmonized" %in% colnames(meta)) meta$condition_harmonized else NA
    emb_df$cell_barcode <- colnames(obj)

    out_path <- file.path(export_dir, paste0(obj_name, "_embedding.csv.gz"))
    message("    Saving ", nrow(emb_df), " cells x ", n_dims, " dims to ", out_path)
    write.csv(emb_df, gzfile(out_path), row.names = FALSE)
    message("    Done: ", obj_name)

    rm(obj, emb, emb_df, meta)
    gc(verbose = FALSE)
  }
}

message("\nEmbedding export complete.")
