#!/usr/bin/env Rscript
# ============================================================================
# Export Seurat RDS to flat files for Python h5ad assembly
#
# Exports: counts matrix (MTX), cell metadata (CSV), gene names (CSV),
# and embeddings (CSV) from Seurat objects. These are reassembled into
# h5ad by the companion Python script (seurat_to_h5ad_assemble.py).
#
# Usage:
#   Rscript scripts/seurat_to_h5ad_bridge.R --workflow cca
#   Rscript scripts/seurat_to_h5ad_bridge.R --workflow cca --object NP
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

# Detect script directory
.get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(dirname(sub("^--file=", "", file_arg[1]))))
  }
  return(normalizePath("."))
}

BASE <- normalizePath(file.path(.get_script_dir(), ".."), mustWork = FALSE)

# Parse args
args <- commandArgs(trailingOnly = TRUE)
workflow <- NULL
object_filter <- NULL
i <- 1
while (i <= length(args)) {
  if (args[i] == "--workflow" && i < length(args)) {
    workflow <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--object" && i < length(args)) {
    object_filter <- args[i + 1]; i <- i + 2
  } else { i <- i + 1 }
}
if (is.null(workflow)) stop("--workflow required (cca or stacas)")

WF_DIR <- file.path(BASE, "data", "integrated", workflow)
OBJECTS <- c("NP", "AF", "CEP", "all_cells")
if (!is.null(object_filter)) OBJECTS <- object_filter

for (obj_name in OBJECTS) {
  rds_path <- file.path(WF_DIR, paste0(obj_name, ".rds"))
  if (!file.exists(rds_path)) {
    message("  SKIP: ", rds_path, " not found")
    next
  }

  export_dir <- file.path(WF_DIR, "bridge_export", obj_name)
  dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)

  message("\n========================================")
  message("Exporting ", workflow, "/", obj_name)
  message("========================================")

  obj <- readRDS(rds_path)
  message("  Cells: ", ncol(obj), ", Genes: ", nrow(obj))
  message("  Assays: ", paste(names(obj@assays), collapse = ", "))
  message("  Reductions: ", paste(names(obj@reductions), collapse = ", "))

  # 1. Export counts matrix (sparse MTX)
  # For CCA with NormalizeData, the RNA assay has layers: counts and data
  assay <- DefaultAssay(obj)
  message("  Default assay: ", assay)

  # Get raw counts
  counts <- NULL
  tryCatch({
    counts <- LayerData(obj, assay = assay, layer = "counts")
    message("  Got counts layer: ", nrow(counts), " x ", ncol(counts))
  }, error = function(e) {
    message("  counts layer not found, trying GetAssayData...")
  })

  if (is.null(counts)) {
    tryCatch({
      counts <- GetAssayData(obj, layer = "counts")
      message("  Got counts via GetAssayData: ", nrow(counts), " x ", ncol(counts))
    }, error = function(e) {
      message("  WARNING: Could not get counts, using data layer")
      counts <<- LayerData(obj, assay = assay, layer = "data")
    })
  }

  # Seurat stores genes x cells; writeMM writes as-is
  # Python side will transpose (cells x genes)
  mtx_path <- file.path(export_dir, "counts.mtx.gz")
  message("  Writing counts matrix...")
  writeMM(counts, file.path(export_dir, "counts.mtx"))
  system2("gzip", c("-f", file.path(export_dir, "counts.mtx")))
  message("  Saved: ", mtx_path)

  # Skip normalized data export — recompute in Python (much faster than exporting dense matrix)
  # Seurat NormalizeData(method="LogNormalize", scale.factor=10000) == sc.pp.normalize_total + sc.pp.log1p
  message("  Skipping normalized data export (will recompute in Python)")

  # 2. Export gene names
  genes <- rownames(counts)
  write.csv(data.frame(gene = genes), file.path(export_dir, "genes.csv"),
            row.names = FALSE)
  message("  Saved ", length(genes), " gene names")

  # 3. Export cell barcodes
  barcodes <- colnames(obj)
  write.csv(data.frame(barcode = barcodes), file.path(export_dir, "barcodes.csv"),
            row.names = FALSE)

  # 4. Export cell metadata
  meta <- obj@meta.data
  meta$barcode <- rownames(meta)
  write.csv(meta, gzfile(file.path(export_dir, "metadata.csv.gz")),
            row.names = FALSE)
  message("  Saved metadata: ", ncol(meta), " columns, ", nrow(meta), " cells")

  # 5. Export embeddings
  for (red_name in names(obj@reductions)) {
    emb <- Embeddings(obj, reduction = red_name)
    n_dims <- ncol(emb)
    out_name <- paste0("embedding_", red_name, ".csv.gz")
    write.csv(emb, gzfile(file.path(export_dir, out_name)), row.names = TRUE)
    message("  Saved embedding: ", red_name, " (", n_dims, " dims)")
  }

  rm(obj, counts)
  gc(verbose = FALSE)
  message("  Export complete for ", obj_name)
}

message("\n\nBridge export complete. Run seurat_to_h5ad_assemble.py next.")
