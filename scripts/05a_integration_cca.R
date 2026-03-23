#!/usr/bin/env Rscript
# ============================================================================
# Module 05 — Workflow A: CCA Integration (PRIMARY)
#
# Integrates cells across studies into four compartment-based objects (NP, AF,
# CEP, all_cells) using Seurat v5 CCA (Canonical Correlation Analysis).
#
# CCA is label-free — finds shared correlation structure without requiring
# prior cell type annotations. This workflow follows the integration approach
# from single_nuclei_r/Sample_QC_Integration.R, updated to use Seurat v5's
# IntegrateLayers API for scalability on large datasets.
#
# Seurat v5 optimizations:
#   - IntegrateLayers(method = CCAIntegration) replaces the v4 pipeline of
#     FindIntegrationAnchors → IntegrateData
#   - Merged object with per-study layers (v5 data model) instead of a list
#     of separate objects
#
# Standard (non-sketch) CCA is used for all objects regardless of size.
#
# Usage:
#   Rscript scripts/05a_integration_cca.R                     # All objects
#   Rscript scripts/05a_integration_cca.R --object NP         # Single object
#   Rscript scripts/05a_integration_cca.R --validate-only      # Validation only
#   Rscript scripts/05a_integration_cca.R --force              # Re-run
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(argparse)
  library(dplyr)
  library(ggplot2)
})

# Allow large objects in future (SCTransform uses future for parallelism)
options(future.globals.maxSize = 16 * 1024^3)  # 16 GB

# ── Paths ──────────────────────────────────────────────────────────────────
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

PROC_DIR   <- file.path(BASE, "data", "processed")
INT_DIR    <- file.path(BASE, "data", "integrated", "cca")
RESULTS_DIR <- file.path(BASE, "results", "integration")

dir.create(INT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Integration parameters ───────────────────────────────────────────────
N_HVG       <- 3000   # Number of highly variable genes
N_DIMS      <- 50     # PCA/CCA dimensions

# ── Study assignments per object ──────────────────────────────────────────
NP_STUDIES <- list(
  list(acc = "GSE160756", comp = "NP"),
  list(acc = "GSE165722", comp = NULL),
  list(acc = "GSE199866", comp = "NP"),
  list(acc = "GSE205535", comp = NULL),
  list(acc = "GSE244889", comp = NULL),
  list(acc = "GSE251686", comp = NULL),
  list(acc = "GSE230809", comp = "NP"),
  list(acc = "CNP0002664", comp = NULL)
)

AF_STUDIES <- list(
  list(acc = "GSE160756", comp = "AF"),
  list(acc = "GSE199866", comp = "AF"),
  list(acc = "GSE230809", comp = "AF")
)

CEP_STUDIES <- list(
  list(acc = "GSE160756", comp = "CEP"),
  list(acc = "GSE255768", comp = NULL),
  list(acc = "GSE242443", comp = NULL)
)

ALL_CELLS_EXTRA <- list(
  list(acc = "GSE189916", comp = NULL)
)

EXCLUDED_SAMPLES <- c("GSE251686_NP3")

# ── Argument parsing ──────────────────────────────────────────────────────
parser <- ArgumentParser(description = "Module 05A: CCA Integration (Seurat v5)")
parser$add_argument("--object", type = "character", default = NULL,
                    choices = c("NP", "AF", "CEP", "all_cells"),
                    help = "Process a single object (default: all)")
parser$add_argument("--validate-only", action = "store_true", default = FALSE,
                    help = "Run validation checks only")
parser$add_argument("--force", action = "store_true", default = FALSE,
                    help = "Re-run even if outputs exist")
args <- parser$parse_args()


# ═══════════════════════════════════════════════════════════════════════════
# HELPER: Get study assignments
# ═══════════════════════════════════════════════════════════════════════════

get_study_assignments <- function(object_name) {
  if (object_name == "NP") return(NP_STUDIES)
  if (object_name == "AF") return(AF_STUDIES)
  if (object_name == "CEP") return(CEP_STUDIES)
  if (object_name == "all_cells") {
    all_studies <- c(NP_STUDIES, AF_STUDIES, CEP_STUDIES, ALL_CELLS_EXTRA)
    seen <- character(0)
    unique_studies <- list()
    for (s in all_studies) {
      key <- paste0(s$acc, "_", ifelse(is.null(s$comp), "ALL", s$comp))
      if (!(key %in% seen)) {
        seen <- c(seen, key)
        unique_studies <- c(unique_studies, list(s))
      }
    }
    return(unique_studies)
  }
  stop(paste("Unknown object:", object_name))
}


# ═══════════════════════════════════════════════════════════════════════════
# HELPER: Load h5ad as Seurat object via Python bridge
# ═══════════════════════════════════════════════════════════════════════════

load_h5ad_as_seurat <- function(h5ad_path) {
  bridge_script <- file.path(BASE, "scripts", "h5ad_to_seurat_bridge.py")
  bridge_dir <- file.path(tempdir(), sub("\\.h5ad$", "", basename(h5ad_path)))

  if (!dir.exists(bridge_dir) || !file.exists(file.path(bridge_dir, "matrix.mtx"))) {
    message("  Converting ", basename(h5ad_path), " via bridge...")
    venv_python <- file.path(BASE, ".venv", "bin", "python3")
    ret <- system2(venv_python, args = c(bridge_script, h5ad_path, bridge_dir),
                   stdout = TRUE, stderr = TRUE)
    if (!file.exists(file.path(bridge_dir, "matrix.mtx"))) {
      stop("Bridge conversion failed: ", paste(ret, collapse = "\n"))
    }
  }

  message("  Loading ", basename(h5ad_path), " from bridge files...")
  counts <- Matrix::readMM(file.path(bridge_dir, "matrix.mtx"))
  barcodes <- read.table(file.path(bridge_dir, "barcodes.tsv"), stringsAsFactors = FALSE)[[1]]
  features <- read.table(file.path(bridge_dir, "features.tsv"), sep = "\t", stringsAsFactors = FALSE)

  colnames(counts) <- barcodes
  rownames(counts) <- features[[2]]
  counts <- as(counts, "CsparseMatrix")

  obj <- CreateSeuratObject(counts = counts, project = sub("\\.h5ad$", "", basename(h5ad_path)))

  meta <- read.csv(file.path(bridge_dir, "metadata.csv"), row.names = 1,
                   stringsAsFactors = FALSE)
  for (col in colnames(meta)) {
    obj@meta.data[[col]] <- meta[colnames(obj), col]
  }

  # Clean up bridge files immediately to save disk space
  unlink(bridge_dir, recursive = TRUE)

  return(obj)
}


# ═══════════════════════════════════════════════════════════════════════════
# LOAD, FILTER, AND MERGE INTO A SINGLE V5 OBJECT
# ═══════════════════════════════════════════════════════════════════════════

load_and_merge <- function(object_name) {
  assignments <- get_study_assignments(object_name)

  message("\n  Loading cells for ", object_name, " object...")
  message("    Studies: ", paste(sapply(assignments, `[[`, "acc"), collapse = ", "))

  seurat_list <- list()

  for (s in assignments) {
    acc <- s$acc
    comp_filter <- s$comp
    h5ad_path <- file.path(PROC_DIR, paste0(acc, ".h5ad"))

    if (!file.exists(h5ad_path)) {
      message("    WARNING: ", h5ad_path, " not found, skipping")
      next
    }

    obj <- tryCatch(
      load_h5ad_as_seurat(h5ad_path),
      error = function(e) {
        message("    WARNING: Failed to load ", acc, ": ", e$message)
        return(NULL)
      }
    )
    if (is.null(obj)) next

    # Apply compartment filter
    if (!is.null(comp_filter) && "compartment" %in% colnames(obj@meta.data)) {
      cells_keep <- which(grepl(comp_filter, obj@meta.data$compartment, ignore.case = TRUE))
      if (length(cells_keep) == 0) {
        message("    WARNING: No cells match compartment filter '", comp_filter, "' for ", acc)
        next
      }
      obj <- subset(obj, cells = colnames(obj)[cells_keep])
    }

    # Apply sample exclusions
    if ("sample_id" %in% colnames(obj@meta.data)) {
      cells_keep <- which(!(obj@meta.data$sample_id %in% EXCLUDED_SAMPLES))
      if (length(cells_keep) < ncol(obj)) {
        obj <- subset(obj, cells = colnames(obj)[cells_keep])
      }
    }

    # Filter to cells with valid cell_class
    if ("cell_class" %in% colnames(obj@meta.data)) {
      valid_classes <- c("mesenchymal", "non_mesenchymal", "unknown")
      cells_keep <- which(obj@meta.data$cell_class %in% valid_classes)
      obj <- subset(obj, cells = colnames(obj)[cells_keep])
    } else {
      message("    WARNING: ", acc, " missing cell_class column, skipping")
      next
    }

    if (ncol(obj) == 0) {
      message("    WARNING: No cells after filtering for ", acc)
      next
    }

    if (!("study" %in% colnames(obj@meta.data))) {
      obj@meta.data$study <- acc
    }

    message("    ", acc, " (", ifelse(is.null(comp_filter), "all", comp_filter),
            "): ", ncol(obj), " cells")
    seurat_list[[paste0(acc, "_", ifelse(is.null(comp_filter), "all", comp_filter))]] <- obj
  }

  if (length(seurat_list) == 0) {
    message("  WARNING: No objects loaded for ", object_name)
    return(NULL)
  }

  # Merge into a single Seurat v5 object
  # In v5, each original object's data becomes a separate layer
  message("    Merging ", length(seurat_list), " objects...")
  if (length(seurat_list) == 1) {
    merged <- seurat_list[[1]]
  } else {
    merged <- merge(seurat_list[[1]], y = seurat_list[-1])
  }
  message("    Merged: ", ncol(merged), " cells, ", nrow(merged), " genes")
  message("    Layers: ", paste(Layers(merged), collapse = ", "))

  rm(seurat_list)
  gc()

  return(merged)
}


# ═══════════════════════════════════════════════════════════════════════════
# CCA INTEGRATION (Seurat v5 — standard CCA for all objects)
# ═══════════════════════════════════════════════════════════════════════════

integrate_cca_standard <- function(obj, object_name) {
  message("  Running standard CCA integration (IntegrateLayers)...")

  # Normalize with NormalizeData (keeps layers split, required for IntegrateLayers)
  # Note: Seurat v5 IntegrateLayers requires split layers. SCTransform joins layers,
  # so we use log-normalization here. The CCA algorithm itself is the same regardless
  # of normalization — it finds shared correlation structure across datasets.
  message("    NormalizeData + FindVariableFeatures + ScaleData (per-layer)...")
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = N_HVG, verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)

  message("    RunPCA (", N_DIMS, " dims)...")
  obj <- RunPCA(obj, npcs = N_DIMS, verbose = FALSE)

  # CCA integration across layers (studies)
  message("    IntegrateLayers (CCA, dims = 1:", N_DIMS, ")...")
  obj <- IntegrateLayers(
    object = obj,
    method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.cca",
    dims = 1:N_DIMS,
    verbose = FALSE
  )

  # Join layers back after integration (required for downstream operations)
  obj <- JoinLayers(obj)

  # Post-integration dimensionality reduction
  message("    UMAP + neighbors on integrated reduction...")
  obj <- RunUMAP(obj, reduction = "integrated.cca", dims = 1:N_DIMS, verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "integrated.cca", dims = 1:N_DIMS, verbose = FALSE)

  message("    Standard CCA complete: ", ncol(obj), " cells")
  return(obj)
}


# ═══════════════════════════════════════════════════════════════════════════
# GENERATE UMAP PLOTS
# ═══════════════════════════════════════════════════════════════════════════

plot_umaps <- function(obj, object_name) {
  color_by <- intersect(c("study", "condition_harmonized", "cell_class", "coarse_label"),
                        colnames(obj@meta.data))
  if (length(color_by) == 0) return(invisible(NULL))

  plots <- lapply(color_by, function(col) {
    DimPlot(obj, group.by = col, reduction = "umap", pt.size = 0.1) +
      ggtitle(paste(object_name, "-", col)) +
      theme(legend.position = "right",
            legend.text = element_text(size = 7))
  })

  p <- patchwork::wrap_plots(plots, ncol = length(plots))
  out_path <- file.path(RESULTS_DIR, paste0("umap_cca_", object_name, ".png"))
  ggsave(out_path, plot = p, width = 5 * length(plots), height = 5, dpi = 150)
  message("  Saved: umap_cca_", object_name, ".png")
}


# ═══════════════════════════════════════════════════════════════════════════
# PROCESS ONE OBJECT
# ═══════════════════════════════════════════════════════════════════════════

process_object <- function(object_name, force = FALSE) {
  output_path <- file.path(INT_DIR, paste0(object_name, ".rds"))

  if (file.exists(output_path) && !force) {
    message("\n=== ", object_name, ": output exists, skipping (use --force) ===")
    return(output_path)
  }

  message("\n", paste(rep("=", 60), collapse = ""))
  message("Processing object: ", object_name)
  message(paste(rep("=", 60), collapse = ""))

  # Load and merge all studies into a single v5 object
  merged <- load_and_merge(object_name)
  if (is.null(merged)) return(NULL)

  message("  Object has ", ncol(merged), " cells — running standard CCA")
  result <- integrate_cca_standard(merged, object_name)

  rm(merged)
  gc()

  if (is.null(result)) return(NULL)

  # Generate UMAPs
  tryCatch(
    plot_umaps(result, object_name),
    error = function(e) message("  WARNING: UMAP plotting failed: ", e$message)
  )

  # Save
  message("  Saving to ", output_path, "...")
  saveRDS(result, output_path)
  message("  ", object_name, " complete: ", ncol(result), " cells")

  rm(result)
  gc()

  return(output_path)
}


# ═══════════════════════════════════════════════════════════════════════════
# COMPUTE INTEGRATION METRICS
# ═══════════════════════════════════════════════════════════════════════════

compute_metrics <- function(object_name) {
  rds_path <- file.path(INT_DIR, paste0(object_name, ".rds"))
  if (!file.exists(rds_path)) return(NULL)

  message("  Computing metrics for CCA/", object_name, "...")
  obj <- readRDS(rds_path)

  metrics <- data.frame(
    object = object_name,
    workflow = "cca",
    n_cells = ncol(obj),
    stringsAsFactors = FALSE
  )

  # Blob check: number of clusters at resolution 0.5
  tryCatch({
    obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
    n_clusters <- length(unique(obj@meta.data$seurat_clusters))
    metrics$n_clusters_res05 <- n_clusters
    message("    Clusters at res=0.5: ", n_clusters)
  }, error = function(e) {
    message("    WARNING: Clustering failed: ", e$message)
    metrics$n_clusters_res05 <<- NA
  })

  rm(obj)
  gc()

  return(metrics)
}


# ═══════════════════════════════════════════════════════════════════════════
# VALIDATION
# ═══════════════════════════════════════════════════════════════════════════

validate_integration <- function() {
  message("\n", paste(rep("=", 60), collapse = ""))
  message("Validating CCA integration outputs")
  message(paste(rep("=", 60), collapse = ""))

  messages <- character(0)
  all_pass <- TRUE

  for (obj_name in c("NP", "AF", "CEP", "all_cells")) {
    rds_path <- file.path(INT_DIR, paste0(obj_name, ".rds"))
    if (file.exists(rds_path)) {
      messages <- c(messages, paste0("PASS: cca/", obj_name, " output exists"))

      obj <- readRDS(rds_path)

      if ("umap" %in% names(obj@reductions)) {
        messages <- c(messages, paste0("PASS: cca/", obj_name, " has UMAP reduction"))
      } else {
        messages <- c(messages, paste0("WARNING: cca/", obj_name, " missing UMAP reduction"))
      }

      # Check for CCA integration reduction
      has_cca <- any(grepl("integrated.cca", names(obj@reductions)))
      if (has_cca) {
        messages <- c(messages, paste0("PASS: cca/", obj_name, " has CCA integration reduction"))
      } else if ("pca" %in% names(obj@reductions)) {
        messages <- c(messages, paste0("PASS: cca/", obj_name, " has PCA reduction"))
      } else {
        messages <- c(messages, paste0("WARNING: cca/", obj_name, " missing integration reduction"))
      }

      rm(obj)
      gc()
    } else {
      messages <- c(messages, paste0("FAIL: cca/", obj_name, " output missing"))
      all_pass <- FALSE
    }
  }

  for (msg in messages) {
    message("  ", msg)
  }

  message("\n", paste(rep("=", 60), collapse = ""))
  message("CCA VALIDATION: ", ifelse(all_pass, "ALL CHECKS PASSED", "SOME CHECKS FAILED"))
  message(paste(rep("=", 60), collapse = ""))

  return(list(passed = all_pass, messages = messages))
}


# ═══════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════

main <- function() {
  if (args$validate_only) {
    result <- validate_integration()
    quit(status = ifelse(result$passed, 0, 1))
  }

  message(paste(rep("=", 60), collapse = ""))
  message("Module 05 — Workflow A: CCA Integration (Seurat v5)")
  message("Seurat version: ", as.character(packageVersion("Seurat")))
  message("HVGs: ", N_HVG, ", Dims: ", N_DIMS)
  message("Integration: standard CCA for all objects (no sketch)")
  message("Started: ", Sys.time())
  message(paste(rep("=", 60), collapse = ""))

  if (!is.null(args$object)) {
    objects_to_process <- args$object
  } else {
    objects_to_process <- c("NP", "AF", "CEP", "all_cells")
  }

  all_metrics <- list()
  for (obj_name in objects_to_process) {
    output <- process_object(obj_name, force = args$force)
    if (!is.null(output)) {
      metrics <- compute_metrics(obj_name)
      if (!is.null(metrics)) {
        all_metrics[[obj_name]] <- metrics
      }
    }
  }

  if (length(all_metrics) > 0) {
    metrics_df <- do.call(rbind, all_metrics)
    metrics_path <- file.path(INT_DIR, "integration_metrics.tsv")
    write.table(metrics_df, metrics_path, sep = "\t", row.names = FALSE, quote = FALSE)
    message("  Saved metrics: ", metrics_path)
  }

  result <- validate_integration()

  message("\nCompleted: ", Sys.time())
  message("Overall: ", ifelse(result$passed, "PASSED", "FAILED"))
}

main()
