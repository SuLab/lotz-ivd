#!/usr/bin/env Rscript
# ============================================================================
# Module 05 — Workflow C: STACAS Integration
#
# Integrates cells across studies into four compartment-based objects (NP, AF,
# CEP, all_cells) using STACAS semi-supervised integration. Uses coarse labels
# from Module 04 (coarse_label / cell_type_coarse) to guide anchor weighting.
#
# STACAS is an anchor-based integration method that uses prior cell type labels
# to weight anchors, reducing overcorrection while preserving biological
# variability. Fully R-native and compatible with Seurat v5.
#
# Usage:
#   Rscript scripts/05c_integration_stacas.R                     # All objects
#   Rscript scripts/05c_integration_stacas.R --object NP         # Single object
#   Rscript scripts/05c_integration_stacas.R --tiered            # Tiered mode
#   Rscript scripts/05c_integration_stacas.R --validate-only     # Validation only
#   Rscript scripts/05c_integration_stacas.R --force             # Re-run
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(STACAS)
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

PROC_DIR    <- file.path(BASE, "data", "processed")
INT_DIR     <- file.path(BASE, "data", "integrated", "stacas")
RESULTS_DIR <- file.path(BASE, "results", "integration")

dir.create(INT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Study assignments per object (same as Workflow A) ─────────────────────

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
parser <- ArgumentParser(description = "Module 05C: STACAS Integration")
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
# HELPER: Load h5ad as Seurat object
# ═══════════════════════════════════════════════════════════════════════════

load_h5ad_as_seurat <- function(h5ad_path) {
  # Use Python bridge to convert h5ad -> mtx + metadata, then load in R
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

  # matrix.mtx is genes x cells (Read10X convention)
  colnames(counts) <- barcodes
  rownames(counts) <- features[[2]]  # gene names
  counts <- as(counts, "dgCMatrix")

  obj <- CreateSeuratObject(counts = counts, project = sub("\\.h5ad$", "", basename(h5ad_path)))

  # Add metadata
  meta <- read.csv(file.path(bridge_dir, "metadata.csv"), row.names = 1,
                   stringsAsFactors = FALSE)
  for (col in colnames(meta)) {
    obj@meta.data[[col]] <- meta[colnames(obj), col]
  }

  return(obj)
}


# ═══════════════════════════════════════════════════════════════════════════
# LOAD AND BUILD OBJECT
# ═══════════════════════════════════════════════════════════════════════════

load_and_build_object <- function(object_name) {
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

    # Ensure study metadata
    if (!("study" %in% colnames(obj@meta.data))) {
      obj@meta.data$study <- acc
    }

    # Ensure coarse_label column exists (check both possible names)
    if (!("coarse_label" %in% colnames(obj@meta.data))) {
      if ("cell_type_coarse" %in% colnames(obj@meta.data)) {
        obj@meta.data$coarse_label <- obj@meta.data$cell_type_coarse
      } else {
        obj@meta.data$coarse_label <- "Unknown"
      }
    }

    message("    ", acc, " (", ifelse(is.null(comp_filter), "all", comp_filter),
            "): ", ncol(obj), " cells")
    seurat_list[[paste0(acc, "_", ifelse(is.null(comp_filter), "all", comp_filter))]] <- obj
  }

  if (length(seurat_list) == 0) {
    message("  WARNING: No objects loaded for ", object_name)
    return(NULL)
  }

  return(seurat_list)
}


# ═══════════════════════════════════════════════════════════════════════════
# STACAS INTEGRATION (FLAT)
# ═══════════════════════════════════════════════════════════════════════════

run_stacas_flat <- function(seurat_list, object_name) {
  n_total <- sum(sapply(seurat_list, ncol))
  message("\n  Running STACAS flat integration for ", object_name,
          " (", n_total, " cells, ", length(seurat_list), " studies)...")

  # For large objects, downsample per study (same approach as CCA v5)
  DOWNSAMPLE_THRESHOLD <- 100000
  CELLS_PER_STUDY <- 2000
  if (n_total > DOWNSAMPLE_THRESHOLD) {
    message("    Downsampling from ", n_total, " to ~",
            CELLS_PER_STUDY * length(seurat_list), " cells...")
    set.seed(42)
    seurat_list <- lapply(seurat_list, function(obj) {
      n_take <- min(CELLS_PER_STUDY, ncol(obj))
      cells <- sample(colnames(obj), n_take)
      subset(obj, cells = cells)
    })
    n_total <- sum(sapply(seurat_list, ncol))
    message("    Downsampled: ", n_total, " cells")
  }

  # Normalize each object with NormalizeData (lighter than SCTransform)
  message("    NormalizeData on each object...")
  seurat_list <- lapply(seurat_list, function(obj) {
    obj <- NormalizeData(obj, verbose = FALSE)
    obj <- FindVariableFeatures(obj, nfeatures = 3000, verbose = FALSE)
    obj
  })

  # Run.STACAS handles anchor finding + integration in one call
  # It uses coarse_label for semi-supervised anchor weighting
  message("    Run.STACAS (dims = 30, anchor.features = 2000)...")
  integrated <- Run.STACAS(
    object.list = seurat_list,
    dims = 30,
    anchor.features = 2000,
    cell.labels = "coarse_label",
    verbose = TRUE
  )

  # Post-integration: PCA on integrated assay, then UMAP
  message("    Running PCA, UMAP, FindNeighbors on integrated assay...")
  DefaultAssay(integrated) <- "integrated"
  integrated <- RunPCA(integrated, npcs = 50, verbose = FALSE)
  integrated <- RunUMAP(integrated, dims = 1:50, verbose = FALSE)
  integrated <- FindNeighbors(integrated, dims = 1:50, verbose = FALSE)

  message("    STACAS flat integration complete: ", ncol(integrated), " cells")
  return(integrated)
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
  out_path <- file.path(RESULTS_DIR, paste0("umap_stacas_", object_name, ".png"))
  ggsave(out_path, plot = p, width = 5 * length(plots), height = 5, dpi = 150)
  message("  Saved: umap_stacas_", object_name, ".png")
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
  message("Processing object: ", object_name, " (STACAS flat)")
  message(paste(rep("=", 60), collapse = ""))

  # Load data
  seurat_list <- load_and_build_object(object_name)
  if (is.null(seurat_list)) return(NULL)

  # Run integration
  result <- run_stacas_flat(seurat_list, object_name)

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

  rm(result, seurat_list)
  gc()

  return(output_path)
}


# ═══════════════════════════════════════════════════════════════════════════
# COMPUTE INTEGRATION METRICS
# ═══════════════════════════════════════════════════════════════════════════

compute_metrics <- function(object_name) {
  rds_path <- file.path(INT_DIR, paste0(object_name, ".rds"))
  if (!file.exists(rds_path)) return(NULL)

  message("  Computing metrics for STACAS/", object_name, "...")
  obj <- readRDS(rds_path)

  metrics <- data.frame(
    object = object_name,
    workflow = "stacas",
    n_cells = ncol(obj),
    stringsAsFactors = FALSE
  )

  # TODO: Compute iLISI using lisi::compute_lisi()
  # TODO: Compute batch-ASW using cluster::silhouette()
  # TODO: Compute condition-ASW

  # Blob check
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
  message("Validating STACAS integration outputs")
  message(paste(rep("=", 60), collapse = ""))

  messages <- character(0)
  all_pass <- TRUE

  for (obj_name in c("NP", "AF", "CEP", "all_cells")) {
    rds_path <- file.path(INT_DIR, paste0(obj_name, ".rds"))
    if (file.exists(rds_path)) {
      messages <- c(messages, paste0("PASS: stacas/", obj_name, " output exists"))

      obj <- readRDS(rds_path)
      if ("umap" %in% names(obj@reductions)) {
        messages <- c(messages, paste0("PASS: stacas/", obj_name, " has UMAP reduction"))
      } else {
        messages <- c(messages, paste0("WARNING: stacas/", obj_name, " missing UMAP reduction"))
      }

      if ("pca" %in% names(obj@reductions)) {
        messages <- c(messages, paste0("PASS: stacas/", obj_name, " has PCA reduction"))
      } else {
        messages <- c(messages, paste0("WARNING: stacas/", obj_name, " missing PCA reduction"))
      }

      rm(obj)
      gc()
    } else {
      messages <- c(messages, paste0("FAIL: stacas/", obj_name, " output missing"))
      all_pass <- FALSE
    }
  }

  for (msg in messages) {
    message("  ", msg)
  }

  message("\n", paste(rep("=", 60), collapse = ""))
  message("STACAS VALIDATION: ", ifelse(all_pass, "ALL CHECKS PASSED", "SOME CHECKS FAILED"))
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
  message("Module 05 — Workflow C: STACAS Integration")
  message("Started: ", Sys.time())
  message(paste(rep("=", 60), collapse = ""))

  # Determine which objects to process
  if (!is.null(args$object)) {
    objects_to_process <- args$object
  } else {
    objects_to_process <- c("NP", "AF", "CEP", "all_cells")
  }

  # Process each object
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

  # Save metrics
  if (length(all_metrics) > 0) {
    metrics_df <- do.call(rbind, all_metrics)
    metrics_path <- file.path(INT_DIR, "integration_metrics.tsv")
    write.table(metrics_df, metrics_path, sep = "\t", row.names = FALSE, quote = FALSE)
    message("  Saved metrics: ", metrics_path)
  }

  # Validation
  result <- validate_integration()

  message("\nCompleted: ", Sys.time())
  message("Overall: ", ifelse(result$passed, "PASSED", "FAILED"))
}

main()
