#!/usr/bin/env Rscript
# ============================================================================
# Module 05 — Workflow A: CCA Integration (PRIMARY)
#
# Integrates cells across studies into four compartment-based objects (NP, AF,
# CEP, all_cells) using Seurat v5 CCA (Canonical Correlation Analysis).
# CCA is label-free — finds shared correlation structure without requiring
# prior cell type annotations.
#
# Usage:
#   Rscript scripts/05a_integration_cca.R                     # All objects
#   Rscript scripts/05a_integration_cca.R --object NP         # Single object
#   Rscript scripts/05a_integration_cca.R --tiered             # Tiered mode
#   Rscript scripts/05a_integration_cca.R --validate-only      # Validation only
#   Rscript scripts/05a_integration_cca.R --force              # Re-run
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(SeuratDisk)
  library(argparse)
  library(dplyr)
  library(ggplot2)
})

# ── Paths ──────────────────────────────────────────────────────────────────
BASE <- normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), ".."),
                      mustWork = FALSE)
# Fallback for non-interactive invocation
if (!dir.exists(BASE)) {
  BASE <- normalizePath(file.path(getwd(), ".."), mustWork = FALSE)
}

PROC_DIR   <- file.path(BASE, "data", "processed")
INT_DIR    <- file.path(BASE, "data", "integrated", "cca")
RESULTS_DIR <- file.path(BASE, "results", "integration")

dir.create(INT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Study assignments per object ──────────────────────────────────────────
# Each entry: list(accession, compartment_filter)
# compartment_filter is NULL for single-compartment studies

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
parser <- ArgumentParser(description = "Module 05A: CCA Integration")
parser$add_argument("--object", type = "character", default = NULL,
                    choices = c("NP", "AF", "CEP", "all_cells"),
                    help = "Process a single object (default: all)")
parser$add_argument("--tiered", action = "store_true", default = FALSE,
                    help = "Use tiered integration (mesenchymal/non-mesenchymal separately)")
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
    # Union of all + extra
    all_studies <- c(NP_STUDIES, AF_STUDIES, CEP_STUDIES, ALL_CELLS_EXTRA)
    # Deduplicate by (acc, comp)
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
  # Convert h5ad to h5seurat, then load
  # TODO: If SeuratDisk is unavailable, use anndata + reticulate as fallback
  h5seurat_path <- sub("\\.h5ad$", ".h5seurat", h5ad_path)

  if (!file.exists(h5seurat_path)) {
    message("  Converting ", basename(h5ad_path), " to h5seurat...")
    Convert(h5ad_path, dest = "h5seurat", overwrite = TRUE)
  }

  message("  Loading ", basename(h5seurat_path), "...")
  obj <- LoadH5Seurat(h5seurat_path, assays = "RNA")
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

    # Load as Seurat
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

    # Filter to cells with valid cell_class (from Module 04)
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
# CCA INTEGRATION (FLAT)
# ═══════════════════════════════════════════════════════════════════════════

run_cca_flat <- function(seurat_list, object_name) {
  message("\n  Running CCA flat integration for ", object_name, "...")

  # Merge all objects
  if (length(seurat_list) == 1) {
    merged <- seurat_list[[1]]
  } else {
    merged <- merge(seurat_list[[1]], y = seurat_list[-1],
                    add.cell.ids = names(seurat_list))
  }
  message("    Merged: ", ncol(merged), " cells x ", nrow(merged), " genes")

  # Split layers by study for integration
  merged[["RNA"]] <- split(merged[["RNA"]], f = merged$study)

  # SCTransform normalization
  message("    Running SCTransform...")
  merged <- SCTransform(merged, vars.to.regress = "pct_counts_mt",
                        verbose = FALSE)

  # CCA integration
  message("    Running CCA integration (dims = 1:50)...")
  merged <- IntegrateLayers(
    object = merged,
    method = CCAIntegration,
    normalization.method = "SCT",
    dims = 1:50,
    verbose = FALSE
  )

  # Rejoin layers after integration
  merged[["RNA"]] <- JoinLayers(merged[["RNA"]])

  # Dimensionality reduction
  message("    Running PCA, UMAP, FindNeighbors...")
  merged <- RunPCA(merged, npcs = 50, verbose = FALSE)
  merged <- RunUMAP(merged, dims = 1:50, verbose = FALSE)
  merged <- FindNeighbors(merged, dims = 1:50, verbose = FALSE)

  message("    CCA flat integration complete: ", ncol(merged), " cells")
  return(merged)
}


# ═══════════════════════════════════════════════════════════════════════════
# CCA INTEGRATION (TIERED)
# ═══════════════════════════════════════════════════════════════════════════

run_cca_tiered <- function(seurat_list, object_name) {
  message("\n  Running CCA tiered integration for ", object_name, "...")

  # Merge first to split by cell_class
  if (length(seurat_list) == 1) {
    merged <- seurat_list[[1]]
  } else {
    merged <- merge(seurat_list[[1]], y = seurat_list[-1],
                    add.cell.ids = names(seurat_list))
  }

  # Split by cell_class
  mes_cells <- which(merged@meta.data$cell_class %in% c("mesenchymal", "unknown"))
  non_mes_cells <- which(merged@meta.data$cell_class == "non_mesenchymal")

  mes_obj <- NULL
  non_mes_obj <- NULL

  # Tier A: Mesenchymal
  if (length(mes_cells) >= 50) {
    message("    Tier A: Mesenchymal (", length(mes_cells), " cells)")
    mes_obj <- subset(merged, cells = colnames(merged)[mes_cells])
    mes_obj[["RNA"]] <- split(mes_obj[["RNA"]], f = mes_obj$study)
    mes_obj <- SCTransform(mes_obj, vars.to.regress = "pct_counts_mt",
                           verbose = FALSE)
    mes_obj <- IntegrateLayers(
      object = mes_obj, method = CCAIntegration,
      normalization.method = "SCT", dims = 1:50, verbose = FALSE
    )
    mes_obj[["RNA"]] <- JoinLayers(mes_obj[["RNA"]])
    mes_obj <- RunPCA(mes_obj, npcs = 50, verbose = FALSE)
    mes_obj <- RunUMAP(mes_obj, dims = 1:50, verbose = FALSE)
    mes_obj <- FindNeighbors(mes_obj, dims = 1:50, verbose = FALSE)
  } else {
    message("    Skipping mesenchymal tier: only ", length(mes_cells), " cells")
  }

  # Tier B: Non-mesenchymal
  if (length(non_mes_cells) >= 200) {
    message("    Tier B: Non-mesenchymal (", length(non_mes_cells), " cells)")
    non_mes_obj <- subset(merged, cells = colnames(merged)[non_mes_cells])
    non_mes_obj[["RNA"]] <- split(non_mes_obj[["RNA"]], f = non_mes_obj$study)
    non_mes_obj <- SCTransform(non_mes_obj, vars.to.regress = "pct_counts_mt",
                               verbose = FALSE)
    non_mes_obj <- IntegrateLayers(
      object = non_mes_obj, method = CCAIntegration,
      normalization.method = "SCT", dims = 1:50, verbose = FALSE
    )
    non_mes_obj[["RNA"]] <- JoinLayers(non_mes_obj[["RNA"]])
    non_mes_obj <- RunPCA(non_mes_obj, npcs = 50, verbose = FALSE)
    non_mes_obj <- RunUMAP(non_mes_obj, dims = 1:50, verbose = FALSE)
    non_mes_obj <- FindNeighbors(non_mes_obj, dims = 1:50, verbose = FALSE)
  } else {
    message("    Skipping non-mesenchymal tier: only ", length(non_mes_cells), " cells")
  }

  # Merge tiers back
  # TODO: Implement tier merging strategy — store tier-specific reductions
  # in separate DimReduc slots and merge the objects
  if (!is.null(mes_obj) && !is.null(non_mes_obj)) {
    # For now, merge and re-run PCA/UMAP on the merged object
    result <- merge(mes_obj, non_mes_obj)
    result <- RunPCA(result, npcs = 50, verbose = FALSE)
    result <- RunUMAP(result, dims = 1:50, verbose = FALSE)
    result <- FindNeighbors(result, dims = 1:50, verbose = FALSE)
  } else if (!is.null(mes_obj)) {
    result <- mes_obj
  } else if (!is.null(non_mes_obj)) {
    result <- non_mes_obj
  } else {
    message("  ERROR: No tiers processed for ", object_name)
    return(NULL)
  }

  message("    CCA tiered integration complete: ", ncol(result), " cells")
  return(result)
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

process_object <- function(object_name, tiered = FALSE, force = FALSE) {
  output_path <- file.path(INT_DIR, paste0(object_name, ".rds"))

  if (file.exists(output_path) && !force) {
    message("\n=== ", object_name, ": output exists, skipping (use --force) ===")
    return(output_path)
  }

  message("\n", paste(rep("=", 60), collapse = ""))
  message("Processing object: ", object_name, " (CCA ",
          ifelse(tiered, "tiered", "flat"), ")")
  message(paste(rep("=", 60), collapse = ""))

  # Load data
  seurat_list <- load_and_build_object(object_name)
  if (is.null(seurat_list)) return(NULL)

  # Run integration
  if (tiered) {
    result <- run_cca_tiered(seurat_list, object_name)
  } else {
    result <- run_cca_flat(seurat_list, object_name)
  }

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

  # Clean up
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

  message("  Computing metrics for CCA/", object_name, "...")
  obj <- readRDS(rds_path)

  metrics <- data.frame(
    object = object_name,
    workflow = "cca",
    n_cells = ncol(obj),
    stringsAsFactors = FALSE
  )

  # TODO: Compute iLISI using lisi::compute_lisi()
  # TODO: Compute batch-ASW using cluster::silhouette()
  # TODO: Compute condition-ASW

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

      # Check for UMAP
      if ("umap" %in% names(obj@reductions)) {
        messages <- c(messages, paste0("PASS: cca/", obj_name, " has UMAP reduction"))
      } else {
        messages <- c(messages, paste0("WARNING: cca/", obj_name, " missing UMAP reduction"))
      }

      # Check for PCA
      if ("pca" %in% names(obj@reductions)) {
        messages <- c(messages, paste0("PASS: cca/", obj_name, " has PCA reduction"))
      } else {
        messages <- c(messages, paste0("WARNING: cca/", obj_name, " missing PCA reduction"))
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
  message("Module 05 — Workflow A: CCA Integration (PRIMARY)")
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
    output <- process_object(obj_name, tiered = args$tiered, force = args$force)
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
