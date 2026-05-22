#!/usr/bin/env Rscript
# ============================================================================
# Module 05k — Tiered v4 Integration for AF / CEP / all_cells
#
# Extends the tiered_v4 strategy validated for NP (see
# `scripts/05g_np_experiment.R` and `docs/np_switch_to_tiered_v4_plan.md`)
# to the AF, CEP, and all_cells compartments. Pipeline matches the NP
# tiered_v4 path: SCTransform per study + Seurat v4 FindIntegrationAnchors
# (CCA) + IntegrateData, with mesenchymal / non_mesenchymal tier split.
#
# Outputs:
#   data/integrated/<compartment>_experiment/tiered_v4/
#     mesenchymal.rds                  (full Seurat v4 integrated obj)
#     non_mesenchymal.rds              (full Seurat v4 integrated obj)
#     mesenchymal/                     (bridge export for h5ad assembly)
#     non_mesenchymal/                 (bridge export)
#   data/integrated/<compartment>_experiment/_anchor_cache/
#     <comp>_tiered_v4_<tier>_anchors.rds
#
# Anchor checkpointing — IntegrateData is the OOM-prone step; if it dies,
# rerunning the script resumes from the saved anchorset and skips the
# expensive SCT + FindIntegrationAnchors recompute.
#
# Usage:
#   Rscript scripts/05k_tiered_v4_compartments.R --compartment AF
#   Rscript scripts/05k_tiered_v4_compartments.R --compartment CEP
#   Rscript scripts/05k_tiered_v4_compartments.R --compartment all_cells
#   Rscript scripts/05k_tiered_v4_compartments.R --compartment AF --skip-rds
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(argparse)
  library(dplyr)
  library(ggplot2)
  library(future)
  library(Matrix)
})

# ── Parallelism ────────────────────────────────────────────────────────────
# v4 SCT pipeline pinned to sequential plan (memory-bound on full-cell objects).
N_WORKERS <- 1
options(future.globals.maxSize = 200 * 1024^3)
message("  Parallelism: workers = ", N_WORKERS,
        " (sequential), BLAS = ", sessionInfo()$BLAS)

# ── Paths ──────────────────────────────────────────────────────────────────
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
PROC_DIR <- file.path(BASE, "data", "processed")

# ── Integration parameters ───────────────────────────────────────────────
N_HVG_DEFAULT <- 3000   # Seurat CCA vignette default
N_DIMS        <- 50

# Adaptive feature selection guards against int32 indexing overflow in
# Seurat v4 IntegrateData. The failing intermediate matrix scales with
# (n_cells × n_features × n_objects). Empirical calibration from
# tiered_v4 integration logs (2026-04 to 2026-05):
#
#   Tier             cells   objects  features  C×F×O     Outcome
#   AF mes          84,568   3        3000      7.6e8     OK
#   CEP mes         50,769   3        3000      4.6e8     OK
#   NP mes         262,951   8        3000      6.3e9     OK
#   all_cells mes  407,179   15       3000      1.83e10   FAIL
#   all_cells mes  407,179   15       1000      6.1e9     OK
#
# Failure boundary lies in (6.3e9, 1.83e10). Default ceiling 1e10 sits
# between the data points with margin on both sides. Override via
# --max-cfo-product (or --features-to-integrate for a fixed value).
INTEGRATION_CFO_CEILING_DEFAULT <- 1e10
N_HVG_FLOOR <- 500     # Never fall below this number of features

# Resolves the actual feature count to use for a given tier.
# Returns the number of integration features after applying the
# adaptive ceiling and any CLI overrides.
adaptive_features <- function(n_cells, n_objects,
                              hvg_cap = N_HVG_DEFAULT,
                              ceiling = INTEGRATION_CFO_CEILING_DEFAULT) {
  budget <- floor(ceiling / (as.numeric(n_cells) * as.numeric(n_objects)))
  chosen <- min(hvg_cap, max(budget, N_HVG_FLOOR))
  message(sprintf(
    "    [adaptive_features] %d cells x %d objects -> budget %d ft; using %d (cap=%d, ceiling=%.1e)",
    n_cells, n_objects, budget, chosen, hvg_cap, ceiling))
  chosen
}

# ── Stage timing helpers ─────────────────────────────────────────────────
.STAGE_T0 <- new.env()
.tic <- function(stage) {
  assign(stage, Sys.time(), envir = .STAGE_T0)
  message("    [t] ", stage, " start: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
}
.toc <- function(stage) {
  if (!exists(stage, envir = .STAGE_T0)) return(invisible(NULL))
  dt <- as.numeric(difftime(Sys.time(), get(stage, envir = .STAGE_T0), units = "secs"))
  message("    [t] ", stage, " end:   ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
          "  (elapsed ", sprintf("%.1f", dt), " s = ", sprintf("%.2f", dt / 60), " min)")
}

# ── Study assignments per compartment ─────────────────────────────────────
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

get_study_assignments <- function(compartment) {
  if (compartment == "NP")        return(NP_STUDIES)
  if (compartment == "AF")        return(AF_STUDIES)
  if (compartment == "CEP")       return(CEP_STUDIES)
  if (compartment == "all_cells") {
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
  stop("Unknown compartment: ", compartment)
}

# ── Argument parsing ─────────────────────────────────────────────────────
parser <- ArgumentParser(description = "Module 05k: Tiered v4 integration for AF/CEP/all_cells")
parser$add_argument("--compartment", type = "character", required = TRUE,
                    choices = c("AF", "CEP", "all_cells", "NP"),
                    help = "Compartment to integrate")
parser$add_argument("--force", action = "store_true", default = FALSE,
                    help = "Re-run even if outputs exist")
parser$add_argument("--skip-rds", action = "store_true", default = FALSE,
                    help = "Skip saveRDS of full integrated object (saves disk space; bridge export still produced)")
parser$add_argument("--features-to-integrate", type = "integer", default = NULL,
                    help = "Limit IntegrateData to top N anchor features (default: all). Use to avoid R Matrix 2^31 nnz overflow on very large objects (e.g. all_cells: 2000).")
parser$add_argument("--max-cfo-product", type = "double", default = INTEGRATION_CFO_CEILING_DEFAULT,
                    help = sprintf("Adaptive ceiling for cells x features x objects (default: %.1e). Override only if you have new calibration data.", INTEGRATION_CFO_CEILING_DEFAULT))
parser$add_argument("--hvg-cap", type = "integer", default = N_HVG_DEFAULT,
                    help = sprintf("Maximum HVGs (Seurat CCA vignette default: %d). Adaptive selection reduces below this when cell or object count is large.", N_HVG_DEFAULT))
args <- parser$parse_args()

EXP_DIR <- file.path(BASE, "data", "integrated", paste0(tolower(args$compartment), "_experiment"))
RESULTS_DIR <- file.path(BASE, "results", "integration", paste0(tolower(args$compartment), "_experiment"))
ANCHOR_CACHE_DIR <- file.path(EXP_DIR, "_anchor_cache")
MODE_DIR <- file.path(EXP_DIR, "tiered_v4")

dir.create(EXP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ANCHOR_CACHE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MODE_DIR, recursive = TRUE, showWarnings = FALSE)


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

  unlink(bridge_dir, recursive = TRUE)
  return(obj)
}


# ═══════════════════════════════════════════════════════════════════════════
# LOAD CELLS FOR THE REQUESTED COMPARTMENT
# ═══════════════════════════════════════════════════════════════════════════

load_compartment_cells <- function(compartment) {
  assignments <- get_study_assignments(compartment)
  message("\n  Loading cells for ", compartment, " (",
          length(assignments), " study/compartment entries)...")

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

    # Sample exclusions
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

    list_key <- paste0(acc, "_", ifelse(is.null(comp_filter), "all", comp_filter))
    message("    ", list_key, ": ", ncol(obj), " cells")
    seurat_list[[list_key]] <- obj
  }

  total <- sum(sapply(seurat_list, ncol))
  message("  Loaded ", length(seurat_list), " study/compartment objects, ", total, " total cells")
  return(seurat_list)
}


# ═══════════════════════════════════════════════════════════════════════════
# SPLIT INTO MESENCHYMAL / NON-MESENCHYMAL TIERS
# ═══════════════════════════════════════════════════════════════════════════

split_tiers <- function(seurat_list) {
  message("\n  Splitting into mesenchymal / non_mesenchymal tiers...")

  mes_list <- list()
  nonmes_list <- list()

  for (key in names(seurat_list)) {
    obj <- seurat_list[[key]]
    cc <- obj@meta.data$cell_class

    mes_cells <- which(cc %in% c("mesenchymal", "unknown"))
    if (length(mes_cells) > 0) {
      mes_list[[key]] <- subset(obj, cells = colnames(obj)[mes_cells])
    }

    nonmes_cells <- which(cc == "non_mesenchymal")
    if (length(nonmes_cells) > 0) {
      nonmes_list[[key]] <- subset(obj, cells = colnames(obj)[nonmes_cells])
    }
  }

  mes_total <- if (length(mes_list) > 0) sum(sapply(mes_list, ncol)) else 0
  nonmes_total <- if (length(nonmes_list) > 0) sum(sapply(nonmes_list, ncol)) else 0

  message("    Mesenchymal tier: ", length(mes_list), " objects, ", mes_total, " cells")
  message("    Non-mesenchymal tier: ", length(nonmes_list), " objects, ", nonmes_total, " cells")

  return(list(mesenchymal = mes_list, non_mesenchymal = nonmes_list))
}


# ═══════════════════════════════════════════════════════════════════════════
# INTEGRATION: Seurat v4 SCT + FindIntegrationAnchors + IntegrateData
# ═══════════════════════════════════════════════════════════════════════════

integrate_v4_cca <- function(seurat_list, label, force = FALSE) {
  message("\n  [", label, "] Running v4 CCA (SCT + FindIntegrationAnchors + IntegrateData)...")
  study_sizes <- sapply(seurat_list, ncol)
  total_cells <- sum(study_sizes)
  smallest <- min(study_sizes)
  k_filter <- min(200, max(5, smallest - 5))
  k_weight <- min(100, max(5, smallest - 5))
  message("    Input: ", length(seurat_list), " objects, ", total_cells,
          " cells (smallest = ", smallest, ")")

  anchor_cache <- file.path(ANCHOR_CACHE_DIR, paste0(label, "_anchors.rds"))

  if (file.exists(anchor_cache) && !force) {
    message("    [cache] Loading anchorset from ", anchor_cache,
            " — skipping SCT + FindIntegrationAnchors")
    .tic("v4_load_anchor_cache")
    anchors <- readRDS(anchor_cache)
    .toc("v4_load_anchor_cache")

    # Trim per-object Seurats inside the anchorset to reduce IntegrateData
    # memory peak and (when --features-to-integrate is set) work around the
    # 2^31 nnz overflow in Seurat v4's sparse difference matrix construction.
    #
    # Two trims:
    #   (1) Drop the RNA assay entirely — IntegrateData(normalization.method="SCT")
    #       reads SCT residuals, not RNA counts/data, so RNA is dead weight
    #       (~3-4 GB per object as dense layers).
    #   (2) If args$features_to_integrate is set, also subset each object's SCT
    #       data/scale.data and anchors@anchor.features to the top N anchor
    #       features. This forces Seurat down the default IntegrateData code
    #       path (no features.to.integrate kwarg, which triggers a much slower
    #       branch) while still capping the integrated feature count.
    if (!is.null(anchors@object.list)) {
      n_features_target <- if (!is.null(args$features_to_integrate) &&
                               args$features_to_integrate > 0)
                            args$features_to_integrate else NULL
      keep_features <- NULL
      if (!is.null(n_features_target)) {
        all_anchor_features <- anchors@anchor.features
        n_use <- min(n_features_target, length(all_anchor_features))
        keep_features <- all_anchor_features[seq_len(n_use)]
        message("    [trim] Will subset SCT data to top ", n_use, " of ",
                length(all_anchor_features), " anchor features (avoids 2^31 overflow)")
      }

      message("    [trim] Stripping RNA assay + reductions from ",
              length(anchors@object.list), " per-object Seurats...")
      .tic("v4_trim_anchors")
      for (i in seq_along(anchors@object.list)) {
        obj <- anchors@object.list[[i]]
        nm <- names(anchors@object.list)[i]

        tryCatch({
          if ("SCT" %in% names(obj@assays)) {
            DefaultAssay(obj) <- "SCT"
          }
        }, error = function(e) {
          message("      (set SCT default failed on ", nm, ": ",
                  conditionMessage(e), ")")
        })

        if ("RNA" %in% names(obj@assays)) {
          tryCatch({
            obj[["RNA"]] <- NULL
          }, error = function(e) {
            message("      (drop RNA failed on ", nm, ": ",
                    conditionMessage(e), ")")
          })
        }

        # Subset SCT to the target feature set so IntegrateData stays under
        # the 2^31 nnz limit when constructing the merge difference matrix.
        if (!is.null(keep_features) && "SCT" %in% names(obj@assays)) {
          tryCatch({
            sct <- obj[["SCT"]]
            present <- intersect(keep_features, rownames(sct))
            if (length(present) < length(keep_features)) {
              message("      (", nm, ": ", length(present), " of ",
                      length(keep_features), " target features present in SCT)")
            }
            obj[["SCT"]] <- subset(sct, features = present)
          }, error = function(e) {
            message("      (SCT feature subset failed on ", nm, ": ",
                    conditionMessage(e), ")")
          })
        }

        obj@reductions <- list()
        obj@graphs <- list()
        obj@neighbors <- list()
        obj@commands <- list()
        anchors@object.list[[i]] <- obj
      }

      # Update the anchorset's anchor.features so IntegrateData picks up the
      # reduced feature set as the default.
      if (!is.null(keep_features)) {
        present_in_all <- Reduce(
          intersect,
          lapply(anchors@object.list, function(o) {
            if ("SCT" %in% names(o@assays)) rownames(o[["SCT"]]) else character(0)
          })
        )
        new_anchor_features <- intersect(keep_features, present_in_all)
        message("    [trim] Setting anchors@anchor.features = ",
                length(new_anchor_features), " features (intersection across all objects)")
        anchors@anchor.features <- new_anchor_features
      }

      gc(verbose = FALSE)
      .toc("v4_trim_anchors")
    }
  } else {
    message("    SCTransform per study/object...")
    .tic("v4_sctransform_total")
    for (key in names(seurat_list)) {
      n <- ncol(seurat_list[[key]])
      message("      ", key, " (", n, " cells)...")
      .tic(paste0("v4_sct_", key))
      seurat_list[[key]] <- SCTransform(
        seurat_list[[key]],
        vars.to.regress = "pct_counts_mt",
        verbose = FALSE
      )
      .toc(paste0("v4_sct_", key))
    }
    .toc("v4_sctransform_total")

    n_cells_tier <- sum(sapply(seurat_list, ncol))
    n_obj_tier <- length(seurat_list)
    n_features_tier <- adaptive_features(
      n_cells_tier, n_obj_tier,
      hvg_cap = args$hvg_cap, ceiling = args$max_cfo_product
    )
    message("    SelectIntegrationFeatures (", n_features_tier, " features)...")
    .tic("v4_select_features")
    features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = n_features_tier)
    .toc("v4_select_features")

    message("    PrepSCTIntegration...")
    .tic("v4_prep_sct")
    seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features)
    .toc("v4_prep_sct")

    message("    FindIntegrationAnchors (CCA, dims = 1:", N_DIMS,
            ", k.filter = ", k_filter, ", sequential)...")
    .tic("v4_find_anchors")
    plan("sequential")
    anchors <- FindIntegrationAnchors(
      object.list = seurat_list,
      normalization.method = "SCT",
      anchor.features = features,
      reduction = "cca",
      dims = 1:N_DIMS,
      k.filter = k_filter,
      verbose = FALSE
    )
    .toc("v4_find_anchors")

    message("    [cache] Saving anchorset to ", anchor_cache, "...")
    .tic("v4_save_anchor_cache")
    saveRDS(anchors, anchor_cache)
    .toc("v4_save_anchor_cache")

    rm(seurat_list)
    gc(verbose = FALSE)
  }

  # If --features-to-integrate was set, the trim step above has already
  # subset each per-object SCT and anchors@anchor.features to the target
  # feature count. IntegrateData uses anchor.features by default.
  message("    IntegrateData (k.weight = ", k_weight,
          ", n anchor features = ", length(anchors@anchor.features), ")...")
  .tic("v4_integrate_data")
  integrated <- IntegrateData(
    anchorset = anchors,
    normalization.method = "SCT",
    dims = 1:N_DIMS,
    k.weight = k_weight,
    verbose = FALSE
  )
  .toc("v4_integrate_data")

  rm(anchors)
  gc(verbose = FALSE)

  DefaultAssay(integrated) <- "integrated"
  message("    ScaleData + RunPCA + UMAP + Neighbors on integrated assay...")
  .tic("v4_scale_pca_umap_neighbors")
  integrated <- ScaleData(integrated, verbose = FALSE)
  integrated <- RunPCA(integrated, npcs = N_DIMS, verbose = FALSE)
  integrated <- RunUMAP(integrated, reduction = "pca", dims = 1:N_DIMS, verbose = FALSE)
  integrated <- FindNeighbors(integrated, reduction = "pca", dims = 1:N_DIMS, verbose = FALSE)
  .toc("v4_scale_pca_umap_neighbors")

  message("    v4 CCA complete: ", ncol(integrated), " cells")
  return(integrated)
}


# ═══════════════════════════════════════════════════════════════════════════
# INTEGRATION: Lightweight merge for small non-mesenchymal tier
# ═══════════════════════════════════════════════════════════════════════════

integrate_simple <- function(seurat_list, label) {
  message("\n  [", label, "] Simple merge + PCA (too few objects for full integration)...")
  total_cells <- sum(sapply(seurat_list, ncol))
  message("    Input: ", length(seurat_list), " objects, ", total_cells, " cells")

  if (length(seurat_list) == 1) {
    merged <- seurat_list[[1]]
  } else {
    merged <- merge(seurat_list[[1]], y = seurat_list[-1])
  }

  merged <- JoinLayers(merged)
  merged <- NormalizeData(merged, verbose = FALSE)
  n_features_simple <- adaptive_features(
    total_cells, length(seurat_list),
    hvg_cap = args$hvg_cap, ceiling = args$max_cfo_product
  )
  merged <- FindVariableFeatures(merged, nfeatures = n_features_simple, verbose = FALSE)
  merged <- ScaleData(merged, verbose = FALSE)

  dims_use <- min(N_DIMS, ncol(merged) - 1, nrow(merged) - 1)
  merged <- RunPCA(merged, npcs = dims_use, verbose = FALSE)
  merged <- RunUMAP(merged, reduction = "pca", dims = 1:dims_use, verbose = FALSE)
  merged <- FindNeighbors(merged, reduction = "pca", dims = 1:dims_use, verbose = FALSE)

  message("    Simple merge complete: ", ncol(merged), " cells")
  return(merged)
}


# ═══════════════════════════════════════════════════════════════════════════
# RAW COUNTS CACHE: merge per-study raw RNA counts before integration so we
# can write actual counts in export_bridge regardless of what survives in
# the integrated Seurat object's assay slots.
#
# IntegrateData(normalization.method = "SCT") returns an object whose RNA
# assay state is unreliable across the trim/cache code path (the trim block
# drops RNA from anchors@object.list to control RAM). Falling back to the
# integrated assay's `data` layer would silently write CCA-corrected dense
# values — wrong for downstream pseudobulk DE.
# ═══════════════════════════════════════════════════════════════════════════

stash_raw_counts <- function(seurat_list, cache_path, force = FALSE) {
  if (file.exists(cache_path) && !force) {
    message("    [raw counts cache] exists, skipping build: ", cache_path)
    return(invisible(NULL))
  }
  if (file.exists(cache_path) && force) {
    message("    [raw counts cache] --force: removing stale cache and rebuilding")
    file.remove(cache_path)
  }
  message("    [raw counts cache] Building merged raw counts → ", cache_path)
  .tic("raw_counts_stash")

  all_genes <- Reduce(union, lapply(seurat_list, rownames))
  message("    [raw counts cache] Gene union: ", length(all_genes), " genes across ",
          length(seurat_list), " objects")

  mats <- vector("list", length(seurat_list))
  names(mats) <- names(seurat_list)
  for (k in names(seurat_list)) {
    m <- LayerData(seurat_list[[k]], assay = "RNA", layer = "counts")
    missing <- setdiff(all_genes, rownames(m))
    if (length(missing) > 0) {
      filler <- Matrix(0, nrow = length(missing), ncol = ncol(m), sparse = TRUE)
      rownames(filler) <- missing
      colnames(filler) <- colnames(m)
      m <- rbind(m, filler)
    }
    mats[[k]] <- m[all_genes, , drop = FALSE]
    message("    [raw counts cache] ", k, ": ", ncol(m), " cells")
  }

  merged <- do.call(cbind, mats)
  rm(mats); gc(verbose = FALSE)
  merged <- as(merged, "CsparseMatrix")
  message("    [raw counts cache] Merged: ", nrow(merged), " genes × ",
          ncol(merged), " cells, nnz = ", length(merged@x))

  saveRDS(merged, cache_path, compress = TRUE)
  rm(merged); gc(verbose = FALSE)
  .toc("raw_counts_stash")
  message("    [raw counts cache] Saved")
}


# ═══════════════════════════════════════════════════════════════════════════
# BRIDGE EXPORT: Seurat → flat files for Python metrics / h5ad assembly
# ═══════════════════════════════════════════════════════════════════════════

export_bridge <- function(obj, export_dir, raw_counts_cache = NULL) {
  dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)
  message("  Exporting bridge files to ", export_dir, "...")
  message("    Cells: ", ncol(obj), ", Genes: ", nrow(obj))
  message("    Assays: ", paste(names(obj@assays), collapse = ", "))
  message("    Reductions: ", paste(names(obj@reductions), collapse = ", "))

  counts <- NULL

  # Preferred path: raw counts cache built before integration. Subset & reorder
  # columns to match the integrated object's cell order.
  if (!is.null(raw_counts_cache) && file.exists(raw_counts_cache)) {
    message("    Loading raw counts from cache: ", raw_counts_cache)
    counts <- readRDS(raw_counts_cache)
    cells_int <- colnames(obj)
    missing_cells <- setdiff(cells_int, colnames(counts))
    if (length(missing_cells) > 0) {
      stop("Raw counts cache is missing ", length(missing_cells),
           " cells from the integrated object (e.g. ",
           paste(head(missing_cells, 3), collapse = ", "), ")")
    }
    counts <- counts[, cells_int, drop = FALSE]
    message("    Raw counts subset to integrated cells: ",
            nrow(counts), " genes × ", ncol(counts), " cells")
  } else {
    # Fallback path (legacy AF/CEP behavior): pull from the surviving RNA assay.
    if ("RNA" %in% names(obj@assays)) {
      tryCatch({
        obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
      }, error = function(e) {
        message("    (JoinLayers skipped: ", conditionMessage(e), ")")
      })
    }
    tryCatch({
      counts <- LayerData(obj, assay = "RNA", layer = "counts")
    }, error = function(e) {
      message("    counts layer not found in RNA, trying GetAssayData...")
    })
    if (is.null(counts)) {
      tryCatch({
        counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
      }, error = function(e) {
        # No silent fallback to integrated `data` — that's CCA-corrected
        # dense values, not raw counts. Caller must provide raw_counts_cache
        # when the RNA assay does not survive integration.
        stop("Could not obtain raw counts from RNA assay and no ",
             "raw_counts_cache was provided. Refusing to write ",
             "CCA-corrected integrated data as 'counts.mtx' — that ",
             "would silently break downstream pseudobulk DE.")
      })
    }
  }

  message("    Writing counts matrix...")
  mtx_tmp <- file.path(export_dir, "counts.mtx")
  writeMM(counts, mtx_tmp)
  system2("gzip", c("-f", mtx_tmp))

  genes <- rownames(counts)
  write.csv(data.frame(gene = genes), file.path(export_dir, "genes.csv"), row.names = FALSE)

  barcodes <- colnames(obj)
  write.csv(data.frame(barcode = barcodes), file.path(export_dir, "barcodes.csv"), row.names = FALSE)

  meta <- obj@meta.data
  meta$barcode <- rownames(meta)
  write.csv(meta, gzfile(file.path(export_dir, "metadata.csv.gz")), row.names = FALSE)
  message("    Saved metadata: ", ncol(meta), " columns")

  for (red_name in names(obj@reductions)) {
    emb <- Embeddings(obj, reduction = red_name)
    out_name <- paste0("embedding_", red_name, ".csv.gz")
    write.csv(emb, gzfile(file.path(export_dir, out_name)), row.names = TRUE)
    message("    Saved embedding: ", red_name, " (", ncol(emb), " dims)")
  }

  message("    Bridge export complete")
}


# ═══════════════════════════════════════════════════════════════════════════
# UMAP PLOTS
# ═══════════════════════════════════════════════════════════════════════════

plot_umaps <- function(obj, run_label) {
  color_by <- intersect(c("study", "condition_harmonized", "cell_class", "coarse_label", "compartment"),
                        colnames(obj@meta.data))
  if (length(color_by) == 0) return(invisible(NULL))

  n_cells <- ncol(obj)
  use_raster <- n_cells > 100000
  pt <- if (use_raster) 2 else 0.3
  rdpi <- c(2048, 2048)

  plots <- lapply(color_by, function(col) {
    DimPlot(obj, group.by = col, reduction = "umap",
            pt.size = pt, raster = use_raster, raster.dpi = rdpi,
            shuffle = TRUE, alpha = 1) +
      ggtitle(paste(run_label, "-", col)) +
      theme(legend.position = "right",
            legend.text = element_text(size = 7))
  })

  p <- patchwork::wrap_plots(plots, ncol = min(length(plots), 2))
  n_rows <- ceiling(length(plots) / 2)
  out_path <- file.path(RESULTS_DIR, paste0("umap_", run_label, ".png"))
  ggsave(out_path, plot = p, width = 12, height = 5 * n_rows, dpi = 150)
  message("  Saved: ", out_path)
}


# ═══════════════════════════════════════════════════════════════════════════
# PROCESS NON-MESENCHYMAL TIER
# ═══════════════════════════════════════════════════════════════════════════

process_nonmes_tier <- function(nonmes_list, mode_prefix, force = FALSE) {
  # Per-study minimum cell count for the non-mesenchymal tier.
  # Set to 50 on 2026-05-21 so every surviving study clears CCA's
  # FindIntegrationAnchors dims=1:50 threshold — the prior 5-cell threshold
  # left objects too small for CCA, forcing a simple-merge fallback that
  # produced study-segregated UMAP lobes (no batch correction). Studies
  # below 50 non-mes cells are dropped from the non-mes tier rather than
  # mixing uncorrected.
  NONMES_MIN_CELLS_PER_STUDY <- 50
  study_sizes <- sapply(nonmes_list, ncol)
  keep <- names(study_sizes)[study_sizes >= NONMES_MIN_CELLS_PER_STUDY]
  dropped <- names(study_sizes)[study_sizes < NONMES_MIN_CELLS_PER_STUDY]

  if (length(dropped) > 0) {
    message("    Dropping ", length(dropped),
            " objects with <", NONMES_MIN_CELLS_PER_STUDY, " non-mes cells: ",
            paste(dropped, " (", study_sizes[dropped], ")", sep = "", collapse = ", "))
    nonmes_list <- nonmes_list[keep]
  }

  if (length(nonmes_list) == 0) {
    message("    WARNING: No objects with >=", NONMES_MIN_CELLS_PER_STUDY,
            " non-mesenchymal cells, skipping tier")
    return(NULL)
  }

  # Post-drop: every surviving object has >= 50 cells (CCA's per-object
  # dims=1:50 threshold). CCA needs >= 2 objects; with a single survivor
  # there is no batch to integrate, so we fall back to simple merge.
  if (length(nonmes_list) < 2) {
    message("    Only ", length(nonmes_list),
            " object(s) survived the >=", NONMES_MIN_CELLS_PER_STUDY,
            "-cell drop — no integration possible, using simple merge")
    result <- integrate_simple(nonmes_list, paste0(mode_prefix, "_non_mesenchymal"))
  } else {
    result <- integrate_v4_cca(nonmes_list, paste0(mode_prefix, "_non_mesenchymal"), force = force)
  }

  return(result)
}


# ═══════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════

run_tiered_v4 <- function(seurat_list, force = FALSE, skip_rds = FALSE) {
  mode_prefix <- paste0(tolower(args$compartment), "_tiered_v4")
  mes_rds <- file.path(MODE_DIR, "mesenchymal.rds")
  nonmes_rds <- file.path(MODE_DIR, "non_mesenchymal.rds")
  mes_bridge <- file.path(MODE_DIR, "mesenchymal")
  nonmes_bridge <- file.path(MODE_DIR, "non_mesenchymal")
  mes_counts_cache <- file.path(MODE_DIR, "_raw_counts_mesenchymal.rds")
  nonmes_counts_cache <- file.path(MODE_DIR, "_raw_counts_non_mesenchymal.rds")

  mes_done <- (skip_rds && dir.exists(mes_bridge) &&
               file.exists(file.path(mes_bridge, "counts.mtx.gz"))) ||
              file.exists(mes_rds)
  nonmes_done <- (skip_rds && dir.exists(nonmes_bridge) &&
                  file.exists(file.path(nonmes_bridge, "counts.mtx.gz"))) ||
                 file.exists(nonmes_rds)

  if (mes_done && nonmes_done && !force) {
    message("\n=== ", args$compartment, " tiered_v4: outputs exist, skipping (use --force) ===")
    return(invisible(NULL))
  }

  message("\n", paste(rep("=", 60), collapse = ""))
  message("MODE: tiered_v4 — Seurat v4 SCT + CCA, mes/non-mes split")
  message("Compartment: ", args$compartment)
  message(paste(rep("=", 60), collapse = ""))

  tiers <- split_tiers(seurat_list)

  # Build raw counts caches BEFORE integration, while per-study RNA assays
  # are still intact. These are the source of truth for the bridge .mtx
  # write — IntegrateData(SCT) does not preserve RNA reliably across the
  # anchor-cache trim path.
  if (!mes_done || force) {
    if (length(tiers$mesenchymal) > 0) {
      stash_raw_counts(tiers$mesenchymal, mes_counts_cache, force = force)
    }
  }
  if (!nonmes_done || force) {
    # Match the per-study threshold used in process_nonmes_tier (≥50 cells).
    NONMES_MIN_CELLS_PER_STUDY <- 50
    nonmes_keep <- tiers$non_mesenchymal[
      sapply(tiers$non_mesenchymal, ncol) >= NONMES_MIN_CELLS_PER_STUDY
    ]
    if (length(nonmes_keep) > 0) {
      stash_raw_counts(nonmes_keep, nonmes_counts_cache, force = force)
    }
  }

  # Mesenchymal tier
  if (mes_done && !force) {
    message("\n  [", mode_prefix, "_mesenchymal] Resuming: outputs already exist, skipping")
  } else {
    mes_result <- integrate_v4_cca(tiers$mesenchymal, paste0(mode_prefix, "_mesenchymal"), force = force)
    plot_umaps(mes_result, paste0(mode_prefix, "_mesenchymal"))
    export_bridge(mes_result, mes_bridge, raw_counts_cache = mes_counts_cache)
    if (!skip_rds) {
      message("  Saving ", mes_rds, "...")
      saveRDS(mes_result, mes_rds)
    } else {
      message("  --skip-rds: not saving mesenchymal.rds")
    }
    rm(mes_result); gc(verbose = FALSE)
  }

  # Non-mesenchymal tier
  if (nonmes_done && !force) {
    message("\n  [", mode_prefix, "_non_mesenchymal] Resuming: outputs already exist, skipping")
  } else {
    nonmes_result <- process_nonmes_tier(tiers$non_mesenchymal, mode_prefix, force = force)
    if (!is.null(nonmes_result)) {
      plot_umaps(nonmes_result, paste0(mode_prefix, "_non_mesenchymal"))
      export_bridge(nonmes_result, nonmes_bridge, raw_counts_cache = nonmes_counts_cache)
      if (!skip_rds) {
        saveRDS(nonmes_result, nonmes_rds)
      } else {
        message("  --skip-rds: not saving non_mesenchymal.rds")
      }
      rm(nonmes_result); gc(verbose = FALSE)
    }
  }

  message("\n  ", args$compartment, " tiered_v4 complete")
}

main <- function() {
  message(paste(rep("=", 60), collapse = ""))
  message("Module 05k — Tiered v4 Integration (", args$compartment, ")")
  message("Seurat version: ", as.character(packageVersion("Seurat")))
  message("HVG cap: ", args$hvg_cap,
          ", adaptive ceiling (C×F×O): ", sprintf("%.1e", args$max_cfo_product),
          ", Dims: ", N_DIMS)
  if (!is.null(args$features_to_integrate)) {
    message("features-to-integrate override: ", args$features_to_integrate,
            " (trims SCT data post-anchor-finding)")
  }
  message("BLAS: ", sessionInfo()$BLAS)
  message("CPUs available: ", parallel::detectCores())
  message("skip_rds: ", args$skip_rds, ", force: ", args$force)
  message("Started: ", Sys.time())
  message(paste(rep("=", 60), collapse = ""))

  seurat_list <- load_compartment_cells(args$compartment)
  if (length(seurat_list) == 0) {
    stop("No cells loaded for ", args$compartment)
  }

  # Deep copy: SCTransform mutates objects; pristine copies isolate state across runs.
  sl_copy <- lapply(seurat_list, function(obj) {
    CreateSeuratObject(
      counts = LayerData(obj, layer = "counts"),
      meta.data = obj@meta.data
    )
  })
  rm(seurat_list)
  gc(verbose = FALSE)

  run_tiered_v4(sl_copy, force = args$force, skip_rds = args$skip_rds)

  message("\nCompleted: ", Sys.time())
}

main()
