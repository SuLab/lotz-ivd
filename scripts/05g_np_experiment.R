#!/usr/bin/env Rscript
# ============================================================================
# NP Integration Quality Experiment
#
# Compares three CCA integration strategies for NP to address over-integration:
#   tiered_v5  — Seurat v5 IntegrateLayers with mesenchymal/non-mesenchymal split
#   flat_v4    — Seurat v4 SCT + FindIntegrationAnchors + IntegrateData (no split)
#   tiered_v4  — Seurat v4 SCT + FindIntegrationAnchors + IntegrateData (with split)
#
# Usage:
#   Rscript scripts/05g_np_experiment.R --mode tiered_v5
#   Rscript scripts/05g_np_experiment.R --mode flat_v4
#   Rscript scripts/05g_np_experiment.R --mode tiered_v4
#   Rscript scripts/05g_np_experiment.R --mode all
#   Rscript scripts/05g_np_experiment.R --mode all --force
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
N_WORKERS <- min(parallel::detectCores(), 2)
options(future.globals.maxSize = 200 * 1024^3)
message("  Parallelism: workers = ", N_WORKERS,
        " (activated during integration), BLAS = ", sessionInfo()$BLAS)

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

PROC_DIR    <- file.path(BASE, "data", "processed")
EXP_DIR     <- file.path(BASE, "data", "integrated", "np_experiment")
RESULTS_DIR <- file.path(BASE, "results", "integration", "np_experiment")

dir.create(EXP_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Integration parameters ───────────────────────────────────────────────
N_HVG  <- 3000
N_DIMS <- 50

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

# ── NP study assignments ─────────────────────────────────────────────────
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

EXCLUDED_SAMPLES <- c("GSE251686_NP3")

# ── Argument parsing ─────────────────────────────────────────────────────
parser <- ArgumentParser(description = "NP Integration Quality Experiment")
parser$add_argument("--mode", type = "character", required = TRUE,
                    choices = c("tiered_v5", "flat_v4", "tiered_v4", "all"),
                    help = "Integration mode to run")
parser$add_argument("--force", action = "store_true", default = FALSE,
                    help = "Re-run even if outputs exist")
args <- parser$parse_args()


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
# LOAD ALL NP CELLS AS LIST OF PER-STUDY SEURAT OBJECTS
# ═══════════════════════════════════════════════════════════════════════════

load_np_cells <- function() {
  message("\n  Loading NP cells from ", length(NP_STUDIES), " studies...")

  seurat_list <- list()

  for (s in NP_STUDIES) {
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

    message("    ", acc, ": ", ncol(obj), " cells")
    seurat_list[[acc]] <- obj
  }

  total <- sum(sapply(seurat_list, ncol))
  message("  Loaded ", length(seurat_list), " studies, ", total, " total cells")
  return(seurat_list)
}


# ═══════════════════════════════════════════════════════════════════════════
# SPLIT INTO MESENCHYMAL / NON-MESENCHYMAL TIERS
# ═══════════════════════════════════════════════════════════════════════════

split_tiers <- function(seurat_list) {
  message("\n  Splitting into mesenchymal / non-mesenchymal tiers...")

  mes_list <- list()
  nonmes_list <- list()

  for (acc in names(seurat_list)) {
    obj <- seurat_list[[acc]]
    cc <- obj@meta.data$cell_class

    # Mesenchymal tier: mesenchymal + unknown
    mes_cells <- which(cc %in% c("mesenchymal", "unknown"))
    if (length(mes_cells) > 0) {
      mes_list[[acc]] <- subset(obj, cells = colnames(obj)[mes_cells])
    }

    # Non-mesenchymal tier
    nonmes_cells <- which(cc == "non_mesenchymal")
    if (length(nonmes_cells) > 0) {
      nonmes_list[[acc]] <- subset(obj, cells = colnames(obj)[nonmes_cells])
    }
  }

  mes_total <- sum(sapply(mes_list, ncol))
  nonmes_total <- if (length(nonmes_list) > 0) sum(sapply(nonmes_list, ncol)) else 0

  message("    Mesenchymal tier: ", length(mes_list), " studies, ", mes_total, " cells")
  message("    Non-mesenchymal tier: ", length(nonmes_list), " studies, ", nonmes_total, " cells")

  return(list(mesenchymal = mes_list, non_mesenchymal = nonmes_list))
}


# ═══════════════════════════════════════════════════════════════════════════
# INTEGRATION: Seurat v5 IntegrateLayers (CCA)
# ═══════════════════════════════════════════════════════════════════════════

integrate_v5_cca <- function(seurat_list, label) {
  message("\n  [", label, "] Running v5 CCA integration (IntegrateLayers)...")
  study_sizes <- sapply(seurat_list, ncol)
  total_cells <- sum(study_sizes)
  message("    Input: ", length(seurat_list), " studies, ", total_cells, " cells")

  # Merge into single v5 layered object
  if (length(seurat_list) == 1) {
    merged <- seurat_list[[1]]
  } else {
    merged <- merge(seurat_list[[1]], y = seurat_list[-1])
  }
  message("    Merged: ", ncol(merged), " cells, ", nrow(merged), " genes")

  rm(seurat_list)
  gc(verbose = FALSE)

  # Normalize + HVGs + Scale (per-layer, preserves split structure)
  message("    NormalizeData + FindVariableFeatures + ScaleData...")
  .tic("v5_norm_hvg_scale")
  merged <- NormalizeData(merged, verbose = FALSE)
  merged <- FindVariableFeatures(merged, nfeatures = N_HVG, verbose = FALSE)
  merged <- ScaleData(merged, verbose = FALSE)
  .toc("v5_norm_hvg_scale")

  message("    RunPCA (", N_DIMS, " dims)...")
  .tic("v5_pca")
  merged <- RunPCA(merged, npcs = N_DIMS, verbose = FALSE)
  .toc("v5_pca")

  # CCA integration — explicit HVG restriction to bound per-worker memory
  hvgs <- VariableFeatures(merged)
  # k.weight must be < smallest study size; Seurat default is 100
  k_weight <- min(100, max(5, min(study_sizes) - 5))
  message("    IntegrateLayers (CCA, dims = 1:", N_DIMS, ", features = ", length(hvgs),
          " HVGs, k.weight = ", k_weight, ", workers = ", N_WORKERS, ")...")
  .tic("v5_integrate_layers")
  plan("multicore", workers = N_WORKERS)
  merged <- IntegrateLayers(
    object = merged,
    method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.cca",
    features = hvgs,
    dims = 1:N_DIMS,
    k.weight = k_weight,
    verbose = FALSE
  )
  plan("sequential")
  .toc("v5_integrate_layers")

  merged <- JoinLayers(merged)

  # Post-integration
  message("    UMAP + neighbors...")
  .tic("v5_umap_neighbors")
  merged <- RunUMAP(merged, reduction = "integrated.cca", dims = 1:N_DIMS, verbose = FALSE)
  merged <- FindNeighbors(merged, reduction = "integrated.cca", dims = 1:N_DIMS, verbose = FALSE)
  .toc("v5_umap_neighbors")

  message("    v5 CCA complete: ", ncol(merged), " cells")
  return(merged)
}


# ═══════════════════════════════════════════════════════════════════════════
# INTEGRATION: Seurat v4 SCT + FindIntegrationAnchors + IntegrateData
# ═══════════════════════════════════════════════════════════════════════════

integrate_v4_cca <- function(seurat_list, label) {
  message("\n  [", label, "] Running v4 CCA integration (SCT + FindIntegrationAnchors)...")
  study_sizes <- sapply(seurat_list, ncol)
  total_cells <- sum(study_sizes)
  smallest <- min(study_sizes)
  # Seurat defaults: k.filter=200 (FindIntegrationAnchors), k.weight=100 (IntegrateData)
  k_filter <- min(200, max(5, smallest - 5))
  k_weight <- min(100, max(5, smallest - 5))
  message("    Input: ", length(seurat_list), " studies, ", total_cells, " cells (smallest = ", smallest, ")")

  # SCTransform each study independently
  message("    SCTransform per study...")
  .tic("v4_sctransform_total")
  for (acc in names(seurat_list)) {
    n <- ncol(seurat_list[[acc]])
    message("      ", acc, " (", n, " cells)...")
    .tic(paste0("v4_sct_", acc))
    seurat_list[[acc]] <- SCTransform(
      seurat_list[[acc]],
      vars.to.regress = "pct_counts_mt",
      verbose = FALSE
    )
    .toc(paste0("v4_sct_", acc))
  }
  .toc("v4_sctransform_total")

  # Select integration features
  message("    SelectIntegrationFeatures (", N_HVG, " features)...")
  .tic("v4_select_features")
  features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = N_HVG)
  .toc("v4_select_features")

  # Prep SCT integration
  message("    PrepSCTIntegration...")
  .tic("v4_prep_sct")
  seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features)
  .toc("v4_prep_sct")

  # Find integration anchors (CCA) — SCT residuals are ~30% denser than log-norm,
  # forcing sequential plan on v4: 2 workers OOM-kill on ~260K cells / 123 GB RAM.
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

  # Integrate data
  message("    IntegrateData (k.weight = ", k_weight, ")...")
  .tic("v4_integrate_data")
  integrated <- IntegrateData(
    anchorset = anchors,
    normalization.method = "SCT",
    dims = 1:N_DIMS,
    k.weight = k_weight,
    verbose = FALSE
  )
  .toc("v4_integrate_data")

  rm(seurat_list, anchors)
  gc(verbose = FALSE)

  # Post-integration on "integrated" assay
  DefaultAssay(integrated) <- "integrated"
  message("    ScaleData + RunPCA on integrated assay...")
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
  message("\n  [", label, "] Simple merge + PCA (too few cells/studies for full integration)...")
  total_cells <- sum(sapply(seurat_list, ncol))
  message("    Input: ", length(seurat_list), " studies, ", total_cells, " cells")

  if (length(seurat_list) == 1) {
    merged <- seurat_list[[1]]
  } else {
    merged <- merge(seurat_list[[1]], y = seurat_list[-1])
  }

  merged <- JoinLayers(merged)
  merged <- NormalizeData(merged, verbose = FALSE)
  merged <- FindVariableFeatures(merged, nfeatures = N_HVG, verbose = FALSE)
  merged <- ScaleData(merged, verbose = FALSE)

  dims_use <- min(N_DIMS, ncol(merged) - 1, nrow(merged) - 1)
  merged <- RunPCA(merged, npcs = dims_use, verbose = FALSE)
  merged <- RunUMAP(merged, reduction = "pca", dims = 1:dims_use, verbose = FALSE)
  merged <- FindNeighbors(merged, reduction = "pca", dims = 1:dims_use, verbose = FALSE)

  message("    Simple merge complete: ", ncol(merged), " cells")
  return(merged)
}


# ═══════════════════════════════════════════════════════════════════════════
# BRIDGE EXPORT: Seurat → flat files for Python metrics
# ═══════════════════════════════════════════════════════════════════════════

export_bridge <- function(obj, export_dir) {
  dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)
  message("  Exporting bridge files to ", export_dir, "...")
  message("    Cells: ", ncol(obj), ", Genes: ", nrow(obj))
  message("    Assays: ", paste(names(obj@assays), collapse = ", "))
  message("    Reductions: ", paste(names(obj@reductions), collapse = ", "))

  # 1. Export raw counts from RNA assay
  counts <- NULL
  tryCatch({
    counts <- LayerData(obj, assay = "RNA", layer = "counts")
  }, error = function(e) {
    message("    counts layer not found in RNA, trying GetAssayData...")
  })
  if (is.null(counts)) {
    tryCatch({
      counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
    }, error = function(e) {
      message("    WARNING: Falling back to default assay data layer")
      counts <<- LayerData(obj, layer = "data")
    })
  }

  message("    Writing counts matrix...")
  mtx_tmp <- file.path(export_dir, "counts.mtx")
  writeMM(counts, mtx_tmp)
  system2("gzip", c("-f", mtx_tmp))

  # 2. Gene names
  genes <- rownames(counts)
  write.csv(data.frame(gene = genes), file.path(export_dir, "genes.csv"),
            row.names = FALSE)

  # 3. Cell barcodes
  barcodes <- colnames(obj)
  write.csv(data.frame(barcode = barcodes), file.path(export_dir, "barcodes.csv"),
            row.names = FALSE)

  # 4. Metadata
  meta <- obj@meta.data
  meta$barcode <- rownames(meta)
  write.csv(meta, gzfile(file.path(export_dir, "metadata.csv.gz")),
            row.names = FALSE)
  message("    Saved metadata: ", ncol(meta), " columns")

  # 5. Embeddings
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
  color_by <- intersect(c("study", "condition_harmonized", "cell_class", "coarse_label"),
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
# PROCESS NON-MESENCHYMAL TIER (shared logic for tiered runs)
# ═══════════════════════════════════════════════════════════════════════════

process_nonmes_tier <- function(nonmes_list, mode_prefix, integrate_fn) {
  # Safeguard: exclude studies with <50 non-mesenchymal cells
  study_sizes <- sapply(nonmes_list, ncol)
  keep <- names(study_sizes)[study_sizes >= 50]
  dropped <- names(study_sizes)[study_sizes < 50]

  if (length(dropped) > 0) {
    message("    Dropping ", length(dropped), " studies with <50 non-mes cells: ",
            paste(dropped, " (", study_sizes[dropped], ")", sep = "", collapse = ", "))
    nonmes_list <- nonmes_list[keep]
  }

  if (length(nonmes_list) == 0) {
    message("    WARNING: No studies with >=50 non-mesenchymal cells, skipping tier")
    return(NULL)
  }

  # If <3 studies remain, use simple merge instead of full integration
  if (length(nonmes_list) < 3) {
    message("    Only ", length(nonmes_list), " studies with >=50 cells — using simple merge")
    result <- integrate_simple(nonmes_list, paste0(mode_prefix, "_non_mesenchymal"))
  } else {
    result <- integrate_fn(nonmes_list, paste0(mode_prefix, "_non_mesenchymal"))
  }

  return(result)
}


# ═══════════════════════════════════════════════════════════════════════════
# RUN MODES
# ═══════════════════════════════════════════════════════════════════════════

run_tiered_v5 <- function(seurat_list, force = FALSE) {
  mode_dir <- file.path(EXP_DIR, "tiered_v5")
  dir.create(mode_dir, recursive = TRUE, showWarnings = FALSE)

  mes_rds <- file.path(mode_dir, "mesenchymal.rds")
  nonmes_rds <- file.path(mode_dir, "non_mesenchymal.rds")

  if (file.exists(mes_rds) && file.exists(nonmes_rds) && !force) {
    message("\n=== tiered_v5: outputs exist, skipping (use --force) ===")
    return(invisible(NULL))
  }

  message("\n", paste(rep("=", 60), collapse = ""))
  message("MODE: tiered_v5 — Seurat v5 IntegrateLayers with tier split")
  message(paste(rep("=", 60), collapse = ""))

  tiers <- split_tiers(seurat_list)

  # Mesenchymal tier — skip if already saved
  if (file.exists(mes_rds) && !force) {
    message("\n  [tiered_v5_mesenchymal] Resuming: mesenchymal.rds already exists, skipping")
  } else {
    mes_result <- integrate_v5_cca(tiers$mesenchymal, "tiered_v5_mesenchymal")
    plot_umaps(mes_result, "tiered_v5_mesenchymal")
    export_bridge(mes_result, file.path(mode_dir, "mesenchymal"))
    message("  Saving ", mes_rds, "...")
    saveRDS(mes_result, mes_rds)
    rm(mes_result); gc(verbose = FALSE)
  }

  # Non-mesenchymal tier
  if (file.exists(nonmes_rds) && !force) {
    message("\n  [tiered_v5_non_mesenchymal] Resuming: non_mesenchymal.rds already exists, skipping")
  } else {
    nonmes_result <- process_nonmes_tier(tiers$non_mesenchymal, "tiered_v5", integrate_v5_cca)
    if (!is.null(nonmes_result)) {
      plot_umaps(nonmes_result, "tiered_v5_non_mesenchymal")
      export_bridge(nonmes_result, file.path(mode_dir, "non_mesenchymal"))
      saveRDS(nonmes_result, nonmes_rds)
      rm(nonmes_result); gc(verbose = FALSE)
    }
  }

  message("\n  tiered_v5 complete")
}


run_flat_v4 <- function(seurat_list, force = FALSE) {
  mode_dir <- file.path(EXP_DIR, "flat_v4")
  dir.create(mode_dir, recursive = TRUE, showWarnings = FALSE)

  rds_path <- file.path(mode_dir, "NP.rds")

  if (file.exists(rds_path) && !force) {
    message("\n=== flat_v4: outputs exist, skipping (use --force) ===")
    return(invisible(NULL))
  }

  message("\n", paste(rep("=", 60), collapse = ""))
  message("MODE: flat_v4 — Seurat v4 SCT + CCA (no tier split)")
  message(paste(rep("=", 60), collapse = ""))

  result <- integrate_v4_cca(seurat_list, "flat_v4")
  plot_umaps(result, "flat_v4")
  export_bridge(result, file.path(mode_dir, "all"))
  message("  Saving ", rds_path, "...")
  saveRDS(result, rds_path)
  rm(result); gc(verbose = FALSE)

  message("\n  flat_v4 complete")
}


run_tiered_v4 <- function(seurat_list, force = FALSE) {
  mode_dir <- file.path(EXP_DIR, "tiered_v4")
  dir.create(mode_dir, recursive = TRUE, showWarnings = FALSE)

  mes_rds <- file.path(mode_dir, "mesenchymal.rds")
  nonmes_rds <- file.path(mode_dir, "non_mesenchymal.rds")

  if (file.exists(mes_rds) && file.exists(nonmes_rds) && !force) {
    message("\n=== tiered_v4: outputs exist, skipping (use --force) ===")
    return(invisible(NULL))
  }

  message("\n", paste(rep("=", 60), collapse = ""))
  message("MODE: tiered_v4 — Seurat v4 SCT + CCA with tier split")
  message(paste(rep("=", 60), collapse = ""))

  tiers <- split_tiers(seurat_list)

  # Mesenchymal tier — skip if already saved
  if (file.exists(mes_rds) && !force) {
    message("\n  [tiered_v4_mesenchymal] Resuming: mesenchymal.rds already exists, skipping")
  } else {
    mes_result <- integrate_v4_cca(tiers$mesenchymal, "tiered_v4_mesenchymal")
    plot_umaps(mes_result, "tiered_v4_mesenchymal")
    export_bridge(mes_result, file.path(mode_dir, "mesenchymal"))
    message("  Saving ", mes_rds, "...")
    saveRDS(mes_result, mes_rds)
    rm(mes_result); gc(verbose = FALSE)
  }

  # Non-mesenchymal tier
  if (file.exists(nonmes_rds) && !force) {
    message("\n  [tiered_v4_non_mesenchymal] Resuming: non_mesenchymal.rds already exists, skipping")
  } else {
    nonmes_result <- process_nonmes_tier(tiers$non_mesenchymal, "tiered_v4", integrate_v4_cca)
    if (!is.null(nonmes_result)) {
      plot_umaps(nonmes_result, "tiered_v4_non_mesenchymal")
      export_bridge(nonmes_result, file.path(mode_dir, "non_mesenchymal"))
      saveRDS(nonmes_result, nonmes_rds)
      rm(nonmes_result); gc(verbose = FALSE)
    }
  }

  message("\n  tiered_v4 complete")
}


# ═══════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════

main <- function() {
  message(paste(rep("=", 60), collapse = ""))
  message("NP Integration Quality Experiment")
  message("Seurat version: ", as.character(packageVersion("Seurat")))
  message("HVGs: ", N_HVG, ", Dims: ", N_DIMS)
  message("Mode: ", args$mode)
  message("BLAS: ", sessionInfo()$BLAS)
  message("CPUs available: ", parallel::detectCores())
  message("Started: ", Sys.time())
  message(paste(rep("=", 60), collapse = ""))

  modes <- if (args$mode == "all") c("tiered_v5", "tiered_v4", "flat_v4") else args$mode

  # Load NP cells once (shared across modes)
  seurat_list <- load_np_cells()
  if (length(seurat_list) == 0) {
    stop("No NP cells loaded")
  }

  for (mode in modes) {
    # Deep-copy the list for each mode so mutations don't carry over
    # (SCTransform modifies objects in-place)
    sl_copy <- lapply(seurat_list, function(obj) {
      CreateSeuratObject(
        counts = LayerData(obj, layer = "counts"),
        meta.data = obj@meta.data
      )
    })

    if (mode == "tiered_v5") {
      run_tiered_v5(sl_copy, force = args$force)
    } else if (mode == "flat_v4") {
      run_flat_v4(sl_copy, force = args$force)
    } else if (mode == "tiered_v4") {
      run_tiered_v4(sl_copy, force = args$force)
    }

    rm(sl_copy)
    gc(verbose = FALSE)
  }

  message("\nCompleted: ", Sys.time())
}

main()
