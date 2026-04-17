#!/usr/bin/env Rscript
# One-off repair for v4 bridge exports: counts.mtx.gz was written from the
# first study's RNA layer only (Seurat v5 split-assay quirk). This loads each
# saved .rds, calls JoinLayers on RNA, and rewrites counts.mtx.gz + genes.csv.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

BASE <- "/home/ubuntu/lotz-ivd"
EXP  <- file.path(BASE, "data/integrated/np_experiment")

repair <- function(rds_path, export_dir) {
  t0 <- Sys.time()
  message("\n=== ", rds_path, " ===")
  message("  [", format(Sys.time(), "%H:%M:%S"), "] Loading rds (",
          round(file.info(rds_path)$size / 1e9, 1), " GB on disk)...")
  obj <- readRDS(rds_path)
  message("  [", format(Sys.time(), "%H:%M:%S"), "] Loaded. ncol=",
          ncol(obj), " nrow=", nrow(obj),
          " assays=", paste(names(obj@assays), collapse = ","))

  message("  [", format(Sys.time(), "%H:%M:%S"), "] JoinLayers(RNA)...")
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

  counts <- LayerData(obj, assay = "RNA", layer = "counts")
  message("  [", format(Sys.time(), "%H:%M:%S"), "] Joined counts: ",
          nrow(counts), " genes x ", ncol(counts), " cells")

  mtx_tmp <- file.path(export_dir, "counts.mtx")
  message("  [", format(Sys.time(), "%H:%M:%S"), "] Writing ", mtx_tmp, "...")
  writeMM(counts, mtx_tmp)
  message("  [", format(Sys.time(), "%H:%M:%S"), "] gzip...")
  system2("gzip", c("-f", mtx_tmp))

  write.csv(data.frame(gene = rownames(counts)),
            file.path(export_dir, "genes.csv"), row.names = FALSE)

  rm(obj, counts); gc(verbose = FALSE)

  message("  [", format(Sys.time(), "%H:%M:%S"), "] done in ",
          round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
}

message("Starting repair at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))

repair(file.path(EXP, "flat_v4", "NP.rds"),
       file.path(EXP, "flat_v4", "all"))

repair(file.path(EXP, "tiered_v4", "mesenchymal.rds"),
       file.path(EXP, "tiered_v4", "mesenchymal"))

repair(file.path(EXP, "tiered_v4", "non_mesenchymal.rds"),
       file.path(EXP, "tiered_v4", "non_mesenchymal"))

message("\nRepair complete at ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
