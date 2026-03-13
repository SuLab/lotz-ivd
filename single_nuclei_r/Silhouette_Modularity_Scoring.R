library(Matrix)
library(dplyr)
library(png)
library(htmltools)
library(glue)
library(Seurat)
library(vctrs)
library(GenomeInfoDb)
library(Signac)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(biovizBase)
library(cowplot)
library(patchwork)
library(JASPAR2020)
library(ggseqlogo)
library(hdf5r)
library(RColorBrewer)
library(SeuratWrappers)
library(clustree)
library(tidyverse)
library(stringr)
library(DESeq2)
library(tidyr)
library(presto)
library(ComplexHeatmap)
library(CellChat)
library(NMF)
library(ggalluvial)
library(scales)
library(cluster)
library(dplyr)
library(ggplot2)
library(igraph)
library(clustree)


#########################################################
# HOW TO CHOOSE THE BEST CLUSTERING RESOLUTION
# Final pipeline: subsampled silhouette + modularity
# handles sparse SNN matrices by converting to igraph
# determine the optimal clustering resolution using silhouette and modularity scores
#########################################################

###############################################
### 1. Parameters
###############################################

#In this code: REJOIN.Synovium.RNA is the name of the Seurat object
#Be sure to process through FindNeighbors function but stop before FindClusters

load("REJOIN.Synovium.RNA.Rdata")

DefaultAssay(object = REJOIN.Synovium.RNA) <- "integrated"

REJOIN.Synovium.RNA <- RunPCA(REJOIN.Synovium.RNA, npcs = 50, nfeatures.print = 10)

REJOIN.Synovium.RNA <- RunUMAP(REJOIN.Synovium.RNA, reduction = "pca", dims = 1:50)

REJOIN.Synovium.RNA <- FindNeighbors(REJOIN.Synovium.RNA, reduction = "pca", dims = 1:50)


### 0. PARAMETERS
res.vec <- seq(0.1, 2.0, by = 0.1)   # resolutions to evaluate
set.seed(123)                        # reproducibility for subsampling

# how many cells to subsample for silhouette calculation
n_sub <- 8000

# PCA dims to try (use up to 50 or fewer if PCA has fewer)
max_pcs <- ncol(Embeddings(REJOIN.Synovium.RNA, "pca"))
pca_dims <- 1:min(50, max_pcs)

#clustering at different resolutions
for (r in res.vec) {
  REJOIN.Synovium.RNA <- FindClusters(REJOIN.Synovium.RNA, resolution = r, verbose = FALSE)
  
  # Store the cluster assignments under stable column names
  colname <- paste0("res.", format(r, nsmall = 1))
  REJOIN.Synovium.RNA[[colname]] <- Idents(REJOIN.Synovium.RNA)
}


#shows how many clusters each resolution produced
a<- sapply(paste0("res.", format(seq(0.1, 2.0, by=0.1), nsmall = 1)), function(x){
  length(unique(REJOIN.Synovium.RNA@meta.data[[x]]))
})

a <- data.frame(a)

### 1. SUBSAMPLE & PREPARE PCA / META
all_cells <- colnames(REJOIN.Synovium.RNA)
n_sub <- min(n_sub, length(all_cells))
sub_cells <- sample(all_cells, n_sub)

# PCA for subsample (rows = cells)
pca_mat <- Embeddings(REJOIN.Synovium.RNA, "pca")[sub_cells, pca_dims]
meta_sub <- REJOIN.Synovium.RNA@meta.data[sub_cells, , drop = FALSE]


### 2. SAFE SILHOUETTE FUNCTION (works on the subsample)
mean_silhouette <- function(labels, embedding) {
  labels <- as.vector(labels)
  if (length(unique(labels)) < 2) return(NA_real_)
  labs_num <- as.numeric(as.factor(labels))
  d <- dist(embedding)
  sil <- silhouette(labs_num, d)
  mean(sil[, 3])
}


### 3. LOCATE SNN GRAPH TO USE FOR MODULARITY
graph_names <- names(REJOIN.Synovium.RNA@graphs)
if (length(graph_names) == 0) {
  stop("No graphs found in REJOIN.Synovium.RNA@graphs. Run FindNeighbors() first.")
}

# prefer integrated_snn if present; otherwise pick the first _snn
if ("integrated_snn" %in% graph_names) {
  snn_graph_name <- "integrated_snn"
} else {
  snn_graphs <- grep("_snn$", graph_names, value = TRUE)
  if (length(snn_graphs) == 0) stop("No SNN graph (ending in '_snn') found in REJOIN.Synovium.RNA@graphs.")
  snn_graph_name <- snn_graphs[1]
}

message("Using SNN graph: ", snn_graph_name)


### 4. PREPARE ADJACENCY MATRIX AND ALIGNMENT
adj <- REJOIN.Synovium.RNA@graphs[[snn_graph_name]]
# check for rownames on the adjacency matrix (these should be cell names)
cell_names_graph <- rownames(adj)
if (is.null(cell_names_graph)) {
  # fallback to Seurat colnames
  cell_names_graph <- colnames(REJOIN.Synovium.RNA)
  warning("SNN adjacency matrix has no rownames. Assuming Seurat colnames order.")
}

# convert sparse adjacency matrix to igraph (weighted undirected)
g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE, diag = FALSE)


### 5. METRICS DATAFRAME
silhouette_scores <- data.frame(
  resolution = res.vec,
  silhouette = NA_real_,
  modularity = NA_real_,
  stringsAsFactors = FALSE
)


### 6. LOOP OVER RESOLUTIONS: compute silhouette (subsample) + modularity (full graph)
for (i in seq_along(res.vec)) {
  res <- res.vec[i]
  colname <- paste0("res.", format(res, nsmall = 1))
  
  # skip if the clustering column is missing
  if (!colname %in% colnames(REJOIN.Synovium.RNA@meta.data)) {
    warning("Skipping missing clustering column: ", colname)
    next
  }
  
  # --- SILHOUETTE on the subsample
  clusters_sub <- meta_sub[[colname]]
  silhouette_scores$silhouette[i] <- mean_silhouette(clusters_sub, pca_mat)
  
  # --- MODULARITY on full graph
  # align cluster labels to the adjacency matrix ordering (very important)
  clusters_full <- REJOIN.Synovium.RNA@meta.data[cell_names_graph, colname, drop = TRUE]
  # convert to numeric factors
  clusters_full_num <- as.numeric(as.factor(clusters_full))
  
  # modularity expects membership vector of length = number of vertices in g
  if (length(clusters_full_num) != vcount(g)) {
    stop("Length mismatch between graph vertices and cluster labels. Check adjacency matrix rownames.")
  }
  
  modularity_score <- modularity(g, clusters_full_num)
  silhouette_scores$modularity[i] <- modularity_score
}


### 7. NORMALIZE & COMBINE SCORES (silhouette + modularity)
score_df <- silhouette_scores %>%
  mutate(
    sil_norm = (silhouette - min(silhouette, na.rm = TRUE)) /
      (max(silhouette, na.rm = TRUE) - min(silhouette, na.rm = TRUE)),
    mod_norm = (modularity - min(modularity, na.rm = TRUE)) /
      (max(modularity, na.rm = TRUE) - min(modularity, na.rm = TRUE)),
    combined = (sil_norm + mod_norm) / 2
  ) %>%
  arrange(desc(combined))

# pick best resolution (highest combined score). If ties/NA, highest silhouette wins.
best_res <- score_df$resolution[1]
message("Best resolution (combined score): ", best_res)

print(score_df)


### 8. PLOT: silhouette vs resolution (best highlighted) and modularity curve
p_sil <- ggplot(silhouette_scores, aes(x = resolution, y = silhouette)) +
  geom_line() +
  geom_point(size = 2) +
  geom_point(data = subset(silhouette_scores, resolution == best_res),
             color = "red", size = 4) +
  theme_classic() +
  labs(title = "Mean Silhouette (subsampled) vs Resolution",
       x = "Resolution", y = "Mean Silhouette (subsampled)")

p_mod <- ggplot(silhouette_scores, aes(x = resolution, y = modularity)) +
  geom_line() +
  geom_point(size = 2) +
  theme_classic() +
  labs(title = "Modularity vs Resolution",
       x = "Resolution", y = "Modularity (full graph)")

print(p_sil)
print(p_mod)


### 9. SET IDS TO BEST RESOLUTION (optional) 
best_colname <- paste0("res.", format(best_res, nsmall = 1))
if (best_colname %in% colnames(REJOIN.Synovium.RNA@meta.data)) {
  Idents(REJOIN.Synovium.RNA) <- REJOIN.Synovium.RNA@meta.data[[best_colname]]
  message("Active identities set to best resolution: ", best_colname)
} else {
  warning("Best resolution column not found in metadata: ", best_colname)
}

#create a data frame containing everything calculated in this code and write the csv file 
df <- score_df[order(score_df$resolution), ]
df <- cbind(a, df)
colnames(df)[1] <- "# of clusters"

df <- df %>% 
  relocate(resolution)

df <- df[order(-df$combined), ]


write.csv(df, file = "Synovium_RNA_Cluster_Silhouette_Modularity.csv")


#create an overlapping single figure

df_plot <- silhouette_scores %>%
  mutate(
    sil_norm = rescale(silhouette),
    mod_norm = rescale(modularity)
  )

p_overlap <- ggplot(df_plot, aes(x = resolution)) +
  geom_line(aes(y = sil_norm), size = 1.2) +
  geom_point(aes(y = sil_norm), size = 2) +
  geom_line(aes(y = mod_norm), linetype = "dashed", size = 1.2) +
  geom_point(aes(y = mod_norm), size = 2) +
  geom_vline(xintercept = best_res, color = "red", linetype = "dashed", linewidth = 1) +
  geom_point(
    data = df_plot[df_plot$resolution == best_res, ], 
    aes(y = sil_norm), 
    color = "red", size = 4
  ) +
  geom_point(
    data = df_plot[df_plot$resolution == best_res, ], 
    aes(y = mod_norm), 
    color = "red", size = 4
  ) +
  scale_x_continuous(breaks = seq(0.1, 2.0, by = 0.1), limits = c(0.1, 2)) +
  theme_classic() +
  labs(
    title = "Silhouette (solid) + Modularity (dashed)\nBest Resolution Highlighted",
    y = "Normalized Score",
    x = "Resolution"
  )

p_overlap


# create an over-under figure

# ---- Add marker to Silhouette ----
p_sil_marked <- p_sil +
  geom_vline(xintercept = best_res, color = "red", linetype = "dashed", linewidth = 1) +
  geom_point(
    data = silhouette_scores[silhouette_scores$resolution == best_res, ],
    aes(x = resolution, y = silhouette),
    color = "red", size = 4
  )

# ---- Add marker to Modularity ----
p_mod_marked <- p_mod +
  geom_vline(xintercept = best_res, color = "red", linetype = "dashed", linewidth = 1) +
  geom_point(
    data = silhouette_scores[silhouette_scores$resolution == best_res, ],
    aes(x = resolution, y = modularity),
    color = "red", size = 4
  )

# ---- Vertical combination ----
p_combined_vert <- p_sil_marked / p_mod_marked

p_combined_vert


############################ Text Summary

#Silhouette Score (Subsampled, PCA Space)
# Measures how well-separated clusters are in PCA space.
# Higher silhouette = clusters are tighter and farther apart.
# Because silhouette uses distance matrices (memory-heavy), it is calculated on a subsample of cells.
# Prevents computational overload while still representing global structure.

#Modularity (Full Graph, SNN)
# Measures how strongly connected cells are within the same cluster vs. between clusters.
# Calculated using the full Seurat SNN graph (integrated_snn if available).
# Higher modularity = clusters form strong graph communities, usually reflecting better biological separation.

#Normalized Combined Score
# For each resolution: combined = silhouette(norm) + modularity(norm)
#                                 -----------------------------------
#                                                 2


# This balances: 1) global PCA structure and 2) local graph community structure
# This combined method is more robust than silhouette or modularity alone

#Interpretation Guides
# Silhouette ↑ Modularity ↑ - Excellent Structure, stable clusters - Choose this resolution
# Silhouette ↑ Modularity ↓ - Over-separation in PCA, poor biological cohesion - Try lower resolution
# Silhouette ↓ Modularity ↑ - Good graph structure but PCA unclear - Check for batch effects 
# Both ↓ - Resolution too low or too high - Avoid
