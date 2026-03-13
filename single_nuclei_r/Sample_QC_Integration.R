
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
library(ggplot2)
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
library(UCell)

####################### synNorm3
counts <- Read10X_h5(filename = "synnorm3_filtered_feature_bc_matrix.h5")
Synovium_Norm3_RNA <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# filter out low quality cells
#DEFAULT
Synovium_Norm3_RNA <- subset(
  x = Synovium_Norm3_RNA,
  subset =
    nCount_RNA < 25000 &
    nCount_RNA > 1000)

Synovium_Norm3_RNA@meta.data$tissue <- "Synovium"
Synovium_Norm3_RNA@meta.data$condition <- "Normal"
Synovium_Norm3_RNA@meta.data$sample <- "Normal_Synovium_3"
Synovium_Norm3_RNA@meta.data$technology <- "10X_Multiome"
Synovium_Norm3_RNA@meta.data$sex <- "Male"


DefaultAssay(Synovium_Norm3_RNA) <- "RNA"
Synovium_Norm3_RNA[["percent.mt"]] <- PercentageFeatureSet(Synovium_Norm3_RNA, pattern = "^MT-")
head(Synovium_Norm3_RNA[[]])
VlnPlot(Synovium_Norm3_RNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Synovium_Norm3_RNA <- subset(Synovium_Norm3_RNA, subset = percent.mt < 5)
Synovium_Norm3_RNA <- SCTransform(Synovium_Norm3_RNA, vars.to.regress = c("percent.mt"))
Synovium_Norm3_RNA <- RunPCA(Synovium_Norm3_RNA)
Synovium_Norm3_RNA<- RunUMAP(
  object = Synovium_Norm3_RNA,
  reduction.name="rna.umap",
  reduction="pca",
  assay = "SCT",
  verbose = TRUE,
  dims=1:50
)

Synovium_Norm3_RNA

save(Synovium_Norm3_RNA, file = "Synovium_Norm3_RNA_Object.RData")


####################### synNorm5
counts <- Read10X_h5(filename = "synnorm5_filtered_feature_bc_matrix.h5")
Synovium_Norm5_RNA <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# filter out low quality cells
#DEFAULT
Synovium_Norm5_RNA <- subset(
  x = Synovium_Norm5_RNA,
  subset = 
    nCount_RNA < 25000 &
    nCount_RNA > 1000)

Synovium_Norm5_RNA


Synovium_Norm5_RNA@meta.data$tissue <- "Synovium"
Synovium_Norm5_RNA@meta.data$condition <- "Normal"
Synovium_Norm5_RNA@meta.data$sample <- "Normal_Synovium_5"
Synovium_Norm5_RNA@meta.data$technology <- "10X_Multiome"
Synovium_Norm5_RNA@meta.data$sex <- "Female"


DefaultAssay(Synovium_Norm5_RNA) <- "RNA"
Synovium_Norm5_RNA[["percent.mt"]] <- PercentageFeatureSet(Synovium_Norm5_RNA, pattern = "^MT-")
head(Synovium_Norm5_RNA[[]])
VlnPlot(Synovium_Norm5_RNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Synovium_Norm5_RNA <- subset(Synovium_Norm5_RNA, subset = percent.mt < 5)
Synovium_Norm5_RNA <- SCTransform(Synovium_Norm5_RNA, vars.to.regress = c("percent.mt"))
Synovium_Norm5_RNA <- RunPCA(Synovium_Norm5_RNA)
Synovium_Norm5_RNA<- RunUMAP(
  object = Synovium_Norm5_RNA,
  reduction.name="rna.umap",
  reduction="pca",
  assay = "SCT",
  verbose = TRUE,
  dims=1:50
)

Synovium_Norm5_RNA

save(Synovium_Norm5_RNA, file = "Synovium_Norm5_RNA_Object.RData")

####################### synNorm6
counts <- Read10X_h5(filename = "synnorm6_filtered_feature_bc_matrix.h5")
Synovium_Norm6_RNA <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# filter out low quality cells
#DEFAULT
Synovium_Norm6_RNA <- subset(
  x = Synovium_Norm6_RNA,
  subset =
    nCount_RNA < 25000 &
    nCount_RNA > 1000)

Synovium_Norm6_RNA

Synovium_Norm6_RNA@meta.data$tissue <- "Synovium"
Synovium_Norm6_RNA@meta.data$condition <- "Normal"
Synovium_Norm6_RNA@meta.data$sample <- "Normal_Synovium_6"
Synovium_Norm6_RNA@meta.data$technology <- "10X_Multiome"
Synovium_Norm6_RNA@meta.data$sex <- "Female"


DefaultAssay(Synovium_Norm6_RNA) <- "RNA"
Synovium_Norm6_RNA[["percent.mt"]] <- PercentageFeatureSet(Synovium_Norm6_RNA, pattern = "^MT-")
head(Synovium_Norm6_RNA[[]])
VlnPlot(Synovium_Norm6_RNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Synovium_Norm6_RNA <- subset(Synovium_Norm6_RNA, subset = percent.mt < 5)
Synovium_Norm6_RNA <- SCTransform(Synovium_Norm6_RNA, vars.to.regress = c("percent.mt"))
Synovium_Norm6_RNA <- RunPCA(Synovium_Norm6_RNA)
Synovium_Norm6_RNA<- RunUMAP(
  object = Synovium_Norm6_RNA,
  reduction.name="rna.umap",
  reduction="pca",
  assay = "SCT",
  verbose = TRUE,
  dims=1:50
)

Synovium_Norm6_RNA


save(Synovium_Norm6_RNA, file = "Synovium_Norm6_RNA_Object.RData")


####################### synNorm10
counts <- Read10X_h5(filename = "synnorm10_filtered_feature_bc_matrix.h5")
Synovium_Norm10_RNA <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# filter out low quality cells
#DEFAULT
Synovium_Norm10_RNA <- subset(
  x = Synovium_Norm10_RNA,
  subset = 
    nCount_RNA < 25000 &
    nCount_RNA > 1000)

Synovium_Norm10_RNA

Synovium_Norm10_RNA@meta.data$tissue <- "Synovium"
Synovium_Norm10_RNA@meta.data$condition <- "Normal"
Synovium_Norm10_RNA@meta.data$sample <- "Normal_Synovium_10"
Synovium_Norm10_RNA@meta.data$technology <- "10X_Multiome"
Synovium_Norm10_RNA@meta.data$sex <- "Male"


DefaultAssay(Synovium_Norm10_RNA) <- "RNA"
Synovium_Norm10_RNA[["percent.mt"]] <- PercentageFeatureSet(Synovium_Norm10_RNA, pattern = "^MT-")
head(Synovium_Norm10_RNA[[]])
VlnPlot(Synovium_Norm10_RNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Synovium_Norm10_RNA <- subset(Synovium_Norm10_RNA, subset = percent.mt < 5)
Synovium_Norm10_RNA <- SCTransform(Synovium_Norm10_RNA, vars.to.regress = c("percent.mt"))
Synovium_Norm10_RNA <- RunPCA(Synovium_Norm10_RNA)
Synovium_Norm10_RNA<- RunUMAP(
  object = Synovium_Norm10_RNA,
  reduction.name="rna.umap",
  reduction="pca",
  assay = "SCT",
  verbose = TRUE,
  dims=1:50
)

Synovium_Norm10_RNA

save(Synovium_Norm10_RNA, file = "Synovium_Norm10_RNA_Object.RData")

####################### synNorm11
counts <- Read10X_h5(filename = "synnorm11_filtered_feature_bc_matrix.h5")
Synovium_Norm11_RNA <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)


# filter out low quality cells
#DEFAULT
Synovium_Norm11_RNA <- subset(
  x = Synovium_Norm11_RNA,
  subset = 
    nCount_RNA < 25000 &
    nCount_RNA > 1000)

Synovium_Norm11_RNA

Synovium_Norm11_RNA@meta.data$tissue <- "Synovium"
Synovium_Norm11_RNA@meta.data$condition <- "Normal"
Synovium_Norm11_RNA@meta.data$sample <- "Normal_Synovium_11"
Synovium_Norm11_RNA@meta.data$technology <- "10X_Multiome"
Synovium_Norm11_RNA@meta.data$sex <- "Male"


DefaultAssay(Synovium_Norm11_RNA) <- "RNA"
Synovium_Norm11_RNA[["percent.mt"]] <- PercentageFeatureSet(Synovium_Norm11_RNA, pattern = "^MT-")
head(Synovium_Norm11_RNA[[]])
VlnPlot(Synovium_Norm11_RNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Synovium_Norm11_RNA <- subset(Synovium_Norm11_RNA, subset = percent.mt < 5)
Synovium_Norm11_RNA <- SCTransform(Synovium_Norm11_RNA, vars.to.regress = c("percent.mt"))
Synovium_Norm11_RNA <- RunPCA(Synovium_Norm11_RNA)
Synovium_Norm11_RNA<- RunUMAP(
  object = Synovium_Norm11_RNA,
  reduction.name="rna.umap",
  reduction="pca",
  assay = "SCT",
  verbose = TRUE,
  dims=1:50
)

Synovium_Norm11_RNA

save(Synovium_Norm11_RNA, file = "Synovium_Norm11_RNA_Object.RData")

####################### synNorm14
counts <- Read10X_h5(filename = "synnorm14_filtered_feature_bc_matrix.h5")
Synovium_Norm14_RNA <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)


# filter out low quality cells
#DEFAULT
Synovium_Norm14_RNA <- subset(
  x = Synovium_Norm14_RNA,
  subset = 
    nCount_RNA < 25000 &
    nCount_RNA > 1000)

Synovium_Norm14_RNA

Synovium_Norm14_RNA@meta.data$tissue <- "Synovium"
Synovium_Norm14_RNA@meta.data$condition <- "Normal"
Synovium_Norm14_RNA@meta.data$sample <- "Normal_Synovium_14"
Synovium_Norm14_RNA@meta.data$technology <- "10X_Multiome"
Synovium_Norm14_RNA@meta.data$sex <- "Female"


DefaultAssay(Synovium_Norm14_RNA) <- "RNA"
Synovium_Norm14_RNA[["percent.mt"]] <- PercentageFeatureSet(Synovium_Norm14_RNA, pattern = "^MT-")
head(Synovium_Norm14_RNA[[]])
VlnPlot(Synovium_Norm14_RNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Synovium_Norm14_RNA <- subset(Synovium_Norm14_RNA, subset = percent.mt < 5)
Synovium_Norm14_RNA <- SCTransform(Synovium_Norm14_RNA, vars.to.regress = c("percent.mt"))
Synovium_Norm14_RNA <- RunPCA(Synovium_Norm14_RNA)
Synovium_Norm14_RNA<- RunUMAP(
  object = Synovium_Norm14_RNA,
  reduction.name="rna.umap",
  reduction="pca",
  assay = "SCT",
  verbose = TRUE,
  dims=1:50
)

Synovium_Norm14_RNA

save(Synovium_Norm14_RNA, file = "Synovium_Norm14_RNA_Object.RData")


####################### synAged13
counts <- Read10X_h5(filename = "synaged13_filtered_feature_bc_matrix.h5")
Synovium_Aged13_RNA <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# filter out low quality cells
#DEFAULT
Synovium_Aged13_RNA <- subset(
  x = Synovium_Aged13_RNA,
  subset = 
    nCount_RNA < 25000 &
    nCount_RNA > 1000)

Synovium_Aged13_RNA

Synovium_Aged13_RNA@meta.data$tissue <- "Synovium"
Synovium_Aged13_RNA@meta.data$condition <- "Aged"
Synovium_Aged13_RNA@meta.data$sample <- "Aged_Synovium_13"
Synovium_Aged13_RNA@meta.data$technology <- "10X_Multiome"
Synovium_Aged13_RNA@meta.data$sex <- "Female"

DefaultAssay(Synovium_Aged13_RNA) <- "RNA"
Synovium_Aged13_RNA[["percent.mt"]] <- PercentageFeatureSet(Synovium_Aged13_RNA, pattern = "^MT-")
head(Synovium_Aged13_RNA[[]])
VlnPlot(Synovium_Aged13_RNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Synovium_Aged13_RNA <- subset(Synovium_Aged13_RNA, subset = percent.mt < 5)
Synovium_Aged13_RNA <- SCTransform(Synovium_Aged13_RNA, vars.to.regress = c("percent.mt"))
Synovium_Aged13_RNA <- RunPCA(Synovium_Aged13_RNA)
Synovium_Aged13_RNA<- RunUMAP(
  object = Synovium_Aged13_RNA,
  reduction.name="rna.umap",
  reduction="pca",
  assay = "SCT",
  verbose = TRUE,
  dims=1:50
)

Synovium_Aged13_RNA

save(Synovium_Aged13_RNA, file = "Synovium_Aged13_RNA_Object.RData")

####################### synAged14
counts <- Read10X_h5(filename = "synaged14_filtered_feature_bc_matrix.h5")
Synovium_Aged14_RNA <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)


# filter out low quality cells
#DEFAULT
Synovium_Aged14_RNA <- subset(
  x = Synovium_Aged14_RNA,
  subset = 
    nCount_RNA < 25000 &
    nCount_RNA > 1000)

Synovium_Aged14_RNA

Synovium_Aged14_RNA@meta.data$tissue <- "Synovium"
Synovium_Aged14_RNA@meta.data$condition <- "Aged"
Synovium_Aged14_RNA@meta.data$sample <- "Aged_Synovium_14"
Synovium_Aged14_RNA@meta.data$technology <- "10X_Multiome"
Synovium_Aged14_RNA@meta.data$sex <- "Male"

DefaultAssay(Synovium_Aged14_RNA) <- "RNA"
Synovium_Aged14_RNA[["percent.mt"]] <- PercentageFeatureSet(Synovium_Aged14_RNA, pattern = "^MT-")
head(Synovium_Aged14_RNA[[]])
VlnPlot(Synovium_Aged14_RNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Synovium_Aged14_RNA <- subset(Synovium_Aged14_RNA, subset = percent.mt < 5)
Synovium_Aged14_RNA <- SCTransform(Synovium_Aged14_RNA, vars.to.regress = c("percent.mt"))
Synovium_Aged14_RNA <- RunPCA(Synovium_Aged14_RNA)
Synovium_Aged14_RNA<- RunUMAP(
  object = Synovium_Aged14_RNA,
  reduction.name="rna.umap",
  reduction="pca",
  assay = "SCT",
  verbose = TRUE,
  dims=1:50
)

Synovium_Aged14_RNA

save(Synovium_Aged14_RNA, file = "Synovium_Aged14_RNA_Object.RData")

####################### synOA2
counts <- Read10X_h5(filename = "synoa2_filtered_feature_bc_matrix.h5")
Synovium_OA2_RNA <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# filter out low quality cells
#DEFAULT
Synovium_OA2_RNA <- subset(
  x = Synovium_OA2_RNA,
  subset = 
    nCount_RNA < 25000 &
    nCount_RNA > 1000)

Synovium_OA2_RNA

Synovium_OA2_RNA@meta.data$tissue <- "Synovium"
Synovium_OA2_RNA@meta.data$condition <- "OA"
Synovium_OA2_RNA@meta.data$sample <- "OA_Synovium_2"
Synovium_OA2_RNA@meta.data$technology <- "10X_Multiome"
Synovium_OA2_RNA@meta.data$sex <- "Male"

DefaultAssay(Synovium_OA2_RNA) <- "RNA"
Synovium_OA2_RNA[["percent.mt"]] <- PercentageFeatureSet(Synovium_OA2_RNA, pattern = "^MT-")
head(Synovium_OA2_RNA[[]])
VlnPlot(Synovium_OA2_RNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Synovium_OA2_RNA <- subset(Synovium_OA2_RNA, subset = percent.mt < 5)
Synovium_OA2_RNA <- SCTransform(Synovium_OA2_RNA, vars.to.regress = c("percent.mt"))
Synovium_OA2_RNA <- RunPCA(Synovium_OA2_RNA)
Synovium_OA2_RNA<- RunUMAP(
  object = Synovium_OA2_RNA,
  reduction.name="rna.umap",
  reduction="pca",
  assay = "SCT",
  verbose = TRUE,
  dims=1:50
)

Synovium_OA2_RNA

save(Synovium_OA2_RNA, file = "Synovium_OA2_RNA_Object.RData")

####################### synOA3
counts <- Read10X_h5(filename = "synoa3_filtered_feature_bc_matrix.h5")
Synovium_OA3_RNA <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)


# filter out low quality cells
#DEFAULT
Synovium_OA3_RNA <- subset(
  x = Synovium_OA3_RNA,
  subset = 
    nCount_RNA < 25000 &
    nCount_RNA > 1000)

Synovium_OA3_RNA

Synovium_OA3_RNA@meta.data$tissue <- "Synovium"
Synovium_OA3_RNA@meta.data$condition <- "OA"
Synovium_OA3_RNA@meta.data$sample <- "OA_Synovium_3"
Synovium_OA3_RNA@meta.data$technology <- "10X_Multiome"
Synovium_OA3_RNA@meta.data$sex <- "Male"

DefaultAssay(Synovium_OA3_RNA) <- "RNA"
Synovium_OA3_RNA[["percent.mt"]] <- PercentageFeatureSet(Synovium_OA3_RNA, pattern = "^MT-")
head(Synovium_OA3_RNA[[]])
VlnPlot(Synovium_OA3_RNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Synovium_OA3_RNA <- subset(Synovium_OA3_RNA, subset = percent.mt < 5)
Synovium_OA3_RNA <- SCTransform(Synovium_OA3_RNA, vars.to.regress = c("percent.mt"))
Synovium_OA3_RNA <- RunPCA(Synovium_OA3_RNA)
Synovium_OA3_RNA<- RunUMAP(
  object = Synovium_OA3_RNA,
  reduction.name="rna.umap",
  reduction="pca",
  assay = "SCT",
  verbose = TRUE,
  dims=1:50
)

Synovium_OA3_RNA

save(Synovium_OA3_RNA, file = "Synovium_OA3_RNA_Object.RData")

####################### synTKA1
counts <- Read10X_h5(filename = "syntka1_filtered_feature_bc_matrix.h5")
Synovium_TKA1_RNA <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# filter out low quality cells
#DEFAULT
Synovium_TKA1_RNA <- subset(
  x = Synovium_TKA1_RNA,
  subset = 
    nCount_RNA < 25000 &
    nCount_RNA > 1000)

Synovium_TKA1_RNA

Synovium_TKA1_RNA@meta.data$tissue <- "Synovium"
Synovium_TKA1_RNA@meta.data$condition <- "TKA"
Synovium_TKA1_RNA@meta.data$sample <- "TKA_Synovium_1"
Synovium_TKA1_RNA@meta.data$technology <- "10X_Multiome"
Synovium_TKA1_RNA@meta.data$sex <- "Female"


DefaultAssay(Synovium_TKA1_RNA) <- "RNA"
Synovium_TKA1_RNA[["percent.mt"]] <- PercentageFeatureSet(Synovium_TKA1_RNA, pattern = "^MT-")
head(Synovium_TKA1_RNA[[]])
VlnPlot(Synovium_TKA1_RNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Synovium_TKA1_RNA <- subset(Synovium_TKA1_RNA, subset = percent.mt < 5)
Synovium_TKA1_RNA <- SCTransform(Synovium_TKA1_RNA, vars.to.regress = c("percent.mt"))
Synovium_TKA1_RNA <- RunPCA(Synovium_TKA1_RNA)
Synovium_TKA1_RNA<- RunUMAP(
  object = Synovium_TKA1_RNA,
  reduction.name="rna.umap",
  reduction="pca",
  assay = "SCT",
  verbose = TRUE,
  dims=1:50
)

Synovium_TKA1_RNA

save(Synovium_TKA1_RNA, file = "Synovium_TKA1_RNA_Object.RData")

####################### synTKA2
counts <- Read10X_h5(filename = "syntka2_filtered_feature_bc_matrix.h5")
Synovium_TKA2_RNA <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)


# filter out low quality cells
#DEFAULT
Synovium_TKA2_RNA <- subset(
  x = Synovium_TKA2_RNA,
  subset = 
    nCount_RNA < 25000 &
    nCount_RNA > 1000)

Synovium_TKA2_RNA

Synovium_TKA2_RNA@meta.data$tissue <- "Synovium"
Synovium_TKA2_RNA@meta.data$condition <- "TKA"
Synovium_TKA2_RNA@meta.data$sample <- "TKA_Synovium_2"
Synovium_TKA2_RNA@meta.data$technology <- "10X_Multiome"
Synovium_TKA2_RNA@meta.data$sex <- "Male"


DefaultAssay(Synovium_TKA2_RNA) <- "RNA"
Synovium_TKA2_RNA[["percent.mt"]] <- PercentageFeatureSet(Synovium_TKA2_RNA, pattern = "^MT-")
head(Synovium_TKA2_RNA[[]])
VlnPlot(Synovium_TKA2_RNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Synovium_TKA2_RNA <- subset(Synovium_TKA2_RNA, subset = percent.mt < 5)
Synovium_TKA2_RNA <- SCTransform(Synovium_TKA2_RNA, vars.to.regress = c("percent.mt"))
Synovium_TKA2_RNA <- RunPCA(Synovium_TKA2_RNA)
Synovium_TKA2_RNA<- RunUMAP(
  object = Synovium_TKA2_RNA,
  reduction.name="rna.umap",
  reduction="pca",
  assay = "SCT",
  verbose = TRUE,
  dims=1:50
)

Synovium_TKA2_RNA

save(Synovium_TKA2_RNA, file = "Synovium_TKA2_RNA_Object.RData")


####################### synTKA3
counts <- Read10X_h5(filename = "syntka3_filtered_feature_bc_matrix.h5")
Synovium_TKA3_RNA <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)


# filter out low quality cells
#DEFAULT
Synovium_TKA3_RNA <- subset(
  x = Synovium_TKA3_RNA,
  subset = 
    nCount_RNA < 25000 &
    nCount_RNA > 1000)

Synovium_TKA3_RNA

Synovium_TKA3_RNA@meta.data$tissue <- "Synovium"
Synovium_TKA3_RNA@meta.data$condition <- "TKA"
Synovium_TKA3_RNA@meta.data$sample <- "TKA_Synovium_3"
Synovium_TKA3_RNA@meta.data$technology <- "10X_Multiome"
Synovium_TKA3_RNA@meta.data$sex <- "Male"


DefaultAssay(Synovium_TKA3_RNA) <- "RNA"
Synovium_TKA3_RNA[["percent.mt"]] <- PercentageFeatureSet(Synovium_TKA3_RNA, pattern = "^MT-")
head(Synovium_TKA3_RNA[[]])
VlnPlot(Synovium_TKA3_RNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Synovium_TKA3_RNA <- subset(Synovium_TKA3_RNA, subset = percent.mt < 5)
Synovium_TKA3_RNA <- SCTransform(Synovium_TKA3_RNA, vars.to.regress = c("percent.mt"))
Synovium_TKA3_RNA <- RunPCA(Synovium_TKA3_RNA)
Synovium_TKA3_RNA<- RunUMAP(
  object = Synovium_TKA3_RNA,
  reduction.name="rna.umap",
  reduction="pca",
  assay = "SCT",
  verbose = TRUE,
  dims=1:50
)

Synovium_TKA3_RNA


save(Synovium_TKA3_RNA, file = "Synovium_TKA3_RNA_Object.RData")

####################### synTKA5
counts <- Read10X_h5(filename = "syntka5_filtered_feature_bc_matrix.h5")
Synovium_TKA5_RNA <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

# filter out low quality cells
#DEFAULT
Synovium_TKA5_RNA <- subset(
  x = Synovium_TKA5_RNA,
  subset = 
    nCount_RNA < 25000 &
    nCount_RNA > 1000)

Synovium_TKA5_RNA

Synovium_TKA5_RNA@meta.data$tissue <- "Synovium"
Synovium_TKA5_RNA@meta.data$condition <- "TKA"
Synovium_TKA5_RNA@meta.data$sample <- "TKA_Synovium_5"
Synovium_TKA5_RNA@meta.data$technology <- "10X_Multiome"
Synovium_TKA5_RNA@meta.data$sex <- "Female"


DefaultAssay(Synovium_TKA5_RNA) <- "RNA"
Synovium_TKA5_RNA[["percent.mt"]] <- PercentageFeatureSet(Synovium_TKA5_RNA, pattern = "^MT-")
head(Synovium_TKA5_RNA[[]])
VlnPlot(Synovium_TKA5_RNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Synovium_TKA5_RNA <- subset(Synovium_TKA5_RNA, subset = percent.mt < 5)
Synovium_TKA5_RNA <- SCTransform(Synovium_TKA5_RNA, vars.to.regress = c("percent.mt"))
Synovium_TKA5_RNA <- RunPCA(Synovium_TKA5_RNA)
Synovium_TKA5_RNA<- RunUMAP(
  object = Synovium_TKA5_RNA,
  reduction.name="rna.umap",
  reduction="pca",
  assay = "SCT",
  verbose = TRUE,
  dims=1:50
)

Synovium_TKA5_RNA


save(Synovium_TKA5_RNA, file = "Synovium_TKA5_RNA_Object.RData")

####################### synTKA11
counts <- Read10X_h5(filename = "syntka11_filtered_feature_bc_matrix.h5")
Synovium_TKA11_RNA <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)


# filter out low quality cells
#DEFAULT
Synovium_TKA11_RNA <- subset(
  x = Synovium_TKA11_RNA,
  subset = 
    nCount_RNA < 25000 &
    nCount_RNA > 1000)

Synovium_TKA11_RNA

Synovium_TKA11_RNA@meta.data$tissue <- "Synovium"
Synovium_TKA11_RNA@meta.data$condition <- "TKA"
Synovium_TKA11_RNA@meta.data$sample <- "TKA_Synovium_11"
Synovium_TKA11_RNA@meta.data$technology <- "10X_Multiome"
Synovium_TKA11_RNA@meta.data$sex <- "Female"


DefaultAssay(Synovium_TKA11_RNA) <- "RNA"
Synovium_TKA11_RNA[["percent.mt"]] <- PercentageFeatureSet(Synovium_TKA11_RNA, pattern = "^MT-")
head(Synovium_TKA11_RNA[[]])
VlnPlot(Synovium_TKA11_RNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Synovium_TKA11_RNA <- subset(Synovium_TKA11_RNA, subset = percent.mt < 5)
Synovium_TKA11_RNA <- SCTransform(Synovium_TKA11_RNA, vars.to.regress = c("percent.mt"))
Synovium_TKA11_RNA <- RunPCA(Synovium_TKA11_RNA)
Synovium_TKA11_RNA<- RunUMAP(
  object = Synovium_TKA11_RNA,
  reduction.name="rna.umap",
  reduction="pca",
  assay = "SCT",
  verbose = TRUE,
  dims=1:50
)

Synovium_TKA11_RNA

save(Synovium_TKA11_RNA, file = "Synovium_TKA11_RNA_Object.RData")

####################### synTKA12
counts <- Read10X_h5(filename = "syntka12_filtered_feature_bc_matrix.h5")
Synovium_TKA12_RNA <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)


# filter out low quality cells
#DEFAULT
Synovium_TKA12_RNA <- subset(
  x = Synovium_TKA12_RNA,
  subset = 
    nCount_RNA < 25000 &
    nCount_RNA > 1000)

Synovium_TKA12_RNA

Synovium_TKA12_RNA@meta.data$tissue <- "Synovium"
Synovium_TKA12_RNA@meta.data$condition <- "TKA"
Synovium_TKA12_RNA@meta.data$sample <- "TKA_Synovium_12"
Synovium_TKA12_RNA@meta.data$technology <- "10X_Multiome"
Synovium_TKA12_RNA@meta.data$sex <- "Female"


DefaultAssay(Synovium_TKA12_RNA) <- "RNA"
Synovium_TKA12_RNA[["percent.mt"]] <- PercentageFeatureSet(Synovium_TKA12_RNA, pattern = "^MT-")
head(Synovium_TKA12_RNA[[]])
VlnPlot(Synovium_TKA12_RNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Synovium_TKA12_RNA <- subset(Synovium_TKA12_RNA, subset = percent.mt < 5)
Synovium_TKA12_RNA <- SCTransform(Synovium_TKA12_RNA, vars.to.regress = c("percent.mt"))
Synovium_TKA12_RNA <- RunPCA(Synovium_TKA12_RNA)
Synovium_TKA12_RNA<- RunUMAP(
  object = Synovium_TKA12_RNA,
  reduction.name="rna.umap",
  reduction="pca",
  assay = "SCT",
  verbose = TRUE,
  dims=1:50
)

Synovium_TKA12_RNA


save(Synovium_TKA12_RNA, file = "Synovium_TKA12_RNA_Object.RData")


#################################################################################

#integration of the RNA only object


load("Synovium_Norm3_RNA_Object.RData")
load("Synovium_Norm5_RNA_Object.RData")
load("Synovium_Norm6_RNA_Object.RData")
load("Synovium_Norm10_RNA_Object.RData")
load("Synovium_Norm11_RNA_Object.RData")
load("Synovium_Norm14_RNA_Object.RData")
load("Synovium_Aged13_RNA_Object.RData")
load("Synovium_Aged14_RNA_Object.RData")
load("Synovium_OA2_RNA_Object.RData")
load("Synovium_OA3_RNA_Object.RData")
load("Synovium_TKA1_RNA_Object.RData")
load("Synovium_TKA2_RNA_Object.RData")
load("Synovium_TKA3_RNA_Object.RData")
load("Synovium_TKA5_RNA_Object.RData")
load("Synovium_TKA11_RNA_Object.RData")
load("Synovium_TKA12_RNA_Object.RData")


DefaultAssay(Synovium_Norm3_RNA) <- "SCT"
DefaultAssay(Synovium_Norm5_RNA) <- "SCT"
DefaultAssay(Synovium_Norm6_RNA) <- "SCT"
DefaultAssay(Synovium_Norm10_RNA) <- "SCT"
DefaultAssay(Synovium_Norm11_RNA) <- "SCT"
DefaultAssay(Synovium_Norm14_RNA) <- "SCT"
DefaultAssay(Synovium_Aged13_RNA) <- "SCT"
DefaultAssay(Synovium_Aged14_RNA) <- "SCT"
DefaultAssay(Synovium_OA2_RNA) <- "SCT"
DefaultAssay(Synovium_OA3_RNA) <- "SCT"
DefaultAssay(Synovium_TKA1_RNA) <- "SCT"
DefaultAssay(Synovium_TKA2_RNA) <- "SCT"
DefaultAssay(Synovium_TKA3_RNA) <- "SCT"
DefaultAssay(Synovium_TKA5_RNA) <- "SCT"
DefaultAssay(Synovium_TKA11_RNA) <- "SCT"
DefaultAssay(Synovium_TKA12_RNA) <- "SCT"


data.list.RNA = list(Synovium_Norm3_RNA, Synovium_Norm5_RNA, Synovium_Norm6_RNA,
                     Synovium_Norm10_RNA, Synovium_Norm11_RNA, Synovium_Norm14_RNA,
                     Synovium_Aged13_RNA, Synovium_Aged14_RNA,
                     Synovium_OA2_RNA, Synovium_OA3_RNA,
                     Synovium_TKA1_RNA, Synovium_TKA2_RNA, Synovium_TKA3_RNA,
                     Synovium_TKA5_RNA, Synovium_TKA11_RNA, Synovium_TKA12_RNA)


save(data.list.RNA, file = "RNA_REJOIN_Normal_Disease_Synovium_Seurat_Objects.Rdata")
load("RNA_REJOIN_Normal_Disease_Synovium_Seurat_Objects.Rdata")

Synovium.features.RNA <- SelectIntegrationFeatures(object.list = data.list.RNA, nfeatures = 3000)

data.list.RNA <- PrepSCTIntegration(object.list = data.list.RNA, anchor.features = Synovium.features.RNA)


#Canonical Correlation Analysis (CCA) for batch effect correction
REJOIN_Synovium_Anchors_RNA <- FindIntegrationAnchors(object.list = data.list.RNA, dims = 1:50, normalization.method = "SCT", 
                                                      anchor.features = Synovium.features.RNA) 
save(REJOIN_Synovium_Anchors_RNA, file = "REJOIN_Synovium_Anchors_RNA.Rdata")
load("REJOIN_Synovium_Anchors_RNA.Rdata")

REJOIN.Synovium.RNA <- IntegrateData(anchorset = REJOIN_Synovium_Anchors_RNA, dims = 1:50, normalization.method = "SCT")
save(REJOIN.Synovium.RNA, file = "REJOIN.Synovium.RNA.Rdata")
load(file = "REJOIN.Synovium.RNA.Rdata")

#################################################################################
#process the integrated RNA only object

load("REJOIN.Synovium.RNA.Rdata")

DefaultAssay(object = REJOIN.Synovium.RNA) <- "integrated"

REJOIN.Synovium.RNA <- RunPCA(REJOIN.Synovium.RNA, npcs = 50, nfeatures.print = 10)

REJOIN.Synovium.RNA <- RunUMAP(REJOIN.Synovium.RNA, reduction = "pca", dims = 1:50)

REJOIN.Synovium.RNA <- FindNeighbors(REJOIN.Synovium.RNA, reduction = "pca", dims = 1:50)

#prior to FindClusters, perform silhouette and modularity scoring to determine optimal resolution, confirm with clustree
REJOIN.Synovium.RNA <- FindClusters(REJOIN.Synovium.RNA, resolution = 0.6, verbose = FALSE, algorithm = 1)

DimPlot(REJOIN.Synovium.RNA, group.by = "seurat_clusters", label = TRUE, repel = TRUE)

###################################################################################

#markers
DefaultAssay(REJOIN.Synovium.RNA) <- "SCT"

REJOIN.Synovium.RNA <- PrepSCTFindMarkers(REJOIN.Synovium.RNA)

markers <- FindAllMarkers(REJOIN.Synovium.RNA, logfc.threshold = 0.25, min.pct = 0.1, 
                          only.pos = TRUE, verbose = TRUE, recorrect_umi = FALSE)

markers <- markers[markers$p_val_adj < 0.05, ]

write.csv(markers, file = "Synovium_RNA_Markers_Harmonized.csv")

