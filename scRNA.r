library(dplyr)
library(Seurat)
library(patchwork)
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "./1.data/sc_RNA/")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
## 12541 features across 2845 samples within 1 assay
## Active assay: RNA (12541 features, 0 variable features)


pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# Visualize QC metrics as a violin plot
pdf("violen.pdf")
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
pdf("featherScatter.pdf")
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
dev.off()
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

pbmc <- NormalizeData(pbmc)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pdf("features.pdf")
plot1 + plot2
dev.off()
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
PC_ 1
Positive:  Gm26917, Gm42418, Lars2, Trav14-1, AL732506.1
Negative:  Rps15a, Rpl30, Rpl39, Rps7, Rpl32
PC_ 2
Positive:  Il7r, Ltb, Cxcr3, Ly6e, Ctla2a
Negative:  Gzma, Zeb2, Lgals1, Gzmb, Klrg1
PC_ 3
Positive:  Hspa1b, Cx3cr1, Gzma, Rps8, Lgals1
Negative:  Tnfaip3, H3f3b, Nfkbia, Hspa5, Pim1
PC_ 4
Positive:  Hspa1b, Dnajb1, Hsp90aa1, Jun, Hspa1a
Negative:  Vps37b, Dennd4a, Shisa5, S100a6, Cd28
PC_ 5
Positive:  Gm42418, Hist1h1c, Gm26917, Hist1h1e, Lars2
Negative:  Dnaja1, Ctla2a, Hsph1, Hsp90aa1, Ifng

pdf("VizDimLoadings.pdf")
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
dev.off()
pdf("pca_plot.pdf")
DimPlot(pbmc, reduction = "pca")
dev.off()
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
pdf("dimHeatmap.pdf")
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
pdf("JackStrawPlot.pdf")
JackStrawPlot(pbmc, dims = 1:15)
dev.off()
ElbowPlot(pbmc)

pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
## 
## Number of nodes: 2638
## Number of edges: 95965
## 
## Running Louvain algorithm...
## Maximum modularity in 10 random starts: 0.8723
## Number of communities: 9
## Elapsed time: 0 seconds
# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)
## AAACATACAACCAC-1 AAACATTGAGCTAC-1 AAACATTGATCAGC-1 AAACCGTGCTTCCG-1 
##                2                3                2                1 
## AAACCGTGTATGCG-1 
##                6 
## Levels: 0 1 2 3 4 5 6 7 8

# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
pbmc <- RunUMAP(pbmc, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
pdf("Umap.pdf")
DimPlot(pbmc, reduction = "umap")
dev.off()
saveRDS(pbmc, file = "./pbmc_tutorial.rds")

# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 5)
##             p_val avg_log2FC pct.1 pct.2    p_val_adj
## IL32 2.593535e-91  1.2154360 0.949 0.466 3.556774e-87
## LTB  7.994465e-87  1.2828597 0.981 0.644 1.096361e-82
## CD3D 3.922451e-70  0.9359210 0.922 0.433 5.379250e-66
## IL7R 1.130870e-66  1.1776027 0.748 0.327 1.550876e-62
## LDHB 4.082189e-65  0.8837324 0.953 0.614 5.598314e-61

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
##                       p_val avg_log2FC pct.1 pct.2     p_val_adj
## FCGR3A        2.150929e-209   4.267579 0.975 0.039 2.949784e-205
## IFITM3        6.103366e-199   3.877105 0.975 0.048 8.370156e-195
## CFD           8.891428e-198   3.411039 0.938 0.037 1.219370e-193
## CD68          2.374425e-194   3.014535 0.926 0.035 3.256286e-190
## RP11-290F20.3 9.308287e-191   2.722684 0.840 0.016 1.276538e-186
# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
## # A tibble: 18 x 7
## # Groups:   cluster [9]
 1 2.62e- 89      1.34  0.956 0.683 3.28e- 85 0       Ccl4
 2 2.33e- 65      1.18  0.642 0.278 2.92e- 61 0       Gzma
 3 5.37e-194      1.12  1     0.881 6.74e-190 1       Rpl12
 4 3.27e- 89      1.17  0.279 0.027 4.10e- 85 1       Id3
 5 3.11e- 88      1.39  0.988 0.931 3.90e- 84 2       Hspa1b
 6 4.82e- 55      1.41  0.649 0.301 6.05e- 51 2       Gzma
 7 4.21e- 62      0.609 1     0.848 5.28e- 58 3       Rpl29
 8 4.70e- 46      0.629 0.829 0.454 5.89e- 42 3       Ly6e
 9 2.99e-141      4.81  0.971 0.954 3.75e-137 4       Gm42418
10 7.20e- 92      5.86  0.823 0.627 9.03e- 88 4       Gm26917
11 3.98e- 64      1.81  1     1     4.99e- 60 5       mt-Co3
12 4.39e- 22      1.88  0.736 0.752 5.50e- 18 5       Cd28
13 3.13e-  8      0.744 0.815 0.705 3.93e-  4 6       Vim
14 7.04e-  5      0.816 0.411 0.292 8.83e-  1 6       Lgals3


VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)
pdf("featurePlot.pdf")
FeaturePlot(pbmc, features = c("Ccl4", "Gzma", "Rpl12", "Id3", "Hspa1b", "Gzma", "Rpl29", "Ly6e", 
    "Gm42418","Gm26917","mt-Co3","Cd28","Vim","Lgals3"))
dev.off()
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", 
    "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(pbmc,"./pbmc3k_final.rds")

