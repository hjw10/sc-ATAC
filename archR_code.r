###################################################################################################
#																								  #
#	                         Don’t do this at this time, please!!!                                #
#                                                                                                 #
###################################################################################################
###ArchR install
###Don’t do this at this time, please!!!
conda create -n ArchR
conda activate ArchR
conda install r-base=4.0.3

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())

library(ArchR)
ArchR::installExtraPackages()
addArchRGenome("hg38")

###################################################################################################
#																								  #
#	                              scATAC-seq analysis class-1                                     #
#                                Start your analysis from here!                                   #
#																								  #
###################################################################################################

###激活ArchR_2023环境
conda activate ArchR_2023

###进入R
R

###Data analysis
library(ArchR)
addArchRGenome("mm10") # hg38, mm9, mm10
addArchRThreads(threads = 8) ###You can modify this value;

setwd("C:/Users/penguin/Desktop/学习资料/高通量测序02/大作业") ###Set up your own workspace;

###输入文件
inputFiles <- getInputFiles(paths = "./3.downsteam")

###生成arrowfiles
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 4, #The minimum numeric transcription start site (TSS) enrichment score required for a cell to pass filtering for use in downstream analyses.
  filterFrags = 1000, #The minimum number of mapped ATAC-seq fragments required per cell to pass filtering for use in downstream analyses.
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

ArrowFiles

###ArchR推断scATAC-seq doublets
doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.The number of cells neighboring a simulated doublet to be considered as putative doublets.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.The name of the dimensionality reduction method to be used for k-nearest neighbors calculation. Possible values are "UMAP" or "LSI".
    LSIMethod = 1 #The list of parameters to pass to the IterativeLSI() function. See IterativeLSI()
)

###创建一个ArchRProject
projPBMC_1 <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "PBMC",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

###检查当前的ArchRProject中存放了哪些矩阵数据，这些数据可以在下游分析中使用。
getAvailableMatrices(projPBMC_1)
# "GeneScoreMatrix" "TileMatrix"

load("71.Rdata")
setwd("./out")
library(ArchR)
###重新绘制质控信息和TSS富集得分的对比图
df <- getCellColData(projPBMC_1, select = c("log10(nFrags)", "TSSEnrichment"))
pdf("TSS-vs-Frags.pdf")
p <- ggPoint(
    x = df[,1],
    y = df[,2],
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
    ) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
p
dev.off()

###为ArchRProject中的样本统计绘图
###根据TSS富集得分为每个样本绘制山脊图
pdf("ridge plot.pdf")
p1 <- plotGroups(
    ArchRProj = projPBMC_1,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "TSSEnrichment",
    plotAs = "ridges"
   )
p1
dev.off()

###根据TSS富集得分为每个样本绘制小提琴图
pdf("violin plot.pdf")
p2 <- plotGroups(
    ArchRProj = projPBMC_1,
    groupBy = "Sample",
    colorBy = "cellColData",
    name = "TSSEnrichment",
    plotAs = "violin",
    alpha = 0.4,
    addBoxPlot = TRUE
   )
p2
dev.off()

###绘制样本的TSS富集谱和Fragment大小分布
pdf("Fragment size.pdf")
p1 <- plotFragmentSizes(ArchRProj = projPBMC_1)
p1
dev.off()

pdf("TSS enrichment.pdf")
p2 <- plotTSSEnrichment(ArchRProj = projPBMC_1)
p2
dev.off()

###过滤doublets
projPBMC_2 <- filterDoublets(projPBMC_1)

###ArchR降维分析
projPBMC_2 <- addIterativeLSI(
    ArchRProj = projPBMC_2,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( 
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30
)

###矫正批次效应
projPBMC_2 <- addHarmony(
    ArchRProj = projPBMC_2,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample"
)

###使用Seurat的FindClusters()函数进行聚类
projPBMC_2 <- addClusters(
    input = projPBMC_2,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8
)

###使用$符号来获取聚类信息，输出每个细胞对应的cluster
head(projPBMC_2$Clusters)
###统计每个cluster的细胞数
table(projPBMC_2$Clusters)

###统计每个样本在不同的cluster的分布情况
cM <- confusionMatrix(paste0(projPBMC_2$Clusters), paste0(projPBMC_2$Sample))
cM
library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
pdf("Cluster heatmap.pdf")
p <- pheatmap::pheatmap(
    mat = as.matrix(cM), 
    color = paletteContinuous("whiteBlue"), 
    border_color = "black"
)
p
dev.off()

###Uniform Manifold Approximation and Projection (UMAP)
projPBMC_2 <- addUMAP(
    ArchRProj = projPBMC_2, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine" ###A number that determines how distance is computed in the reducedDims to compute a UMAP. This argument is passed to metric in uwot::umap().
)

###绘制UMAP图
###根据sample上色
p1 <- plotEmbedding(ArchRProj = projPBMC_2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
pdf("UMAP sampe.pdf")
p1
dev.off()
###根据Cluster上色
p2 <- plotEmbedding(ArchRProj = projPBMC_2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
pdf("UMAP cluster.pdf")
p2
dev.off()
plotPDF(p1,p2, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = projPBMC_2, addDOC = FALSE, width = 5, height = 5)


###t-SNE图
projPBMC_3 <- addTSNE(
    ArchRProj = projPBMC_2, 
    reducedDims = "IterativeLSI", 
    name = "TSNE", 
    perplexity = 30
)
p1 <- plotEmbedding(ArchRProj = projPBMC_3, colorBy = "cellColData", name = "Sample", embedding = "TSNE")
pdf("t-SNE sampe.pdf")
p1
dev.off()
p2 <- plotEmbedding(ArchRProj = projPBMC_3, colorBy = "cellColData", name = "Clusters", embedding = "TSNE")
pdf("t-SNE cluster.pdf")
p2
dev.off()
ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-TSNE-Sample-Clusters.pdf", ArchRProj = projPBMC_3, addDOC = FALSE, width = 5, height = 5)

###Harmony后降维
projPBMC_2 <- addUMAP(
    ArchRProj = projPBMC_2, 
    reducedDims = "Harmony", 
    name = "UMAPHarmony", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)
pdf("UMAP Harmony.pdf")
p3 <- plotEmbedding(ArchRProj = projPBMC_2, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
p4 <- plotEmbedding(ArchRProj = projPBMC_2, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
ggAlignPlots(p3, p4, type = "h")
dev.off()
pdf("Plot-UMAP2Harmony-Sample-Clusters.pdf")
ggAlignPlots(p1, p2, type = "h")
ggAlignPlots(p3, p4, type = "h")
dev.off()


###鉴定标记基因
markersGS <- getMarkerFeatures(
    ArchRProj = projPBMC_2, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25") ###cutoff值，可根据需要修改；
markerList$C2

###可视化特定的基因，可替换，寻找感兴趣的基因or使用database中的Marker genes;
markerGenes  <- c(
    "Ccr9",  
    "Rgs9", 
    "Me1", "Trerf1", "Prdm1", "Trps1",
    "Snap47", "Cd69", "S1pr1",
    "Cxcr4", "Ccr2"
)

markerGenes2  <- c(
    "CD34", #Early Progenitor
    "GATA1", #Erythroid
    "PAX5", "MS4A1", "EBF1", "MME", #B-Cell Trajectory
    "CD14", "CEBPB", "MPO", #Monocytes
    "IRF8", 
    "CD3D", "CD8A", "TBX21", "IL7R" #TCells
  )


heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = projPBMC_2, addDOC = FALSE)

###Or use C1 cluster expression genes (inferred)
markerGenes3  <- c(
    "Tox", "Tox2", "Pdcd1", "Lag3"
  )


p <- plotEmbedding(
    ArchRProj = projPBMC_2, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes, 
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL
)

p2 <- plotEmbedding(
    ArchRProj = projPBMC_2, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes2, 
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL
)

p3 <- plotEmbedding(
    ArchRProj = projPBMC_2, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes3, 
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL
)


###CD34，可替换，markerGenes中任意基因均可；
pdf("Ccr9.pdf")
p$Ccr9
dev.off()


###绘制全部基因
p12 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})
pdf("all genes.pdf")
do.call(cowplot::plot_grid, c(list(ncol = 3),p12))
dev.off()

p22 <- lapply(p2, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})
pdf("all genes2.pdf")
do.call(cowplot::plot_grid, c(list(ncol = 3),p22))
dev.off()

p32 <- lapply(p3, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})
pdf("all genes3.pdf")
do.call(cowplot::plot_grid, c(list(ncol = 3),p32))
dev.off()

