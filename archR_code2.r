###################################################################################################
#																								  #
#	                              scATAC-seq analysis class-2                                     #
#                                                                                                 #
###################################################################################################
###使用ArchRBrowser绘制Track，可替换，寻找感兴趣的基因or使用database中的Marker genes;
######使用文献中的Marker Gene；
markerGenes4  <- c(
    "PAX5",  
    "MS4A1", 
    "CD14", "CEBPB", "CD3G", 
    "CD4", "KLRB1", 
    "TBX21", "IL7R", "IL6", "IL1B"
)

###解决服务器报错；
pdf("Marker gene.pdf")

p4 <- plotEmbedding(
    ArchRProj = projPBMC_2, 
    colorBy = "GeneScoreMatrix", 
    name = markerGenes4, 
    embedding = "UMAPHarmony",  ######Harmony后的结果；
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL
)

###输出单个基因分布情况；
pdf("CD14.pdf")
p4$CD14
dev.off()

###输出全部Marker Gene分布清苦；
pdf("all genes4.pdf",width = 20)

p42 <- lapply(p4, function(x){
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

do.call(cowplot::plot_grid, c(list(ncol = 3),p42))
dev.off()

###Marker Gene附近Peak分布情况；
markerGenes5 <- c(
    "Tox", "Tcf7", "Pdcd1","Ccr7", "Klrg1"
)
p5 <- plotBrowserTrack(
    ArchRProj = projPBMC_2, 
    groupBy = "Clusters", 
    geneSymbol = markerGenes5, 
    upstream = 50000,
    downstream = 50000
)
p4 <- plotBrowserTrack(
    ArchRProj = projPBMC_2, 
    groupBy = "Clusters", 
    geneSymbol = markerGenes4, 
    upstream = 50000,
    downstream = 50000
)

###CD14基因Peak分布情况；
grid::grid.newpage()
pdf("CD14 track.pdf")
grid::grid.draw(p4$CD14)
dev.off()

###IL6基因Peak分布情况
pdf("IL6 track.pdf")
grid::grid.draw(p4$IL6)
dev.off()

###输出全部marker gene track;
pdf("all gene tracks.pdf",onefile=TRUE)
for(i in 1:length(p5)){
	grid::grid.newpage()
	grid::grid.draw(p5[[i]])	
}
dev.off()

### scATAC-seq细胞和scRNA-seq细胞跨平台连接
###输入自己rds文件路径；
seRNA <- readRDS("../1.data/scRNA-Hematopoiesis-Granja-2019.rds")
seRNA

###查看scRNA-seq信息
colnames(colData(seRNA))

######scATAC-seq与scRNA-seq联合分析；
###无约束整合
projPBMC_2 <- addGeneIntegrationMatrix(
    ArchRProj = projPBMC_2, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = FALSE,
    groupRNA = "BioClassification",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

###保存projPBMC_2
saveArchRProject(ArchRProj = projPBMC_2, outputDirectory = "Save-ProjPBMC_2", load = FALSE)

###查看整合效果
###获取scRNA-seq中的细胞类型信息，并赋颜色信息；
pal <- paletteDiscrete(values = colData(seRNA)$BioClassification)
pdf("Unconstrained intergration.pdf")
p1 <- plotEmbedding(
    projPBMC_2, 
    colorBy = "cellColData", 
    name = "predictedGroup_Un", 
    embedding = "UMAPHarmony",
    pal = pal
)
p1
dev.off()

################
#   约束整合    #
################
###构建混淆矩阵，比较scRNA-seq与scATAC-seq中预测所得的Cluster
cM <- as.matrix(confusionMatrix(projPBMC_2$Clusters, projPBMC_2$predictedGroup_Un))

###筛选出scRNA-seq cluster与scATAC-seq相似的Cluster
preClust <- colnames(cM)[apply(cM,1,which.max)]
###建立对应关系
cbind(preClust, rownames(cM))


 [1,] "12_CD14.Mono.2" "C8"
 [2,] "12_CD14.Mono.2" "C9"
 [3,] "12_CD14.Mono.2" "C7"
 [4,] "12_CD14.Mono.2" "C6"
 [5,] "12_CD14.Mono.2" "C10"
 [6,] "12_CD14.Mono.2" "C5"
 [7,] "05_CMP.LMPP"    "C3"
 [8,] "22_CD4.M"       "C1"
 [9,] "12_CD14.Mono.2" "C4"
[10,] "12_CD14.Mono.2" "C2"


###查看无约束整合中出现的细胞类型；
unique(unique(projPBMC_2$predictedGroup_Un))
###输出
#> unique(unique(projPBMC_2$predictedGroup_Un))
 [1] "12_CD14.Mono.2" "01_HSC"         "22_CD4.M"       "05_CMP.LMPP"
 [5] "24_CD8.CM"      "11_CD14.Mono.1" "21_CD4.N2"      "02_Early.Eryth"
 [9] "20_CD4.N1"      "13_CD16.Mono"   "07_GMP"         "25_NK"
[13] "14_Unk"         "23_CD8.EM"      "17_B"           "08_GMP.Neut"
[17] "19_CD8.N"



######根据先验知识得知，PBMC中存在T细胞、NK细胞。因此在scRNA-seq和scATAC-seq中应该均能够获得这两类细胞的信息。

###在scRNA-seq数据中提取T、NK细胞信息；
###根据上一条命令输出，在scRNA-seq数据中，T细胞和NK细胞的cluster为：19 ~ 25；
###新建变量，存贮这个信息；
cTNK <- paste0(paste0(19:25),collapse="|")
###剩余的cluster为非T、NK细胞。
cNonTNK <- paste0(c(paste0("0", c(1:2, 5 , 7:8)), 11:13, 17), collapse="|")

####在scATAC-seq数据中提取T、NK细胞信息；
###匹配scATAC-seq中T、NK细胞对应的cluster;
clustTNK <- rownames(cM)[grep(cTNK, preClust)]
###匹配scATAC-seq中非T、NK细胞对应的clusters;
clustNonTNK <- rownames(cM)[grep(cNonTNK,preClust)]

####在scRNA-seq数据中提取T、NK细胞信息；
rnaTNK <- colnames(seRNA)[grep(cTNK, colData(seRNA)$BioClassification)]
###筛选scRNA-seq中的非T、NK细胞
rnaNonTNK <- colnames(seRNA)[grep(cNonTNK,colData(seRNA)$BioClassification)]

####建立嵌套list，即将scATAC-seq和scRNA-seq属于同种细胞的簇对应起来；
groupList <- SimpleList(
  TNK = SimpleList(
    ATAC = projPBMC_2$cellNames[projPBMC_2$Clusters %in% clustTNK],
    RNA = rnaTNK 
  ),
  NonTNK = SimpleList(
    ATAC = projPBMC_2$cellNames[projPBMC_2$Clusters %in% clustNonTNK],
    RNA = rnaNonTNK 
  )
)

###将该列表传递给addGeneIntegrationMatrix()中的groupList参数；
projPBMC_2 <- addGeneIntegrationMatrix(
    ArchRProj = projPBMC_2, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = FALSE, 
    groupList = groupList,
    groupRNA = "BioClassification",
    nameCell = "predictedCell_Co",
    nameGroup = "predictedGroup_Co",
    nameScore = "predictedScore_Co"
)

####对比约束整合和无约束整合结果
pal <- paletteDiscrete(values = colData(seRNA)$BioClassification)
pdf("Unconstrained integration.pdf")
p1 <- plotEmbedding(
    projPBMC_2, 
    colorBy = "cellColData", 
    name = "predictedGroup_Un", 
    pal = pal
)
p1
dev.off()

pdf("Constrained integration.pdf")
p2 <- plotEmbedding(
    projPBMC_2, 
    colorBy = "cellColData", 
    name = "predictedGroup_Co", 
    embedding = "UMAPHarmony",
    pal = pal
)
p2
dev.off()

######存储projPBMC_2
saveArchRProject(ArchRProj = projPBMC_2, outputDirectory = "Save-projPBMC_2", load = FALSE)

####为每个scATAC-seq细胞增加拟scRNA-seq谱
projPBMC_3 <- addGeneIntegrationMatrix(
    ArchRProj = projPBMC_2, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = TRUE,
    force= TRUE,
    groupList = groupList,
    groupRNA = "BioClassification",
    nameCell = "predictedCell",
    nameGroup = "predictedGroup",
    nameScore = "predictedScore"
)

###检查GeneIntegrationMatrix矩阵是否已经被添加到Arrow文件中；
getAvailableMatrices(projPBMC_3)

###在新的GeneIntegrationMatrix中，可以比较连接的基因表达量和根据基因得分推断的基因表达量
###需要先确保在项目中加入了填充权重值(impute weights)
projPBMC_3 <- addImputeWeights(projPBMC_3)

###选取进行比较的Marker Gene;
markerGenes  <- c(
    "Tox" ,"Pdcd1", "Tcf7" ,
    "Gzma", "Gzmb", "Klrcl",
    "CEBPB", "TCF3", "TCF12", "SNAI",
    "Ccr9", "S1prl", "Cd69",
    "Foxp1", "Id3", "Fyn", "Ptger4", "Btg1", "Rgs1"
  )

markerGenes  <- c(
    "CD34", #Early Progenitor
    "GATA1", #Erythroid
    "PAX5", "MS4A1", #B-Cell Trajectory
    "CD14", #Monocytes
    "CD3D", "CD8A", "TBX21", "IL7R" #TCells
  )
###基于scRNA-seq获得的基因表达值绘图；
p1 <- plotEmbedding(
    ArchRProj = projPBMC_3, 
    colorBy = "GeneIntegrationMatrix", 
    name = markerGenes, 
    continuousSet = "horizonExtra",
    embedding = "UMAPHarmony",
    imputeWeights = getImputeWeights(projPBMC_3)
)

###基于scATAC-seq获得的基因表达值绘图；
p2 <- plotEmbedding(
    ArchRProj = projPBMC_3, 
    colorBy = "GeneScoreMatrix", 
    continuousSet = "horizonExtra",
    name = markerGenes, 
    embedding = "UMAPHarmony",
    imputeWeights = getImputeWeights(projPBMC_3)
)

###绘制标记基因UMAP图
p1c <- lapply(p1, function(x){
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

p2c <- lapply(p2, function(x){
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

###Marker Gene表达情况输出PDF;
pdf("Marker Gene Value scRNA-seq.pdf")
do.call(cowplot::plot_grid, c(list(ncol = 3), p1c))
dev.off()

pdf("Marker Gene Value scATAC-seq.pdf")
do.call(cowplot::plot_grid, c(list(ncol = 3), p2c))
dev.off()

#######使用scRNA-seq信息标记scATAC-seq聚类，定义细胞类型；
###在scATAC-seq和整合分析得到predictedGroup之间构建一个混淆矩阵，对应两者之间簇的关系；
cM <- confusionMatrix(projPBMC_3$Clusters,projPBMC_3$predictedGroup)
###捕获scATAC-seq cluster
labelOld <- rownames(cM)
labelOld

###鉴定每个scATAC-seq cluster在scRNA-seq中对应的cluster(细胞类型)，对应细胞数目最多的Cluster;
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
labelNew

###查看scRNA-seq中的Cluster名称
unique(seRNA$BioClassification)
 [1] "05_CMP.LMPP"    "08_GMP.Neut"    "01_HSC"         "06_CLP.1"
 [5] "15_CLP.2"       "02_Early.Eryth" "07_GMP"         "09_pDC"
 [9] "04_Early.Baso"  "03_Late.Eryth"  "17_B"           "12_CD14.Mono.2"
[13] "16_Pre.B"       "10_cDC"         "11_CD14.Mono.1" "25_NK"
[17] "21_CD4.N2"      "22_CD4.M"       "23_CD8.EM"      "19_CD8.N"
[21] "24_CD8.CM"      "26_Unk"         "20_CD4.N1"      "14_Unk"
[25] "13_CD16.Mono"   "18_Plasma"

###对scRNA-seq中的cluster重新命名，对其进行简化
remapClust <- c(
    "01_HSC" = "Progenitor",
    "02_Early.Eryth" = "Erythroid",
    "03_Late.Eryth" = "Erythroid",
    "04_Early.Baso" = "Basophil",
    "05_CMP.LMPP" = "Progenitor",
    "06_CLP.1" = "CLP",
    "07_GMP" = "GMP",
    "08_GMP.Neut" = "GMP",
    "09_pDC" = "pDC",
    "10_cDC" = "cDC",
    "11_CD14.Mono.1" = "Mono",
    "12_CD14.Mono.2" = "Mono",
    "13_CD16.Mono" = "Mono",
    "15_CLP.2" = "CLP",
    "16_Pre.B" = "PreB",
    "17_B" = "B",
    "18_Plasma" = "Plasma",
    "19_CD8.N" = "CD8.N",
    "20_CD4.N1" = "CD4.N",
    "21_CD4.N2" = "CD4.N",
    "22_CD4.M" = "CD4.M",
    "23_CD8.EM" = "CD8.EM",
    "24_CD8.CM" = "CD8.CM",
    "25_NK" = "NK"
)

###提取scATAC-seq中出现的Cluster名称
remapClust <- remapClust[names(remapClust) %in% labelNew]

###将原Cluster名称简化
labelNew2 <- mapLabels(labelNew, oldLabels = names(remapClust), newLabels = remapClust)

###利用简化后的scRNA-seq cluster名称对scATAC-seq cluster进行命名，在projPBMC_3中增加重命名后的Clusters信息；
projPBMC_3$Clusters2 <- mapLabels(projPBMC_3$Clusters, newLabels = labelNew2, oldLabels = labelOld)

###绘制UMAP图
pdf("Rename Clusters UMAP.pdf")
p1 <- plotEmbedding(projPBMC_3, colorBy = "cellColData", name = "Clusters2", embedding = "UMAPHarmony",)
p1
dev.off()

###保存projPBMC_3
saveArchRProject(ArchRProj = projPBMC_3, outputDirectory = "Save-ProjPBMC_3", load = FALSE)

#########模拟混池重复，进行Peak分析
###利用cluster定义scATAC-seq分组信息
library(BSgenome.Mmusculus.UCSC.mm10)
projPBMC_4 <- addGroupCoverages(ArchRProj = projPBMC_3, groupBy = "Clusters2")

###使用MACS2鉴定peak
###确认MACS2安装成功
pathToMacs2 <- findMacs2()
###成功信息
## Searching For MACS2..
## Found with $path!

###利用macs2对不同clusters进行peak calling；
pdf("macs2.pdf") ######由于服务器无法调用图形界面，增加该命令解决报错问题；
projPBMC_4 <- addReproduciblePeakSet(
    ArchRProj = projPBMC_4, 
    groupBy = "Clusters2", 
    pathToMacs2 = pathToMacs2
)

###查看peak信息
getPeakSet(projPBMC_4)

###保存projPBMC_4，该project中包含Macs2得到的peak信息；
saveArchRProject(ArchRProj = projPBMC_4, outputDirectory = "Save-ProjPBMC_4", load = FALSE)

###增加Peak矩阵,
projPBMC_5 <- addPeakMatrix(projPBMC_4)
getAvailableMatrices(projPBMC_5)
## [1] "GeneIntegrationMatrix" "GeneScoreMatrix" "PeakMatrix"
## [4] "TileMatrix"

#####使用ArchR鉴定标记Peak
###可先查看project中的细胞类型，及细胞数量；
table(projPBMC_5$Clusters2)	

###通过getMarkerFeatures获得细胞类型特异的peak；
markerPeaks <- getMarkerFeatures(
    ArchRProj = projPBMC_5, 
    useMatrix = "PeakMatrix", 
    groupBy = "Clusters2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerPeaks


###提取细胞类型特异的peak中，感兴趣的部分，用于后续研究；
markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")	###通过设置cutoff值进行筛选；
markerList

###在ArchR中绘制Marker Peaks
###Marker Peak MA和火山图
###MA plot
pdf("MA.pdf")
pma <- markerPlot(seMarker = markerPeaks, name = "Mono", cutOff = "FDR <= 0.1 & Log2FC >= 1 | Log2FC <= -1", plotAs = "MA")
pma
dev.off()

###Volcano plot
pdf("Volcano plot.pdf")
pv <- markerPlot(seMarker = markerPeaks, name = "Mono", cutOff = "Log2FC >= 1 | Log2FC <= -1", plotAs = "Volcano")
pv
dev.off()

pdf("Volcano_CD4.M_plot.pdf")
pv <- markerPlot(seMarker = markerPeaks, name = "CD4.M", cutOff = "Log2FC >= 1 | Log2FC <= -1", plotAs = "Volcano")
pv
dev.off()

###Browser Tracks的Marker Peak
p <- plotBrowserTrack(
    ArchRProj = projPBMC_5, 
    groupBy = "Clusters2", 
    geneSymbol = c("Klrg1"),	###可替换其他marker基因名称；
    features =  getMarkers(markerPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)["Mono"],
    upstream = 50000,
    downstream = 50000
)
grid::grid.newpage()
pdf("Mono Klrg1 track.pdf")
grid::grid.draw(p$Klrg1)
dev.off()

###不同细胞类型之间配对检验
markerTest <- getMarkerFeatures(
  ArchRProj = projPBMC_5, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters2",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Mono",
  bgdGroups = "CD4.M"
)

###MA plot
pma <- markerPlot(seMarker = markerTest, name = "Mono", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
pdf("Mono vs CD4.M cell Markers MA.pdf")
pma
dev.off()

###Volcano plot
pv <- markerPlot(seMarker = markerTest, name = "Mono", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
pdf("Mono vs CD4.M cell Volcano plot.pdf")
pv
dev.off()

###保存projPBMC_5
saveArchRProject(ArchRProj = projPBMC_5, outputDirectory = "Save-projPBMC_5", load = FALSE)

###差异peak中的motif富集
projPBMC_5 <- addMotifAnnotations(ArchRProj = projPBMC_5, motifSet = "cisbp", name = "Motif")

###筛选NK与B cell开放程度存在差异的染色质区域；
motifsUp <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = projPBMC_5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"	###通过改变Log2FC的正负，能够筛选开放程度升高or降低的区域；
  )
save.image("2_509.Rdata")  
###ggplot2绘图
df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

###查看差异peak中的转录因子；
head(df)

###NK细胞中开放程度高的染色质区域绘图；
ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

pdf("TFs in difference open stage chromatin region up.pdf")
ggUp
dev.off()

### ArchR的footprinting分析
###获取Motif的位置
motifPositions <- getPositions(projPBMC_5)
motifPositions

###提取感兴趣的TF
motifs <- c("Smad5_883", " Smad1_863", "Sp3_802"," Wt1_872"," Klf12_236") ###注意，不是所有的TF都会结合；
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs

###根据提示信息判断，是否需要计算分组覆盖深度；
###前序步骤，scATAC-seq模拟混池分析中已经计算过该信息；
projPBMC_5 <- addGroupCoverages(ArchRProj = projPBMC_5, groupBy = "Clusters2")
projPBMC_5 <- addGroupCoverages(ArchRProj = projPBMC_5, groupBy = "Clusters2",force = TRUE)

###获取motifs的footprint信息
seFoot <- getFootprints(
  ArchRProj = projPBMC_5, 
  positions = motifPositions[markerMotifs], 
  groupBy = "Clusters2"
)

###对footprint信息进行可视化；
###Tn5片段化过程中，会产生偏移，可选择不同的校正方法；
###除以Tn5偏好信号
plotFootprints(
  seFoot = seFoot,
  ArchRProj = projPBMC_5, 
  normMethod = "Divide",
  plotName = "Footprints-Divide-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)

###减去Tn5偏好信号
plotFootprints(
  seFoot = seFoot,
  ArchRProj = projPBMC_5, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)

###不进行校正；
plotFootprints(
  seFoot = seFoot,
  ArchRProj = projPBMC_5, 
  normMethod = "None",
  plotName = "Footprints-No-Normalization",
  addDOC = FALSE,
  smoothWindow = 5
)

###TSS上下游特征
###设置范围为TSS上下游各1 kb
seTSS <- getFootprints(
  ArchRProj = projPBMC_5, 
  positions = GRangesList(TSS = getTSS(projPBMC_5)), 
  groupBy = "Clusters2",
  flank = 1000
)

###输出结果
plotFootprints(
  seFoot = seTSS,
  ArchRProj = projPBMC_5, 
  normMethod = "None",
  plotName = "TSS-No-Normalization",
  addDOC = FALSE,
  flank = 1000,
  flankNorm = 100
)

###scATAC-seq共开放分析
###计算scATAC-seq中的共开放信息；
projPBMC_5 <- addCoAccessibility(
    ArchRProj = projPBMC_5,
    reducedDims = "IterativeLSI"
)

###提取project中的共开放信息；
cA <- getCoAccessibility(
    ArchRProj = projPBMC_5,
    corCutOff = 0.5,
    ###resolution：loop的碱基对分辨率。当resolution=1时，输出连接每个peak中心的loop
    resolution = 1, 
    returnLoops = FALSE
)

###在browser track中绘制共开放
###设置需要绘制的Marker Gene；
markerGenes  <- c(
    "CD34", #Early Progenitor
    "GATA1", #Erythroid
    "PAX5", "MS4A1", #B-Cell Trajectory
    "CD14", #Monocytes
    "CD3D", "CD8A", "TBX21", "IL7R", #T Cells
    "CCL5", "IL1B"
)
markerGenes  <- c(
    "Kcnip4" ,"Rln1", "Klra14-ps", "Chn2",
    "Gzma", "Gzmb", "Klrcl", "Cx3cr1",
    "CEBPB", "TCF3", "Sp5", " Ptpn5",
    "Ccr7", "Otx1", "Agap1", "Tnf",
    "Klra3", "Gja1", "Bcl6", "Klra21"
)


p <- plotBrowserTrack(
    ArchRProj = projPBMC_5, 
    groupBy = "Clusters2", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000,
    loops = getCoAccessibility(projPBMC_5)
)

plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", 
    ArchRProj = projPBMC_5, 
    addDOC = FALSE, width = 5, height = 5
)
plotPDF(plotList = p, 
    name = "Plot2-Tracks-Marker-Genes-with-CoAccessibility.pdf", 
    ArchRProj = projPBMC_5, 
    addDOC = FALSE, width = 5, height = 5
)

###scATAC-seq peak关联gene分析；
###整合scRNA-seq信息关联peak & gene;
projPBMC_5 <- addPeak2GeneLinks(
    ArchRProj = projPBMC_5,
    reducedDims = "IterativeLSI"
)

###提取peak & gene links
p2g <- getPeak2GeneLinks(
    ArchRProj = projPBMC_5,
    corCutOff = 0.45,
    resolution = 1,
    returnLoops = FALSE
)

###在browser track中展示peak & gene links
p <- plotBrowserTrack(
    ArchRProj = projPBMC_5, 
    groupBy = "Clusters2", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000,
    loops = getPeak2GeneLinks(projPBMC_5)
)

###输出Marker Gene track
plotPDF(plotList = p, 
    name = "Plot2-Tracks-Marker-Genes-with-Peak2GeneLinks.pdf", 
    ArchRProj = projPBMC_5, 
    addDOC = FALSE, width = 5, height = 5
)



##轨迹分析
p1 <- plotEmbedding(ArchRProj = projPBMC_5, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = projPBMC_5, colorBy = "cellColData", name = "Clusters2", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")

trajectory <- c("C8","C9","C10")

projPBMC_5 <- addTrajectory(
    ArchRProj = projPBMC_5, 
    name = "MyeloidU", 
    groupBy = "Clusters",
    trajectory = trajectory, 
    embedding = "UMAP", 
    force = TRUE
)

head(projPBMC_5$MyeloidU[!is.na(projPBMC_5$MyeloidU)])
# [1] 46.07479 53.75027 44.82834 43.18828 47.49617 43.21015

p <- plotTrajectory(projPBMC_5, trajectory = "MyeloidU", colorBy = "cellColData", name = "MyeloidU")
p[[1]]

plotPDF(p, name = "Plot-MyeloidU-Traj-UMAP.pdf", ArchRProj = projPBMC_5, addDOC = FALSE, width = 5, height = 5)
save.image("2_730.Rdata")
projPBMC_5 <- addImputeWeights(projPBMC_5)

p1 <- plotTrajectory(projPBMC_5, trajectory = "MyeloidU", colorBy = "GeneScoreMatrix", name = "Zeb1", continuousSet = "horizonExtra")

p2 <- plotTrajectory(projPBMC_5, trajectory = "MyeloidU", colorBy = "GeneIntegrationMatrix", name = "Zeb1", continuousSet = "blueYellow")
ggAlignPlots(p1[[1]], p1[[2]], type = "h")

##拟时间热图
trajMM  <- getTrajectory(ArchRProj = projPBMC_5, name = "MyeloidU", useMatrix = "MotifMatrix", log2Norm = FALSE)
p1 <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))

trajGSM <- getTrajectory(ArchRProj = projPBMC_5, name = "MyeloidU", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- trajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))

trajGIM <- getTrajectory(ArchRProj = projPBMC_5, name = "MyeloidU", useMatrix = "GeneIntegrationMatrix", log2Norm = FALSE)
p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"))

trajPM  <- getTrajectory(ArchRProj = projPBMC_5, name = "MyeloidU", useMatrix = "PeakMatrix", log2Norm = TRUE)
p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"))

plotPDF(p1, p2, p3, p4, name = "Plot-MyeloidU-Traj-Heatmaps.pdf", ArchRProj = projPBMC_5, addDOC = FALSE, width = 6, height = 8)

corGSM_MM <- correlateTrajectories(trajGSM, trajMM)

trajGSM2 <- trajGSM[corGSM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGSM_MM[[1]]$name2, ]

trajCombined <- trajGSM2
assay(trajCombined) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))

combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))

ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)

ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)

ht1 + ht2


