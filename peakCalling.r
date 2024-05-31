projPBMC_3 <- addGroupCoverages(ArchRProj = projPBMC_2, groupBy = "Clusters")
projPBMC_3 <- addReproduciblePeakSet(
    ArchRProj = projPBMC_3, 
    groupBy = "Clusters", 
    pathToMacs2 = pathToMacs2
)
getPeakSet(projPBMC_3)
###使用saveArchRProject()函数保存原始图像。这ArchRProject包含MACS2衍生的合并峰集。
saveArchRProject(ArchRProj = projPBMC_3, outputDirectory = "Save-ProjPBMC_2_3", load = FALSE)
projPBMC_4 <- addPeakMatrix(projPBMC_3)

###正在使用的细胞类型projPBMC_4及其相对比例
table(projPBMC_4$Clusters)
###标记识别峰
markersPeaks <- getMarkerFeatures(
    ArchRProj = projPBMC_4, 
    useMatrix = "PeakMatrix", 
    groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markersPeaks
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList
markerList$C3

library(devtools)
install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)
pdf("peak heatmap.pdf")
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
dev.off()
plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = projPBMC_4, addDOC = FALSE)
###绘制MA图
pdf("MA plot.pdf")
pma <- markerPlot(seMarker = markersPeaks, name = "C3", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "MA")
pma
dev.off()
###绘制火山图
pdf("Volcano plot.pdf")
pv <- markerPlot(seMarker = markersPeaks, name = "C3", cutOff = "Log2FC >= 1 | Log2FC <= -1", plotAs = "Volcano")
pv
dev.off()
###Browser Tracks的Marker Peak
p <- plotBrowserTrack(
    ArchRProj = projPBMC_4, 
    groupBy = "Clusters", 
    geneSymbol = c("Klra3"), 
    features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)["C1"],
    upstream = 50000,
    downstream = 50000
)
grid::grid.newpage()
pdf("C3 Klra3 track.pdf")
grid::grid.draw(p$Klra3)
dev.off()



