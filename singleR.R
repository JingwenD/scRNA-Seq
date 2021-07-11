library(SingleR)
library(celldex)
library(Seurat)
# T Cell only
ref <- DatabaseImmuneCellExpressionData() 
# Decent immune cell granularity, does not contain finer monocyte and dendritic cell subtypes.
ref <- BlueprintEncodeData()
#Most of the labels refer to blood subpopulations but cell types from other tissues are also available.
ref <- HumanPrimaryCellAtlasData()
# Best fit for PBMCs, came from Singaporean-Chinese individuals 
ref <- MonacoImmuneData() 
library("scRNAseq")
ref <- StoeckiusHashingData(mode='human')
colData(ref)
ref <- ref[,!is.na(ref$label)]

ref$label.main
ref$label.fine


load("~/Project/10_SingleCell/outs/PBMC_merge.RData")
#condition.list<-SplitObject(all.big, split.by = "condition")
sce_for_SingleR = GetAssayData(all.big, slot="data")
labelmain <- SingleR(test=sce_for_SingleR, ref=ref, labels=ref$label.main, assay.type.test = 1)
labelmain_log <- SingleR(test=sce_for_SingleR, ref=ref, labels=ref$label.main, assay.type.test = "logcounts")
table(labelmain$labels)
labelfine <- SingleR(test=sce_for_SingleR, ref=ref, labels=ref$label.fine, assay.type.test = 1)
table(labelfine$labels)

pdf("CellScoreHeatmap.pdf",height = 6, width = 10)
plotScoreHeatmap(labelmain)
dev.off()
pdf("CellScoreDistribution.pdf",height = 8, width = 15)
plotDeltaDistribution(labelmain, ncol = 5)
dev.off()
summary(is.na(labelmain$pruned.labels))

pdf("CellScoreHeatmap1.pdf",height = 7, width = 13)
plotScoreHeatmap(labelfine)
dev.off()
pdf("CellScoreDistribution1.pdf",height = 20, width = 15)
plotDeltaDistribution(labelfine, ncol = 5)
dev.off()
summary(is.na(labelfine$pruned.labels))

labelmain.markers <- metadata(labelmain)$de.genes
labelfine.markers <- metadata(labelfine)$de.genes
all.big[["cell.type1"]] <- labelmain$labels
all.big[["cell.type2"]] <- labelfine$labels
save(all.big, file="PBMC_merge.RData")

pdf("All_umap_labelmain.pdf",height = 6, width = 9)
DimPlot(all.big, group.by = "cell.type1",label = TRUE)
dev.off()
plot1<-DimPlot(subset(all.big, subset = condition=='PsO'),reduction = "umap",group.by="cell.type1",label = TRUE,pt.size = 0.5)
plot2<-DimPlot(subset(all.big, subset = condition=='HC'),reduction = "umap",group.by="cell.type1",label = TRUE,pt.size = 0.5)
pdf("labelmain_PsOvsHC.pdf",height = 7, width = 19)
CombinePlots(plots = list(plot1, plot2))
dev.off()

pdf("All_umap_labelfine.pdf",height = 9, width = 15)
DimPlot(all.big, group.by = "cell.type2",label = TRUE)
dev.off()
plot1<-DimPlot(subset(all.big, subset = condition=='PsO'),reduction = "umap",group.by="cell.type2",label = TRUE,pt.size = 0.5)
plot2<-DimPlot(subset(all.big, subset = condition=='HC'),reduction = "umap",group.by="cell.type2",label = TRUE,pt.size = 0.5)
pdf("labelfine_PsOvsHC.pdf",height = 9, width = 30)
CombinePlots(plots = list(plot1, plot2))
dev.off()

##比较2组数据中特定cell.type的差异基因
all.big@meta.data$sample_type <- paste(all.big@meta.data$condition, all.big@meta.data$cell.type2, sep = "_")
celltype<-unique(all.big@meta.data$cell.type2)
marker<-FindMarkers(all.big, group.by="sample_type",ident.1 = "PsO_Th1 cells", ident.2 = "HC_Th1 cells", logfc.threshold = 0, min.pct = 0.1)
marker$cluster<-"Th1 cells"


for (i in celltype[2:length(celltype)]){
  ident.1= paste0("PsO_",i)
  ident.2= paste0("HC_",i)
  markerx<-FindMarkers(all.big, group.by="sample_type",ident.1 = ident.1, ident.2 = ident.2, logfc.threshold = 0, min.pct = 0.1)
  markerx$cluster<-i
  marker<-rbind(marker,markerx)
}
write.csv(marker,file = "DEGs_PsOvsHC_celltype2.csv")

###############Cluster-level annotation################
labelfine_cluster <- SingleR(test=sce_for_SingleR, ref=ref, labels=ref$label.fine, clusters=all.big@meta.data$seurat_clusters)
labelfine_cluster$labels
