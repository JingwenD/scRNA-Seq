load("~/Project/10_SingleCell/outs/PBMC_merge.RData")
library(Seurat)
plot1<-DimPlot(subset(all.big, subset = condition=='PsO'),reduction = "umap",label = TRUE,pt.size = 1.5)
plot2<-DimPlot(subset(all.big, subset = condition=='HC'),reduction = "umap",label = TRUE,pt.size = 1.5)
pdf("umap_PsOvsHC.pdf",height = 7, width = 15)
CombinePlots(plots = list(plot1, plot2))
dev.off()



##比较2组数据中特定cluster的差异基因
all.big@meta.data$sample_type <- paste(all.big@meta.data$condition, all.big@meta.data$SCT_snn_res.0.8, sep = "_")
marker<-FindMarkers(all.big, group.by="sample_type",ident.1 = "PsO_0", ident.2 = "HC_0", logfc.threshold = 0, min.pct = 0.1)
marker$cluster<-0

cluster<- 28
for (i in 1:cluster){
  ident.1= paste0("PsO_",i)
  ident.2= paste0("HC_",i)
  markerx<-FindMarkers(all.big, group.by="sample_type",ident.1 = ident.1, ident.2 = ident.2, logfc.threshold = 0, min.pct = 0.1)
  markerx$cluster<-i
  marker<-rbind(marker,markerx)
}

write.csv(marker,file = "DEGs_PsOvsHC.csv")

##细胞周期归类
all.big <- CellCycleScoring(object = all.big, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
head(x = all.big@meta.data)
pdf("CellCycle_all.pdf",height = 6, width = 9)
DimPlot(all.big,reduction = "umap",label = TRUE,group.by="Phase",pt.size = 1.5, cols =c("#FF8000FF","blue","#4D9221"))
dev.off()

plot1<-DimPlot(subset(all.big, subset = condition=='PsO'),reduction = "umap",group.by="Phase",label = TRUE,pt.size = 1.5)
plot2<-DimPlot(subset(all.big, subset = condition=='HC'),reduction = "umap",group.by="Phase",label = TRUE,pt.size = 1.5)
pdf("CellCycle_PsOvsHC.pdf",height = 7, width = 15)
CombinePlots(plots = list(plot1, plot2))
dev.off()

save(all.big, file="PBMC_merge.RData")

##寻找每一cluster的marker，其中ident.1参数设置待分析的细胞类别，min.pct表示该基因表达数目占该类细胞总数的比例
merged.markers <- FindAllMarkers(all.big, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)#如何选方向
merged.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

##存储marker
write.csv(merged.markers,file="merged_allmarker.csv")


##各种绘图
##绘制Marker 基因的umap图
FeaturePlot(all.big, features = c("CD4", "CD8A","CD3E","IL2RA"),cols = c("gray", "red"))

##绘制Marker 基因的小提琴图
VlnPlot(all.big, features = c("RPL38", "RPL39"),split.by = "condition")

##绘制分cluster的热图
top10 <- merged.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(all.big, features = top10$gene) + NoLegend()

##各种绘图
##绘制Marker 基因的umap图
FeaturePlot(all.big, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"),cols = c("gray", "red"))

##绘制Marker 基因的小提琴图
VlnPlot(all.big, features = c("MS4A1", "CD79A"))

VlnPlot(all.big, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

VlnPlot(all.big, features = c("CD4", "CD8A"))



##绘制气泡图
DotPlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"),cols = c("blue", "red"))

##绘制RidgePlot
RidgePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14"), ncol = 2)

############取出一部分cluster用来亚分群
sub_pbmc<-subset(pbmc, idents = '0')

sub_pbmc <- FindVariableFeatures(sub_pbmc, selection.method = "vst", nfeatures = 2000)

sub_pbmc <- ScaleData(sub_pbmc, features = all.genes)
sub_pbmc <- ScaleData(sub_pbmc, vars.to.regress = "percent.mt")

sub_pbmc <- RunPCA(sub_pbmc, features = VariableFeatures(object = pbmc))

sub_pbmc <- FindNeighbors(sub_pbmc, dims = 1:10)

sub_pbmc <- FindClusters(sub_pbmc, resolution = 0.5)

sub_pbmc <- RunTSNE(sub_pbmc, dims = 1:10)

DimPlot(sub_pbmc, reduction = "tsne",label = TRUE, pt.size = 1.5)


##每个cluster按照marker规定细胞类型
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 1.5) + NoLegend()

write.table(Idents(pbmc),"D:/cell-type.xls",sep="\t")

##存储细胞归类结果
saveRDS(pbmc, file = "D:/pbmc3k_final.rds")

##按照细胞类型绘制Marker基因的小提琴图
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

##按照细胞类型绘制分cluster的热图
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

##按照细胞类型绘制气泡图
DotPlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"),cols = c("blue", "red"))

##按照细胞类型绘制RidgePlot
RidgePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14"), ncol = 2)