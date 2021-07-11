##检查PCA分群结果, 这里只展示前5个PC,每个PC只显示5个基因；
print(all.big[["pca"]], dims = 1:20, nfeatures = 10)

##展示主成分基因分值
pdf("All_pca_loading.pdf",height = 6, width = 9)
VizDimLoadings(all.big, dims = 1:2, reduction = "pca")
dev.off()
##绘制pca散点图
pdf("All_pca_disease.pdf",height = 6, width = 9)
DimPlot(all.big, reduction = "pca",group.by="condition")
dev.off()

##画第1个或15个主成分的热图；
pdf("All_pca_loadingHeatmap.pdf",height = 30, width = 15)
DimHeatmap(all.big, dims = 1:30, cells = 500, balanced = TRUE)
dev.off()

pdf("All_ElbowPlot.pdf",height = 6, width = 9)
ElbowPlot(all.big, ndims = 50)
dev.off()

all.big <- RunPCA(all.big, verbose = FALSE,group.by="condition")
all.big <- RunUMAP(all.big, dims = 1:22, verbose = FALSE)
all.big <- FindNeighbors(all.big, dims = 1:22, verbose = FALSE)
#先执行不同resolution 下的分群
all.big1 <- FindClusters(object = all.big, resolution = c(seq(.4,1.6,.2)))
library(clustree)
pdf("All_clustree.pdf",height = 6, width = 12)
clustree(all.big1@meta.data, prefix = "SCT_snn_res.")
dev.off()
# Choose optimal resolution
all.big <- FindClusters(all.big, verbose = FALSE) #resolution = 0.5

head(Idents(all.big), 5)
pdf("All_umap_cluster_22PC.pdf",height = 6, width = 9)
DimPlot(all.big, label = TRUE)
dev.off()
pdf("All_umap_patient_22PC.pdf",height = 6, width = 9)
DimPlot(all.big, group.by = "orig.ident",label = TRUE)
dev.off()
pdf("All_umap_condition_22PC.pdf",height = 6, width = 9)
DimPlot(all.big, group.by = "condition",label = TRUE)
dev.off()
write.csv(all.big@meta.data, file = "PBMC_metaData.csv")
write.csv(all.big@meta.data, file = "PBMC_metaData_22PC.csv")

