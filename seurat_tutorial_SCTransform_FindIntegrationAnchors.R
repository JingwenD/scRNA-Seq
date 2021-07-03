###该代码同时消除样本内及样本间的异质性，样本内用SCTransform，样本间用锚点方法。

##载入seurat包
rm(list = ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
options(future.globals.maxSize = 4000 * 1024^2)


##读入pbmc数据
pbmc.data <- Read10X(data.dir = "D:/BC11_TUMOR1/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "TUMOR1", min.cells = 3, min.features = 500)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & percent.mt < 5)
pbmc<-SCTransform(pbmc, verbose = FALSE)

##读入pbmc1数据
pbmc1.data <- Read10X(data.dir = "D:/BC11_TUMOR2/")
pbmc1 <- CreateSeuratObject(counts = pbmc1.data, project = "TUMOR2", min.cells = 3, min.features = 500)
pbmc1[["percent.mt"]] <- PercentageFeatureSet(pbmc1, pattern = "^MT-")
pbmc1 <- subset(pbmc1, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & percent.mt < 5)
pbmc1<-SCTransform(pbmc1, verbose = FALSE)

features <- SelectIntegrationFeatures(object.list = list(pbmc, pbmc1), nfeatures = 3000)
list <- PrepSCTIntegration(object.list = list(pbmc, pbmc1), anchor.features = features, verbose = FALSE)


##然后使用FindIntegrationAnchors函数识别锚点，该函数将Seurat对象列表作为输入，并使用这些锚点将两个数据集集成在一起IntegrateData。
merged.anchors <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT", anchor.features = features, verbose = FALSE)
merged <- IntegrateData(anchorset = merged.anchors, normalization.method = "SCT",verbose = FALSE)


##线性降维PCA；
merged <- RunPCA(object = merged, verbose = FALSE)

##Seurat提供了几种非线性降维的方法进行数据可视化（在低维空间把相似的细胞聚在一起）,比如UMAP和t-SNE,运行UMAP需要先安装'umap-learn'包,这里不做介绍。
merged <- RunTSNE(merged, dims = 1:30)
merged <- FindNeighbors(merged, reduction = "pca", dims = 1:30)
merged <- FindClusters(merged, resolution = 0.5)


##用DimPlot()函数绘制散点图,reduction = "tsne",指定绘制类型；如果不指定,默认先从搜索 umap,然后 tsne, 再然后 pca；也可以直接使用这3个函数PCAPlot()、TSNEPlot()、UMAPPlot()； cols,pt.size分别调整分组颜色和点的大小；

DimPlot(merged,reduction = "tsne",label = TRUE,pt.size = 1.5)
DimPlot(merged,reduction = "tsne",label = TRUE,group.by="orig.ident",pt.size = 1.5)
