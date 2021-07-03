###该代码主要是为了消除样本之间的异质性，即聚类会按照样本的情况下用该代码。

##载入seurat包
library(dplyr)
library(Seurat)

##读入pbmc数据
pbmc.data <- Read10X(data.dir = "D:/BC11_TUMOR1/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "TUMOR1", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
pbmc <- NormalizeData(pbmc, verbose = FALSE)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000,verbose = FALSE)


##读入pbmc1数据
pbmc1.data <- Read10X(data.dir = "D:/BC9_TUMOR1/")
pbmc1 <- CreateSeuratObject(counts = pbmc1.data, project = "TUMOR2", min.cells = 3, min.features = 200)
pbmc1[["percent.mt"]] <- PercentageFeatureSet(pbmc1, pattern = "^MT-")
pbmc1 <- subset(pbmc1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
pbmc1 <- NormalizeData(pbmc1, verbose = FALSE)
pbmc1 <- FindVariableFeatures(pbmc1, selection.method = "vst", nfeatures = 2000,verbose = FALSE)

##然后使用FindIntegrationAnchors函数识别锚点，该函数将Seurat对象列表作为输入，并使用这些锚点将两个数据集集成在一起IntegrateData。这里选的维度是30，作者建议可以在10-50间调试
merged.anchors <- FindIntegrationAnchors(object.list = list(pbmc, pbmc1), dims = 1:30)
# 进行数据集整合
# 已经整合后的表达矩阵存储在Assay中，未处理的表达矩阵在RNA对象中
merged <- IntegrateData(anchorset = merged.anchors, dims = 1:30)


DefaultAssay(merged) <- "integrated"
##标准化数据
merged <- ScaleData(merged,verbose = FALSE)

##线性降维（PCA）,默认用高变基因集,但也可通过features参数自己指定；
merged <- RunPCA(merged, npcs = 30, verbose = FALSE)

##Seurat提供了几种非线性降维的方法进行数据可视化（在低维空间把相似的细胞聚在一起）,比如UMAP和t-SNE,运行UMAP需要先安装'umap-learn'包,这里不做介绍。
merged <- RunTSNE(merged, reduction = "pca", dims = 1:30)
merged <- FindNeighbors(merged, reduction = "pca", dims = 1:30)
merged <- FindClusters(merged, resolution = 0.5)


##用DimPlot()函数绘制散点图,reduction = "tsne",指定绘制类型；如果不指定,默认先从搜索 umap,然后 tsne, 再然后 pca；也可以直接使用这3个函数PCAPlot()、TSNEPlot()、UMAPPlot()； cols,pt.size分别调整分组颜色和点的大小；

DimPlot(merged,reduction = "tsne",label = TRUE,pt.size = 1.5)
DimPlot(merged,reduction = "tsne",label = TRUE,group.by="orig.ident",pt.size = 1.5)
