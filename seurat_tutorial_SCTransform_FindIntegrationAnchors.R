###�ô���ͬʱ���������ڼ�������������ԣ���������SCTransform����������ê�㷽����

##����seurat��
rm(list = ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
options(future.globals.maxSize = 4000 * 1024^2)


##����pbmc����
pbmc.data <- Read10X(data.dir = "D:/BC11_TUMOR1/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "TUMOR1", min.cells = 3, min.features = 500)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & percent.mt < 5)
pbmc<-SCTransform(pbmc, verbose = FALSE)

##����pbmc1����
pbmc1.data <- Read10X(data.dir = "D:/BC11_TUMOR2/")
pbmc1 <- CreateSeuratObject(counts = pbmc1.data, project = "TUMOR2", min.cells = 3, min.features = 500)
pbmc1[["percent.mt"]] <- PercentageFeatureSet(pbmc1, pattern = "^MT-")
pbmc1 <- subset(pbmc1, subset = nFeature_RNA > 500 & nFeature_RNA < 3000 & percent.mt < 5)
pbmc1<-SCTransform(pbmc1, verbose = FALSE)

features <- SelectIntegrationFeatures(object.list = list(pbmc, pbmc1), nfeatures = 3000)
list <- PrepSCTIntegration(object.list = list(pbmc, pbmc1), anchor.features = features, verbose = FALSE)


##Ȼ��ʹ��FindIntegrationAnchors����ʶ��ê�㣬�ú�����Seurat�����б���Ϊ���룬��ʹ����Щê�㽫�������ݼ�������һ��IntegrateData��
merged.anchors <- FindIntegrationAnchors(object.list = list, normalization.method = "SCT", anchor.features = features, verbose = FALSE)
merged <- IntegrateData(anchorset = merged.anchors, normalization.method = "SCT",verbose = FALSE)


##���Խ�άPCA��
merged <- RunPCA(object = merged, verbose = FALSE)

##Seurat�ṩ�˼��ַ����Խ�ά�ķ����������ݿ��ӻ����ڵ�ά�ռ�����Ƶ�ϸ������һ��,����UMAP��t-SNE,����UMAP��Ҫ�Ȱ�װ'umap-learn'��,���ﲻ�����ܡ�
merged <- RunTSNE(merged, dims = 1:30)
merged <- FindNeighbors(merged, reduction = "pca", dims = 1:30)
merged <- FindClusters(merged, resolution = 0.5)


##��DimPlot()��������ɢ��ͼ,reduction = "tsne",ָ���������ͣ������ָ��,Ĭ���ȴ����� umap,Ȼ�� tsne, ��Ȼ�� pca��Ҳ����ֱ��ʹ����3������PCAPlot()��TSNEPlot()��UMAPPlot()�� cols,pt.size�ֱ����������ɫ�͵�Ĵ�С��

DimPlot(merged,reduction = "tsne",label = TRUE,pt.size = 1.5)
DimPlot(merged,reduction = "tsne",label = TRUE,group.by="orig.ident",pt.size = 1.5)