###�ô�����Ҫ��Ϊ����������֮��������ԣ�������ᰴ��������������øô��롣

##����seurat��
library(dplyr)
library(Seurat)

##����pbmc����
pbmc.data <- Read10X(data.dir = "D:/BC11_TUMOR1/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "TUMOR1", min.cells = 3, min.features = 200)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
pbmc <- NormalizeData(pbmc, verbose = FALSE)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000,verbose = FALSE)


##����pbmc1����
pbmc1.data <- Read10X(data.dir = "D:/BC9_TUMOR1/")
pbmc1 <- CreateSeuratObject(counts = pbmc1.data, project = "TUMOR2", min.cells = 3, min.features = 200)
pbmc1[["percent.mt"]] <- PercentageFeatureSet(pbmc1, pattern = "^MT-")
pbmc1 <- subset(pbmc1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
pbmc1 <- NormalizeData(pbmc1, verbose = FALSE)
pbmc1 <- FindVariableFeatures(pbmc1, selection.method = "vst", nfeatures = 2000,verbose = FALSE)

##Ȼ��ʹ��FindIntegrationAnchors����ʶ��ê�㣬�ú�����Seurat�����б���Ϊ���룬��ʹ����Щê�㽫�������ݼ�������һ��IntegrateData������ѡ��ά����30�����߽��������10-50�����
merged.anchors <- FindIntegrationAnchors(object.list = list(pbmc, pbmc1), dims = 1:30)
# �������ݼ�����
# �Ѿ����Ϻ�ı������洢��Assay�У�δ�����ı��������RNA������
merged <- IntegrateData(anchorset = merged.anchors, dims = 1:30)


DefaultAssay(merged) <- "integrated"
##��׼������
merged <- ScaleData(merged,verbose = FALSE)

##���Խ�ά��PCA��,Ĭ���ø߱����,��Ҳ��ͨ��features�����Լ�ָ����
merged <- RunPCA(merged, npcs = 30, verbose = FALSE)

##Seurat�ṩ�˼��ַ����Խ�ά�ķ����������ݿ��ӻ����ڵ�ά�ռ�����Ƶ�ϸ������һ��,����UMAP��t-SNE,����UMAP��Ҫ�Ȱ�װ'umap-learn'��,���ﲻ�����ܡ�
merged <- RunTSNE(merged, reduction = "pca", dims = 1:30)
merged <- FindNeighbors(merged, reduction = "pca", dims = 1:30)
merged <- FindClusters(merged, resolution = 0.5)


##��DimPlot()��������ɢ��ͼ,reduction = "tsne",ָ���������ͣ������ָ��,Ĭ���ȴ����� umap,Ȼ�� tsne, ��Ȼ�� pca��Ҳ����ֱ��ʹ����3������PCAPlot()��TSNEPlot()��UMAPPlot()�� cols,pt.size�ֱ����������ɫ�͵�Ĵ�С��

DimPlot(merged,reduction = "tsne",label = TRUE,pt.size = 1.5)
DimPlot(merged,reduction = "tsne",label = TRUE,group.by="orig.ident",pt.size = 1.5)