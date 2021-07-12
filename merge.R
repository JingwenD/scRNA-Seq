library(dplyr)
library(Seurat)
library(patchwork)
# Load the PBMC dataset
Y1.data <- Read10X(data.dir = "./PBMC/Y1/filtered_feature_bc_matrix/")
# Initialize the Seurat object with the raw (non-normalized data).
Y1 <- CreateSeuratObject(counts = Y1.data, project = "PsO",assay = "Y1", min.cells = 3, min.features = 200)
Y1
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
Y1[["percent.mt"]] <- PercentageFeatureSet(Y1, pattern = "^MT-")
Y1[["percent.rp"]] <- PercentageFeatureSet(Y1, pattern = "^RP[SL]")
VlnPlot(Y1, features = c("nFeature_Y1", "nCount_Y1", "percent.mt","percent.rp"), ncol = 4)

plot1 <- FeatureScatter(Y1, feature1 = "nCount_Y1", feature2 = "percent.mt")
plot2 <- FeatureScatter(Y1, feature1 = "nCount_Y1", feature2 = "nFeature_Y1")
plot1 + plot2

Y1 <- subset(Y1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

folder<-"./PBMC/HC"
sample_name = list.files(folder)            
dir = paste0(folder,"/",sample_name,"/filtered_feature_bc_matrix/")   
n = length(dir)
#QC
for (sample in sample_name){
  data <- Read10X(paste0(folder,"/",sample,"/filtered_feature_bc_matrix/"))
  # Initialize the Seurat object with the raw (non-normalized data).
  data <- CreateSeuratObject(counts = data, project = "PsO", min.cells = 3, min.features = 200)
  data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
  data[["percent.rp"]] <- PercentageFeatureSet(data, pattern = "^RP[SL]")
  pdf(paste0(folder,"/QC/",sample,"_QC1.pdf"),height = 5, width = 10)
  VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rp"), ncol = 4)
  dev.off()
  pdf(paste0(folder,"/QC/",sample,"_QC2.pdf"),height = 5, width = 8)
  plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  plot1 + plot2
  dev.off()
}
folder<-"./PBMC/PsO"
folder<-"./PBMC/HC"
sample_name1 = list.files(folder)
# read data and perform
for (sample in sample_name){
  seurat_data <- Read10X(paste0(folder,"/",sample,"/filtered_feature_bc_matrix/"))
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.cells = 3, 
                                   min.features = 200, 
                                   project = sample)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
  seurat_obj[["percent.rp"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  seurat_obj<- SCTransform(seurat_obj, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
  assign(sample, seurat_obj)
  rm(seurat_data)
  rm(seurat_obj)
}

pso.big <- merge(Y1, y = c(Y10,Y11,Y2,Y6,Y7,Y8,Y9), add.cell.ids = sample_name, project = "PsO")
# These are now standard steps in the Seurat workflow for visualization and clustering

obj.features <- SelectIntegrationFeatures(object.list = c(Y1,Y10,Y11,Y2,Y6,Y7,Y8,Y9), nfeatures = 2000)
VariableFeatures(pso.big[["SCT"]]) <- obj.features
pso.big <- RunPCA(pso.big, verbose = FALSE)
pso.big <- RunUMAP(pso.big, dims = 1:30, verbose = FALSE)
pso.big <- FindNeighbors(pso.big, dims = 1:30, verbose = FALSE)
pso.big <- FindClusters(pso.big, verbose = FALSE)
pdf("PsO_umap_cluster.pdf",height = 6, width = 9)
DimPlot(pso.big, label = TRUE)
dev.off()
pdf("PsO_umap_patient.pdf",height = 6, width = 9)
DimPlot(pso.big, group.by = "orig.ident",label = TRUE)
dev.off()

hc.big <- merge(H1, y = c(H2,H3,H4,H6,H7,M1,M2), add.cell.ids = sample_name, project = "HC")
# These are now standard steps in the Seurat workflow for visualization and clustering
obj.features1 <- SelectIntegrationFeatures(object.list = c(H1,H2,H3,H4,H6,H7,M1,M2), nfeatures = 2000)
VariableFeatures(hc.big[["SCT"]]) <- obj.features1
hc.big <- RunPCA(hc.big, verbose = FALSE)
hc.big <- RunUMAP(hc.big, dims = 1:30, verbose = FALSE)
hc.big <- FindNeighbors(hc.big, dims = 1:30, verbose = FALSE)
hc.big <- FindClusters(hc.big, verbose = FALSE)
pdf("HC_umap_cluster.pdf",height = 6, width = 9)
DimPlot(hc.big, label = TRUE) #cols,pt.size分别调整分组颜色和点的大小
dev.off()
pdf("HC_umap_patient.pdf",height = 6, width = 9)
DimPlot(hc.big, group.by = "orig.ident",label = TRUE)
dev.off()

pdf("umap_overlap_PsO.pdf",height = 6, width = 9)
DimPlot(all.big,reduction = "umap",label = TRUE,group.by="condition",pt.size = 0.5, cols =c("gray","red"))
dev.off()

all.big <- merge(H1, y = c(H2,H3,H4,H6,H7,M1,M2,Y1,Y10,Y11,Y2,Y6,Y7,Y8,Y9), add.cell.ids = c(sample_name1,sample_name), project = "ALL")
# These are now standard steps in the Seurat workflow for visualization and clustering

all.big[["batch"]]<-"B"
all.big@meta.data$batch[all.big@meta.data$orig.ident %in% c("Y1","Y2")]<-"A"
all.big@meta.data$batch[all.big@meta.data$orig.ident %in% c("Y9","Y10","Y11","H4","H6","H7")]<-"C"
all.big@meta.data$batch[all.big@meta.data$orig.ident %in% c("M1","M2")]<-"D"
table(all.big@meta.data$batch)

# # 2. 特征提取
# all.big = FindVariableFeatures(object = all.big,selection.method = "vst", nfeatures = 2000)
# # 所有特征
# VariableFeatures(sce)
# top20=head(VariableFeatures(sce),20)# 提取差异最大的 top20 基因
# plot1= VariableFeaturePlot(sce)
# plot2=LabelPoints(plot = plot1, points = top20, repel = TRUE)
# plot2


obj.features <- SelectIntegrationFeatures(object.list = c(H1,H2,H3,H4,H6,H7,M1,M2,Y1,Y10,Y11,Y2,Y6,Y7,Y8,Y9), nfeatures = 2000)
VariableFeatures(all.big[["SCT"]]) <- obj.features

all.big <- AddMetaData(all.big, substr(all.big@meta.data$orig.ident,1,1), col.name = "condition")
all.big@meta.data$condition<-gsub("H","HC",all.big@meta.data$condition)
all.big@meta.data$condition<-gsub("M","HC",all.big@meta.data$condition)
all.big@meta.data$condition<-gsub("Y","PsO",all.big@meta.data$condition)


library(tidyverse)
all.big@meta.data%>% 
  ggplot(aes(x=condition, fill=condition)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
ggsave("NCells_condition.pdf",height = 5, width = 5)

all.big@meta.data%>% 
  ggplot(aes(x=orig.ident, fill=orig.ident)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")
ggsave("NCells_sample.pdf",height = 5, width = 10)

all.big@meta.data%>% 
  ggplot(aes(color=condition, x=nFeature_RNA, fill= condition)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10()  #+ geom_vline(xintercept = 300)
ggsave("Feature_per_cell.pdf",height = 5, width = 7)

all.big@meta.data%>% 
  ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10()  #+ geom_vline(xintercept = 300)
ggsave("Feature_per_cell1.pdf",height = 5, width = 7)


all.big@meta.data%>% 
  ggplot(aes(x=orig.ident, y=log10(nFeature_RNA), fill=orig.ident)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) #+ggtitle("NCells vs NGenes")
ggsave("Feature_per_cell2.pdf",height = 5, width = 7)


all.big@meta.data%>% 
  ggplot(aes(color=condition, x=percent.mt, fill= condition)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10()  #+ geom_vline(xintercept = 300)
ggsave("Mt_per_cell.pdf",height = 5, width = 7)

pdf("features_petients.pdf",height = 6, width = 9)
VlnPlot(all.big, features = 'nFeature_RNA', group.by = 'orig.ident')
dev.off()
pdf("features_disease.pdf",height = 6, width = 6)
VlnPlot(all.big, features = 'nFeature_RNA', group.by = 'condition',pt.size=0.001)
dev.off()
save(all.big, file="PBMC_merge.RData")
