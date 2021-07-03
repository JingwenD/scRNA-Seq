##��װseurat
#install.packages('Seurat')

setwd("C:/work/data/PBMCsc")
##����seurat��
library(dplyr)
library(Seurat)

##����pbmc����
y1_5p <- Read10X(data.dir = "C:/work/data/PBMCsc/Y1/filtered_feature_bc_matrix")
y2_5p <- Read10X(data.dir = "C:/work/data/PBMCsc/Y2/filtered_feature_bc_matrix")

##import TCR/BCR data
add_clonotype <- function(tcr_prefix, seurat_obj, type="t"){    
  tcr <- read.csv(paste(tcr_prefix,"filtered_contig_annotations.csv", sep=""))    
  # Remove the -1 at the end of each barcode.��ע�⣬�˲���������ʹ�ò�ͬ��barcode��������˸�-1,����ʹ�� tcr$barcode <- gsub("-1", "", tcr$barcode)������ȡ��    
  # Subsets so only the first line of each barcode is kept,   
  # as each entry for given barcode will have same clonotype.    
  tcr <- tcr[!duplicated(tcr$barcode), ]    # Only keep the barcode and clonotype columns.    
  # We'll get additional clonotype info from the clonotype table.    
  tcr <- tcr[,c("barcode", "raw_clonotype_id")]    
  names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"    
  # Clonotype-centric info.    
  clono <- read.csv(paste(tcr_prefix,"clonotypes.csv", sep=""))    
  # Slap the AA sequences onto our original table by clonotype_id.    
  tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa")])    
  names(tcr)[names(tcr) == "cdr3s_aa"] <- "cdr3s_aa"    
  # Reorder so barcodes are first column and set them as rownames.    
  tcr <- tcr[, c(2,1,3)]    
  rownames(tcr) <- tcr[,1]    
  tcr[,1] <- NULL    
  colnames(tcr) <- paste(type, colnames(tcr), sep="_")    
  # Add to the Seurat object's metadata.    
  clono_seurat <- AddMetaData(object=seurat_obj, metadata=tcr)    
  return(clono_seurat)
  }

s_balbc_pbmc <- add_clonotype("vdj_v1_mm_balbc_pbmc/vdj_v1_mm_balbc_pbmc_t_", s_balbc_pbmc, "t")s_balbc_pbmc <- add_clonotype("vdj_v1_mm_balbc_pbmc/vdj_v1_mm_balbc_pbmc_b_", s_balbc_pbmc, "b")head(s_balbc_pbmc[[]])

##����Seurat���������ݹ���
y1_5p <- CreateSeuratObject(counts = y1_5p, project = "Y1", min.cells = 3, min.features = 200) #��������200�������ϸ��������3��ϸ�����ǵĻ���
y2_5p <- CreateSeuratObject(counts = y2_5p, project = "Y2", min.cells = 3, min.features = 200)

merged<-merge(y1_5p,y2_5p, add.cell.ids=c("Y1","Y2"))#Ϊ�˷�ֹ���ݼ�֮���barcodes�ص�����ϸ������ID��

##����ÿ��ϸ�������������ת¼�����İٷֱȣ�%��,ʹ��[[ ]] ��������ŵ�metadata�У�MT-�����������ʱ����Сд
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")
merged[["percent.ribo"]] <- PercentageFeatureSet(merged, pattern = "^RP[SL]")

##չʾ����������ٷֱ�
VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 2)

plot1 <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(merged, feature1= "nCount_RNA", feature2="percent.ribo")#feature1������"nCount_RNA"����"nFeature_RNA"
CombinePlots(plots = list(plot1, plot2))

##����ϸ��������gene������200С��5000��ϸ����Ŀ����ȥ����GEMs��1��GEMs����2������ϸ�������ݣ�����������������ת¼��������10%��ϸ��,Ϊ�˹��˵���ϸ���ȵ�������ϸ�����ݡ�
merged <- subset(merged, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

##���������ݱ�׼��,LogNormalize���㷨��A = log( 1 + ( UMIA �� UMITotal ) �� 10000 ,��Ϊcpm�����scale factorΪ1000000����ΪTPM?)
merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000)
#merged <- NormalizeData(merged) ������Ĭ�ϵ�

# ���ÿ��������ϸ������������洢��metadata��orig.ident�У����Զ�����Ϊactive ident��
merged
table(Idents(merged))
head(alldata@meta.data)

### Identify the 10 most highly variable gene������ϸ����������߱�Ļ���(2000����,�������η���,��PCA��
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)

##��ȡ��������仯��ߵ�10������
top20 <- head(VariableFeatures(merged), 20)
top20

plot1 <- VariableFeaturePlot(merged)
plot2 <- LabelPoints(plot = plot1, points = top20)
CombinePlots(plots = list(plot1, plot2))

##PCA������
#PCA��������׼��,ʹ��ScaleData()�������ݹ�һ����Ĭ��ֻ�Ǳ�׼���߱����2000����,�ٶȸ���,��Ӱ��PCA�ͷ�Ⱥ,��Ӱ����ͼ�Ļ��ơ�
#merged <- ScaleData(merged,vars.to.regress ="percent.mt")
 
##Scaling the data��scale֮������������Ϊ10-50���������л�����б�׼���ķ������£�
all.genes <- rownames(merged)
# merged <- ScaleData(merged, features = all.genes)
merged <- ScaleData(merged, features = all.genes, vars.to.regress = "percent.mt") #vars.to.regress ="nUMI"/"percent.mt".

##Perform linear dimensional reduction�����Խ�ά��PCA��,Ĭ���ø߱����,��Ҳ��ͨ��features�����Լ�ָ����
merged <- RunPCA(merged, features = VariableFeatures(object = merged))
 
##���PCA��Ⱥ���, ����ֻչʾǰ5��PC,ÿ��PCֻ��ʾ5������
print(merged[["pca"]], dims = 1:12, nfeatures = 10)

##չʾ���ɷֻ����ֵ
VizDimLoadings(merged, dims = 1:2, reduction = "pca")

##����pcaɢ��ͼ
DimPlot(merged, reduction = "pca")

##����1����15�����ɷֵ���ͼ��
DimHeatmap(merged, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(merged, dims = 1:15, cells = 500, balanced = TRUE)

##ȷ�����ݼ��ķ�Ⱥ���� ?
#����1��Jackstraw�û������㷨���ظ�ȡ����ԭ���ݵ�1%��,����PCA,����p-value��С��PC�����㡮null distribution��(����������ʱ)ʱ�Ļ���scores��
#merged <- JackStraw(merged, num.replicate = 100)
#merged <- ScoreJackStraw(merged, dims = 1:20)
#JackStrawPlot(merged, dims = 1:15)

#����2��ѡ������ɷֵ���Ŀ�����ں���ϸ�����ࡣ�ⲿͼ����ʯͼ��,����ÿ�����ɷֶԷ�������ʵ�������
ElbowPlot(merged, ndims = 50)

##��Ⱥ��������ѡ��12,���鳢��ѡ�������ɷָ��������η���,������Ӱ�첻����ѡ��˲���ʱ,����ѡ��ƫ�ߵ����֣�Ϊ�˻�ȡ�����ϡ�з�Ⱥ,��������ȱ��������Щ��Ⱥ�ܺ���,���û������֪ʶ,���ѽ����ִ�С�����ݼ��뱳���������ֿ�����

##�����Խ�ά��UMAP/tSNE)����PCA�ռ��е�ŷ�Ͼ������nearest neighbor graph,�Ż���������ϸ����ľ���Ȩ�أ�������һ���õ���PCά������
merged <- FindNeighbors(merged, dims = 1:12)

##�����Ż�ģ��,resolution�����������ξ�������õ��ķ�Ⱥ��,����3K���ҵ�ϸ��,��Ϊ0.4-1.2 �ܵõ��ϺõĽ��(�ٷ�˵��)���������������,�ò���ҲӦ���ʵ�����
merged <- FindClusters(merged, resolution = 1.0)

##ʹ��Idents���������ɲ鿴��ͬϸ���ķ�Ⱥ��
head(Idents(merged), 5)

##Seurat�ṩ�˼��ַ����Խ�ά�ķ����������ݿ��ӻ����ڵ�ά�ռ�����Ƶ�ϸ������һ��,����UMAP��t-SNE,����UMAP��Ҫ�Ȱ�װ'umap-learn'��,���ﲻ�����ܡ�
merged <- RunTSNE(merged, dims = 1:12) #�����duplicates����������tsne.method = "Rtsne",check_duplicates = FALSE

##��DimPlot()��������ɢ��ͼ,reduction = "tsne",ָ���������ͣ������ָ��,Ĭ���ȴ����� umap,Ȼ�� tsne, ��Ȼ�� pca��Ҳ����ֱ��ʹ����3������PCAPlot()��TSNEPlot()��UMAPPlot()�� cols,pt.size�ֱ����������ɫ�͵�Ĵ�С��

DimPlot(merged,reduction = "tsne",label = TRUE,pt.size = 1.5)

plot1<-DimPlot(subset(merged, subset = orig.ident=='Y1'),reduction = "tsne",label = TRUE,pt.size = 1.5)
plot2<-DimPlot(subset(merged, subset = orig.ident=='Y2'),reduction = "tsne",label = TRUE,pt.size = 1.5)
CombinePlots(plots = list(plot1, plot2))

DimPlot(merged,reduction = "tsne",label = TRUE,group.by="orig.ident",pt.size = 1.5)

merged <- RunUMAP(merged, dims = 1:12)# # note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
plot1<-DimPlot(subset(merged, subset = orig.ident=='Y1'),reduction = "umap",label = TRUE,pt.size = 1.5)
plot2<-DimPlot(subset(merged, subset = orig.ident=='Y2'),reduction = "umap",label = TRUE,pt.size = 1.5)
CombinePlots(plots = list(plot1, plot2))



##�Ƚ�2���������ض�cluster�Ĳ������
merged@meta.data$sample_type <- paste(merged@meta.data$orig.ident, merged@meta.data$RNA_snn_res.0.5, sep = "_")
marker<-FindMarkers(merged, group.by="sample_type",ident.1 = "Y1_6", ident.2 = "Y2_6", logfc.threshold = 0, min.pct = 0.1)
marker

##ϸ�����ڹ���
merged<- CellCycleScoring(object = merged, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
head(x = merged@meta.data)
DimPlot(merged,reduction = "tsne",label = TRUE,group.by="Phase",pt.size = 1.5, cols =c("#FF8000FF","blue","#4D9221"))


##�洢���
saveRDS(merged, file = "D:/merged_tutorial.rds")
save(merged,file="D:/merged_res0.5.Robj") 

##Ѱ��cluster 1��marker
cluster1.markers <- FindMarkers(merged, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

##Ѱ��ÿһcluster��marker������ident.1�������ô�������ϸ�����min.pct��ʾ�û��������Ŀռ����ϸ�������ı���
merged.markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)#���ѡ����
merged.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

##�洢marker
write.table(merged.markers,file="D:/merged_allmarker.txt")


##���ֻ�ͼ
##����Marker �����tsneͼ
FeaturePlot(merged, features = c("CD4", "CD8A","CD3E","IL2RA"),cols = c("gray", "red"))

##����Marker �����С����ͼ
VlnPlot(merged, features = c("RPL38", "RPL39"),split.by = "orig.ident")

##���Ʒ�cluster����ͼ
top10 <- merged.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(merged, features = top10$gene) + NoLegend()

##���ֻ�ͼ
##����Marker �����tsneͼ
FeaturePlot(merged, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"),cols = c("gray", "red"))

##����Marker �����С����ͼ
VlnPlot(merged, features = c("MS4A1", "CD79A"))

VlnPlot(merged, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

VlnPlot(merged, features = c("CD4", "CD8A"))

##���Ʒ�cluster����ͼ
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

##��������ͼ
DotPlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"),cols = c("blue", "red"))

##����RidgePlot
RidgePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14"), ncol = 2)

############ȡ��һ����cluster�����Ƿ�Ⱥ
sub_pbmc<-subset(pbmc, idents = '0')

sub_pbmc <- FindVariableFeatures(sub_pbmc, selection.method = "vst", nfeatures = 2000)

sub_pbmc <- ScaleData(sub_pbmc, features = all.genes)
sub_pbmc <- ScaleData(sub_pbmc, vars.to.regress = "percent.mt")

sub_pbmc <- RunPCA(sub_pbmc, features = VariableFeatures(object = pbmc))

sub_pbmc <- FindNeighbors(sub_pbmc, dims = 1:10)

sub_pbmc <- FindClusters(sub_pbmc, resolution = 0.5)

sub_pbmc <- RunTSNE(sub_pbmc, dims = 1:10)

DimPlot(sub_pbmc, reduction = "tsne",label = TRUE, pt.size = 1.5)


##ÿ��cluster����marker�涨ϸ������
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 1.5) + NoLegend()

write.table(Idents(pbmc),"D:/cell-type.xls",sep="\t")

##�洢ϸ��������
saveRDS(pbmc, file = "D:/pbmc3k_final.rds")

##����ϸ�����ͻ���Marker�����С����ͼ
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

##����ϸ�����ͻ��Ʒ�cluster����ͼ
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

##����ϸ�����ͻ�������ͼ
DotPlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"),cols = c("blue", "red"))

##����ϸ�����ͻ���RidgePlot
RidgePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14"), ncol = 2)