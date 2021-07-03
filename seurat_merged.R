##安装seurat
#install.packages('Seurat')

setwd("C:/work/data/PBMCsc")
##载入seurat包
library(dplyr)
library(Seurat)

##读入pbmc数据
y1_5p <- Read10X(data.dir = "C:/work/data/PBMCsc/Y1/filtered_feature_bc_matrix")
y2_5p <- Read10X(data.dir = "C:/work/data/PBMCsc/Y2/filtered_feature_bc_matrix")

##import TCR/BCR data
add_clonotype <- function(tcr_prefix, seurat_obj, type="t"){    
  tcr <- read.csv(paste(tcr_prefix,"filtered_contig_annotations.csv", sep=""))    
  # Remove the -1 at the end of each barcode.（注意，此步骤如果标记使用不同的barcode，比如多了个-1,可以使用 tcr$barcode <- gsub("-1", "", tcr$barcode)进行提取）    
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

##创建Seurat对象与数据过滤
y1_5p <- CreateSeuratObject(counts = y1_5p, project = "Y1", min.cells = 3, min.features = 200) #过滤少于200个基因的细胞和少于3个细胞覆盖的基因
y2_5p <- CreateSeuratObject(counts = y2_5p, project = "Y2", min.cells = 3, min.features = 200)

merged<-merge(y1_5p,y2_5p, add.cell.ids=c("Y1","Y2"))#为了防止数据集之间的barcodes重叠，给细胞添加ID。

##计算每个细胞的线粒体基因转录本数的百分比（%）,使用[[ ]] 操作符存放到metadata中，MT-线粒体基因，有时候是小写
merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")
merged[["percent.ribo"]] <- PercentageFeatureSet(merged, pattern = "^RP[SL]")

##展示基因及线粒体百分比
VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo"), ncol = 2)

plot1 <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(merged, feature1= "nCount_RNA", feature2="percent.ribo")#feature1到底是"nCount_RNA"还是"nFeature_RNA"
CombinePlots(plots = list(plot1, plot2))

##过滤细胞：保留gene数大于200小于5000的细胞；目的是去掉空GEMs和1个GEMs包含2个以上细胞的数据；而保留线粒体基因的转录本数低于10%的细胞,为了过滤掉死细胞等低质量的细胞数据。
merged <- subset(merged, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

##表达量数据标准化,LogNormalize的算法：A = log( 1 + ( UMIA ÷ UMITotal ) × 10000 ,此为cpm，如果scale factor为1000000，则为TPM?)
merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000)
#merged <- NormalizeData(merged) 或者用默认的

# 检查每个样本的细胞数量，将其存储在metadata的orig.ident中，并自动设置为active ident。
merged
table(Idents(merged))
head(alldata@meta.data)

### Identify the 10 most highly variable gene，鉴定细胞间表达量高变的基因(2000个）,用于下游分析,如PCA；
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)

##提取表达量变变化最高的10个基因；
top20 <- head(VariableFeatures(merged), 20)
top20

plot1 <- VariableFeaturePlot(merged)
plot2 <- LabelPoints(plot = plot1, points = top20)
CombinePlots(plots = list(plot1, plot2))

##PCA分析：
#PCA分析数据准备,使用ScaleData()进行数据归一化；默认只是标准化高变基因（2000个）,速度更快,不影响PCA和分群,但影响热图的绘制。
#merged <- ScaleData(merged,vars.to.regress ="percent.mt")
 
##Scaling the data，scale之后基因表达区间为10-50？而对所有基因进行标准化的方法如下：
all.genes <- rownames(merged)
# merged <- ScaleData(merged, features = all.genes)
merged <- ScaleData(merged, features = all.genes, vars.to.regress = "percent.mt") #vars.to.regress ="nUMI"/"percent.mt".

##Perform linear dimensional reduction，线性降维（PCA）,默认用高变基因集,但也可通过features参数自己指定；
merged <- RunPCA(merged, features = VariableFeatures(object = merged))
 
##检查PCA分群结果, 这里只展示前5个PC,每个PC只显示5个基因；
print(merged[["pca"]], dims = 1:12, nfeatures = 10)

##展示主成分基因分值
VizDimLoadings(merged, dims = 1:2, reduction = "pca")

##绘制pca散点图
DimPlot(merged, reduction = "pca")

##画第1个或15个主成分的热图；
DimHeatmap(merged, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(merged, dims = 1:15, cells = 500, balanced = TRUE)

##确定数据集的分群个数 ?
#方法1：Jackstraw置换检验算法；重复取样（原数据的1%）,重跑PCA,鉴定p-value较小的PC；计算‘null distribution’(即零假设成立时)时的基因scores。
#merged <- JackStraw(merged, num.replicate = 100)
#merged <- ScoreJackStraw(merged, dims = 1:20)
#JackStrawPlot(merged, dims = 1:15)

#方法2：选择出主成分的数目，用于后续细胞分类。肘部图（碎石图）,基于每个主成分对方差解释率的排名。
ElbowPlot(merged, ndims = 50)

##分群个数这里选择12,建议尝试选择多个主成分个数做下游分析,对整体影响不大；在选择此参数时,建议选择偏高的数字（为了获取更多的稀有分群,“宁滥勿缺”）；有些亚群很罕见,如果没有先验知识,很难将这种大小的数据集与背景噪声区分开来。

##非线性降维（UMAP/tSNE)基于PCA空间中的欧氏距离计算nearest neighbor graph,优化任意两个细胞间的距离权重（输入上一步得到的PC维数）。
merged <- FindNeighbors(merged, dims = 1:12)

##接着优化模型,resolution参数决定下游聚类分析得到的分群数,对于3K左右的细胞,设为0.4-1.2 能得到较好的结果(官方说明)；如果数据量增大,该参数也应该适当增大。
merged <- FindClusters(merged, resolution = 1.0)

##使用Idents（）函数可查看不同细胞的分群；
head(Idents(merged), 5)

##Seurat提供了几种非线性降维的方法进行数据可视化（在低维空间把相似的细胞聚在一起）,比如UMAP和t-SNE,运行UMAP需要先安装'umap-learn'包,这里不做介绍。
merged <- RunTSNE(merged, dims = 1:12) #如果有duplicates，可以添加tsne.method = "Rtsne",check_duplicates = FALSE

##用DimPlot()函数绘制散点图,reduction = "tsne",指定绘制类型；如果不指定,默认先从搜索 umap,然后 tsne, 再然后 pca；也可以直接使用这3个函数PCAPlot()、TSNEPlot()、UMAPPlot()； cols,pt.size分别调整分组颜色和点的大小；

DimPlot(merged,reduction = "tsne",label = TRUE,pt.size = 1.5)

plot1<-DimPlot(subset(merged, subset = orig.ident=='Y1'),reduction = "tsne",label = TRUE,pt.size = 1.5)
plot2<-DimPlot(subset(merged, subset = orig.ident=='Y2'),reduction = "tsne",label = TRUE,pt.size = 1.5)
CombinePlots(plots = list(plot1, plot2))

DimPlot(merged,reduction = "tsne",label = TRUE,group.by="orig.ident",pt.size = 1.5)

merged <- RunUMAP(merged, dims = 1:12)# # note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
plot1<-DimPlot(subset(merged, subset = orig.ident=='Y1'),reduction = "umap",label = TRUE,pt.size = 1.5)
plot2<-DimPlot(subset(merged, subset = orig.ident=='Y2'),reduction = "umap",label = TRUE,pt.size = 1.5)
CombinePlots(plots = list(plot1, plot2))



##比较2组数据中特定cluster的差异基因
merged@meta.data$sample_type <- paste(merged@meta.data$orig.ident, merged@meta.data$RNA_snn_res.0.5, sep = "_")
marker<-FindMarkers(merged, group.by="sample_type",ident.1 = "Y1_6", ident.2 = "Y2_6", logfc.threshold = 0, min.pct = 0.1)
marker

##细胞周期归类
merged<- CellCycleScoring(object = merged, g2m.features = cc.genes$g2m.genes, s.features = cc.genes$s.genes)
head(x = merged@meta.data)
DimPlot(merged,reduction = "tsne",label = TRUE,group.by="Phase",pt.size = 1.5, cols =c("#FF8000FF","blue","#4D9221"))


##存储结果
saveRDS(merged, file = "D:/merged_tutorial.rds")
save(merged,file="D:/merged_res0.5.Robj") 

##寻找cluster 1的marker
cluster1.markers <- FindMarkers(merged, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

##寻找每一cluster的marker，其中ident.1参数设置待分析的细胞类别，min.pct表示该基因表达数目占该类细胞总数的比例
merged.markers <- FindAllMarkers(merged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)#如何选方向
merged.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

##存储marker
write.table(merged.markers,file="D:/merged_allmarker.txt")


##各种绘图
##绘制Marker 基因的tsne图
FeaturePlot(merged, features = c("CD4", "CD8A","CD3E","IL2RA"),cols = c("gray", "red"))

##绘制Marker 基因的小提琴图
VlnPlot(merged, features = c("RPL38", "RPL39"),split.by = "orig.ident")

##绘制分cluster的热图
top10 <- merged.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(merged, features = top10$gene) + NoLegend()

##各种绘图
##绘制Marker 基因的tsne图
FeaturePlot(merged, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"),cols = c("gray", "red"))

##绘制Marker 基因的小提琴图
VlnPlot(merged, features = c("MS4A1", "CD79A"))

VlnPlot(merged, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

VlnPlot(merged, features = c("CD4", "CD8A"))

##绘制分cluster的热图
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

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
