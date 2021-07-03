##安装monocle
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("monocle")


##载入monocle包
library(Seurat)
library(dplyr)
library(Matrix)
library(monocle)


##读入pbmc数据
load('D:/res0.5.Robj')

##取需要做轨迹的2个亚群
sub_pbmc <- subset(pbmc, idents=c(3,4))

markers <- FindAllMarkers(sub_pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

save(sub_pbmc,file="D:/sub.Robj") 


##往monocle加载Seurat对象
seurat_object <- sub_pbmc
data <- as(as.matrix(seurat_object@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = seurat_object@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)


#Construct monocle cds
cds <- newCellDataSet(data,
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit = 0.5,
                         expressionFamily = negbinomial.size());


cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
expressed_genes <- row.names(subset(fData(cds)))
save(cds,file="D:/cds_normal.Robj")

###用提供的marker建轨迹
ordering_genes<-as.matrix(top10)
cds <- setOrderingFilter(cds, ordering_genes = ordering_genes)
cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree',auto_param_selection = F) # take long time
cds <- orderCells(cds)
save(cds,file="D:/orderCells.Robj")
write.table(pData(cds),file="D:/my_pseudotime.txt")

plot_cell_trajectory(cds, color_by = "State")
plot_cell_trajectory(cds, color_by = "Pseudotime")

pData(cds)

plot_cell_trajectory(cds, color_by = "RNA_snn_res.0.5")


###展示感兴趣基因在轨迹上表达量
plot_cell_trajectory(cds, markers="CD79A", cell_size=0.5, use_color_gradient=T) + scale_color_gradient2(low="gray",mid="yellow",high="red")


###展示感兴趣基因热图(一个方向走)
load('D:/orderCells.Robj')
my_pseudotime_de <- differentialGeneTest(cds,fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 5)

save(my_pseudotime_de,file="D:/my_pseudotime_de.Robj")

write.table(my_pseudotime_de,file="D:/my_pseudotime_de.txt")

c<-subset(my_pseudotime_de, qval < 10^-20)

write.table(c[order(c[,4]),],file="D:/my_pseudotime_de1.txt")

sig_gene_names <- row.names(subset(my_pseudotime_de, qval < 10^-20))

plot_pseudotime_heatmap(cds[sig_gene_names,], num_clusters = 4,cores = 5,use_gene_short_name = TRUE,show_rownames = TRUE)


###展示感兴趣基因拟时间表达变化(一个方向走)
to_be_tested <- row.names(subset(fData(cds),gene_short_name %in% c("CD79A", "NKG7", "HLA-DRA", "CCL5")))
cds_subset <- cds[to_be_tested,]

diff_test_res <- differentialGeneTest(cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)")

diff_test_res[,c("gene_short_name", "pval", "qval")]


plot_genes_in_pseudotime(cds_subset,color_by="RNA_snn_res.0.5")


###展示感兴趣基因热图(两个方向走)
plot_genes_branched_heatmap(cds[sig_gene_names,], num_clusters = 4,cores = 5,use_gene_short_name = TRUE,show_rownames = TRUE)


###展示感兴趣基因拟时间表达变化(两个方向走)
plot_genes_branched_pseudotime(cds_subset, branch_point = 2,color_by="RNA_snn_res.0.5")
