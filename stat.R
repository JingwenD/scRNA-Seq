library(tidyverse)
celltype2inCluster<-all.big@meta.data %>%
  group_by(cell.type2,seurat_clusters) %>%
  summarise(n=n())

Clusterincelltype2<-all.big@meta.data %>%
  group_by(seurat_clusters,cell.type2) %>%
  summarise(n=n()) 

Clusterincelltype2<-all.big@meta.data %>%
  group_by(seurat_clusters,cell.type2) %>%
  summarise(n=n()) %>%
  spread(seurat_clusters,n)
Clusterincelltype2<- as.data.frame(Clusterincelltype2)
rownames(Clusterincelltype2)<-Clusterincelltype2$cell.type2
Clusterincelltype2<- Clusterincelltype2[,2:30]
Clusterincelltype2<- apply(Clusterincelltype2,2,unlist)

library(corrplot)
pdf("ClustervsCelltype2.pdf",height = 10, width = 12)
corrplot(Clusterincelltype2, method = "circle", 
         type = "full", add = FALSE, col = "lightblue", 
         bg = "white", title = "", na.label = " ",addCoef.col = "black",is.corr = F, diag = T)
dev.off()

all.big@meta.data %>%
  group_by(condition,cell.type2,orig.ident) %>%
  summarise(n=n()) %>%
  ggplot(aes(x=n, y=condition)) +
  geom_boxplot(width = 0.5,position=position_dodge(0.9))+ #绘制箱线图
  geom_jitter(aes(fill= condition),width =0.01,shape = 21,size=2.5)+
  facet_wrap(~cell.type2, scales="free_y", ncol=5) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("celltype2_condition.pdf",height = 11, width = 12)
  
all.big@meta.data %>%
  group_by(condition,seurat_clusters,orig.ident) %>%
  summarise(n=n()) %>%
  ggplot(aes(x=n, y=condition)) +
  geom_boxplot(width = 0.5,position=position_dodge(0.9))+ #绘制箱线图
  geom_jitter(aes(fill= condition),width =0.01,shape = 21,size=2.5)+
  facet_wrap(~seurat_clusters, scales="free_y", ncol=5) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("cluster_condition.pdf",height = 11, width = 12)  


library(rstatix)
all.big@meta.data %>%
  group_by(condition,cell.type2,orig.ident) %>%
  summarise(n=n()) %>%
  group_by(cell.type2) %>% 
  rstatix::t_test(n ~ condition, data = .) %>%
  adjust_pvalue() %>%
  add_significance("p.adj") %>%
  write_csv("celltype2_PsOvsHC.csv")

all.big@meta.data %>%
  group_by(condition,seurat_clusters,orig.ident) %>%
  summarise(n=n()) %>%
  group_by(seurat_clusters) %>% 
  rstatix::t_test(n ~ condition, data = .) %>%
  adjust_pvalue() %>%
  add_significance("p.adj") %>%
  write_csv("clusterP_PsOvsHC.csv")
