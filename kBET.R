library(kBET)
library(ggplot2)
#data: a matrix (rows: samples, columns: features (genes))
#batch: vector or factor with batch label of each cell 
load("~/Project/10_SingleCell/outs/PBMC_merge.RData")
data<-t(as.matrix(all.big@assays$RNA@data))
batch<-all.big@meta.data$batch

subset_size <- 0.1 #subsample to 10% of the data
subset_id <- sample.int(n = length(batch), size = floor(subset_size * length(batch)), replace=FALSE)
batch.estimate <- kBET(data[subset_id,], batch[subset_id], plot=FALSE)

plot.data <- data.frame(class=rep(c('observed', 'expected'), 
                                  each=length(batch.estimate$stats$kBET.observed)), 
                        data =  c(batch.estimate$stats$kBET.observed,
                                  batch.estimate$stats$kBET.expected))
ggplot(plot.data, aes(class, data)) + geom_boxplot() + 
  labs(x='Test', y='Rejection rate',title='kBET test results') +
  theme_bw() +  
  scale_y_continuous(limits=c(0,1))
ggsave("Batch_raw01.pdf",height = 5, width = 7)
write.csv(batch.estimate, file = "batchEstimate.csv")