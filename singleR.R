library(Seurat)
library(SingleR)
library(celldex)
browseVignettes("celldex")
#Most of the labels refer to blood subpopulations but cell types from other tissues are also available.
ref <- HumanPrimaryCellAtlasData()
# Decent immune cell granularity, does not contain finer monocyte and dendritic cell subtypes.
ref <- BlueprintEncodeData()