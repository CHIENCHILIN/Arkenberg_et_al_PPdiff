library(Seurat)
library(tidyverse)
library(ggplot2)

# Path
Path <- "C:/ExampleFolder"


# Load 10X data in Seurat objects
CTR <- Read10X(paste0(Path, "/2D/Data")) %>% CreateSeuratObject(); CTR$Cond <- "CTR"
GEL <- Read10X(paste0(Path, "/3D_Gel/Data/")) %>% CreateSeuratObject(); GEL$Cond <- "GEL"


#### Filtering functions ####
# Filters for (log-transformed) mitochondrial percentage genes
filterPlot <- function(obj) {
  obj[["mt.per"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  mt.per.vec <- obj$mt.per
  obj <- subset(obj, mt.per < 15)
  return(obj)
}

# Filtering for housekeeping gene expression
# Appends metadata column "hk" that stores the log normalized expression of a selected housekeeping gene
# if subset = TRUE, will subset the seurat object discarding cells ~std.dev~
# standard deviations away from the mean 
housekeeping.filter <- function(object, housekeeping.gene = NULL, 
                                subset = F, std.dev = 2) {
  #Sanitize input
  housekeeping.gene <- as.character(housekeeping.gene)
  # Creates a vector containing the reads found in each cell and log normalizing
  expression.vector <- object@assays$RNA@counts[housekeeping.gene, ] %>% log1p()
  #Add metadata column to the Seurat object to facilitate subsetting
  object <- AddMetaData(object, metadata = expression.vector,
                        col.name = "hk")
  # Determines whether or not the object should be subset such that
  # cells greater or less than ~std.dev~ standard deviations are discarded
  if (subset){
    object <- subset(object, hk > mean(object$hk) - std.dev*sd(object$hk) &
                       hk < mean(object$hk) + std.dev*sd(object$hk))
  }
  return(object)
}


#### Filtering ####
# GroupDatasets into a list #
list.all <- list(CTR, GEL)

# Filters each dataset in list.all according to mitochondrial percentage,
# gives number of cells removed for each dataset
counts.list.pre <- sapply(list.all, ncol)
list.all <- lapply(list.all, filterPlot)
counts.list.post <- sapply(list.all, ncol)
print(counts.list.pre - counts.list.post)

# Filtering based on housekeeping gene expression; RPL27 was used here
list.all <- lapply(list.all, housekeeping.filter, housekeeping.gene = "RPL27", subset = T)
counts.list.post.2 <- sapply(list.all, ncol)
print(counts.list.post - counts.list.post.2)


#### Merge datasets ####
MGD <- merge(list.all[[1]], c(list.all[[2]]))
saveRDS(MGD, file = paste0(Path, "/Output/Merged.rds"), compress = F)


#### SCTransform, PCA ####
MGD <- SCTransform(MGD, ncells = 15000, vars.to.regress = "mt.per")
MGD <- RunPCA(MGD)

tiff(paste0(Path, "/Output/DimPlot_Merged.tiff"), width = 6, height = 4, units = "in", res = 600)
  DimPlot(MGD)
dev.off()

tiff(paste0(Path, "/Output/ElbowPlot_Merged.tiff"), width = 6, height = 4, units = "in", res = 600)
  ElbowPlot(MGD, ndims = 50, reduction = "pca")
dev.off()


#### PC determination ####
pct <- MGD[["pca"]]@stdev / sum(MGD[["pca"]]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)


#### Construct SNN graph and UMAP ####
dim = 12
res = 0.8
MGD <- FindNeighbors(MGD, k.param = 160, nn.method = "rann", dims = 1:dim) %>% FindClusters(resolution = res)
MGD <- RunUMAP(MGD, dims = 1:dim, n.neighbors = 10, min.dist = 0.25)

saveRDS(MGD, paste0(Path, "/Output/dim", dim, "_res", res, ".rds"))
# MGD <- readRDS(paste0(Path, Proj, "/Output/dim", dim, "_res", res, ".rds"))

tiff(paste0(Path, "/Output/UMAP_Merged.tiff"), width = 6, height = 4, units = "in", res = 600)
 DimPlot(MGD, pt.size = 0.4, split.by = "Cond")
dev.off()


#### Find Markers ####
Wilcox_tab <- FindAllMarkers(MGD, logfc.threshold = 0.25, test.use = "wilcox", min.pct = 0.1,
                             verbose = T, only.pos = T, return.thresh = 0.01) %>%
  select(c("cluster", "gene", "pct.1", "pct.2", "p_val_adj", "avg_log2FC")) %>%
  group_by(cluster) %>% top_n(30) %>% arrange(cluster, desc(avg_log2FC))
write.csv(Wilcox_tab, paste0(Path, "/Output/dim", dim, "_res", res, ".csv"))



#### Manuscript Outs ####

### Cluster IDs ###
PP <- c("GATA4", "NKX6-1", "PDX1", "PTF1A")
Endocrine <- c("CHGA", "INS", "NEUROG3", "NKX2-2")
AF <- c("IRX3", "IRX5", "SOX2", "SOX21")
Duct <- c("AKAP12", "AREG", "SERPINA1", "SPP1")
Intestinal <- c("CDX2", "HOXA9", "HOXB6", "HOXB9")
Hepatic <- c("ALB", "CYP1A1", "FGB", "TF")
Endothelial <- c("CD34", "CDH5", "ESM1", "KDR")
Fibroblast <- c("ACTA2", "CD248", "PDGFRA", "LUM")
Neural <- c("HES5", "PANTR1", "S100B", "STMN4")
EE <- c("GATA3", "TFAP2A", "TFAP2C", "KRT7")
EE2 <- c("GATA2", "GATA3", "TFAP2A", "TFAP2C", "KRT7", "TP63", "CDH5", "HLA-G", 
         "MMP9", "PLAU")

tiff(paste0(Path, "/Output/PP-dot.tiff"), width = 6, height = 4, units = "in", res = 600)
  DotPlot(MGD, features = PP)
dev.off()

tiff(paste0(Path, "/Output/Endocrine-dot.tiff"), width = 6, height = 4, units = "in", res = 600)
  DotPlot(MGD, features = Endocrine)
dev.off()

tiff(paste0(Path, "/Output/Foregut-dot.tiff"), width = 6, height = 4, units = "in", res = 600)
  DotPlot(MGD, features = AF)
dev.off()

tiff(paste0(Path, "/Output/Duct-dot.tiff"), width = 6, height = 4, units = "in", res = 600)
  DotPlot(MGD, features = Duct)
dev.off()

tiff(paste0(Path, "/Output/Intestinal-dot.tiff"), width = 6, height = 4, units = "in", res = 600)
  DotPlot(MGD, features = Intestinal)
dev.off()

tiff(paste0(Path, "/Output/Hepatic-dot.tiff"), width = 6, height = 4, units = "in", res = 600)
  DotPlot(MGD, features = Hepatic)
dev.off()

tiff(paste0(Path, "/Output/Endothelial-dot.tiff"), width = 6, height = 4, units = "in", res = 600)
  DotPlot(MGD, features = Endothelial)
dev.off()

tiff(paste0(Path, "/Output/Fibroblast-dot.tiff"), width = 6, height = 4, units = "in", res = 600)
  DotPlot(MGD, features = Fibroblast)
dev.off()

tiff(paste0(Path, "/Output/Neural-dot.tiff"), width = 6, height = 4, units = "in", res = 600)
  DotPlot(MGD, features = Neural)
dev.off()

tiff(paste0(Path, "/Output/EE-dot.tiff"), width = 6, height = 4, units = "in", res = 600)
  DotPlot(MGD, features = EE)
dev.off()

tiff(paste0(Path, "/Output/EE2-dot.tiff"), width = 9, height = 4, units = "in", res = 600)
  DotPlot(MGD, features = EE2, group.by = "Cond") + RotatedAxis()
dev.off()

#### UMAP Plots ####
tiff(paste0(Path, "/Output/UMAP_Merged.tiff"), width = 8, height = 4, units = "in", res = 600)
  DimPlot(MGD, pt.size = 0.4, split.by = "Cond")
dev.off()

#### Feature Plots ####
tiff(paste0(Path, "/Output/SPP1-Feat.tiff"), width = 8, height = 4, units = "in", res = 600)
  FeaturePlot(MGD, features = "SPP1", pt.size = 0.4, split.by = "Cond") & theme(legend.position = c(1, 0.2), legend.key.size = unit(0.15, 'in'))
dev.off()

tiff(paste0(Path, "/Output/AKAP12-Feat.tiff"), width = 8, height = 4, units = "in", res = 600)
  FeaturePlot(MGD, features = "AKAP12", pt.size = 0.4, split.by = "Cond") & theme(legend.position = c(1, 0.2), legend.key.size = unit(0.15, 'in'))
dev.off()

tiff(paste0(Path, "/Output/AREG-Feat.tiff"), width = 8, height = 4, units = "in", res = 600)
  FeaturePlot(MGD, features = "AREG", pt.size = 0.4, split.by = "Cond") & theme(legend.position = c(1, 0.2), legend.key.size = unit(0.15, 'in'))
dev.off()

tiff(paste0(Path, "/Output/MMP2-Feat.tiff"), width = 8, height = 4, units = "in", res = 600)
  FeaturePlot(MGD, features = "MMP2", pt.size = 0.4, split.by = "Cond") & theme(legend.position = c(1, 0.2), legend.key.size = unit(0.15, 'in'))
dev.off()

tiff(paste0(Path, "/Output/MMP9-Feat.tiff"), width = 8, height = 4, units = "in", res = 600)
  FeaturePlot(MGD, features = "MMP9", pt.size = 0.4, split.by = "Cond") & theme(legend.position = c(1, 0.2), legend.key.size = unit(0.15, 'in'))
dev.off()

tiff(paste0(Path, "/Output/GATA3-Feat.tiff"), width = 8, height = 4, units = "in", res = 600)
  FeaturePlot(MGD, features = "GATA3", pt.size = 0.4, split.by = "Cond") & theme(legend.position = c(1, 0.2), legend.key.size = unit(0.15, 'in'))
dev.off()

tiff(paste0(Path, "/Output/TFAP2A-Feat.tiff"), width = 8, height = 4, units = "in", res = 600)
  FeaturePlot(MGD, features = "TFAP2A", pt.size = 0.4, split.by = "Cond") & theme(legend.position = c(1, 0.2), legend.key.size = unit(0.15, 'in'))
dev.off()

tiff(paste0(Path, "/Output/KRT7-Feat.tiff"), width = 8, height = 4, units = "in", res = 600)
  FeaturePlot(MGD, features = "KRT7", pt.size = 0.4, split.by = "Cond") & theme(legend.position = c(1, 0.2), legend.key.size = unit(0.15, 'in'))
dev.off()

tiff(paste0(Path, "/Output/WNT3-Feat.tiff"), width = 8, height = 4, units = "in", res = 600)
  FeaturePlot(MGD, features = "WNT3", pt.size = 0.4, split.by = "Cond") & theme(legend.position = c(1, 0.2), legend.key.size = unit(0.15, 'in'))
dev.off()

tiff(paste0(Path, "/Output/WNT5A-Feat.tiff"), width = 8, height = 4, units = "in", res = 600)
  FeaturePlot(MGD, features = "WNT5A", pt.size = 0.4, split.by = "Cond") & theme(legend.position = c(1, 0.2), legend.key.size = unit(0.15, 'in'))
dev.off()

tiff(paste0(Path, "/Output/WNT6-Feat.tiff"), width = 8, height = 4, units = "in", res = 600)
  FeaturePlot(MGD, features = "WNT6", pt.size = 0.4, split.by = "Cond") & theme(legend.position = c(1, 0.2), legend.key.size = unit(0.15, 'in'))
dev.off()

tiff(paste0(Path, "/Output/WNT9B-Feat.tiff"), width = 8, height = 4, units = "in", res = 600)
  FeaturePlot(MGD, features = "WNT9B", pt.size = 0.4, split.by = "Cond") & theme(legend.position = c(1, 0.2), legend.key.size = unit(0.15, 'in'))
dev.off()

tiff(paste0(Path, "/Output/BMP4-Feat.tiff"), width = 8, height = 4, units = "in", res = 600)
  FeaturePlot(MGD, features = "BMP4", pt.size = 0.4, split.by = "Cond") & theme(legend.position = c(1, 0.2), legend.key.size = unit(0.15, 'in'))
dev.off()

tiff(paste0(Path, "/Output/BMP5-Feat.tiff"), width = 8, height = 4, units = "in", res = 600)
  FeaturePlot(MGD, features = "BMP5", pt.size = 0.4, split.by = "Cond") & theme(legend.position = c(1, 0.2), legend.key.size = unit(0.15, 'in'))
dev.off()

tiff(paste0(Path, "/Output/BMP7-Feat.tiff"), width = 8, height = 4, units = "in", res = 600)
  FeaturePlot(MGD, features = "BMP7", pt.size = 0.4, split.by = "Cond") & theme(legend.position = c(1, 0.2), legend.key.size = unit(0.15, 'in'))
dev.off()

tiff(paste0(Path, "/Output/GDF6-Feat.tiff"), width = 8, height = 4, units = "in", res = 600)
  FeaturePlot(MGD, features = "GDF6", pt.size = 0.4, split.by = "Cond") & theme(legend.position = c(1, 0.2), legend.key.size = unit(0.15, 'in'))
dev.off()

#### Violin Plots ####
PP2 <- c("PDX1", "NKX6-1", "GATA4", "SOX9")
PP3 <- c("CTRB2", "PTF1A", "CPA1", "CPA2")

tiff(paste0(Path, "/Output/PP2-Violin.tiff"), width = 8, height = 8, units = "in", res = 600)
  VlnPlot(MGD, features = PP2, pt.size = 0, group.by = "Cond", ncol = 2) 
dev.off()

tiff(paste0(Path, "/Output/PP3-Violin.tiff"), width = 8, height = 8, units = "in", res = 600)
  VlnPlot(MGD, features = PP3, pt.size = 0, group.by = "Cond", ncol = 2, idents = 0) 
dev.off()


#### End of Code ####
