library(Seurat, quietly = TRUE)
library(SeuratObject, quietly = TRUE)

Path <- "C:/ExampleFolder"
data <- readRDS(paste0(Path, "/dim12_res0.8.rds"))

#### for cluster 2, comparison of CTR vs. GEL
MGD <- subset(data, ident = c("2"))#, subset = (Cond == "CTR" | Cond == "GEL"))

#%%%%%%%%%%%%%%%%%%%%%%%%%% DEA with DESeq2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("scRNAseq")

install.packages("remotes")
remotes::install_github("statOmics/zingeR")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

library(tidyverse, quietly = TRUE)
library(scRNAseq, quietly = TRUE)
library(matrixStats, quietly = TRUE)
library(zingeR, quietly = TRUE)
library(DESeq2, quietly = TRUE)
library(BiocParallel, quietly = TRUE)

# register(MulticoreParam(8)) # only used when working on linux/OSX/hpc
# print("Cores = Registered")

#Create Summarized Experiment object
MGD_sce <- SummarizedExperiment(as.matrix(MGD@assays$RNA@counts),
                                colData = MGD@meta.data)

rm(data, MGD); gc(); gc() 
#print("SCE created, Seurat removed")

#Filter to only look at highly expressed genes
filter <- rowSums(assay(MGD_sce)) > 5
MGD_sce <- MGD_sce[filter,]
vars <- assay(MGD_sce) %>% log1p() %>% rowVars()
names(vars) <- rownames(MGD_sce)
vars <- sort(vars, decreasing = T)

#Subset to only look at the top 15,000 genes
MGD_sce <- MGD_sce[names(vars[1:15000]),]
assayNames(MGD_sce)[1] <- "counts"

#Construct Inputs for DESeq2 Test
MGD.colData <- data.frame(Cond = colData(MGD_sce)$Cond)
iCounts <- assay(MGD_sce, i = "counts")

#Set factor levels for variables of interest
MGD.colData$Cond <- as.factor(MGD.colData$Cond)

#Ensure that the control group is the first factor level
MGD.colData$Cond <- relevel(MGD.colData$Cond, ref = "CTR")

design <- model.matrix(~MGD.colData$Cond)

dse <- DESeqDataSetFromMatrix(countData = iCounts, colData = MGD.colData, design = ~ Cond)

rm(MGD_sce); gc(); gc()

weights <- zingeR::zeroWeightsLS(counts = iCounts, 
                                design = design,
                                maxit = 500, normalization = "DESeq2_poscounts",
                                colData = MGD.colData,
                                designFormula = ~ Cond, 
                                verbose = TRUE)

assays(dse)[["weights"]] <- weights

dse <- DESeq2::DESeq(dse,
                     test = "Wald", 
                     sfType = "poscounts", 
                     useT = TRUE, 
                     betaPrior = TRUE, 
                     modelMatrixType="expanded",
                     parallel = FALSE) # Can be parallelized on OSX/Linux
# dse <- DESeq2::estimateSizeFactors(dse, type="poscounts")
# dse <- estimateDispersions(dse)
# dse <- nbinomWaldTest(dse, modelMatrixType="standard", betaPrior=TRUE, useT=TRUE, df=rowSums(weights)-(length(unique(d20.colData$Condition))-1) )

#Save DESeq results
full.results <- results(dse, parallel = FALSE)
full.results.df <- data.frame(gene = row.names(full.results),
                              log2FoldChange = full.results$log2FoldChange,
                              pvalue = full.results$pvalue,
                              padj = full.results$padj,
                              baseMean = full.results$baseMean,
                              logFC_SE = full.results$lfcSE)
write.csv(full.results.df, file = "C:/Users/markenbe/OneDrive - Indiana University/7- R/Data/Pancreatic/Output/dim12_res0.8_cl2_CTRGEL_MGD_DESeq2.csv")


#%%%%%%%%%%%%%%%%%%%%%%%%%% Volcano %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(Seurat)
library(tidyverse)
library(EnhancedVolcano)

#Load the data
data <- readRDS(paste0(Path, "Data/dim12_res0.8.rds"))

# Select cluster #2
data <- subset(data, idents = "2")
FeaturePlot(data, "MMP9", pt.size = 1, max.cutoff = 3)
FeaturePlot(data, "MMP9", pt.size = 1, max.cutoff = 3, split.by = "Cond")

DimPlot(data, pt.size = 1)    # checking loaded data
VlnPlot(data, "MMP9")         # checking loaded data

# MMP9(+) = MMP9 > 0
Idents(data, cells = names(Idents(subset(data, subset = ( MMP9 > 0 ))))) <- "2.1"

DimPlot(data, pt.size = 1)    # checking cluster separation
VlnPlot(data, "MMP9")         # checking cluster separation

# # MMP9(+) = MMP9 >= 1
# Idents(data, cells = names(Idents(subset(data, subset = ( MMP9 >= 1 ))))) <- "2.1"
# 
# DimPlot(data, pt.size = 1)    # checking cluster separation
# VlnPlot(data, "MMP9")         # checking cluster separation


#Function for outputting a volcano plot
VOLC <- function(object, interest = NULL, control = NULL, subset = NULL, group = NULL, tiff = F) {
  Wilcox_tab <<- NULL
  Wilcox_tab <<- FindMarkers(object,
                             group.by = group,
                             ident.1 = interest,
                             ident.2 = control,
                             subset.ident = subset,
                             logfc.threshold = 0, test.use = "wilcox", min.pct = 0.1,
                             verbose = T, only.pos = F, return.thresh = 0.01) %>%
    select(c("pct.1", "pct.2", "p_val_adj", "avg_log2FC")) %>%
    arrange(desc(avg_log2FC)) %>%
    mutate(pval = -log10(p_val_adj))
  
  if (tiff == T) {
    tiff(paste0(Path, "Output/volc_", interest, "_vs_", control, "_in_cluster_", subset, ".tiff"), 
         width = 10, height = 8, units = "in", res = 300)
  }
  print(
    EnhancedVolcano(Wilcox_tab,
                    lab = rownames(Wilcox_tab),
                    x = 'avg_log2FC',
                    y = 'p_val_adj',
                    xlim = c(-3, 3),
                    title = "", subtitle = "",
                    pointSize = 4,
                    pCutoff = 1e-50,
                    colAlpha = 0.4,
                    col=c('grey60', 'grey60', 'grey60', 'red2'),
                    legendPosition = "none",
                    # legendLabels = c("NS", expression(Log[2] ~ FC), "p-value", expression(p - value ~ and ~ log[2] ~ FC)),
                    # legendLabSize = 10,
                    # legendIconSize = 5,
                    drawConnectors = T,
                    widthConnectors = 0.2,
                    lengthConnectors = unit(0.01, "npc"),
                    colConnectors = "grey10",
                    directionConnectors = "both",
                    arrowheads = F,
                    maxoverlapsConnectors = 80)#+ coord_flip()
  )
  if (tiff == T) {
    dev.off()
  }
}

VOLC(data, interest = "2.1", control = "2",
     # group = "Cond",
     # subset = "0",
     tiff = F)
