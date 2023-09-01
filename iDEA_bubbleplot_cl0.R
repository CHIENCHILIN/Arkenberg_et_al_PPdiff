#### iDEA ####
library(tidyverse)
library(scRNAseq)
library(matrixStats)
library(iDEA)
library(zingeR)
library(DESeq2)
library(BiocParallel)
library(data.table)

message("iDEA Library loaded")

res <- fread(file = "C:/Users/.../Output/dim12_res0.8_cl0_CTRGEL_MGD_DESeq2.csv")
GS <- readRDS(file = "C:/Users/.../Output/CompleteSets.rds")

message("Data loaded")

de.summary <- data.frame(log2FoldChange = res$log2FoldChange,
                         lfcSE2 = res$logFC_SE^2,
                         row.names = res$gene)

de.summary.up <- de.summary[de.summary$log2FoldChange > 0, ]

idea <- CreateiDEAObject(de.summary.up, GS, num_core = 1); message("CreateiDEAObject done")
idea <- iDEA.fit(idea, modelVariant = F); message("iDEA.fit done")
idea <- iDEA.louis(idea); message("iDEA.louis done")

write.csv(idea@gsea, "C:/Users/.../Output/iDEA_cl0_GEL-CTR(r).csv")

message("Data saved")

##-------------------------------------------------------------
## Simulation Plot: Bubble plot
##-------------------------------------------------------------
library(data.table)
library(iDEA)
library(tidyverse)
library(stringr)
library(RColorBrewer)

Path = "C:/ExampleFolder"

# Building the GSEA + Gene Set info dataframe:

customGeneSetsInfo <- readRDS(file = "C:/Users/.../Output/customGeneSetsInfo.rds")

### The initial DESeq2 analysis seems to switched the up- and down- regulated genes, which results in switched GSEA results.
# The following script converts the "up" GSEA results into "down" files while adding gene set info, and vice versa.

plotdata<- fread(paste0(Path, "/Output/iDEA_cl0_GEL-CTR(r).csv"))
# plotdata<- fread("C:/Users/.../Output/GSEA/iDEA_dn.csv")
bp.data <- merge.data.table(plotdata, customGeneSetsInfo, by.x = "annot_id", by.y = "gset") #I think this is easier in SQL

bp.data$Category <- droplevels(bp.data$gsetBioName)

#### Bubble plot 

includedCats <- c("GO BIOLOGICAL PROCESS",
                  "GO MOLECULAR FUNCTION",
                  "GO CELLULAR COMPONENT",
                  "REACTOME",
                  "KEGG",
                  "PID",
                  "POSITIONAL",
                  "TRANSCRIPTION FACTORS")

caseCats <- c("GO biological process",
              "GO molecular function",
              "GO cellular component",
              "Reactome",
              "KEGG",
              "PID",
              "Positional",
              "Transcription factors")

bp.data <- bp.data[toupper(bp.data$Category) %in% includedCats, ] #subset data, removing all geneset categories not in includedCats
bp.data$Category <- droplevels(bp.data$Category) %>% factor(caseCats) # Drop unused levels, reorder factor so that plot legend matches x position
bp.data <- bp.data[order(match(toupper(bp.data$Category), includedCats)), ] #Organizes the data so that categories are grouped on the x axis
bp.data$IDNum <- row.names(bp.data) %>% as.integer() #IDNum will define x position on bubble plot
bp.data$Log10_Pvalue_Louis <- -1*log10(bp.data$pvalue_louis)

write.csv(bp.data, paste0(Path, "/Output/iDEA_cl0_GEL-CTR(r)_iDEAup.full.csv"))


##-------------------------------------------------------------
## Configurable Variables
##-------------------------------------------------------------


# Colors
bp.colors <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5")

labeled.genesets <- c("REACTOME_TRANSLATION",
                      "REACTOME_CELL_CYCLE_MITOTIC",
                      "HALLMARK_G2M_CHECKPOINT",
                      "REACTOME_CELL_CYCLE",
                      "REACTOME_MITOTIC_M_M_G1_PHASES",
                      "GO_DNA_REPLICATION",
                      "GO_MITOTIC_NUCLEAR_DIVISION",
                      "GO_CELL_CYCLE_PHASE_TRANSITION",
                      "GO_CELL_DIVISION",
                      "KEGG_CELL_CYCLE",
                      "KEGG_DNA_REPLICATION",
                      "TEAD2_TARGET_GENES",
                      "GO_REGULATION_OF_CELL_DIVISION",
                      "GO_REGULATION_OF_EPITHELIAL_CELL_MIGRATION",
                      "GO_NEGATIVE_REGULATION_OF_WNT_SIGNALING_PATHWAY",
                      "GO_REGULATION_OF_CANONICAL_WNT_SIGNALING_PATHWAY",
                      "GO_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY",
                      "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                      "GO_POSITIVE_REGULATION_OF_CANONICAL_WNT_SIGNALING_PATHWAY",
                      "PID_BETA_CATENIN_NUC_PATHWAY",
                      "GO_POSITIVE_REGULATION_OF_CELL_DIVISION",
                      "GO_POSITIVE_REGULATION_OF_WNT_SIGNALING_PATHWAY",
                      "KEGG_ECM_RECEPTOR_INTERACTION",
                      "REACTOME_SIGNALING_BY_WNT",
                      "REACTOME_INTEGRIN_CELL_SURFACE_INTERACTIONS",
                      "HALLMARK_WNT_BETA_CATENIN_SIGNALING")


included.categories <- c("GO biological process",
                         "GO molecular function",
                         "GO cellular component",
                         "Reactome",
                         "KEGG",
                         "PID",
                         "Positional",
                         "Transcription factors")


# The following step should not be performed with the current bp.data structure, as the first column is the annot_id column.

# You have to index when loading the data in (e.g. bp.data[,-1]). 
# This removes an extra column at the start of the dataframe
# If it's not removed, you will see this error "Error: `data` must be uniquely named but has duplicate columns"

#data.to.plot <- bp.data[ , -1]

data.to.plot <- bp.data

#Write title for the bubble plot
# plot.title <- expression(paste("Upregulated genes in ", italic("HC"), " cluster"))


##-------------------------------------------------------------
## Create data frame for for plot labels
##-------------------------------------------------------------

#### Select the biologically significant gene set in the bubbleplot to show
Sig <- data.to.plot[which(data.to.plot$annot_id %in% labeled.genesets),]


#### You can modify your text label here (here we're converting the geneset names to lowercase)
Sig$Term <- tolower(Sig$annot_id)

##-------------------------------------------------------------
## ggplot function for bubble plot
##-------------------------------------------------------------
library(ggplot2)
library(ggrepel)


bp <- ggplot(data.to.plot, aes(x = IDNum, y = Log10_Pvalue_Louis, color = Category)) +
  geom_point(shape = 19, alpha=1, size = 4) + 
  labs(x = "",
       y = expression(paste(bold(-log[10]), bold("("), bolditalic(p), bold("-value)")))) +
  ggtitle(label = "Cluster 0; upregulated in GEL (vs. CTR)") + 
  theme(plot.margin = margin(1, 1, 1, 1, "cm"),
        plot.title = element_text(size = 30, face="bold"),
        axis.text.y = element_text(size = 30),
        axis.text.x = element_blank(),
        axis.line = element_line(colour = 'black'),
        axis.ticks = element_line(colour = 'grey80'),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 40, face = 'bold'),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
        legend.title=element_text(size=16,face = 'bold'),
        legend.text=element_text(size=16),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = 'white'),
        plot.background = element_blank(),
        legend.key = element_rect(color = "transparent", fill = "transparent")) +
  geom_hline(yintercept = 1.82, col = 'black', linetype = 2, size=2) +
  scale_color_manual(values = bp.colors) +
  scale_fill_manual(values = bp.colors) +
  theme(legend.direction = "vertical") +
  theme(legend.position = c(0.8, 0.95)) +
  theme(legend.box = "horizontal") +
  theme(legend.title.align = 0) +
  geom_label_repel(
    data = Sig,
    aes(label = Term),
    col = 'black',
    size = 6,
    nudge_y = 6,
    max.overlaps = 20)


setEPS(reset = T)
setEPS(width = 20, height = 20)
postscript(paste0(Path, "/Output/iDEA-bubble_cl0_GEL-CTR(r).eps"))
bp
dev.off()

ggsave(filename = paste0(Path, "/Output/iDEA-bubble_cl0_GEL-CTR(r).eps"),
       plot = bp,
       width = 20,
       height = 20/1.66667,
       units = "in",
       dpi = 300,
       device = "eps")



