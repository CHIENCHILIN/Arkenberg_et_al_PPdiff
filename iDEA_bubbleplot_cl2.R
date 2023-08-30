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

customGeneSetsInfo <- readRDS(file = "C:/Users/yosueda/Documents/RStudio/SOX2 project/Combine/20220112/Output/customGeneSetsInfo.rds")

### The initial DESeq2 analysis seems to switched the up- and down- regulated genes, which results in switched GSEA results.
# The following script converts the "up" GSEA results into "down" files while adding gene set info, and vice versa.

plotdata<- fread(paste0(Path, "20230112/Output/iDEA_cl2_GEL-CTR(r).csv"))
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

write.csv(bp.data, paste0(Path, "20230112/Output/iDEA_cl2_GEL-CTR(r)_iDEAup.full.csv"))


##-------------------------------------------------------------
## Configurable Variables
##-------------------------------------------------------------


# Colors
bp.colors <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5")

labeled.genesets <- c("HALLMARK_INTERFERON_ALPHA_RESPONSE",
                      "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
                      "GO_EMBRYONIC_ORGAN_MORPHOGENESIS",
                      "GO_MORPHOGENESIS_OF_A_BRANCHING_STRUCTURE",
                      "GO_HEART_MORPHOGENESIS",
                      "GO_REGULATION_OF_EPITHELIAL_CELL_PROLIFERATION",
                      "GO_TUBE_MORPHOGENESIS",
                      "GO_MORPHOGENESIS_OF_AN_EPITHELIUM",
                      "GO_BRANCHING_MORPHOGENESIS_OF_AN_EPITHELIAL_TUBE",
                      "GO_REGULATION_OF_EPITHELIAL_CELL_MIGRATION",
                      "GO_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_PROLIFERATION",
                      "GO_REGULATION_OF_BMP_SIGNALING_PATHWAY",
                      "GO_POSITIVE_REGULATION_OF_CELL_MORPHOGENESIS_INVOLVED_IN_DIFFERENTIATION",
                      "GO_REGULATION_OF_NOTCH_SIGNALING_PATHWAY",
                      "GO_POSITIVE_REGULATION_OF_CELL_ADHESION",
                      "GO_DIGESTIVE_SYSTEM_DEVELOPMENT",
                      "GO_REGULATION_OF_CANONICAL_WNT_SIGNALING_PATHWAY",
                      "GO_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY",
                      "GO_RESPONSE_TO_TRANSFORMING_GROWTH_FACTOR_BETA",
                      "GO_RESPONSE_TO_BMP",
                      "GO_POSITIVE_REGULATION_OF_WNT_SIGNALING_PATHWAY",
                      "PID_BMP_PATHWAY",
                      "KEGG_HEDGEHOG_SIGNALING_PATHWAY",
                      "GO_REGULATION_OF_WNT_SIGNALING_PATHWAY",
                      "GO_SMAD_PROTEIN_SIGNAL_TRANSDUCTION",
                      "HALLMARK_HEDGEHOG_SIGNALING",
                      "KEGG_TGF_BETA_SIGNALING_PATHWAY",
                      "GO_INTEGRIN_MEDIATED_SIGNALING_PATHWAY",
                      "GO_NEGATIVE_REGULATION_OF_BMP_SIGNALING_PATHWAY",
                      "GO_POSITIVE_REGULATION_OF_CANONICAL_WNT_SIGNALING_PATHWAY",
                      "PID_AVB3_INTEGRIN_PATHWAY",
                      "GO_CANONICAL_WNT_SIGNALING_PATHWAY",
                      "GO_POSITIVE_REGULATION_OF_DEVELOPMENTAL_GROWTH",
                      "GO_WNT_SIGNALING_PATHWAY",
                      "GO_POSITIVE_REGULATION_OF_EPITHELIAL_CELL_DIFFERENTIATION",
                      "PID_WNT_SIGNALING_PATHWAY",
                      "WNT_SIGNALING")



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
  ggtitle(label = "Cluster 2; upregulated in GEL (vs. CTR)") + 
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
postscript(paste0(Path, "20230112/Output/iDEA-bubble_cl2_GEL-CTR(r).eps"))
bp
dev.off()

ggsave(filename = paste0(Path, "20230112/Output/iDEA-bubble_cl2_GEL-CTR(r).eps"),
       plot = bp,
       width = 20,
       height = 20/1.66667,
       units = "in",
       dpi = 300,
       device = "eps")
