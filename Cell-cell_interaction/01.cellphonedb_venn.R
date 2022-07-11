rm(list=ls())
options(stringsAsFactors = F)

library(Seurat)
library(patchwork)
library(harmony)
library(reshape2)
library(RColorBrewer)
library(scales)
library(ggplot2)
library(ggpubr)
library(ggplotify)
library(pheatmap)
library(reshape2)
library(stringr)
library(dplyr)
library(cowplot)
library(ggrepel)
library(magrittr)
library(VennDiagram)

setwd('/home/cellphonedb/out_yoursample/venn')
set.seed(123456789)
dpi =300

###preB3 for example
preB3_ITP<-read.table("preB3_ITP.txt",sep='\t')
preB3_HC<-read.table("preB3_HC.txt",sep='\t')

preB3_ITP<-preB3_ITP[,1]
preB3_HC<-preB3_HC[,1]

library(VennDiagram)
venn_list <- list("preB3_ITP" = preB3_ITP, "preB3_HC" = preB3_HC)
venn.diagram(venn_list,
             filename = 'venn_ITP&HC_preB3.png', imagetype = "png",
             col = "black",
             fill = c("#DE5678", "#4F7FDC"),
             alpha = c(0.8,0.6),
             cex = 0.8,
             cat.col = 'black',
             cat.cex = 0.8,
             cat.fontface = "bold",
             margin = 0.05,
             main = "Complex Venn Diagram",
             main.cex = 1.2
)
inter <- get.venn.partitions(venn_list)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
inter <- apply(inter,2,as.character)
write.table(inter, 'venn_ITP&HC_preB3_inter.txt', row.names = FALSE, sep = '\t', quote = FALSE)

venn.plot <-venn.diagram(venn_list,
             filename = NULL, 
             col = "black",
             fill = c("#DE5678", "#4F7FDC"),
             alpha = c(0.8,0.6),
             cex = 0.8,
             cat.col = 'black',
             cat.cex = 0.8,
             cat.fontface = "bold",
             margin = 0.05,
             main = "venn_ITP&HC_preB3",
             main.cex = 1.2
)
pdf(file="venn_ITP&HC_preB3.pdf")
grid.draw(venn.plot)
dev.off()
