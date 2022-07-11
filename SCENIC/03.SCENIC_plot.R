rm(list=ls())
library(SCENIC)
library(AUCell)
library(ComplexHeatmap)
setwd('/home/SCENIC/results/')
regulonAUC <- importAUCfromText("yoursample.auc.csv")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
SO<-readRDS("/home/yoursample.rds") 
macro.cellinfo$seurat_clusters <-macro.cellinfo$cell.types
regulonActivity_byCellType <- sapply(split(rownames(macro.cellinfo), macro.cellinfo$seurat_clusters),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType<-t(regulonActivity_byCellType)
regulonActivity_byCellType_Scaled <- t(scale(regulonActivity_byCellType, center = T, scale=T))
rownames(regulonActivity_byCellType_Scaled)<-gsub("\\(","",rownames(regulonActivity_byCellType_Scaled))
rownames(regulonActivity_byCellType_Scaled)<-gsub("\\)","",rownames(regulonActivity_byCellType_Scaled))
rownames(regulonActivity_byCellType_Scaled)<-gsub("\\+","",rownames(regulonActivity_byCellType_Scaled))
sel<-c()#select topreg
regulonActivity_byCellType_Scaled <-regulonActivity_byCellType_Scaled[sel,]
pdf(file="SCENIC-heatmap.pdf", width = 8, height = 5) 
pheatmap::pheatmap(mat_cluster, fontsize_row = 5, fontsize_col = 12, 
                   color=colorRampPalette(c("#00a9ff","white","#F8766D"))(100), breaks=seq(-3, 3, length.out = 100),
                   cluster_rows = T,cluster_cols=F, border_color=NA)
dev.off()

##preB3 for example
regulonAUC <- importAUCfromText("yoursample.auc.csv") 
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
SO<-readRDS("/home/yoursample.rds") 
macro.cellinfo<-SO@meta.data
macro.cellinfo$seurat_clusters <-macro.cellinfo$cell.types
cellInfo2<-subset(macro.cellinfo,macro.cellinfo$seurat_clusters=="preB3")
regulonActivity_byCellType <- sapply(split(rownames(cellInfo2), cellInfo2$orig.ident),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- regulonActivity_byCellType
regulonActivity_byCellType_Scaled <- t(regulonActivity_byCellType_Scaled)
colnames(regulonActivity_byCellType_Scaled)<-gsub("\\(","",colnames(regulonActivity_byCellType_Scaled))
colnames(regulonActivity_byCellType_Scaled)<-gsub("\\)","",colnames(regulonActivity_byCellType_Scaled))
colnames(regulonActivity_byCellType_Scaled)<-gsub("\\+","",colnames(regulonActivity_byCellType_Scaled))
rownames(regulonActivity_byCellType_Scaled)<-c("HC_1" , "HC_2",  "HC_3" , "HC_4" , "ITP_1", "ITP_2", "ITP_3" ,"ITP_4")
annotation_row = data.frame(group = c("HC","HC","HC","HC","ITP", "ITP","ITP", "ITP"))
rownames(annotation_row)<-rownames(regulonActivity_byCellType_Scaled)
ann_colors =  list(group = c("HC" = "#619cff", "ITP" = "#f8766d"))
##NKT
#selTF<-c("PKNOX1", "GMEB1", "MTF1", "BHLHE41","ZFX", "ZNF879")
##preB1
#selTF<-c("IRF3", "ZNF219", "ZNF787", "ARID3A","TAL1", "HOXA5","IRF5", "TFEC", "GATA2", "HOXA3",  "CREB3L2",  "GATA1", "MITF", "HOXB5", "CEBPE", "STAT4", "ZNF781","FOXJ1","MTF1",  "TP53",  "ZNF22","BACH2", "ATF4", "ZFX")
##preB3
selTF<-c("CREB3L2","CEBPD", "FOSL2", "CEBPE", "CEBPB",  "MITF","KLF5", "JDP2", "CEBPA", "ASCL2", "HSF1","NFKB2", "FOSB", "KLF1",  "ZNF781","NFE2", "ZNF787","BHLHE41",  "PRDM1",  "HES6", "E2F4",  "NR1I2", "FOXJ1", "IRF4",  "IRF3", "ARID3A", "NR1H2", "EGR3", "MXD1", "FOSL1",  "MTF1", "STAT4", "HOXA5",  "TFDP1",  "IRF5",  "BATF", "GATA1", "TFEC","GATA2","E2F2", "TCF3", "ELK3", "FOXO1", "PAX5", "LEF1", "TCF4", "VEZF1", "EBF1", "ZFX", "NR3C1", "BACH2", "CREB1", "SP3", "BCL11A", "REST", "MEF2A", "ARID5B", "ZNF423", "MLXIP", "SP1", "MGA", "MYBL1","ZNF217","TCF12", "IKZF1", "RREB1", "STAT5B", "KLF12", "MEF2D",  "POU2F1",  "HIVEP2", "ZBTB14", "NFYA", "ARID5A")

regulonActivity_byCellType_Scaled<-regulonActivity_byCellType_Scaled[,selTF]
pdf(file="heatmap-yoursample-ITPversusHC-preB3.pdf", width =15, height =4)
pheatmap::pheatmap(regulonActivity_byCellType_Scaled, #fontsize_row=3, 
                   color=colorRampPalette(c("#00a9ff","white","#F8766D"))(100), breaks=seq(-3, 3, length.out = 100),
                   border_color=NA,
                   cluster_rows = F,cluster_cols=F,
                   annotation_row = annotation_row,
                   annotation_colors = ann_colors,
                   scale="column",
                   main = "preB3")
dev.off()
write.table(t(regulonActivity_byCellType_Scaled),file='heatmap-yoursample-ITPversusHC-preB3.csv',row.names = TRUE,quote = FALSE,sep='\t')

##Transpose
regulonActivity_byCellType_Scaled<-regulonActivity_byCellType_Scaled[,selTF]
regulonActivity_byCellType_Scaled<-t(regulonActivity_byCellType_Scaled)
annotation_col = data.frame(group = c("HC","HC","HC","HC","ITP", "ITP","ITP", "ITP"))
rownames(annotation_col)<-colnames(regulonActivity_byCellType_Scaled)
ann_colors =  list(group = c("HC" = "#619cff", "ITP" = "#f8766d"))
pdf(file="heatmap-yoursample-ITPversusHC-preB3_Transpose.pdf", width =4, height =3)
pheatmap::pheatmap(regulonActivity_byCellType_Scaled, #fontsize_row=3, 
                   color=colorRampPalette(c("#00a9ff","white","#F8766D"))(100), breaks=seq(-3, 3, length.out = 100),
                   border_color=NA,
                   cluster_rows = F,cluster_cols=F,
                   annotation_col = annotation_col,
                   annotation_colors = ann_colors,
                   scale="row",
                   main = "preB3")
dev.off()
