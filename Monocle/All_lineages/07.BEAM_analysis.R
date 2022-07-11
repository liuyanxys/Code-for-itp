rm(list=ls())
###BEAM 
library(monocle)
library(Seurat)
library(dplyr)
setwd("/home/monocle2/")

load('my_cds_subset_yoursample_1100.Rdata')
expressed_genes <- row.names(subset(fData(my_cds_subset),
                                    num_cells_expressed >= 10))
BEAM_res <- BEAM(my_cds_subset[expressed_genes, ], branch_point =1, cores = 4)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
BEAM_res <-subset(BEAM_res,qval < 1e-4)
pData(my_cds_subset)<-pData(my_cds_subset)[order(pData(my_cds_subset)$State,decreasing=F),]
my_cds_subset$State<-factor(my_cds_subset$State,levels=c("1","2","3"), ordered=TRUE)
BEAMgene<-row.names(BEAM_res)[1:1000]
pdf("BEAM-yoursample-branch1-6clusters.pdf", width = 10, height = 12)
plot<-plot_genes_branched_heatmap(my_cds_subset[BEAMgene,],
                                          branch_point = 1,
                                          num_clusters = 6,
                                          cores = 4,
                                          use_gene_short_name = T,
                                          return_heatmap  =  T, 
                                          #branch_colors = c("#979797", "#F05662", "#7990C8"), #pre-branch, Cell fate 1, Cell fate 2 --color
                                          show_rownames = F)
dev.off()
save(plot,file="BEAM-yoursample-branch1-6clusters.Rdata")
annotation_row<-plot$annotation_row
write.table(annotation_row, "BEAM-yoursample-branch1-6clusters-annotation_row.csv", row.names=T, sep='\t', append=F, quote=F)
 
###GO enrichment for BEAM
annotation_row$Gene <- rownames(annotation_row)
gene_group=annotation_row
library(clusterProfiler)
library(org.Hs.eg.db)
allcluster_go=data.frame()
for (i in unique(gene_group$Cluster)) {
  small_gene_group=filter(gene_group,gene_group$Cluster==i)
  df_name=bitr(small_gene_group$Gene, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
  go <- enrichGO(gene         = unique(df_name$ENTREZID),
                 OrgDb         = org.Hs.eg.db,
                 keyType       = 'ENTREZID',
                 ont           = "ALL",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)
  go_res=go@result
  if (dim(go_res)[1] != 0) {
    go_res$cluster=i
    allcluster_go=rbind(allcluster_go,go_res)
  }
}
head(allcluster_go[,c("ID","Description","qvalue","cluster")])
go<-allcluster_go[,c("ID","Description","qvalue","cluster")]
head(subset(go,go$cluster=="1"))
write.table(go, "monocle2BEAM-yoursample-branch1-6clusters-annotation_row-GO.csv", row.names=T, sep='\t', append=T, quote=F)
