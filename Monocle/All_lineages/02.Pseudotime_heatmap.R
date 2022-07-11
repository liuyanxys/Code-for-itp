rm(list=ls())
###pseudotime heatmap
library(monocle)
library(Seurat)
library(dplyr)
setwd("/home/monocle2/")

load('my_cds_subset_yoursample_1100.Rdata')
head(pData(my_cds_subset))
my_pseudotime_de <- differentialGeneTest(my_cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 22)
write.table(my_pseudotime_de,file='my_pseudotime_de_yoursample_1100.csv', row.names=T, sep='\t', append=T, quote=F)
my_pseudotime_de<-read.csv('my_pseudotime_de_yoursample_1100.csv', sep='\t')
my_pseudotime_de %>% arrange(qval) %>% head()
my_pseudotime_de %>% arrange(qval) %>% head() %>% dplyr::select(gene_short_name) -> my_pseudotime_gene
#sig_gene_names <- row.names(subset(my_pseudotime_de, qval < 0.1))
my_pseudotime_de <- my_pseudotime_de[order(my_pseudotime_de$qval),]
sig_gene_names <- row.names(my_pseudotime_de[1:2000,])
pdf("pseudotime_heatmap_2000_yoursample-6cluster.pdf", width = 12, height = 14)
plot<-plot_pseudotime_heatmap(my_cds_subset[sig_gene_names,],
                num_clusters = 6,
                cores = 4,
                show_rownames = F,return_heatmap =T )
dev.off()
save(plot,file="pseudotime_heatmap_2000_yoursample-6cluster-plot.Rdata")

t <- as.data.frame(cutree(plot$tree_row, k=6))
colnames(t) <- "Cluster"
t$Gene <- rownames(t)
write.table(t, "pseudotime-yoursample-6clusters-annotation_row.csv", row.names=F, sep='\t', append=F, quote=F)

###GO enrichment for pseudotime heatmap
load("pseudotime_heatmap_2000_yoursample-6cluster-plot.Rdata")
t <- as.data.frame(cutree(plot$tree_row, k=6))
colnames(t) <- "Cluster"
t$Gene <- rownames(t)
gene_group=t
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
write.table(go, "pseudotime-yoursample-6clusters-annotation_row-GO.csv", row.names=T, sep='\t', append=T, quote=F)
