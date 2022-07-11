rm(list=ls())
library(monocle)
library(Seurat)
library(dplyr)
library(sampling)
setwd("/home/monocle2/")
SO<-readRDS("yoursample.rds")

###randomly select half of the cells
cluster0 <- SO@meta.data[SO@meta.data$cell.types=="CLP",]
cluster1 <- SO@meta.data[SO@meta.data$cell.types=="EBMP",]
cluster2 <- SO@meta.data[SO@meta.data$cell.types=="EryP",]
cluster3 <- SO@meta.data[SO@meta.data$cell.types=="GMP",]
cluster4 <- SO@meta.data[SO@meta.data$cell.types=="HSC",]
cluster5 <- SO@meta.data[SO@meta.data$cell.types=="MDP",]
cluster6 <- SO@meta.data[SO@meta.data$cell.types=="MEP",]
cluster7 <- SO@meta.data[SO@meta.data$cell.types=="MkP1",]
cluster8 <- SO@meta.data[SO@meta.data$cell.types=="MkP2",]
cluster9 <- SO@meta.data[SO@meta.data$cell.types=="MPP",]
cluster10 <- SO@meta.data[SO@meta.data$cell.types=="NeuP",]
cluster11 <- SO@meta.data[SO@meta.data$cell.types=="NK/Tp",]
cluster12 <- SO@meta.data[SO@meta.data$cell.types=="preB1",]
cluster13 <- SO@meta.data[SO@meta.data$cell.types=="preB2",]
cluster14 <- SO@meta.data[SO@meta.data$cell.types=="preB3",]

cluster0HC <-cluster0[cluster0$group=="HC",]
cluster1HC <-cluster1[cluster1$group=="HC",]
cluster2HC <-cluster2[cluster2$group=="HC",]
cluster3HC <-cluster3[cluster3$group=="HC",]
cluster4HC <-cluster4[cluster4$group=="HC",]
cluster5HC <-cluster5[cluster5$group=="HC",]
cluster6HC <-cluster6[cluster6$group=="HC",]
cluster7HC <-cluster7[cluster7$group=="HC",]
cluster8HC <-cluster8[cluster8$group=="HC",]
cluster9HC <-cluster9[cluster9$group=="HC",]
cluster10HC <-cluster10[cluster10$group=="HC",]
cluster11HC <-cluster11[cluster11$group=="HC",]
cluster12HC <-cluster12[cluster12$group=="HC",]
cluster13HC <-cluster13[cluster13$group=="HC",]
cluster14HC <-cluster14[cluster14$group=="HC",]

cluster0ITP <-cluster0[cluster0$group=="ITP",]
cluster1ITP <-cluster1[cluster1$group=="ITP",]
cluster2ITP <-cluster2[cluster2$group=="ITP",]
cluster3ITP <-cluster3[cluster3$group=="ITP",]
cluster4ITP <-cluster4[cluster4$group=="ITP",]
cluster5ITP <-cluster5[cluster5$group=="ITP",]
cluster6ITP <-cluster6[cluster6$group=="ITP",]
cluster7ITP <-cluster7[cluster7$group=="ITP",]
cluster8ITP <-cluster8[cluster8$group=="ITP",]
cluster9ITP <-cluster9[cluster9$group=="ITP",]
cluster10ITP <-cluster10[cluster10$group=="ITP",]
cluster11ITP <-cluster11[cluster11$group=="ITP",]
cluster12ITP <-cluster12[cluster12$group=="ITP",]
cluster13ITP <-cluster13[cluster13$group=="ITP",]
cluster14ITP <-cluster14[cluster14$group=="ITP",]

cell_id0HC<- sample(rownames(cluster0HC),(length(rownames(cluster0HC))/2))
cell_id1HC<- sample(rownames(cluster1HC),(length(rownames(cluster1HC))/2))
cell_id2HC<- sample(rownames(cluster2HC),(length(rownames(cluster2HC))/2))
cell_id3HC<- sample(rownames(cluster3HC),(length(rownames(cluster3HC))/2))
cell_id4HC<- sample(rownames(cluster4HC),(length(rownames(cluster4HC))/2))
cell_id5HC<- sample(rownames(cluster5HC),(length(rownames(cluster5HC))/2))
cell_id6HC<- sample(rownames(cluster6HC),(length(rownames(cluster6HC))/2))
cell_id7HC<- sample(rownames(cluster7HC),(length(rownames(cluster7HC))/2))
cell_id8HC<- sample(rownames(cluster8HC),(length(rownames(cluster8HC))/2))
cell_id9HC<- sample(rownames(cluster9HC),(length(rownames(cluster9HC))/2))
cell_id10HC<- sample(rownames(cluster10HC),(length(rownames(cluster10HC))/2))
cell_id11HC<- sample(rownames(cluster11HC),(length(rownames(cluster11HC))/2))
cell_id12HC<- sample(rownames(cluster12HC),(length(rownames(cluster12HC))/2))
cell_id13HC<- sample(rownames(cluster13HC),(length(rownames(cluster13HC))/2))
cell_id14HC<- sample(rownames(cluster14HC),(length(rownames(cluster14HC))/2))

cell_id0ITP<- sample(rownames(cluster0ITP),(length(rownames(cluster0ITP))/2))
cell_id1ITP<- sample(rownames(cluster1ITP),(length(rownames(cluster1ITP))/2))
cell_id2ITP<- sample(rownames(cluster2ITP),(length(rownames(cluster2ITP))/2))
cell_id3ITP<- sample(rownames(cluster3ITP),(length(rownames(cluster3ITP))/2))
cell_id4ITP<- sample(rownames(cluster4ITP),(length(rownames(cluster4ITP))/2))
cell_id5ITP<- sample(rownames(cluster5ITP),(length(rownames(cluster5ITP))/2))
cell_id6ITP<- sample(rownames(cluster6ITP),(length(rownames(cluster6ITP))/2))
cell_id7ITP<- sample(rownames(cluster7ITP),(length(rownames(cluster7ITP))/2))
cell_id8ITP<- sample(rownames(cluster8ITP),(length(rownames(cluster8ITP))/2))
cell_id9ITP<- sample(rownames(cluster9ITP),(length(rownames(cluster9ITP))/2))
cell_id10ITP<- sample(rownames(cluster10ITP),(length(rownames(cluster10ITP))/2))
cell_id11ITP<- sample(rownames(cluster11ITP),(length(rownames(cluster11ITP))/2))
cell_id12ITP<- sample(rownames(cluster12ITP),(length(rownames(cluster12ITP))/2))
cell_id13ITP<- sample(rownames(cluster13ITP),(length(rownames(cluster13ITP))/2))
cell_id14ITP<-sample(rownames(cluster14ITP),(length(rownames(cluster14ITP))/2))

cell_id<-c(cell_id0HC,cell_id1HC,cell_id2HC,cell_id3HC,cell_id4HC,cell_id5HC,cell_id6HC,cell_id7HC,cell_id8HC,cell_id9HC,cell_id10HC,cell_id11HC,cell_id12HC,cell_id13HC,cell_id14HC,
           cell_id0ITP,cell_id1ITP,cell_id2ITP,cell_id3ITP,cell_id4ITP,cell_id5ITP,cell_id6ITP,cell_id7ITP,cell_id8ITP,cell_id9ITP,cell_id10ITP,cell_id11ITP,cell_id12ITP,cell_id13ITP,cell_id14ITP)

mtx_sml <- SO@assays$RNA@counts[, cell_id]
cell_sml <- SO@meta.data[cell_id, ]

data <- as(as.matrix(mtx_sml), 'sparseMatrix') 
pd <- new('AnnotatedDataFrame', data =cell_sml)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
###Construct monocle cds
cds <- newCellDataSet(data,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
list_cluster <-cell_sml$cell.types
names(list_cluster) <- SO@assays[["RNA"]]@data[, cell_id]@Dimnames[[2]]
pData(cds)$cell.types <- list_cluster

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
my_cds <- detectGenes(cds, min_expr = 0.1)
rm(cds)
head(fData(my_cds))
head(pData(my_cds))
pData(my_cds)$UMI <- Matrix::colSums(exprs(my_cds))
disp_table <- dispersionTable(my_cds)
head(disp_table)
table(disp_table$mean_expression>=0.1)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
my_cds <- setOrderingFilter(my_cds, unsup_clustering_genes$gene_id)

expressed_genes <- row.names(subset(fData(my_cds), num_cells_expressed >= 10)) 
my_cds_subset <- my_cds[expressed_genes, ]
rm(my_cds)
my_cds_subset
head(pData(my_cds_subset))
my_cds_subset <- detectGenes(my_cds_subset, min_expr = 0.1)
fData(my_cds_subset)$use_for_ordering <- fData(my_cds_subset)$num_cells_expressed > 0.05 * ncol(my_cds_subset)
table(fData(my_cds_subset)$use_for_ordering)

my_cds_subset <- reduceDimension(my_cds_subset,max_components = 2,norm_method = 'log',num_dim = 10,reduction_method = 'tSNE',verbose = TRUE)

my_cds_subset <- clusterCells(my_cds_subset,rho_threshold = 2,delta_threshold = 10,skip_rho_sigma = T,verbose = FALSE,cores = 20)
table(pData(my_cds_subset)$Cluster)

head(pData(my_cds_subset))
clustering_DEG_genes <- differentialGeneTest(my_cds_subset,fullModelFormulaStr = '~cell.types',cores = 5)
dim(clustering_DEG_genes)

library(dplyr)
clustering_DEG_genes %>% arrange(qval) %>% head()
my_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1100]
my_cds_subset <- setOrderingFilter(my_cds_subset, ordering_genes = my_ordering_genes)
deg <- subset(clustering_DEG_genes, qval < 0.01) 
deg <- deg[order(deg$qval,decreasing=F),]
head(deg)
write.table(deg,file="train.monocle.DEG.csv",col.names=T,row.names=F,sep="\t",quote=F)
my_ordering_genes <-  rownames(deg) [1:1100]
my_cds_subset <- setOrderingFilter(my_cds_subset, ordering_genes = my_ordering_genes)

my_cds_subset <- reduceDimension(my_cds_subset, method = 'DDRTree') 
GM_state <- function(my_cds_subset){
  if (length(unique(pData(my_cds_subset)$State)) > 1){
    T0_counts <- table(pData(my_cds_subset)$State, pData(my_cds_subset)$cell.types)[,"HSC"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

pData(my_cds_subset)<-pData(my_cds_subset)[order(pData(my_cds_subset)$Pseudotime),]
my_cds_subset <- orderCells(my_cds_subset)
my_cds_subset <- orderCells(my_cds_subset,root_state = GM_state(my_cds_subset))

save(my_cds_subset,file='my_cds_subset_yoursample_1100.Rdata')

cols<-c("HSC"="#D9534F","MPP"="#96CEB4", "GMP"="#CBE86B", "NeuP"="#EDE574","MDP"="#0099CC","EBMP"="#66CCB9","CLP"="#FF9999", "preB1"="#99CCFF","preB2"="#88A6AD","preB3"="#996666","NK/Tp"="#3AC569","MEP"="#D071A9","MkP1"="#9ED47D","MkP2"="#F57F2C","EryP"="#CCCCFF")

pdf("plot_trajectory1.pdf", width = 20, height = 20)
plot1<-plot_cell_trajectory(my_cds_subset, color_by = "State")
print(plot1)
dev.off()

pdf("plot_trajectory2.pdf", width = 20, height = 20)
plot2<-plot_cell_trajectory(my_cds_subset, color_by = "Cluster")
print(plot2)
dev.off()

pdf("plot_trajectory3_1100.pdf", width = 20, height = 20)
plot3<-plot_cell_trajectory(my_cds_subset,cell_size = 2.2, color_by = "cell.types")+
       scale_colour_manual(values = cols,breaks = c("HSC","MPP","GMP","NeuP","MDP","EBMP","CLP","preB1","preB2","preB3","NK/Tp","MEP","MkP1","MkP2","EryP"))
print(plot3)
dev.off()

pdf("plot_trajectory4_1100.pdf", width = 30, height = 10)
plot4<-plot_cell_trajectory(my_cds_subset, cell_size = 2.5, color_by = "cell.types") +
       scale_colour_manual(values = cols)+
       facet_wrap(~cell.types, nrow = 3)
print(plot4)
dev.off()

pdf("plot_trajectory_state_1100_yoursample.pdf", width = 10, height = 10)
plot5<-plot_cell_trajectory(my_cds_subset, color_by = "State")
print(plot5)
dev.off()

##HC and ITP display separately
pdf("plot_trajectory_HC_versus_ITP-2.2.pdf", width = 20, height = 10)
plot6<-plot_cell_trajectory(my_cds_subset, cell_size = 0.5, color_by = "cell.types") +
       scale_colour_manual(values = cols, breaks = c("HSC","MPP","GMP","NeuP","MDP","G2M","CLP","preB1","preB2","preB3","NK/Tp","MEP","MkP1","MkP2","EryP"))+
       facet_wrap(~group, nrow = 1)
print(plot6)
dev.off()

head(pData(my_cds_subset))
my_pseudotime_de <- differentialGeneTest(my_cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 22)
#write.table(my_pseudotime_de,file='my_pseudotime_de_yoursample.csv', row.names=T, sep='\t', append=T, quote=F)
my_pseudotime_de %>% arrange(qval) %>% head()
my_pseudotime_de %>% arrange(qval) %>% head() %>% dplyr::select(gene_short_name) -> my_pseudotime_gene

pdf("plot_pseudotime_1100_yoursample.pdf", width = 15, height = 15)
plot7<-plot_cell_trajectory(my_cds_subset, cell_size = 2.5, color_by = "Pseudotime")
print(plot7)
dev.off()

A=c('HLF','NOG',"GATA2","CD38","ITGA2B","MEIS1",'FLI1','PLEK','PBX1')
my_pseudotime_gene <-A
pdf("plot_pseudotime1_1100_yoursample_gene.pdf", width = 20, height = 20)
plot8<-plot_genes_in_pseudotime(my_cds_subset[my_pseudotime_gene,])
print(plot8)
dev.off()

pdf("plot_pseudotime2_1100_yoursample.pdf", width = 15, height = 15)
plot9<-plot_genes_in_pseudotime(my_cds_subset[my_pseudotime_gene,], color_by = "cell.types")
print(plot9)
dev.off()

pdf("plot_pseudotime3.pdf", width = 20, height = 20)
plot10<-plot_genes_in_pseudotime(my_cds_subset[my_pseudotime_gene,], color_by = "Cluster")
print(plot10)
dev.off()

pdf("ITPplot_heatmap.pdf", width = 10, height = 50)
plot11<-plot_pseudotime_heatmap(my_cds_subset[my_ordering_genes,],
                                num_clusters = 4,
                                cores = 1,
                                show_rownames = F,return_heatmap=T)
print(plot11)
dev.off()
