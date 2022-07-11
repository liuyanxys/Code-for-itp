rm(list=ls())
library(monocle)
library(Seurat)
library(dplyr)
setwd("/home/monocle2/MkEry/")

SO<-readRDS("/home/monocle2/yoursample.rds")
HSCcell<-rownames(subset(SO@meta.data,SO@meta.data$cell.types=="HSC"))
MPPcell<-rownames(subset(SO@meta.data,SO@meta.data$cell.types=="MPP"))
MkP1cell<-rownames(subset(SO@meta.data,SO@meta.data$cell.types=="MkP1"))
MkP2cell<-rownames(subset(SO@meta.data,SO@meta.data$cell.types=="MkP2"))
MEPcell<-rownames(subset(SO@meta.data,SO@meta.data$cell.types=="MEP"))
EryPcell<-rownames(subset(SO@meta.data,SO@meta.data$cell.types=="EryP"))
CLPcell<-rownames(subset(SO@meta.data,SO@meta.data$cell.types=="CLP"))
preB1cell<-rownames(subset(SO@meta.data,SO@meta.data$cell.types=="preB1"))
preB2cell<-rownames(subset(SO@meta.data,SO@meta.data$cell.types=="preB2"))
preB3cell<-rownames(subset(SO@meta.data,SO@meta.data$cell.types=="preB3"))
NKTpcell<-rownames(subset(SO@meta.data,SO@meta.data$cell.types=="NK/Tp"))
preB3cell<-rownames(subset(SO@meta.data,SO@meta.data$cell.types=="preB3"))
GMPcell<-rownames(subset(SO@meta.data,SO@meta.data$cell.types=="GMP"))
NeuPcell<-rownames(subset(SO@meta.data,SO@meta.data$cell.types=="NeuP"))
EBMPcell<-rownames(subset(SO@meta.data,SO@meta.data$cell.types=="EBMP"))
MDPcell<-rownames(subset(SO@meta.data,SO@meta.data$cell.types=="MDP"))

MkErycell<-c(HSCcell,MPPcell,MkP1cell,MkP2cell,MEPcell,EryPcell)

mtx_sml <- SO@assays$RNA@counts[, MkErycell]
cell_sml <- SO@meta.data[MkErycell, ]
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
names(list_cluster) <- SO@assays[["RNA"]]@data[,LYnoMPPcell]@Dimnames[[2]]
pData(cds)$cell.types <- list_cluster
rm(SO)
rm(data)
rm(pd)
rm(fData)
rm(fd)
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
my_cds_subset <- reduceDimension(my_cds_subset,max_components = 2,norm_method = 'log',num_dim = 10,reduction_method = 'tSNE',verbose = TRUE)
gc()
my_cds_subset <- clusterCells(my_cds_subset, verbose = FALSE)
table(pData(my_cds_subset)$Cluster)
save(my_cds_subset,file='my_cds_subset_yoursample_MkEry.Rdata')
#load('my_cds_subset_yoursample_MkEry.Rdata')
head(pData(my_cds_subset))
clustering_DEG_genes <- differentialGeneTest(my_cds_subset,fullModelFormulaStr = '~cell.types',cores = 5)
dim(clustering_DEG_genes)
save(clustering_DEG_genes,file='clustering_DEG_genes_yoursample_MkEry.Rdata')
#load('clustering_DEG_genes_yoursample_MkEry.Rdata')
library(dplyr)
deg <- subset(clustering_DEG_genes, qval < 0.01)
deg <- deg[order(deg$qval,decreasing=F),]
head(deg)
#write.table(deg,file="train.monocle.DEG.csv",col.names=T,row.names=F,sep="\t",quote=F)
my_ordering_genes <-  rownames(deg) [1:500]
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

my_cds_subset <- orderCells(my_cds_subset)
gc()
my_cds_subset <- orderCells(my_cds_subset,root_state = GM_state(my_cds_subset))
gc()
save(my_cds_subset,file='my_cds_subset_yoursample_MkEry_500.Rdata')
#load('my_cds_subset_yoursample_MkEry_500.Rdata')

cols<-c("HSC"="#D9534F","MPP"="#96CEB4", "GMP"="#CBE86B", "NeuP"="#EDE574","MDP"="#0099CC","EBMP"="#66CCB9","CLP"="#FF9999","preB1"="#99CCFF","preB2"="#88A6AD","preB3"="#996666","NK/Tp"="#3AC569","MEP"="#D071A9","MkP1"="#9ED47D","MkP2"="#F57F2C","EryP"="#CCCCFF")

pdf("plot_trajectory_yoursample_MkEry_500.pdf", width = 15, height = 15)
plot<-plot_cell_trajectory(my_cds_subset,cell_size = 2.2, color_by = "cell.types")+scale_colour_manual(values = cols)
print(plot)
dev.off()

pdf("plot_trajectory_celltype_yoursample_MkEry_500.pdf", width = 30, height = 20)
plot<-plot_cell_trajectory(my_cds_subset, cell_size = 2.2, color_by = "cell.types") +
     scale_colour_manual(values = cols)+
     facet_wrap(~cell.types, nrow = 2)
print(plot)
dev.off()

pdf("plot_trajectory_state_yoursample_MkEry_500.pdf", width = 15, height = 15)
plot<-plot_cell_trajectory(my_cds_subset, cell_size = 2.2, color_by = "State")
print(plot)
dev.off()

##HC and ITP display separately
pdf("plot_trajectory_HC_versus_ITP_yoursample_MkEry_500.pdf", width = 20, height = 10)
plot<-plot_cell_trajectory(my_cds_subset, cell_size = 2.2, color_by = "cell.types") +
scale_colour_manual(values = cols)+
facet_wrap(~group, nrow = 1)
print(plot)
dev.off()

pdf("plot_trajectory_eachclusters_HC_versus_ITP_yoursample_MkEry_500.pdf", width = 20, height = 45)
plot<-plot_cell_trajectory(my_cds_subset, cell_size = 2.2, color_by = "cell.types") +
scale_colour_manual(values = cols)+
facet_grid(cell.types~group)
print(plot)
dev.off()
