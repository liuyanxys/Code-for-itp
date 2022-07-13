rm(list=ls())
options(stringsAsFactors = F)

library(Seurat)
library(reshape2)
library(stringr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(harmony)

set.seed(123456789)
dpi =300

setwd('/home/erin/Desktop/itp/output')

HC_1 <- readRDS('/home/erin/Desktop/itp/samplesQC/removegenes/reHC_1.rds')
HC_1 <- subset(HC_1, subset = DF.classification == "Singlet")

HC_2 <- readRDS('/home/erin/Desktop/itp/samplesQC/removegenes/reHC_2.rds')
HC_2 <- subset(HC_2, subset = DF.classification == "Singlet")

HC_3 <- readRDS('/home/erin/Desktop/itp/samplesQC/removegenes/reHC_3.rds')
HC_3 <- subset(HC_3, subset = DF.classification == "Singlet")

HC_4 <- readRDS('/home/erin/Desktop/itp/samplesQC/removegenes/reHC_4.rds')
HC_4 <- subset(HC_4, subset = DF.classification == "Singlet")

ITP_1 <- readRDS('/home/erin/Desktop/itp/samplesQC/removegenes/reITP_1.rds')
ITP_1 <- subset(ITP_1, subset = DF.classification == "Singlet")

ITP_2 <- readRDS('/home/erin/Desktop/itp/samplesQC/removegenes/reITP_2.rds')
ITP_2 <- subset(ITP_2, subset = DF.classification == "Singlet")

ITP_3 <- readRDS('/home/erin/Desktop/itp/samplesQC/removegenes/reITP_3.rds')
ITP_3 <- subset(ITP_3, subset = DF.classification == "Singlet")

ITP_4 <- readRDS('/home/erin/Desktop/itp/samplesQC/removegenes/reITP_4.rds')
ITP_4 <- subset(ITP_4, subset = DF.classification == "Singlet")

itp.data <- merge(HC_1, c(HC_2,HC_3,HC_4,ITP_1,ITP_2,ITP_3,ITP_4), 
                  add.cell.ids=c("HC_1","HC_2","HC_3","HC_4","ITP_1","ITP_2","ITP_3","ITP_4"))
head(itp.data@meta.data)

itp.data$pANN_0.25_0.18_445 <- NULL
itp.data$DF.classifications_0.25_0.18_445 <- NULL
itp.data$pANN_0.25_0.01_415 <- NULL
itp.data$DF.classifications_0.25_0.01_415 <- NULL
itp.data$pANN_0.25_0.09_770 <- NULL
itp.data$DF.classifications_0.25_0.09_770 <- NULL
itp.data$pANN_0.25_0.03_449 <- NULL
itp.data$DF.classifications_0.25_0.03_449 <- NULL
itp.data$pANN_0.25_0.07_275 <- NULL
itp.data$DF.classifications_0.25_0.07_275 <- NULL
itp.data$pANN_0.25_0.01_599 <- NULL
itp.data$DF.classifications_0.25_0.01_599 <- NULL
itp.data$pANN_0.25_0.22_760 <- NULL
itp.data$DF.classifications_0.25_0.22_760 <- NULL
itp.data$pANN_0.25_0.01_369 <- NULL
itp.data$DF.classifications_0.25_0.01_369 <- NULL

head(itp.data@meta.data)
dim(itp.data)
table(Idents(itp.data))

itp.data <- NormalizeData(itp.data)
itp.data <- FindVariableFeatures(itp.data, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(itp.data)
itp.data <- ScaleData(itp.data, vars.to.regress = c("nCount_RNA", "percent.mt","percent.rb","percent.hsp"), features = all.genes)
itp.data <- RunPCA(itp.data, pc.genes = itp.data@var.genes, npcs = 50, verbose = FALSE)

saveRDS(itp.data, file = "itp.data")

rm(HC_1)
rm(HC_2)
rm(HC_3)
rm(HC_4)
rm(ITP_1)
rm(ITP_2)
rm(ITP_3)
rm(ITP_4)

options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = itp.data, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = itp.data, features = "PC_1", group.by = "orig.ident", pt.size = .1)
plot_grid(p1,p2)


options(repr.plot.height = 2.5, repr.plot.width = 6)
itp.data <- itp.data %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(itp.data, 'harmony')
harmony_embeddings[1:5, 1:5]
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = itp.data, reduction = "harmony", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = itp.data, features = "harmony_1", group.by = "orig.ident", pt.size = .1)
plot_grid(p1,p2)

itp.data <- itp.data %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution =  1.5) %>% 
  identity() 
head(itp.data@meta.data)

DimPlot(itp.data, reduction = "umap", label = TRUE, pt.size = .1)

#addmeta.data$cell.types
cell.types <- vector("logical",length=ncol(itp.data))
names(cell.types) <- colnames(itp.data)

cell.types[itp.data@meta.data$seurat_clusters=="0"] <- 'preB2'
cell.types[itp.data@meta.data$seurat_clusters=="1"] <- 'MEP'
cell.types[itp.data@meta.data$seurat_clusters=="2"] <- 'MPP'
cell.types[itp.data@meta.data$seurat_clusters=="3"] <- 'HSC'
cell.types[itp.data@meta.data$seurat_clusters=="4"] <- 'GMP'
cell.types[itp.data@meta.data$seurat_clusters=="5"] <- 'HSC'
cell.types[itp.data@meta.data$seurat_clusters=="6"] <- 'NeuP'
cell.types[itp.data@meta.data$seurat_clusters=="7"] <- 'EryP'
cell.types[itp.data@meta.data$seurat_clusters=="8"] <- 'MkP1'
cell.types[itp.data@meta.data$seurat_clusters=="9"] <- 'preB2'
cell.types[itp.data@meta.data$seurat_clusters=="10"] <- 'MPP'
cell.types[itp.data@meta.data$seurat_clusters=="11"] <- 'MEP'
cell.types[itp.data@meta.data$seurat_clusters=="12"] <- 'NeuP'
cell.types[itp.data@meta.data$seurat_clusters=="13"] <- 'preB1'
cell.types[itp.data@meta.data$seurat_clusters=="14"] <- 'preB1'
cell.types[itp.data@meta.data$seurat_clusters=="15"] <- 'CLP'
cell.types[itp.data@meta.data$seurat_clusters=="16"] <- 'CLP'
cell.types[itp.data@meta.data$seurat_clusters=="17"] <- 'MDP'
cell.types[itp.data@meta.data$seurat_clusters=="18"] <- 'EBMP'
cell.types[itp.data@meta.data$seurat_clusters=="19"] <- 'preB1'
cell.types[itp.data@meta.data$seurat_clusters=="20"] <- 'NeuP'
cell.types[itp.data@meta.data$seurat_clusters=="21"] <- 'MkP2'
cell.types[itp.data@meta.data$seurat_clusters=="22"] <- 'EryP'
cell.types[itp.data@meta.data$seurat_clusters=="23"] <- 'MDP'
cell.types[itp.data@meta.data$seurat_clusters=="24"] <- 'EryP'
cell.types[itp.data@meta.data$seurat_clusters=="25"] <- 'EryP'
cell.types[itp.data@meta.data$seurat_clusters=="26"] <- 'EryP'
cell.types[itp.data@meta.data$seurat_clusters=="27"] <- 'NeuP'
cell.types[itp.data@meta.data$seurat_clusters=="28"] <- 'NK/T'
cell.types[itp.data@meta.data$seurat_clusters=="29"] <- 'preB1'
cell.types[itp.data@meta.data$seurat_clusters=="30"] <- 'preB3'


cell.types[itp.data@meta.data$seurat_clusters=="31"] <- 'Unknown1'
cell.types[itp.data@meta.data$seurat_clusters=="32"] <- 'Unknown2'

itp.data[["cell.types"]] <- cell.types
head(itp.data@meta.data)

Idents(itp.data) = 'cell.types'

n_cells <- FetchData(itp.data, 
                     vars = c("ident", "sample.id")) %>%
  dplyr::count(ident, sample.id) %>%
  tidyr::spread(ident, n)

View(n_cells)

write.table(n_cells,file = "celltypes_n_cells.xls",sep = "\t",row.names = F,quote = F)

itp.data$cell.types <- factor(x =itp.data$cell.types, levels = c("HSC","MPP","GMP","NeuP","MDP","EBMP","CLP","preB1","preB2","preB3",
                                                                 "NK/T","MEP","MkP1","MkP2","EryP","Unknown1","Unknown2"))

allcolour=c("#D9534F", "#96CEB4", "#CBE86B", "#EDE574", "#0099CC", "#66CCB9", "#FF9999", "#99CCFF", "#88A6AD", "#996666",
            "#3AC569", "#D071A9", "#9ED47D", "#F57F2C", "#CCCCFF","#888b8d","#888b8d")

itp.data$group <- factor(x =itp.data$group, levels = c("HC","ITP"))
allcolour2=c("#619CFF","#F8766D")

itp.data<-AddMetaData(itp.data,itp.data@reductions$umap@cell.embeddings,col.name = colnames(itp.data@reductions$umap@cell.embeddings))

head(itp.data@meta.data)

p1 <-  ggplot(itp.data@meta.data ,aes(x=UMAP_1,y=UMAP_2))+
  geom_point(aes(color=cell.types), alpha=1, size=0.5)+
  scale_color_manual(values = allcolour)
p1
p2 <- p1 +theme_bw()+ theme(panel.grid=element_blank())
p2
ggplot2::ggsave(file="vectorgraph0.pdf", plot = p2, device = 'pdf',width = 10, height = 9, units = "in",dpi = dpi,limitsize = FALSE)

p1 <-  ggplot(itp.data@meta.data ,aes(x=UMAP_1,y=UMAP_2))+
  geom_point(aes(color=group), alpha=0.5, size=0.2)+
  scale_color_manual(values = allcolour2)
p1
p2 <- p1 +theme_bw()+ theme(panel.grid=element_blank())
p2
ggplot2::ggsave(file="vectorgraph0'.pdf", plot = p2, device = 'pdf',width = 10, height = 9, units = "in",dpi = dpi,limitsize = FALSE)

p4 <- p2 + facet_wrap( ~ sample.id, nrow = 2)
p4

ggplot2::ggsave(file="vectorgraph0''.pdf", plot = p4, device = 'pdf',width = 20, height = 9, units = "in",dpi = dpi,limitsize = FALSE)

#remove unknown cluster
Idents(itp.data) = 'seurat_clusters'
itp.data <- subset(itp.data, idents=c(0:30))

# Create metadata dataframe
metadata <- itp.data@meta.data

# Add cell IDs to metadata
metadata$cells <- rownames(metadata)

# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

# Visualize the number UMIs/transcripts per cell
metadata %>%
  ggplot(aes(color=group, x=nUMI, fill= group)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  scale_fill_manual(values=c("#619CFF","#F8766D"))+
  scale_color_manual(values = c("#619CFF","#F8766D")) +
  ylab("log10 cell density")

ggsave("nUMI.pdf", plot = last_plot(), 
       width = 7, height = 5, units = "in")

# Visualize the distribution of genes detected per cell via histogram
metadata %>%
  ggplot(aes(color=group, x=nGene, fill= group)) +
  geom_density(alpha = 0.2) +
  scale_x_log10() +
  theme_classic() +
  scale_fill_manual(values=c("#619CFF","#F8766D"))+
  scale_color_manual(values = c("#619CFF","#F8766D")) +
  ylab("log10 cell density")

ggsave("nGene.pdf", plot = last_plot(), 
       width = 7, height = 5, units = "in")

metadata %>%
  ggplot(aes(x=nGene, y=nUMI, color=group)) +
  geom_point(alpha = 0.5,size = .4) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  scale_color_manual(values = c("#619CFF","#F8766D")) +
  facet_wrap(~group)

ggsave("correlation_1.pdf", plot = last_plot(), 
       width = 7, height = 7, units = "in")

metadata %>%
  ggplot(aes(x=nGene, y=nUMI, color=cell.types)) +
  geom_point(size = .4) +
  scale_x_log10() +
  scale_y_log10() +
  theme_classic() +
  scale_color_manual(values=c("#D9534F", "#96CEB4", "#CBE86B", "#EDE574", "#0099CC", "#66CCB9", "#FF9999", "#99CCFF", "#88A6AD", "#996666",
                              "#3AC569", "#D071A9", "#9ED47D", "#F57F2C", "#CCCCFF")) +
  facet_wrap(~group)

ggsave("correlation_2.pdf", plot = last_plot(), 
       width = 7, height = 7, units = "in")


# Visualize the distribution of Genes detected per cell via boxplot
metadata %>%
  ggplot(aes(x=cell.types, y=nGene, fill=cell.types)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  scale_fill_manual(values=c("#D9534F", "#96CEB4", "#CBE86B", "#EDE574", "#0099CC", "#66CCB9", "#FF9999", "#99CCFF", "#88A6AD", "#996666",
                             "#3AC569", "#D071A9", "#9ED47D", "#F57F2C", "#CCCCFF"))

ggsave("nGene_2.pdf", plot = last_plot(), 
       width = 8, height = 5, units = "in")

# Visualize the distribution of UMIs detected per cell via boxplot
metadata %>%
  ggplot(aes(x=cell.types, y=log10(nUMI), fill=cell.types)) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  scale_fill_manual(values=c("#D9534F", "#96CEB4", "#CBE86B", "#EDE574", "#0099CC", "#66CCB9", "#FF9999", "#99CCFF", "#88A6AD", "#996666",
                             "#3AC569", "#D071A9", "#9ED47D", "#F57F2C", "#CCCCFF"))

ggsave("nUMI_2.pdf", plot = last_plot(), 
       width = 8, height = 5, units = "in")

#cellcycle
cc.genes
itp.data_cc <- CellCycleScoring(
  object = itp.data,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes
)
VlnPlot(itp.data_cc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","G2M.Score","S.Score"), ncol = 3)
umapem<-itp.data_cc@reductions$umap@cell.embeddings
metada = itp.data_cc@meta.data
dim(umapem);dim(metada)
ccdata<-data.frame(umapem,metada)
plot<-ggplot(ccdata, aes(UMAP_1, UMAP_2,label=Phase))+geom_point(aes(colour = factor(Phase)))+
  labs(x = "", y="") 
plot=plot+
  theme_bw()+theme(panel.grid=element_blank(),legend.title=element_blank(),legend.text = element_text(color="black", size = 10, face = "bold"))
plot<-plot+guides(colour = guide_legend(override.aes = list(size=5))) +theme(plot.title = element_text(hjust = 0.5))
plot

ggplot2::ggsave(file="cellcycle.pdf",plot = plot,device = 'pdf',width = 10, height = 9, units = "in",dpi = dpi,limitsize = FALSE)

itp.data.markers <- FindAllMarkers(itp.data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(itp.data.markers,file = "markers.xls",sep = "\t",row.names = F,quote = F)
write.csv(itp.data.markers, file = "markers.csv")

top10 <- itp.data.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
write.table(top10,file = "top10_markers.xls",sep = "\t",row.names = F,quote = F)

#heat map
hm = DoHeatmap(subset(itp.data,downsample = 180), features = top10$gene,label = T,size = 6, group.bar.height = 0.02,group.by ="cell.types",
                group.colors = c("#D9534F", "#96CEB4", "#CBE86B", "#EDE574", "#0099CC", "#66CCB9", "#FF9999", "#99CCFF", "#88A6AD", "#996666",
                                 "#3AC569", "#D071A9", "#9ED47D", "#F57F2C", "#CCCCFF")) +
  scale_fill_gradientn(colors = c("#00A9FF", "#FFFFFF", "#F8766D"))
ggplot2::ggsave(file="NEW-feature.pdf",plot = hm,device = 'pdf',width = 15, height = 15, units = "in",dpi = dpi,limitsize = FALSE)

#VlnPlot
features <- c('HLF','NOG','PROM1','MPO','AZU1','CD68','SPIB','JCHAIN',
              'CLC','DNTT','VPREB3','PAX5','CD72','CD96','KLRB1',
              'GATA1','GATA2','PLEK','PPBP','PF4','AHSP','GYPA')

p <- VlnPlot(itp.data, features, stack = TRUE, sort = F, fill.by = 'ident',adjust = 1.5,flip = TRUE,
              cols = c("#CCCCFF","#F57F2C", "#9ED47D","#D071A9", "#3AC569","#996666","#88A6AD", "#99CCFF", "#FF9999",  "#66CCB9",
                       "#0099CC", "#EDE574",  "#CBE86B","#96CEB4", "#D9534F")) +
  theme(legend.position = "none") + ggtitle("identities on x-axis")
ggplot2::ggsave(file="MarkerGene.pdf",plot = p4,device = 'pdf',width = 10, height = 15, units = "in",dpi = dpi,limitsize = FALSE)

saveRDS(itp.data, file = "itp.data")
