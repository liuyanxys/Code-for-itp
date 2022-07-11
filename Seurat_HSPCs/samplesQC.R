rm(list=ls())
options(stringsAsFactors = F)

library(Seurat)
library(DoubletFinder)
library(reshape2)
library(stringr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork)

set.seed(123456789)
dpi =300

dir.create('/home/erin/Desktop/itp/samplesQC3/merge')
setwd('/home/erin/Desktop/itp/samplesQC3/merge')

HC_1.data <- Read10X(data.dir = '/home/erin/Desktop/itp/input/H01_ZTJ/filtered_feature_bc_matrix')
HC_1 <- CreateSeuratObject(counts = HC_1.data, project = "HC_1", min.cells = 3, min.features = 200)
HC_1

HC_2.data <- Read10X(data.dir = '/home/erin/Desktop/itp/input/H02_WJQ/filtered_feature_bc_matrix')
HC_2<- CreateSeuratObject(counts = HC_2.data, project = "HC_2", min.cells = 3, min.features = 200)
HC_2

HC_3.data <- Read10X(data.dir = '/home/erin/Desktop/itp/input/H03_ZSQ/filtered_feature_bc_matrix')
HC_3 <- CreateSeuratObject(counts = HC_3.data, project = "HC_3", min.cells = 3, min.features = 200)
HC_3

HC_4.data <- Read10X(data.dir = '/home/erin/Desktop/itp/input/H04_JXM/filtered_feature_bc_matrix')
HC_4<- CreateSeuratObject(counts = HC_4.data, project = "HC_4", min.cells = 3, min.features = 200)
HC_4

ITP_1.data <- Read10X(data.dir = '/home/erin/Desktop/itp/input/P01_SCC/filtered_feature_bc_matrix')
ITP_1 <- CreateSeuratObject(counts = ITP_1.data, project = "ITP_1", min.cells = 3, min.features = 200)
ITP_1

ITP_3.data <- Read10X(data.dir = '/home/erin/Desktop/itp/input/P03_ZX/filtered_feature_bc_matrix')
ITP_3 <- CreateSeuratObject(counts = ITP_3.data, project = "ITP_3", min.cells = 3, min.features = 200)
ITP_3

ITP_4.data <- Read10X(data.dir = '/home/erin/Desktop/itp/input/P04_ZJX/filtered_feature_bc_matrix')
ITP_4 <- CreateSeuratObject(counts = ITP_4.data, project = "ITP_4", min.cells = 3, min.features = 200)
ITP_4

ITP_5.data <- Read10X(data.dir =  '/home/erin/Desktop/itp/input/P05_YRK/filtered_feature_bc_matrix')
ITP_5 <- CreateSeuratObject(counts =ITP_5.data, project = "ITP_5", min.cells = 3, min.features = 200)
ITP_5

itp.data <- merge(HC_1, c(HC_2,HC_3,HC_4,ITP_1,ITP_3,ITP_4,ITP_5), 
                  add.cell.ids=c("HC_1","HC_2","HC_3","HC_4","ITP_1","ITP_3","ITP_4","ITP_5"))

dim(itp.data)
#[1] 23107 68056
table(Idents(itp.data))
#HC_1  HC_2  HC_3  HC_4 ITP_1 ITP_3 ITP_4 ITP_5 
#7729  7743 11101  8590  6331  9136 10271  7155 

head(itp.data@meta.data)

quantile(itp.data@meta.data$nCount_RNA,c(0.025,0.975))
#2.5%   97.5% 
#1061.375 63723.000  
quantile(itp.data@meta.data$nCount_RNA,c(0.050,0.950))
#5%     95% 
#2303 52938 
quantile(itp.data@meta.data$nFeature_RNA,c(0.025,0.975))
#2.5%   97.5% 
#471.000 7190.625
quantile(itp.data@meta.data$nFeature_RNA,c(0.050,0.950))
#5%     95% 
#885.75 6623.00 


itp.data[["percent.mt"]] <- PercentageFeatureSet(itp.data, pattern = "^MT-")
head(itp.data@meta.data)

VlnPlot(itp.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
quantile(itp.data@meta.data$percent.mt,c(0.025,0.975))
#2.5%     97.5% 
#2.305898 42.133762 
quantile(itp.data@meta.data$percent.mt,c(0.025,0.950))
#2.5%     95% 
#2.305898 17.123945 
quantile(itp.data@meta.data$percent.mt,c(0.025,0.925))
#2.5%     92.5% 
#2.305898 10.871715 
quantile(itp.data@meta.data$percent.mt,c(0.025,0.900))
#2.5%     90% 
#2.305898 9.242237  
mt.genes <- rownames(itp.data)[grep("^MT-",rownames(itp.data))]
mt.genes
length(mt.genes)
#[1] 13

itp.data[["percent.rb"]] <- PercentageFeatureSet(itp.data, pattern = "^RP[SL]")
head(itp.data@meta.data)
VlnPlot(itp.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4)
quantile(itp.data@meta.data$percent.rb,c(0.025,0.975))
#2.5%     97.5% 
#5.414278 48.401861     
rb.genes <- rownames(itp.data)[grep("^RP[SL]",rownames(itp.data))]
rb.genes
length(rb.genes)
#[1] 101

itp.data[["percent.hsp"]] <- PercentageFeatureSet(itp.data, pattern = "^HSP")
head(itp.data@meta.data)
VlnPlot(itp.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb","percent.hsp"), ncol = 5)
quantile(itp.data@meta.data$percent.hsp,c(0.025,0.975))
#2.5%     97.5% 
#0.1624072 2.6350887
quantile(itp.data@meta.data$percent.hsp,c(0.025,0.950))
#2.5%       95% 
#0.1624072 2.0470973  
quantile(itp.data@meta.data$percent.hsp,c(0.025,0.925))
#2.5%     92.5% 
#0.1624072 1.7508975 
quantile(itp.data@meta.data$percent.hsp,c(0.025,0.900))
#2.5%       90% 
#0.1624072 1.5390551 
hsp.genes <- rownames(itp.data)[grep("^HSP",rownames(itp.data))]
hsp.genes
length(hsp.genes)
#[1] 28

itp.data$log10GenesPerUMI <- log10(itp.data$nFeature_RNA) / log10(itp.data$nCount_RNA)
head(itp.data@meta.data)
p <- VlnPlot(itp.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb","percent.hsp","log10GenesPerUMI"), ncol = 3)
png(file="qc1.png", width = dpi*12, height = dpi*10, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
quantile(itp.data@meta.data$log10GenesPerUMI,c(0.025,0.975))
#2.5%     97.5% 
#0.7975623 0.9091164 

itp.data <- NormalizeData(itp.data)
itp.data <- ScaleData(itp.data)
itp.data <- FindVariableFeatures(itp.data, selection.method = "vst", nfeatures = 2000)

variablefeatures <- VariableFeatures(itp.data)
mt.genes <- grep("^MT-",variablefeatures,value=T)
mt.genes
#[1] "MT-ND6"  "MT-ND5"  "MT-CO1"  "MT-CO3"  "MT-ND1"  "MT-ND4"  "MT-CO2"  "MT-ND2"  "MT-ATP6" "MT-CYB" 
length(mt.genes)
#[1] 10
rb.genes <- grep("^RP[SL]",variablefeatures,value=T)
rb.genes
#[1] "RPS4Y2"  "RPS4Y1"  "RPL22L1"
length(rb.genes)
#[1] 3
hsp.genes <- grep("^HSP",variablefeatures,value=T)
hsp.genes 
#[1] "HSPA6"    "HSPA1A"   "HSPA1B"   "HSPH1"    "HSPA5"    "HSP90AA1" "HSPA8"    "HSPE1"    "HSPD1"    "HSPA9"    "HSPB1"
length(hsp.genes)
#[1] 11

head(itp.data@meta.data)
p1 <- FeatureScatter(itp.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
png(file="qc2.png", width = dpi*9, height = dpi*9, units = "px",res = dpi,type='cairo')
print(p1)
dev.off()
p2 <- FeatureScatter(itp.data, feature1="nCount_RNA", feature2="percent.mt")
png(file="qc3.png", width = dpi*9, height = dpi*9, units = "px",res = dpi,type='cairo')
print(p2)
dev.off()
p3 <- FeatureScatter(itp.data, feature1="nFeature_RNA", feature2="percent.mt")
png(file="qc4.png", width = dpi*9, height = dpi*9, units = "px",res = dpi,type='cairo')
print(p3)
dev.off()
p4 <- FeatureScatter(itp.data, feature1="nFeature_RNA", feature2="percent.rb")
png(file="qc5.png", width = dpi*9, height = dpi*9, units = "px",res = dpi,type='cairo')
print(p4)
dev.off()
p5 <- FeatureScatter(itp.data, feature1="nFeature_RNA", feature2="percent.hsp")
png(file="qc6.png", width = dpi*9, height = dpi*9, units = "px",res = dpi,type='cairo')
print(p5)
dev.off()

p6 <- FeatureScatter(itp.data, feature1="percent.mt", feature2="percent.rb")
png(file="qc7.png", width = dpi*9, height = dpi*9, units = "px",res = dpi,type='cairo')
print(p6)
dev.off()
p7 <- FeatureScatter(itp.data, feature1="percent.mt", feature2="percent.hsp")
png(file="qc8.png", width = dpi*9, height = dpi*9, units = "px",res = dpi,type='cairo')
print(p7)
dev.off()
p8 <- FeatureScatter(itp.data, feature1="percent.rb", feature2="percent.hsp")
png(file="qc9.png", width = dpi*9, height = dpi*9, units = "px",res = dpi,type='cairo')
print(p8)
dev.off()
##########################################################################################

rm(list=ls())
options(stringsAsFactors = F)

library(Seurat)
library(DoubletFinder)
library(reshape2)
library(stringr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork) 

set.seed(123456789)
dpi =300

#######
dir.create('/home/erin/Desktop/itp/samplesQC3/HC_1')
setwd('/home/erin/Desktop/itp/samplesQC3/HC_1')
HC_1.data <- Read10X(data.dir = '/home/erin/Desktop/itp/input/H01_ZTJ/filtered_feature_bc_matrix')
HC_1 <- CreateSeuratObject(counts = HC_1.data, project = "HC_1", min.cells = 3, min.features = 200)
dim(HC_1)
#[1] 20392  7729

HC_1$group<-"HC"
HC_1$gender<-"M"

HC_1[["percent.mt"]] <- PercentageFeatureSet(HC_1, pattern = "^MT-")
HC_1[["percent.rb"]] <- PercentageFeatureSet(HC_1, pattern = "^RP[SL]")
HC_1[["percent.hsp"]] <- PercentageFeatureSet(HC_1, pattern = "^HSP")

head(HC_1@meta.data)

HC_1 <- subset(HC_1, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 15 & percent.rb < 48.4 & percent.hsp < 2.6)
dim(HC_1)
#[1] 20392  7416
HC_1 <- NormalizeData(HC_1)
HC_1 <- FindVariableFeatures(HC_1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(HC_1)
HC_1 <- ScaleData(HC_1, features = all.genes)

HC_1 <- RunPCA(HC_1, npcs = 50)

HC_1 <- JackStraw(HC_1, dims = 50, num.replicate = 100)
HC_1 <- ScoreJackStraw(HC_1, dims = 1:50)
p1 <- JackStrawPlot(HC_1, dims = 1:50)
png(file="pc_pvalue.png", width = dpi*12, height = dpi*6, units = "px",res = dpi,type='cairo')
print(p1)
dev.off()
p2 <- ElbowPlot(HC_1,ndims =50)
png(file="pc_elbowplot.png", width = dpi*9, height = dpi*6, units = "px",res = dpi,type='cairo')
print(p2)
dev.off()

HC_1  <- RunUMAP(HC_1 , dims = 1:30)
p <- DimPlot(HC_1)
png(file="dimplot1.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(HC_1,
                 features = "CD34")
png(file="dimplot2.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(HC_1,
                 features = "HLF")#hsc
png(file="dimplot3.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(HC_1,
                 features = "NOG")#hsc
png(file="dimplot4.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(HC_1,
                 features = "EBF1")#lymphoid
png(file="dimplot5.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(HC_1,
                 features = "CEBPD")#myeloid
png(file="dimplot6.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(HC_1,
                 features = "GATA1")#erythroid")
png(file="dimplot7.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()

sweep.res.list_HC_1 <- paramSweep_v3(HC_1, PCs = 1:30, sct = FALSE)

sweep.stats_HC_1 <- summarizeSweep(sweep.res.list_HC_1, GT = FALSE)

bcmvn_HC_1 <- find.pK(sweep.stats_HC_1)

mpK<-as.numeric(as.vector(bcmvn_HC_1$pK[which.max(bcmvn_HC_1$BCmetric)]))
mpK

annotations <- HC_1@meta.data$seurat_clusters

homotypic.prop <- modelHomotypic(annotations)   
nExp_poi <- round(0.060*length(HC_1@meta.data$orig.ident))  #7813
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

HC_1 <- doubletFinder_v3(HC_1, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
head(HC_1@meta.data)
table(HC_1@meta.data$ DF.classifications_0.25_0.18_445)
#Doublet Singlet 
#445    6971 
p <- DimPlot(HC_1,group.by = 'DF.classifications_0.25_0.18_445')
png(file="dimplot8.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- DimPlot(HC_1,split.by = 'DF.classifications_0.25_0.18_445')
png(file="dimplot9.png", width = dpi*9, height = dpi*5, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
HC_1@meta.data$DF.classification <- HC_1@meta.data$DF.classifications_0.25_0.18_445
head(HC_1@meta.data)

saveRDS(HC_1, file = 'HC_1.rds')

HC_1<-readRDS('/home/erin/Desktop/itp/samplesQC3/HC_1/HC_1.rds')

dir.create('/home/erin/Desktop/itp/samplesQC3/removegenes')
setwd('/home/erin/Desktop/itp/samplesQC3/removegenes')
mt.genes <- rownames(HC_1)[grep("^MT-",rownames(HC_1))]
rb.genes <- rownames(HC_1)[grep("^RP[SL]",rownames(HC_1))]
hsp.genes <- rownames(HC_1)[grep("^HSP",rownames(HC_1))]
length(c(mt.genes,rb.genes,hsp.genes))
#[1] 140
dim(HC_1)
#[1] 20392  7416
counts <- GetAssayData(HC_1, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c(mt.genes,rb.genes,hsp.genes))),]
HC_1 <- subset(HC_1, features = rownames(counts))
dim(HC_1)
#[1] 20252  7416
saveRDS(HC_1, file = 'reHC_1.rds')

##########
dir.create('/home/erin/Desktop/itp/samplesQC3/HC_2')
setwd('/home/erin/Desktop/itp/samplesQC3/HC_2')
HC_2.data <- Read10X(data.dir = '/home/erin/Desktop/itp/input/H02_WJQ/filtered_feature_bc_matrix')
HC_2<- CreateSeuratObject(counts = HC_2.data, project = "HC_2", min.cells = 3, min.features = 200)
dim(HC_2)
#[1] 20304  7743

HC_2$group<-"HC"
HC_2$gender<-"M"

HC_2[["percent.mt"]] <- PercentageFeatureSet(HC_2, pattern = "^MT-")
HC_2[["percent.rb"]] <- PercentageFeatureSet(HC_2, pattern = "^RP[SL]")
HC_2[["percent.hsp"]] <- PercentageFeatureSet(HC_2, pattern = "^HSP")

head(HC_2@meta.data)

HC_2 <- subset(HC_2, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 15 & percent.rb < 48.4 & percent.hsp < 2.6)
dim(HC_2)
#[1] 20304  6799

HC_2 <- NormalizeData(HC_2)
HC_2 <- FindVariableFeatures(HC_2, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(HC_2)
HC_2 <- ScaleData(HC_2, features = all.genes)

HC_2 <- RunPCA(HC_2, npcs = 50)

HC_2 <- JackStraw(HC_2, dims = 50, num.replicate = 100)
HC_2 <- ScoreJackStraw(HC_2, dims = 1:50)
p1 <- JackStrawPlot(HC_2, dims = 1:50)
png(file="pc_pvalue.png", width = dpi*12, height = dpi*6, units = "px",res = dpi,type='cairo')
print(p1)
dev.off()
p2 <- ElbowPlot(HC_2,ndims =50)
png(file="pc_elbowplot.png", width = dpi*9, height = dpi*6, units = "px",res = dpi,type='cairo')
print(p2)
dev.off()

HC_2  <- RunUMAP(HC_2 , dims = 1:30)
p <- DimPlot(HC_2)
png(file="dimplot1.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(HC_2,
                 features = "CD34")
png(file="dimplot2.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(HC_2,
                 features = "HLF")#hsc
png(file="dimplot3.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(HC_2,
                 features = "NOG")#hsc
png(file="dimplot4.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(HC_2,
                 features = "EBF1")#lymphoid
png(file="dimplot5.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(HC_2,
                 features = "CEBPD")#myeloid
png(file="dimplot6.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(HC_2,
                 features = "GATA1")#erythroid")
png(file="dimplot7.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()

sweep.res.list_HC_2 <- paramSweep_v3(HC_2, PCs = 1:30, sct = FALSE)

sweep.stats_HC_2 <- summarizeSweep(sweep.res.list_HC_2, GT = FALSE)

bcmvn_HC_2 <- find.pK(sweep.stats_HC_2)

mpK<-as.numeric(as.vector(bcmvn_HC_2$pK[which.max(bcmvn_HC_2$BCmetric)]))
mpK

annotations <- HC_2@meta.data$seurat_clusters

homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.061*length(HC_2@meta.data$orig.ident))  #8033
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

HC_2 <- doubletFinder_v3(HC_2, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
head(HC_2@meta.data)
table(HC_2@meta.data$DF.classifications_0.25_0.01_415)
#Doublet Singlet 
#415    6384 
p <- DimPlot(HC_2,group.by = 'DF.classifications_0.25_0.01_415')
png(file="dimplot8.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- DimPlot(HC_2,split.by = 'DF.classifications_0.25_0.01_415')
png(file="dimplot9.png", width = dpi*9, height = dpi*5, units = "px",res = dpi,type='cairo')
print(p)
dev.off()

HC_2@meta.data$DF.classification <- HC_2@meta.data$DF.classifications_0.25_0.01_415
head(HC_2@meta.data)

saveRDS(HC_2, file = 'HC_2.rds')

setwd('/home/erin/Desktop/itp/samplesQC3/removegenes')
mt.genes <- rownames(HC_2)[grep("^MT-",rownames(HC_2))]
rb.genes <- rownames(HC_2)[grep("^RP[SL]",rownames(HC_2))]
hsp.genes <- rownames(HC_2)[grep("^HSP",rownames(HC_2))]
length(c(mt.genes,rb.genes,hsp.genes))
#[1] 139
dim(HC_2)
#[1] 20304  6799
counts <- GetAssayData(HC_2, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c(mt.genes,rb.genes,hsp.genes))),]
HC_2 <- subset(HC_2, features = rownames(counts))
dim(HC_2)
#[1] 20165  6799
saveRDS(HC_2, file = 'reHC_2.rds')

#############
dir.create('/home/erin/Desktop/itp/samplesQC3/HC_3')
setwd('/home/erin/Desktop/itp/samplesQC3/HC_3')
HC_3.data <- Read10X(data.dir = '/home/erin/Desktop/itp/input/H03_ZSQ/filtered_feature_bc_matrix')
HC_3 <- CreateSeuratObject(counts = HC_3.data, project = "HC_3", min.cells = 3, min.features = 200)
dim(HC_3)
#[1] 20878 11101

HC_3$group<-"HC"
HC_3$gender<-"F"

HC_3[["percent.mt"]] <- PercentageFeatureSet(HC_3, pattern = "^MT-")
HC_3[["percent.rb"]] <- PercentageFeatureSet(HC_3, pattern = "^RP[SL]")
HC_3[["percent.hsp"]] <- PercentageFeatureSet(HC_3, pattern = "^HSP")

head(HC_3@meta.data)

HC_3 <- subset(HC_3, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 15 & percent.rb < 48.4 & percent.hsp < 2.6)
dim(HC_3)
#[1] 20878  9059
HC_3 <- NormalizeData(HC_3)
HC_3 <- FindVariableFeatures(HC_3, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(HC_3)
HC_3 <- ScaleData(HC_3, features = all.genes)

HC_3 <- RunPCA(HC_3, npcs = 50)

HC_3 <- JackStraw(HC_3, dims = 50, num.replicate = 100)
HC_3 <- ScoreJackStraw(HC_3, dims = 1:50)
p1 <- JackStrawPlot(HC_3, dims = 1:50)
png(file="pc_pvalue.png", width = dpi*12, height = dpi*6, units = "px",res = dpi,type='cairo')
print(p1)
dev.off()
p2 <- ElbowPlot(HC_3,ndims =50)
png(file="pc_elbowplot.png", width = dpi*9, height = dpi*6, units = "px",res = dpi,type='cairo')
print(p2)
dev.off()

HC_3  <- RunUMAP(HC_3 , dims = 1:30)
p <- DimPlot(HC_3)
png(file="dimplot1.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(HC_3,
                 features = "CD34")
png(file="dimplot2.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(HC_3,
                 features = "HLF")#hsc
png(file="dimplot3.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(HC_3,
                 features = "NOG")#hsc
png(file="dimplot4.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(HC_3,
                 features = "EBF1")#lymphoid
png(file="dimplot5.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(HC_3,
                 features = "CEBPD")#myeloid
png(file="dimplot6.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(HC_3,
                 features = "GATA1")#erythroid")
png(file="dimplot7.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()

sweep.res.list_HC_3 <- paramSweep_v3(HC_3, PCs = 1:30, sct = FALSE)

sweep.stats_HC_3 <- summarizeSweep(sweep.res.list_HC_3, GT = FALSE)

bcmvn_HC_3 <- find.pK(sweep.stats_HC_3)

mpK<-as.numeric(as.vector(bcmvn_HC_3$pK[which.max(bcmvn_HC_3$BCmetric)]))
mpK

annotations <- HC_3@meta.data$seurat_clusters

homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.085*length(HC_3@meta.data$orig.ident))  #11346
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

HC_3 <- doubletFinder_v3(HC_3, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
head(HC_3@meta.data)
table(HC_3@meta.data$DF.classifications_0.25_0.09_770)
#Doublet Singlet 
#770    8289 
p <- DimPlot(HC_3,group.by = 'DF.classifications_0.25_0.09_770')
png(file="dimplot8.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- DimPlot(HC_3,split.by = 'DF.classifications_0.25_0.09_770')
png(file="dimplot9.png", width = dpi*9, height = dpi*5, units = "px",res = dpi,type='cairo')
print(p)
dev.off()

HC_3@meta.data$DF.classification <- HC_3@meta.data$DF.classifications_0.25_0.09_770
head(HC_3@meta.data)

saveRDS(HC_3, file = 'HC_3.rds')

setwd('/home/erin/Desktop/itp/samplesQC3/removegenes')
mt.genes <- rownames(HC_3)[grep("^MT-",rownames(HC_3))]
rb.genes <- rownames(HC_3)[grep("^RP[SL]",rownames(HC_3))]
hsp.genes <- rownames(HC_3)[grep("^HSP",rownames(HC_3))]
length(c(mt.genes,rb.genes,hsp.genes))
#[1] 139
dim(HC_3)
#[1] 20878  9059
counts <- GetAssayData(HC_3, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c(mt.genes,rb.genes,hsp.genes))),]
HC_3 <- subset(HC_3, features = rownames(counts))
dim(HC_3)
#[1] 20739  9059
saveRDS(HC_3, file = 'reHC_3.rds')

##########
dir.create('/home/erin/Desktop/itp/samplesQC3/HC_4')
setwd('/home/erin/Desktop/itp/samplesQC3/HC_4')
HC_4.data <- Read10X(data.dir = '/home/erin/Desktop/itp/input/H04_JXM/filtered_feature_bc_matrix')
HC_4<- CreateSeuratObject(counts = HC_4.data, project = "HC_4", min.cells = 3, min.features = 200)
dim(HC_4)
#[1] 20265  8590

HC_4$group<-"HC"
HC_4$gender<-"F"

HC_4[["percent.mt"]] <- PercentageFeatureSet(HC_4, pattern = "^MT-")
HC_4[["percent.rb"]] <- PercentageFeatureSet(HC_4, pattern = "^RP[SL]")
HC_4[["percent.hsp"]] <- PercentageFeatureSet(HC_4, pattern = "^HSP")

head(HC_4@meta.data)

HC_4 <- subset(HC_4, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 15 & percent.rb < 48.4 & percent.hsp < 2.6)
dim(HC_4)
#[1] 20265  6610
HC_4 <- NormalizeData(HC_4)
HC_4 <- FindVariableFeatures(HC_4, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(HC_4)
HC_4 <- ScaleData(HC_4, features = all.genes)

HC_4 <- RunPCA(HC_4, npcs = 50)

HC_4 <- JackStraw(HC_4, dims = 50, num.replicate = 100)
HC_4 <- ScoreJackStraw(HC_4, dims = 1:50)
p1 <- JackStrawPlot(HC_4, dims = 1:50)
png(file="pc_pvalue.png", width = dpi*12, height = dpi*6, units = "px",res = dpi,type='cairo')
print(p1)
dev.off()
p2 <- ElbowPlot(HC_4,ndims =50)
png(file="pc_elbowplot.png", width = dpi*9, height = dpi*6, units = "px",res = dpi,type='cairo')
print(p2)
dev.off()

HC_4  <- RunUMAP(HC_4 , dims = 1:30)
p <- DimPlot(HC_4)
png(file="dimplot1.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(HC_4,
                 features = "CD34")
png(file="dimplot2.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(HC_4,
                 features = "HLF")#hsc
png(file="dimplot3.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(HC_4,
                 features = "NOG")#hsc
png(file="dimplot4.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(HC_4,
                 features = "EBF1")#lymphoid
png(file="dimplot5.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(HC_4,
                 features = "CEBPD")#myeloid
png(file="dimplot6.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(HC_4,
                 features = "GATA1")#erythroid")
png(file="dimplot7.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()

sweep.res.list_HC_4 <- paramSweep_v3(HC_4, PCs = 1:30, sct = FALSE)

sweep.stats_HC_4 <- summarizeSweep(sweep.res.list_HC_4, GT = FALSE)

bcmvn_HC_4 <- find.pK(sweep.stats_HC_4)

mpK<-as.numeric(as.vector(bcmvn_HC_4$pK[which.max(bcmvn_HC_4$BCmetric)]))
mpK

annotations <- HC_4@meta.data$seurat_clusters

homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.068*length(HC_4@meta.data$orig.ident))  #8737
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

HC_4 <- doubletFinder_v3(HC_4, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
head(HC_4@meta.data)
table(HC_4@meta.data$DF.classifications_0.25_0.03_449)
#Doublet Singlet 
#449    6161 
p <- DimPlot(HC_4,group.by = 'DF.classifications_0.25_0.03_449')
png(file="dimplot8.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- DimPlot(HC_4,split.by = 'DF.classifications_0.25_0.03_449')
png(file="dimplot9.png", width = dpi*9, height = dpi*5, units = "px",res = dpi,type='cairo')
print(p)
dev.off()

HC_4@meta.data$DF.classification <- HC_4@meta.data$DF.classifications_0.25_0.03_449
head(HC_4@meta.data)

saveRDS(HC_4, file = 'HC_4.rds')

setwd('/home/erin/Desktop/itp/samplesQC3/removegenes')
mt.genes <- rownames(HC_4)[grep("^MT-",rownames(HC_4))]
rb.genes <- rownames(HC_4)[grep("^RP[SL]",rownames(HC_4))]
hsp.genes <- rownames(HC_4)[grep("^HSP",rownames(HC_4))]
length(c(mt.genes,rb.genes,hsp.genes))
#[1] 139
dim(HC_4)
#[1] 20265  6610
counts <- GetAssayData(HC_4, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c(mt.genes,rb.genes,hsp.genes))),]
HC_4 <- subset(HC_4, features = rownames(counts))
dim(HC_4)
#[1] 20126  6610
saveRDS(HC_4, file = 'reHC_4.rds')

############
dir.create('/home/erin/Desktop/itp/samplesQC3/ITP_1')
setwd('/home/erin/Desktop/itp/samplesQC3/ITP_1')
ITP_1.data <- Read10X(data.dir = '/home/erin/Desktop/itp/input/P01_SCC/filtered_feature_bc_matrix')
ITP_1 <- CreateSeuratObject(counts = ITP_1.data, project = "ITP_1", min.cells = 3, min.features = 200)
dim(ITP_1)
#[1] 20900  6331

ITP_1$group<-"ITP"
ITP_1$gender<-"F"

ITP_1[["percent.mt"]] <- PercentageFeatureSet(ITP_1, pattern = "^MT-")
ITP_1[["percent.rb"]] <- PercentageFeatureSet(ITP_1, pattern = "^RP[SL]")
ITP_1[["percent.hsp"]] <- PercentageFeatureSet(ITP_1, pattern = "^HSP")

head(ITP_1@meta.data)
ITP_1 <- subset(ITP_1, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 15 & percent.rb < 48.4 & percent.hsp < 2.6)
dim(ITP_1)
#[1] 20900  5726
ITP_1 <- NormalizeData(ITP_1)
ITP_1 <- FindVariableFeatures(ITP_1, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(ITP_1)
ITP_1 <- ScaleData(ITP_1, features = all.genes)

ITP_1 <- RunPCA(ITP_1, npcs = 50)

ITP_1 <- JackStraw(ITP_1, dims = 50, num.replicate = 100)
ITP_1 <- ScoreJackStraw(ITP_1, dims = 1:50)
p1 <- JackStrawPlot(ITP_1, dims = 1:50)
png(file="pc_pvalue.png", width = dpi*12, height = dpi*6, units = "px",res = dpi,type='cairo')
print(p1)
dev.off()
p2 <- ElbowPlot(ITP_1,ndims =50)
png(file="pc_elbowplot.png", width = dpi*9, height = dpi*6, units = "px",res = dpi,type='cairo')
print(p2)
dev.off()

ITP_1  <- RunUMAP(ITP_1 , dims = 1:30)
p <- DimPlot(ITP_1)
png(file="dimplot1.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(ITP_1,
                 features = "CD34")
png(file="dimplot2.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(ITP_1,
                 features = "HLF")#hsc
png(file="dimplot3.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(ITP_1,
                 features = "NOG")#hsc
png(file="dimplot4.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(ITP_1,
                 features = "EBF1")#lymphoid
png(file="dimplot5.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(ITP_1,
                 features = "CEBPD")#myeloid
png(file="dimplot6.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(ITP_1,
                 features = "GATA1")#erythroid")
png(file="dimplot7.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()

sweep.res.list_ITP_1 <- paramSweep_v3(ITP_1, PCs = 1:30, sct = FALSE)

sweep.stats_ITP_1 <- summarizeSweep(sweep.res.list_ITP_1, GT = FALSE)

bcmvn_ITP_1 <- find.pK(sweep.stats_ITP_1)

mpK<-as.numeric(as.vector(bcmvn_ITP_1$pK[which.max(bcmvn_ITP_1$BCmetric)]))
mpK

annotations <- ITP_1@meta.data$seurat_clusters

homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.048*length(ITP_1@meta.data$orig.ident))  #6388
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

ITP_1 <- doubletFinder_v3(ITP_1, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
head(ITP_1@meta.data)
table(ITP_1@meta.data$DF.classifications_0.25_0.07_275)
#Doublet Singlet 
#275    5451
p <- DimPlot(ITP_1,group.by = 'DF.classifications_0.25_0.07_275')
png(file="dimplot8.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- DimPlot(ITP_1,split.by = 'DF.classifications_0.25_0.07_275')
png(file="dimplot9.png", width = dpi*9, height = dpi*5, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
ITP_1@meta.data$DF.classification <- ITP_1@meta.data$DF.classifications_0.25_0.07_275
head(ITP_1@meta.data)

saveRDS(ITP_1, file = 'ITP_1.rds')


setwd('/home/erin/Desktop/itp/samplesQC3/removegenes')
mt.genes <- rownames(ITP_1)[grep("^MT-",rownames(ITP_1))]
rb.genes <- rownames(ITP_1)[grep("^RP[SL]",rownames(ITP_1))]
hsp.genes <- rownames(ITP_1)[grep("^HSP",rownames(ITP_1))]
length(c(mt.genes,rb.genes,hsp.genes))
#[1] 139
dim(ITP_1)
#[1] 20900  5726
counts <- GetAssayData(ITP_1, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c(mt.genes,rb.genes,hsp.genes))),]
ITP_1 <- subset(ITP_1, features = rownames(counts))
dim(ITP_1)
#[1] 20761  5726
saveRDS(ITP_1, file = 'reITP_1.rds')

########
dir.create('/home/erin/Desktop/itp/samplesQC3/ITP_3')
setwd('/home/erin/Desktop/itp/samplesQC3/ITP_3')
ITP_3.data <- Read10X(data.dir = '/home/erin/Desktop/itp/input/P03_ZX/filtered_feature_bc_matrix')
ITP_3 <- CreateSeuratObject(counts = ITP_3.data, project = "ITP_3", min.cells = 3, min.features = 200)
dim(ITP_3)
#[1] 20343  9136

ITP_3$group<-"ITP"
ITP_3$gender<-"M"

ITP_3[["percent.mt"]] <- PercentageFeatureSet(ITP_3, pattern = "^MT-")
ITP_3[["percent.rb"]] <- PercentageFeatureSet(ITP_3, pattern = "^RP[SL]")
ITP_3[["percent.hsp"]] <- PercentageFeatureSet(ITP_3, pattern = "^HSP")

head(ITP_3@meta.data)

ITP_3 <- subset(ITP_3, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 15 & percent.rb < 48.4 & percent.hsp < 2.6)
dim(ITP_3)
#[1] 20343  8443
ITP_3 <- NormalizeData(ITP_3)
ITP_3 <- FindVariableFeatures(ITP_3, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(ITP_3)
ITP_3 <- ScaleData(ITP_3, features = all.genes)

ITP_3 <- RunPCA(ITP_3 , npcs = 50)

ITP_3 <- JackStraw(ITP_3, dims = 50, num.replicate = 100)
ITP_3 <- ScoreJackStraw(ITP_3, dims = 1:50)
p1 <- JackStrawPlot(ITP_3, dims = 1:50)
png(file="pc_pvalue.png", width = dpi*12, height = dpi*6, units = "px",res = dpi,type='cairo')
print(p1)
dev.off()
p2 <- ElbowPlot(ITP_3,ndims =50)
png(file="pc_elbowplot.png", width = dpi*9, height = dpi*6, units = "px",res = dpi,type='cairo')
print(p2)
dev.off()

ITP_3 <- RunUMAP(ITP_3, dims = 1:30)
p <- DimPlot(ITP_3)
png(file="dimplot1.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(ITP_3,
                 features = "CD34")
png(file="dimplot2.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(ITP_3,
                 features = "HLF")#hsc
png(file="dimplot3.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(ITP_3,
                 features = "NOG")#hsc
png(file="dimplot4.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(ITP_3,
                 features = "EBF1")#lymphoid
png(file="dimplot5.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(ITP_3,
                 features = "CEBPD")#myeloid
png(file="dimplot6.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(ITP_3,
                 features = "GATA1")#erythroid")
png(file="dimplot7.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()

sweep.res.list_ITP_3 <- paramSweep_v3(ITP_3, PCs = 1:30, sct = FALSE)

sweep.stats_ITP_3 <- summarizeSweep(sweep.res.list_ITP_3, GT = FALSE)

bcmvn_ITP_3 <- find.pK(sweep.stats_ITP_3)

mpK<-as.numeric(as.vector(bcmvn_ITP_3$pK[which.max(bcmvn_ITP_3$BCmetric)]))
mpK

annotations <- ITP_3@meta.data$seurat_clusters

homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.071*length(ITP_3@meta.data$orig.ident))  #9245
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

ITP_3 <- doubletFinder_v3(ITP_3, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
head(ITP_3@meta.data)
table(ITP_3@meta.data$DF.classifications_0.25_0.01_599)
#Doublet Singlet 
#599    7844  
p <- DimPlot(ITP_3,group.by = 'DF.classifications_0.25_0.01_599')
png(file="dimplot8.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- DimPlot(ITP_3,split.by = 'DF.classifications_0.25_0.01_599')
png(file="dimplot9.png", width = dpi*9, height = dpi*5, units = "px",res = dpi,type='cairo')
print(p)
dev.off()

ITP_3@meta.data$DF.classification <- ITP_3@meta.data$DF.classifications_0.25_0.01_599
head(ITP_3@meta.data)

saveRDS(ITP_3, file = 'ITP_3.rds')

setwd('/home/erin/Desktop/itp/samplesQC3/removegenes')
mt.genes <- rownames(ITP_3)[grep("^MT-",rownames(ITP_3))]
rb.genes <- rownames(ITP_3)[grep("^RP[SL]",rownames(ITP_3))]
hsp.genes <- rownames(ITP_3)[grep("^HSP",rownames(ITP_3))]
length(c(mt.genes,rb.genes,hsp.genes))
#[1] 141
dim(ITP_3)
#[1] 20343  8443
counts <- GetAssayData(ITP_3, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c(mt.genes,rb.genes,hsp.genes))),]
ITP_3 <- subset(ITP_3, features = rownames(counts))
dim(ITP_3)
#[1] 20202  8443
saveRDS(ITP_3, file = 'reITP_3.rds')

#######
dir.create('/home/erin/Desktop/itp/samplesQC3/ITP_4')
setwd('/home/erin/Desktop/itp/samplesQC3/ITP_4')
ITP_4.data <- Read10X(data.dir = '/home/erin/Desktop/itp/input/P04_ZJX/filtered_feature_bc_matrix')
ITP_4 <- CreateSeuratObject(counts = ITP_4.data, project = "ITP_4", min.cells = 3, min.features = 200)
dim(ITP_4)
#[1] 20980 10271

ITP_4$group<-"ITP"
ITP_4$gender<-"F"

ITP_4[["percent.mt"]] <- PercentageFeatureSet(ITP_4, pattern = "^MT-")
ITP_4[["percent.rb"]] <- PercentageFeatureSet(ITP_4, pattern = "^RP[SL]")
ITP_4[["percent.hsp"]] <- PercentageFeatureSet(ITP_4, pattern = "^HSP")

head(ITP_4@meta.data)

ITP_4 <- subset(ITP_4, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 15 & percent.rb < 48.4 & percent.hsp < 2.6)
dim(ITP_4)
#[1] 20980  9626
ITP_4 <- NormalizeData(ITP_4)
ITP_4 <- FindVariableFeatures(ITP_4, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(ITP_4)
ITP_4 <- ScaleData(ITP_4, features = all.genes)

ITP_4 <- RunPCA(ITP_4, npcs = 50)

ITP_4 <- JackStraw(ITP_4, dims = 50, num.replicate = 100)
ITP_4 <- ScoreJackStraw(ITP_4, dims = 1:50)
p1 <- JackStrawPlot(ITP_4, dims = 1:50)
png(file="pc_pvalue.png", width = dpi*12, height = dpi*6, units = "px",res = dpi,type='cairo')
print(p1)
dev.off()
p2 <- ElbowPlot(ITP_4,ndims =50)
png(file="pc_elbowplot.png", width = dpi*9, height = dpi*6, units = "px",res = dpi,type='cairo')
print(p2)
dev.off()

ITP_4 <- RunUMAP(ITP_4 , dims = 1:30)
p <- DimPlot(ITP_4)
png(file="dimplot1.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(ITP_4,
                 features = "CD34")
png(file="dimplot2.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(ITP_4,
                 features = "HLF")#hsc
png(file="dimplot3.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(ITP_4,
                 features = "NOG")#hsc
png(file="dimplot4.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(ITP_4,
                 features = "EBF1")#lymphoid
png(file="dimplot5.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(ITP_4,
                 features = "CEBPD")#myeloid
png(file="dimplot6.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(ITP_4,
                 features = "GATA1")#erythroid")
png(file="dimplot7.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()

sweep.res.list_ITP_4 <- paramSweep_v3(ITP_4, PCs = 1:30, sct = FALSE)

sweep.stats_ITP_4 <- summarizeSweep(sweep.res.list_ITP_4, GT = FALSE)

bcmvn_ITP_4 <- find.pK(sweep.stats_ITP_4)

mpK<-as.numeric(as.vector(bcmvn_ITP_4$pK[which.max(bcmvn_ITP_4$BCmetric)]))
mpK

annotations <- ITP_4@meta.data$seurat_clusters

homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.079*length(ITP_4@meta.data$orig.ident))  #10437
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

ITP_4 <- doubletFinder_v3(ITP_4, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
head(ITP_4@meta.data)
table(ITP_4@meta.data$DF.classifications_0.25_0.22_760)
#Doublet Singlet 
#760    8866 
p <- DimPlot(ITP_4,group.by = 'DF.classifications_0.25_0.22_760')
png(file="dimplot8.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- DimPlot(ITP_4,split.by = 'DF.classifications_0.25_0.22_760')
png(file="dimplot9.png", width = dpi*9, height = dpi*5, units = "px",res = dpi,type='cairo')
print(p)
dev.off()

ITP_4@meta.data$DF.classification <- ITP_4@meta.data$DF.classifications_0.25_0.22_760
head(ITP_4@meta.data)

saveRDS(ITP_4, file = 'ITP_4.rds')

setwd('/home/erin/Desktop/itp/samplesQC3/removegenes')
mt.genes <- rownames(ITP_4)[grep("^MT-",rownames(ITP_4))]
rb.genes <- rownames(ITP_4)[grep("^RP[SL]",rownames(ITP_4))]
hsp.genes <- rownames(ITP_4)[grep("^HSP",rownames(ITP_4))]
length(c(mt.genes,rb.genes,hsp.genes))
#[1] 138
dim(ITP_4)
#[1] 20980  9626
counts <- GetAssayData(ITP_4, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c(mt.genes,rb.genes,hsp.genes))),]
ITP_4 <- subset(ITP_4, features = rownames(counts))
dim(ITP_4)
#[1] 20842  9626
saveRDS(ITP_4, file = 'reITP_4.rds')


##########
dir.create('/home/erin/Desktop/itp/samplesQC3/ITP_5')
setwd('/home/erin/Desktop/itp/samplesQC3/ITP_5')
ITP_5.data <- Read10X(data.dir =  '/home/erin/Desktop/itp/input/P05_YRK/filtered_feature_bc_matrix')
ITP_5 <- CreateSeuratObject(counts =ITP_5.data, project = "ITP_5", min.cells = 3, min.features = 200)
dim(ITP_5)
#[1] 19879  7155

ITP_5$group<-"ITP"
ITP_5$gender<-"M"

ITP_5[["percent.mt"]] <- PercentageFeatureSet(ITP_5, pattern = "^MT-")
ITP_5[["percent.rb"]] <- PercentageFeatureSet(ITP_5, pattern = "^RP[SL]")
ITP_5[["percent.hsp"]] <- PercentageFeatureSet(ITP_5, pattern = "^HSP")

head(ITP_5@meta.data)

ITP_5 <- subset(ITP_5, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 15 & percent.rb < 48.4 & percent.hsp < 2.6)
dim(ITP_5)
#[1] 19879  6715 
ITP_5 <- NormalizeData(ITP_5)
ITP_5 <- FindVariableFeatures(ITP_5, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(ITP_5)
ITP_5 <- ScaleData(ITP_5, features = all.genes)

ITP_5 <- RunPCA(ITP_5, npcs = 50)

ITP_5 <- JackStraw(ITP_5, dims = 50, num.replicate = 100)
ITP_5 <- ScoreJackStraw(ITP_5, dims = 1:50)
p1 <- JackStrawPlot(ITP_5, dims = 1:50)
png(file="pc_pvalue.png", width = dpi*12, height = dpi*6, units = "px",res = dpi,type='cairo')
print(p1)
dev.off()
p2 <- ElbowPlot(ITP_5,ndims =50)
png(file="pc_elbowplot.png", width = dpi*9, height = dpi*6, units = "px",res = dpi,type='cairo')
print(p2)
dev.off()

ITP_5 <- RunUMAP(ITP_5, dims = 1:30)
p <- DimPlot(ITP_5)
png(file="dimplot1.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(ITP_5,
                 features = "CD34")
png(file="dimplot2.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(ITP_5,
                 features = "HLF")#hsc
png(file="dimplot3.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(ITP_5,
                 features = "NOG")#hsc
png(file="dimplot4.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(ITP_5,
                 features = "EBF1")#lymphoid
png(file="dimplot5.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(ITP_5,
                 features = "CEBPD")#myeloid
png(file="dimplot6.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- FeaturePlot(ITP_5,
                 features = "GATA1")#erythroid")
png(file="dimplot7.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()

sweep.res.list_ITP_5 <- paramSweep_v3(ITP_5, PCs = 1:30, sct = FALSE)

sweep.stats_ITP_5 <- summarizeSweep(sweep.res.list_ITP_5, GT = FALSE)

bcmvn_ITP_5 <- find.pK(sweep.stats_ITP_5)

mpK<-as.numeric(as.vector(bcmvn_ITP_5$pK[which.max(bcmvn_ITP_5$BCmetric)]))
mpK

annotations <- ITP_5@meta.data$seurat_clusters

homotypic.prop <- modelHomotypic(annotations)           
nExp_poi <- round(0.055*length(ITP_5@meta.data$orig.ident))  #7179
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

ITP_5 <- doubletFinder_v3(ITP_5, PCs = 1:30, pN = 0.25, pK = mpK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
head(ITP_5@meta.data)
table(ITP_5@meta.data$DF.classifications_0.25_0.01_369)
#Doublet Singlet 
#369    6346
p <- DimPlot(ITP_5,group.by = 'DF.classifications_0.25_0.01_369')
png(file="dimplot8.png", width = dpi*9, height = dpi*8, units = "px",res = dpi,type='cairo')
print(p)
dev.off()
p <- DimPlot(ITP_5,split.by = 'DF.classifications_0.25_0.01_369')
png(file="dimplot9.png", width = dpi*9, height = dpi*5, units = "px",res = dpi,type='cairo')
print(p)
dev.off()

ITP_5@meta.data$DF.classification <- ITP_5@meta.data$DF.classifications_0.25_0.01_369
head(ITP_5@meta.data)

saveRDS(ITP_5, file = 'ITP_5.rds')

setwd('/home/erin/Desktop/itp/samplesQC3/removegenes')
mt.genes <- rownames(ITP_5)[grep("^MT-",rownames(ITP_5))]
rb.genes <- rownames(ITP_5)[grep("^RP[SL]",rownames(ITP_5))]
hsp.genes <- rownames(ITP_5)[grep("^HSP",rownames(ITP_5))]
length(c(mt.genes,rb.genes,hsp.genes))
#[1] 140
dim(ITP_5)
#[1] 19879  6715
counts <- GetAssayData(ITP_5, assay = "RNA")
counts <- counts[-(which(rownames(counts) %in% c(mt.genes,rb.genes,hsp.genes))),]
ITP_5 <- subset(ITP_5, features = rownames(counts))
dim(ITP_5)
#[1] 19739  6715
saveRDS(ITP_5, file = 'reITP_5.rds')
