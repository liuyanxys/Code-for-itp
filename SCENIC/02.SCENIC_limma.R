rm(list=ls())
library(edgeR)
library(limma)
library("impute")
library(SCENIC)
library(AUCell)
library(data.table)
setwd("/home/SCENIC/limma/")
regulonAUC <- importAUCfromText("/home/SCENIC/results/yoursample.auc.csv")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
SO<-readRDS("/home/yoursample.rds")
macro.cellinfo<-SO@meta.data
macro.cellinfo$seurat_clusters <-macro.cellinfo$cell.types 

macro.cellinfo<-subset(macro.cellinfo,macro.cellinfo$seurat_clusters=="preB3")#preB3 for example
targets<-data.table(FileName=rownames(macro.cellinfo),Target=macro.cellinfo$group)
lev<-unique(targets$Target)##Use the unique () function for de-repetition
f <- factor(targets$Target, levels=lev) 
design <- model.matrix(~0+f) #The sample matrix
colnames(design) <- lev #Change the column name to levels name  
eset=getAUC(regulonAUC) #expression matrix
eset=eset[,targets$FileName]
eset<-t(scale(t(eset))) #scale
###The topTable function was used to find the difference features for comparison between pairs
cont.wt <- makeContrasts("ITP-HC",levels=design) 
fit <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, cont.wt) 
fit2 <- eBayes(fit2) 
tT=topTable(fit2, adjust="BH",sort.by="logFC",n=Inf)
tT = subset(tT, select=c("adj.P.Val","P.Value","logFC"))
colnames(tT)=c("FDR","P.Value","logFC")
write.table(tT,file="tT-preB3.csv",sep="\t",quote=F)
logFoldChange<-0.25
adjustP<-0.05
diffSig <- tT[with(tT, (abs(logFC)>logFoldChange & FDR < adjustP )), ]
write.table(diffSig,file="diff-preB3.csv",sep="\t",quote=F)
diffUp <- tT[with(tT, (logFC>logFoldChange & FDR < adjustP )), ]
write.table(diffUp,file="up-preB3.csv",sep="\t",quote=F)
diffDown <- tT[with(tT, (logFC<(-logFoldChange) & FDR < adjustP )), ]
write.table(diffDown,file="down-preB3.csv",sep="\t",quote=F)

##To estimate significant TFs in each ITP sample
macro.cellinfo.ITP1<-subset(macro.cellinfo,macro.cellinfo$orig.ident=="ITP_1")
macro.cellinfo.ITP2<-subset(macro.cellinfo,macro.cellinfo$orig.ident=="ITP_2")
macro.cellinfo.ITP3<-subset(macro.cellinfo,macro.cellinfo$orig.ident=="ITP_3")
macro.cellinfo.ITP4<-subset(macro.cellinfo,macro.cellinfo$orig.ident=="ITP_4")
macro.cellinfo.HC<-subset(macro.cellinfo,macro.cellinfo$group=="HC")

macro.cellinfo.ITP1.HC<-rbind(macro.cellinfo.ITP1,macro.cellinfo.HC)
targets<-data.table(FileName=rownames(macro.cellinfo.ITP1.HC),Target=macro.cellinfo.ITP1.HC$group)
lev<-unique(targets$Target)##Use the unique () function for de-repetition
f <- factor(targets$Target, levels=lev) 
design <- model.matrix(~0+f) #The sample matrix
colnames(design) <- lev #Change the column name to levels name  
eset=getAUC(regulonAUC) #expression matrix
eset=eset[,targets$FileName]
eset<-t(scale(t(eset))) #scale
###The topTable function was used to find the difference features for comparison between pairs
cont.wt <- makeContrasts("ITP-HC",levels=design) 
fit <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, cont.wt) 
fit2 <- eBayes(fit2) 
tT=topTable(fit2, adjust="BH",sort.by="logFC",n=Inf)
tT = subset(tT, select=c("adj.P.Val","P.Value","logFC"))
colnames(tT)=c("FDR","P.Value","logFC")
write.table(tT,file="tT-preB3.ITP1.HC.csv",sep="\t",quote=F)
logFoldChange<-0.25
adjustP<-0.05
diffSig <- tT[with(tT, (abs(logFC)>logFoldChange & FDR < adjustP )), ]
write.table(diffSig,file="diff-preB3.ITP1.HC.csv",sep="\t",quote=F)
diffUp <- tT[with(tT, (logFC>logFoldChange & FDR < adjustP )), ]
write.table(diffUp,file="up-preB3.ITP1.HC.csv",sep="\t",quote=F)
diffDown <- tT[with(tT, (logFC<(-logFoldChange) & FDR < adjustP )), ]
write.table(diffDown,file="down-preB3.ITP1.HC.csv",sep="\t",quote=F)

macro.cellinfo.ITP2.HC<-rbind(macro.cellinfo.ITP2,macro.cellinfo.HC)
targets<-data.table(FileName=rownames(macro.cellinfo.ITP2.HC),Target=macro.cellinfo.ITP2.HC$group)
lev<-unique(targets$Target)##Use the unique () function for de-repetition
f <- factor(targets$Target, levels=lev) 
design <- model.matrix(~0+f) #The sample matrix
colnames(design) <- lev #Change the column name to levels name  
eset=getAUC(regulonAUC) #expression matrix
eset=eset[,targets$FileName]
eset<-t(scale(t(eset))) #scale
###The topTable function was used to find the difference features for comparison between pairs
cont.wt <- makeContrasts("ITP-HC",levels=design) 
fit <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, cont.wt) 
fit2 <- eBayes(fit2) 
tT=topTable(fit2, adjust="BH",sort.by="logFC",n=Inf)
tT = subset(tT, select=c("adj.P.Val","P.Value","logFC"))
colnames(tT)=c("FDR","P.Value","logFC")
write.table(tT,file="tT-preB3.ITP2.HC.csv",sep="\t",quote=F)
logFoldChange<-0.25
adjustP<-0.05
diffSig <- tT[with(tT, (abs(logFC)>logFoldChange & FDR < adjustP )), ]
write.table(diffSig,file="diff-preB3.ITP2.HC.csv",sep="\t",quote=F)
diffUp <- tT[with(tT, (logFC>logFoldChange & FDR < adjustP )), ]
write.table(diffUp,file="up-preB3.ITP2.HC.csv",sep="\t",quote=F)
diffDown <- tT[with(tT, (logFC<(-logFoldChange) & FDR < adjustP )), ]
write.table(diffDown,file="down-preB3.ITP2.HC.csv",sep="\t",quote=F)

macro.cellinfo.ITP3.HC<-rbind(macro.cellinfo.ITP3,macro.cellinfo.HC)
targets<-data.table(FileName=rownames(macro.cellinfo.ITP3.HC),Target=macro.cellinfo.ITP3.HC$group)
lev<-unique(targets$Target)##Use the unique () function for de-repetition
f <- factor(targets$Target, levels=lev) 
design <- model.matrix(~0+f) #The sample matrix
colnames(design) <- lev #Change the column name to levels name  
eset=getAUC(regulonAUC) #expression matrix
eset=eset[,targets$FileName]
eset<-t(scale(t(eset))) #scale
###The topTable function was used to find the difference features for comparison between pairs
cont.wt <- makeContrasts("ITP-HC",levels=design) 
fit <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, cont.wt) 
fit2 <- eBayes(fit2) 
tT=topTable(fit2, adjust="BH",sort.by="logFC",n=Inf)
tT = subset(tT, select=c("adj.P.Val","P.Value","logFC"))
colnames(tT)=c("FDR","P.Value","logFC")
write.table(tT,file="tT-preB3.ITP3.HC.csv",sep="\t",quote=F)
logFoldChange<-0.25
adjustP<-0.05
diffSig <- tT[with(tT, (abs(logFC)>logFoldChange & FDR < adjustP )), ]
write.table(diffSig,file="diff-preB3.ITP3.HC.csv",sep="\t",quote=F)
diffUp <- tT[with(tT, (logFC>logFoldChange & FDR < adjustP )), ]
write.table(diffUp,file="up-preB3.ITP3.HC.csv",sep="\t",quote=F)
diffDown <- tT[with(tT, (logFC<(-logFoldChange) & FDR < adjustP )), ]
write.table(diffDown,file="down-preB3.ITP3.HC.csv",sep="\t",quote=F)

macro.cellinfo.ITP4.HC<-rbind(macro.cellinfo.ITP4,macro.cellinfo.HC)
targets<-data.table(FileName=rownames(macro.cellinfo.ITP4.HC),Target=macro.cellinfo.ITP4.HC$group)
lev<-unique(targets$Target)##Use the unique () function for de-repetition
f <- factor(targets$Target, levels=lev) 
design <- model.matrix(~0+f) #The sample matrix
colnames(design) <- lev #Change the column name to levels name  
eset=getAUC(regulonAUC) #expression matrix
eset=eset[,targets$FileName]
eset<-t(scale(t(eset))) #scale
###The topTable function was used to find the difference features for comparison between pairs
cont.wt <- makeContrasts("ITP-HC",levels=design) 
fit <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, cont.wt) 
fit2 <- eBayes(fit2) 
tT=topTable(fit2, adjust="BH",sort.by="logFC",n=Inf)
tT = subset(tT, select=c("adj.P.Val","P.Value","logFC"))
colnames(tT)=c("FDR","P.Value","logFC")
write.table(tT,file="tT-preB3.ITP4.HC.csv",sep="\t",quote=F)
logFoldChange<-0.25
adjustP<-0.05
diffSig <- tT[with(tT, (abs(logFC)>logFoldChange & FDR < adjustP )), ]
write.table(diffSig,file="diff-preB3.ITP4.HC.csv",sep="\t",quote=F)
diffUp <- tT[with(tT, (logFC>logFoldChange & FDR < adjustP )), ]
write.table(diffUp,file="up-preB3.ITP4.HC.csv",sep="\t",quote=F)
diffDown <- tT[with(tT, (logFC<(-logFoldChange) & FDR < adjustP )), ]
write.table(diffDown,file="down-preB3.ITP4.HC.csv",sep="\t",quote=F)
