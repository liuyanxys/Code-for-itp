rm(list=ls())
###Loess-smoothed curves fitted to the z-scored averaged expression of genes in the modules 1 to 4
library(monocle)
library(Seurat)
setwd("/home/monocle2/")

load('my_cds_subset_yoursample_1100.Rdata')
SO<-readRDS("yoursample.rds")
exp<- exprs(my_cds_subset)
t<-read.table("pseudotime-yoursample-6clusters-annotation_row.csv", sep='\t',header=T)
rownames(t)<-t$Gene
exp<-exp[t$Gene,]
exp0<- SO@assays$RNA@data 
exp0<- exp0[rownames(exp),colnames(exp)] 
cluster1gene<-t[t$Cluster=="6",]$Gene
cluster2gene<-t[t$Cluster=="5",]$Gene
cluster3gene<-t[t$Cluster=="3",]$Gene
cluster4gene<-t[t$Cluster=="2",]$Gene
cluster5gene<-t[t$Cluster=="1",]$Gene
cluster6gene<-t[t$Cluster=="4",]$Gene
exp1<-as.matrix(exp0[cluster1gene,])
exp2<-as.matrix(exp0[cluster2gene,])
exp1<-rbind(exp1,exp2)
exp3<-as.matrix(exp0[cluster3gene,])
exp4<-as.matrix(exp0[cluster4gene,])
exp2<-rbind(exp3,exp4)
exp5<-as.matrix(exp0[cluster5gene,])
exp3<-exp5
exp6<-as.matrix(exp0[cluster6gene,])
exp4<-exp6
standardize <- function(x) 
{
  rowmean <- apply(x, 1, mean)
  rowsd <- apply(x, 1, sd)  
  rv <- sweep(x, 1, rowmean,"-")
  rv <- sweep(rv, 1, rowsd, "/")
  return(rv)
}

exp1<-standardize(exp1)
exp2<-standardize(exp2)
exp3<-standardize(exp3)
exp4<-standardize(exp4)
mean1=apply(exp1,2,mean)
mean2=apply(exp2,2,mean)
mean3=apply(exp3,2,mean)
mean4=apply(exp4,2,mean)
cds_exprs<-rbind(mean1,mean2,mean3,mean4)
cds_exprs<- reshape2::melt(cds_exprs)
colnames(cds_exprs) <- c("f_id", "Cell", "expression")
cds_subset<-my_cds_subset
cds_pData <- pData(cds_subset)
cds_fData <- fData(cds_subset)
cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
cds_exprs$adjusted_expression <- cds_exprs$expression
cds_exprs$f_id <- as.character(cds_exprs$f_id)
new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime)

monocle_theme_opts <- function()
{
    theme(strip.background = element_rect(colour = 'white', fill = 'white')) +
    theme(panel.border = element_blank()) +
    theme(axis.line.x = element_line(size=0.25, color="black")) +
    theme(axis.line.y = element_line(size=0.25, color="black")) +
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
    theme(panel.background = element_rect(fill='white')) +
    theme(legend.key=element_blank())
}

cols<-c("mean1" = "#f29228", "mean2" = "#f2eb3b","mean3"="#72bd56","mean4"="#58c2df")
cds_exprs1<-subset(cds_exprs,f_id=="mean1")
q <- ggplot(aes(Pseudotime, expression), data = cds_exprs1)
q <- q +monocle_theme_opts()+
     geom_smooth(aes(x = Pseudotime, y = expression,color=f_id,fill=f_id),method='loess',se = T)+
     scale_colour_manual(values = cols)+scale_fill_manual(values =cols)+
     ylab("Expression")+ xlab("Pseudotime")+
     scale_x_continuous(breaks=seq(-4,18,4))
ggsave(q, file="zscored_expression_1.pdf", width = 10, height = 5)

cds_exprs2<-subset(cds_exprs,f_id=="mean2")
q <- ggplot(aes(Pseudotime, expression), data = cds_exprs2)
q <- q +monocle_theme_opts()+
     geom_smooth(aes(x = Pseudotime, y = expression,color=f_id,fill=f_id),method='loess',se = T)+
     scale_colour_manual(values = cols)+scale_fill_manual(values =cols)+
     ylab("Expression")+ xlab("Pseudotime")+
     scale_x_continuous(breaks=seq(-4,18,4))
ggsave(q, file="zscored_expression_2.pdf", width = 10, height = 5)

cds_exprs3<-subset(cds_exprs,f_id=="mean3")
q <- ggplot(aes(Pseudotime, expression), data = cds_exprs3)
q <- q +monocle_theme_opts()+
     geom_smooth(aes(x = Pseudotime, y = expression,color=f_id,fill=f_id),method='loess',se = T)+
     scale_colour_manual(values = cols)+scale_fill_manual(values =cols)+
     ylab("Expression")+ xlab("Pseudotime")+
     scale_x_continuous(breaks=seq(-4,18,4))
ggsave(q, file="zscored_expression_3.pdf", width = 10, height = 5)

cds_exprs4<-subset(cds_exprs,f_id=="mean4")
q <- ggplot(aes(Pseudotime, expression), data = cds_exprs4)
q <- q +monocle_theme_opts()+
     geom_smooth(aes(x = Pseudotime, y = expression,color=f_id,fill=f_id),method='loess',se = T)+
     scale_colour_manual(values = cols)+scale_fill_manual(values =cols)+
     ylab("Expression")+ xlab("Pseudotime")+
     scale_x_continuous(breaks=seq(-4,18,4))
ggsave(q, file="zscored_expression_4.pdf", width = 10, height = 5)
