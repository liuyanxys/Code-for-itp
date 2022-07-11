rm(list=ls())
###GSVA score
library(GSVA)
library(GSEABase)
library(GSVAdata)
library(limma)
library(data.table)
library(ggplot2)
setwd("/home/monocle2/")

geneset <- getGmt("c5.gmt")  
load('my_cds_subset_yoursample_1100.Rdata')
SO<-readRDS("yoursample.rds")
exp<- exprs(my_cds_subset)
exp0<- SO@assays$RNA@data 
exp0<- exp0[rownames(exp),colnames(exp)] 
rm(SO)
es <- gsva(as.matrix(exp0), geneset, method="zscore", verbose=TRUE, parallel.sz=10)
#save(es,file="es.c5.Rdata")
#load(file="es.c5.Rdata")
selpathway<-c("GOBP_REGULATION_OF_HUMORAL_IMMUNE_RESPONSE","GOBP_REGULATION_OF_COMPLEMENT_ACTIVATION","HP_ABNORMALITY_OF_COMPLEMENT_SYSTEM","HALLMARK_HEME_METABOLISM","HP_RETICULOCYTOSIS","GOMF_COMPLEMENT_BINDING")
selpathway%in%rownames(es)
cds_exprs<-es[selpathway,]
cds_exprs<- reshape2::melt(cds_exprs)
colnames(cds_exprs) <- c("f_id", "Cell", "expression")
f_id<-as.character(cds_exprs$f_id)
f_id[f_id=="GOBP_REGULATION_OF_HUMORAL_IMMUNE_RESPONSE"]<-"Regulation of humoral immune response" 
f_id[f_id=="GOBP_REGULATION_OF_COMPLEMENT_ACTIVATION"]<-"Regulation of complement activation" 
f_id[f_id=="HP_ABNORMALITY_OF_COMPLEMENT_SYSTEM"]<-"Abnormality of complement system"
f_id[f_id=="HALLMARK_HEME_METABOLISM"]<-"Heme metabolism"
f_id[f_id=="HP_RETICULOCYTOSIS"]<-"Reticulocytosis"
f_id[f_id=="GOMF_COMPLEMENT_BINDING"]<-"Complement binding"
cds_exprs$f_id<-f_id
cds_subset<-my_cds_subset
cds_pData <- pData(cds_subset)
cds_fData <- fData(cds_subset)
cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
cds_exprs$adjusted_expression <- cds_exprs$expression
cds_exprs$f_id <- as.character(cds_exprs$f_id)
new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime)
cds_exprs$group<-factor(cds_exprs$group)
cds_exprs$f_id<-factor(cds_exprs$f_id)
cols<-c("HC" = "#619cff", "ITP" = "#f8766d")
q <- ggplot(aes(Pseudotime, expression,color=group,group=group,fill=group), data = cds_exprs)
q <- q +monocle_theme_opts()+
     geom_smooth(aes(x = Pseudotime, y = expression,color=group,fill=group),method='loess',se = T)+
     scale_colour_manual(values = cols)+scale_fill_manual(values = cols)+
     facet_wrap(~f_id, nrow = 5,ncol = 2, scales = "free_y")+
     ylab("Feaure scores")+ xlab("Pseudotime")+
     scale_x_continuous(breaks=seq(-4,18,4))
ggsave(q, file="plot_genes_in_pseudotime_yoursample_pathway_sel.pdf", width = 16, height = 12)
