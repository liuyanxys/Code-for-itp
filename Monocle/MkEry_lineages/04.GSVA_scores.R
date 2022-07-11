###GSVA score
library(GSVA)
library(GSEABase)
library(GSVAdata)
library(limma)
library(data.table)
library(ggplot2)
setwd("/home/monocle2/MkEry/")

geneset <- getGmt("pathway.gmt")
load('my_cds_subset_yoursample_MkEry_500.Rdata')
SO<-readRDS("/home/monocle2/yoursample.rds")
exp<- exprs(my_cds_subset)
exp0<- SO@assays$RNA@data
exp0<- exp0[rownames(exp),colnames(exp)] 
rm(SO)
es <- gsva(as.matrix(exp0), geneset,method="zscore", verbose=TRUE, parallel.sz=10)
save(es,file="es.pathway.MkEry.Rdata")
#load('es.pathway.MkEry.Rdata')
cds_exprs<-es
cds_exprs<- reshape2::MkErylt(cds_exprs)
colnames(cds_exprs) <- c("f_id", "Cell", "expression")
f_id<-as.character(cds_exprs$f_id)
f_id[f_id=="GOBP_COMPLEMENT_ACTIVATION"]<-"Complement activation" 
f_id[f_id=="GOBP_HUMORAL_IMMUNE_RESPONSE"]<-"Humoral immune response"
f_id[f_id=="GOBP_HUMORAL_IMMUNE_RESPONSE_MEDIATED_BY_CIRCULATING_IMMUNOGLOBULIN"]<-"Humoral immune response mediated by circulating immunoglobulin"
f_id[f_id=="GOBP_REGULATION_OF_COMPLEMENT_ACTIVATION"]<-"Regulation of complement activation"
f_id[f_id=="GOBP_REGULATION_OF_HUMORAL_IMMUNE_RESPONSE"]<-"Regulation of humoral immune response"
f_id[f_id=="HP_ABNORMALITY_OF_COMPLEMENT_SYSTEM"]<-"Abnormality of complement system"
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
q <- ggplot(aes(Pseudotime, expression,color=group,group=group,fill=group), data = cds_exprs0)+monocle_theme_opts()+
     geom_smooth(aes(x = Pseudotime, y = expression,color=group,fill=group),method='loess',se = T, span=0.4)+
     scale_colour_manual(values = cols)+scale_fill_manual(values = cols)+
     facet_wrap(~f_id, nrow = 3,ncol = 2, scales = "free_y")+
     ylab("Feaure scores")+ xlab("Pseudotime")+
     scale_x_continuous(breaks=seq(-5,30,5))
ggsave(q, file="GSVA_pathway_zscore_MkEry_500_ITP_versus_HC.pdf", width = 16, height =12)
