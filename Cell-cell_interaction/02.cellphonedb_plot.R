rm(list=ls())
library(ggplot2)

filename = 'plot.png'
means_path = '/home/cellphonedb/out_yoursample/means.txt'
pvalues_path = '/home/cellphonedb/out_yoursample/pvalues.txt'
means_separator = '\t'
pvalues_separator= '\t'

all_pval = read.table(pvalues_path, header=T, stringsAsFactors = F, sep=means_separator, comment.char = '', check.names=F)
all_means = read.table(means_path, header=T, stringsAsFactors = F, sep=pvalues_separator, comment.char = '', check.names=F)

intr_pairs = all_pval$interacting_pair
all_pval = all_pval[,-c(1:11)]
all_means = all_means[,-c(1:11)]

###preB3 fore example
selected_columns = c("HC_preB3|HC_HSC","HC_preB3|HC_MEP","HC_preB3|HC_MkP1","HC_preB3|HC_MkP2","HC_preB3|HC_EryP","ITP_preB3|ITP_HSC","ITP_preB3|ITP_MEP","ITP_preB3|ITP_MkP1","ITP_preB3|ITP_MkP2","ITP_preB3|ITP_EryP")
   
#select pairs by venn
ITPpair<-c("ADGRG5_FAM3C", "ADRB2_VEGFB", "BMP8B_PLAUR", "BTLA_TNFRSF14", "CCL3_IDE", "CCR10_CCL28", "CD1D_LILRB2", "CD44_HBEGF", "CD48_CD244", "CD55_ADGRE5", "CDH1_a2b1 complex", "CDH1_aEb7 complex", "CDH1_KLRG1", "COL18A1_a2b1 complex", "COL4A4_a2b1 complex", "COL9A2_a2b1 complex", "COL9A3_a2b1 complex", "FBN1_a5b1 complex", "FCER2_aMb2 complex", "FCER2_aVb3 complex", "FCER2_CR2", "FGF2_aVb3 complex", "FGF2_CD44", "FGF2_FGFR1", "FGF2_FGFR3", "FGF2_FGFRL1", "GRN_SORT1", "ICAM1_aMb2 complex", "ICAM1_AREG", 
           "ICAM1_SPN", "ICAM2_aLb2 complex", "ICAM3_aLb2 complex", "IDE_CCL23", "IGF1_a6b4 complex", "IGF1_IGF1R", "LAMC1_a2b1 complex", "LAMC1_a6b1 complex", "LAMC1_aVb3 complex", "MDK_LRP1", "MIF_TNFRSF14", "PLAUR_a4b1 complex", "SEMA4A_PLXND1", "SIRPA_CD47", "TNFRSF10A_TNFSF10", "TNFRSF13C_TNFSF13B", "TNFRSF17_TNFSF13B", "TNFSF12_TNFRSF25", "WNT4_FZD1", "WNT4_SMO")
HCpair<-c("HLA-DPA1_GAL", "COL5A1_a2b1 complex", "LAMP1_FAM3C", "LAMP1_VSTM1", "HLA-DPB1_TNFSF13B", "COPA_SORT1", "PLXNB2_SEMA4C", "NOTCH1_DLK1", "IL7 receptor_IL7", "LGALS9_SORL1", "BMP2_SMO", "PLXNB2_SEMA4D", "CD72_SEMA4D", "GDF11_TGFR_AVR2A", "GDF11_TGFR_AVR2B", "EPHB6_EFNB2", "PLXNC1_SEMA7A")
overlappair<-c("CD40_CD40LG", "CD40_TNFSF13B", "CD74_APP", "CD74_COPA", "CD74_MIF", "COL24A1_a2b1 complex", "FAM3C_CLEC2D", "HLA-C_FAM3C", "HLA-F_LILRB2", "LGALS9_CD44", "MDK_SORL1", "PECAM1_CD38", "TNFRSF13B_TNFSF13B")
 
#select pairs through pvalue<0.05 in ITP  
selected_rows = c("ADGRG5_FAM3C", "ADRB2_VEGFB", "BMP8B_PLAUR", "BTLA_TNFRSF14", "CCL3_IDE", "CCR10_CCL28", "CD1D_LILRB2", "CD44_HBEGF", "CD48_CD244", "CD55_ADGRE5", "CDH1_a2b1 complex", "CDH1_aEb7 complex", "CDH1_KLRG1", "COL18A1_a2b1 complex", "COL4A4_a2b1 complex", "COL9A2_a2b1 complex", "COL9A3_a2b1 complex", "FBN1_a5b1 complex", "FCER2_aMb2 complex", "FCER2_aVb3 complex", "FCER2_CR2", "FGF2_aVb3 complex", "FGF2_CD44", "FGF2_FGFR1", "FGF2_FGFR3", "FGF2_FGFRL1", "GRN_SORT1", "ICAM1_aMb2 complex", "ICAM1_AREG", 
      "ICAM1_SPN", "ICAM2_aLb2 complex", "ICAM3_aLb2 complex", "IDE_CCL23", "IGF1_a6b4 complex", "IGF1_IGF1R", "LAMC1_a2b1 complex", "LAMC1_a6b1 complex", "LAMC1_aVb3 complex", "MDK_LRP1", "MIF_TNFRSF14", "PLAUR_a4b1 complex", "SEMA4A_PLXND1", "SIRPA_CD47", "TNFRSF10A_TNFSF10", "TNFRSF13C_TNFSF13B", "TNFRSF17_TNFSF13B", "TNFSF12_TNFRSF25", "WNT4_FZD1", "WNT4_SMO","CD40_CD40LG", "CD40_TNFSF13B", "CD74_APP", "CD74_COPA", "CD74_MIF", "COL24A1_a2b1 complex", "FAM3C_CLEC2D", "HLA-C_FAM3C", "HLA-F_LILRB2", "LGALS9_CD44", "MDK_SORL1", "PECAM1_CD38", "TNFRSF13B_TNFSF13B")
  
sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]
df_names = expand.grid(selected_rows, selected_columns)
pval = unlist(sel_pval)
pval[pval==0] = 0.0009
plot.data = cbind(df_names,pval)
pr = unlist(as.data.frame(sel_means))
#pr[pr==0] = 1
#plot.data = cbind(plot.data,log2(pr+1))
plot.data = cbind(plot.data,pr)
colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')

plot.data$group<-ifelse(substr(plot.data$clusters,1,2)=="HC","HC","ITP")
plot.data$venn<-plot.data$group

for (i in 1:length(plot.data$pair)){
  if (plot.data$pair[i]%in%ITPpair){
    plot.data$venn[i]<-"ITP"
  } else {
    plot.data$venn[i]<-"Overlap"
  }
}

plot.data$venn<-factor(plot.data$venn,levels=c("ITP","Overlap"), ordered=TRUE)
my_palette <- colorRampPalette(c("blue", "red"), alpha=TRUE)(n=399)
ggplot(plot.data,aes(x=clusters,y=pair)) +
      geom_point(aes(size=-log10(pvalue+ 0.001),color=mean)) +scale_size_continuous(range=c(1,7),breaks = c(0,0.5,1.0,1.5,2.0))+
      scale_color_gradientn('Mean expression', colors=my_palette,limits=c(0,4)) +
      theme_bw() +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.text=element_text(size=10, colour = "black"),
            axis.text.x = element_text(angle = 90, hjust = 1),
            axis.text.y = element_text(size=8, colour = "black"),
            axis.title=element_blank(),
            panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))+
      facet_grid(venn~group,scales = "free",space="free")
ggsave('D://X//ITP//cellphonedb//out_itp.data1.2//preB3&NKT_venn//cellphonedb-preB3.pdf', width = 10, height =18, limitsize=F)
