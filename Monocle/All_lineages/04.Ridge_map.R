rm(list=ls())
###The ridge map
library(monocle)
library(tidyverse)
library(ggridges)
library(RColorBrewer)
library(scales)

load('my_cds_subset_yoursample_1100.Rdata')
plotdf=pData(my_cds_subset)
plotdf$cell.types<-factor(plotdf$cell.types,levels=c("EryP","MkP2","MkP1","MEP","NK/Tp","preB3","preB2","preB1","CLP","G2M","MDP","NeuP","GMP","MPP","HSC"), ordered=TRUE)
cols<-c("HSC"="#D9534F","MPP"="#96CEB4", "GMP"="#CBE86B", "NeuP"="#EDE574","MDP"="#0099CC","G2M"="#66CCB9","CLP"="#FF9999", "preB1"="#99CCFF","preB2"="#88A6AD","preB3"="#996666","NK/Tp"="#3AC569","MEP"="#D071A9","MkP1"="#9ED47D","MkP2"="#F57F2C","EryP"="#CCCCFF")
ggplot(plotdf, aes(x=Pseudotime,y=cell.types,fill=cell.types))+
       geom_density_ridges(scale=1) +
       scale_y_discrete("")+
       theme_minimal()+
       theme(panel.grid = element_blank())+scale_fill_manual(values = cols)+ scale_x_continuous(breaks=seq(-4,18,4))
ggsave("pseudotime-ridgemap-yoursample-1100.pdf",width = 15,height = 6)
