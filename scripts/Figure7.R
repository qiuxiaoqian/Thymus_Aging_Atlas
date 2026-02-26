setwd('D:/OneDrive/Work/Thymus/Manuscript/Figure7')

rm(list=ls())


library(tidyverse)
library(dplyr)
library(reshape2)
library(stringr)
library(Seurat)
library(openxlsx)
library(ggplot2)
library(nichenetr)

#ANGPTI4
dat = readRDS('D:/OneDrive/Work/Thymus/ScSn_integration/anno2/Thymus_scsn_anno2_finalname.rds')
dat = subset(dat,Platform =='snRNA')
dat = subset(dat,Level1 != 'Doublets')

celltype_metadata = read.xlsx('D:/OneDrive/Work/Thymus/ScSn_integration/data/celltype_metadata.xlsx')
dat$Level1 = factor(dat$Level1,levels = rev(unique(celltype_metadata$Level1)))


gene = 'ANGPTL4'
gene_set = 'ANGPTL4'

gene = 'FGF13'
gene_set = 'FGF13'

png(str_glue('dat_umap_{gene_set}.png'),res = 300,height = 5.5,width = 6.2,units = 'in')
FeaturePlot(dat,gene,cols = c("lightgrey", "#ff0000"))+
  theme(plot.title = element_text(face = "italic"))+NoAxes()
dev.off() 

pdf(str_glue('dat_umap_{gene_set}.pdf'),height = 5.5,width = 6.2)
FeaturePlot(dat,gene,cols = c("lightgrey", "#ff0000"))+
  theme(plot.title = element_text(face = "italic"))+NoAxes()
dev.off()  


#Dotplot
pdf(str_glue('dat_{gene_set}_dotplot_Level1.pdf'),height = 4,width = 3)
DotPlot(dat,features = gene,cols = c("lightgrey", "#ff0000"),
        scale = T,assay = 'RNA',group.by  = "Level1")+
  theme(panel.border = element_rect(colour = "black", fill=NA,size = 0.4),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black",size = 9),
        axis.ticks = element_line(linewidth = 0.2),
        plot.title = element_text(hjust = 0.5,face = 'italic',size = 10),
        legend.title = element_text(size = 10),
        egend.text = element_text(size = 10),
        axis.line = element_blank())+
  labs(x = '',y = 'Cell type',title = gene_set)
dev.off()





#NicheNet

input_path = 'D:/OneDrive/Work/Thymus/ScSn_integration/Analysis/NicheNet'


colname.cluster = 'Level1_label'
receiver = 'c07_cTEC'
receiver = 'c08_mTEC'
print(receiver)


platform = 'snRNA'
print(platform)
regulation = 'significant'


# ## Read in the expression data of interacting cells
# dat = readRDS('D:/OneDrive/Work/Thymus/ScSn_integration/anno2/Thymus_scsn_anno2_finalname.rds')
# cat('Total cell of dat:',dim(dat)[2],'\n')
# table(dat$Platform,dat$Tissue)
# 
# 
# dat = subset(dat,Tissue == 'Thymus')
# cat(str_glue('Total cell of Thymus:'),dim(dat)[2],'\n')
# 
# dat$Level1_label = factor(dat$Level1_label,sort(unique(dat$Level1_label)))
# dat$Level2_label = factor(dat$Level2_label,sort(unique(dat$Level2_label)))
# 
# names(dat@meta.data)
# 
# data = subset(dat, Platform == platform)
# cat(str_glue('Total cell of Thymus by {platform}'),dim(data)[2],'\n')#113676
# saveRDS(data,'Thymus_sn_anno2.rds')

data = readRDS('Thymus_sn_anno2.rds')




ligand_activities = read.xlsx(str_glue('{input_path}/{colname.cluster}/{platform}_{receiver}_{regulation}_ligand_activities.xlsx'),rowNames=F)

#1. ligand图   
best_upstream_ligands = ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()

#We can also visualize the ligand activity measure (AUPR) of these top-ranked ligands:
vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)

p_ligand_aupr = vis_ligand_aupr %>% 
  make_heatmap_ggplot( "Prioritized ligands", "Ligand activity",
                       legend_title = "AUPR", color = "darkorange") +
  geom_tile(color = "white", size = 0.05) + 
  theme(axis.text.x.top = element_blank(),
        legend.position = 'right')

p_ligand_aupr
ggsave(str_glue('{platform}_{receiver}_{regulation}_ligand_aupr.pdf'),p_ligand_aupr,width = 3,height = 5)

    
# 2. ligand表达
Idents(data) = colname.cluster; width = 9
pdf(str_glue('{platform}_{receiver}_{regulation}_best_upstream_ligands_dotplot.pdf'),height = 8,width = width )
p=DotPlot(data,features = rev(best_upstream_ligands), scale = T,assay = 'RNA')+
  theme(panel.border = element_rect(colour = "black", fill=NA,size = 1))+
  theme(axis.text.x = element_text(angle = 90,hjust =1,vjust = 0.5))+
  coord_flip()+
  labs(x='',y='')
print(p)
dev.off()


# 3. ligand-target gene
vis_ligand_target = read.xlsx(str_glue('{input_path}/{colname.cluster}/{platform}_{receiver}_{regulation}_ligand_target_matrix.xlsx'),rowNames = T)
vis_ligand_target = as.matrix(vis_ligand_target)    

p_ligand_target_network = vis_ligand_target %>% 
  make_heatmap_ggplot("Prioritized ligands","Predicted target genes", 
                      color = "purple",legend_position = "right", x_axis_position = "top",
                      legend_title = "Regulatory potential")  + 
  theme(axis.text.x = element_text(face = "italic",size = 8,hjust = 0.5),
        plot.title = element_text(hjust = 0.5)) + 
  labs(title =str_glue('{receiver}_{regulation}'))+
  scale_fill_gradient2(low = "#FFC9C9",  high = "#7F0000")#, (low = "whitesmoke",  high = "purple")breaks = c(0,0.0045,0.0090)

p_ligand_target_network 
ggsave(str_glue('{platform}_{receiver}_{regulation}_ligand_target_network.pdf'),p_ligand_target_network,width = 15,height = 5)


    
#②Receptors of top-ranked ligands
#which receptors have the highest interaction potential with the top-ranked ligands.
vis_ligand_receptor_network = read.xlsx(str_glue('{input_path}//{colname.cluster}/{platform}_{receiver}_{regulation}_ligand_receptor_matrix.xlsx'),rowNames=T)
vis_ligand_receptor_network = vis_ligand_receptor_network[rev(best_upstream_ligands),]
vis_ligand_receptor_network = as.matrix(t(vis_ligand_receptor_network))

p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% 
  make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred",
                      legend_position = "right",  #x_axis_position = "top",
                      legend_title = "Prior interaction potential")+
  theme(axis.text.x = element_text(face = "italic",colour = 'black'),
        axis.text.y = element_text(colour = 'black'),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  geom_tile(color = "white", size = 0.1) + 
  labs(title =str_glue('{receiver}_{regulation}'))

p_ligand_receptor_network
ggsave(str_glue('{platform}_{receiver}_{regulation}_ligand_receptor_network.pdf'),p_ligand_receptor_network,width = 8,height = 4.8)




