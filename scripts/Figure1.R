
setwd('D:/OneDrive/Work/Thymus/Manuscript/Figure1/Plot')
setwd('C:/Users/qxq/OneDrive/Work/Thymus/Manuscript/Figure1')


rm(list=ls())
library(dplyr)
library(reshape2)
library(stringr)
library(Seurat)
library(openxlsx)
library(ggplot2)


dat = readRDS('../../../ScSn_integration/anno2/Thymus_scsn_anno2_finalname.rds')

table(dat$Level2_label) %>% data.frame()


# ① 绘制umap
dim(dat) #33400 405395
dat = subset(dat,Major != 'Doublets')
dim(dat)#33400 394124

metadata = dat@meta.data

L1_color = c("#A6CEE3", "#1F78B4", "#EDB869", "#FB9A99", "#B2DF8A", "#84D7F7", 
             "#D8B6F7", "#E3EDC0", "#E78AC3", "#D4EEF6", "#FC8D62", "#8DA0CB", 
             "#66C2A5", "#6F89D3", "#5050FF", "#F9C9BD", "#F0E685", "#33A02C", 
             "#F07590", "#4DBBD5", "#FFD92F", "#E31A1C", "#E16CB7", "#B3E2CD", 
             "#EEBED6", "#7FC28E", "#5DB4F3", "#FFF2AE", "#B589BC", "#CD7560", 
             "#FDBF6F", "#DCDCDC") %>% unique()


png('dat_umap_level1_main.png',res = 300,height = 5.5,width = 6,units = 'in')
DimPlot(dat,label = T,cols = L1_color,raster = F,group.by = 'L1_label',
        label.size = 3,repel = T,alpha = 0.6)+labs(title = '')+
  NoAxes()+theme(legend.position = 'none')
dev.off()



L2_color = c("#F07590", "#A6CEE3", "#F7D39B","#00BFC4","#FF62BC","#1F78B4", 
             "#66C2A5", "#F0E685", "#EFE0E7","#EDB869","#84D7F7","#E76BF3",
             "#00BF7D", "#00B0F6", "#D8B6F7","#D89000","#CCEBC5","#FF7F00", 
             "#A3A500", "#A032CB", "#619CFF","#F8766D","#E78AC3","#FFF2AE", 
             "#D4EEF6", "#EEBED6", "#CD7560","#00C08B","#917393",
             "#FC8D62", "#77A1EC", "#EA94CD","#8DA0CB","#6F89D3","#5050FF",# 35 Mast
             "#F9C9BD", "#FEF6AE", "#CCECFF","#FB9A99","#39B600","#FFC107", 
             "#4DBBD5", "#7FC28E", "#98CBE4","#E31A1C","#DEDBEE","#FED9A6",# 
             "#33A02C", "#F5BCC9", "#B2DF8A","#FC7440","#C76DA8","#BAB7FF",#53
             "#B3E2CD", "#E16CB7", "#F48F89","#F7F4B7", #57
             "#5394C3", "#BEAED4", "#FFD92F","#F4CAE4","#D4B595", "#F7E3DB",
             "#CBD5E8", "#5DB4F3", "#B099B5","#AF88BB", "#D59CBD", "#91CDC1") %>% unique()



options(ggrepel.max.overlaps = Inf)
png('dat_umap_level2_label.png',res = 300,height = 5.5,width = 8,units = 'in')
DimPlot(dat,label = T,cols = L2_color,raster = F,group.by = 'L2_label',
        label.size = 3,repel = T)+labs(title = '')+#,alpha = 0.6
  guides(color = guide_legend(override.aes = list(size = 4),ncol = 3))
  #+NoAxes()+theme(legend.position = 'none')
dev.off()


#pdf(str_glue('dat_umap_level2_label.pdf'),height = 5.5,width = 6)
png('dat_umap_level2_main.png',res = 300,height = 5.5,width = 6,units = 'in')
DimPlot(dat,label = F,cols = L2_color,raster = F,group.by = 'L2_label',
        label.size = 3,repel = T,alpha = 0.6)+labs(title = '')+
  NoAxes()+theme(legend.position = 'none')
dev.off()



#label
obj = dat
meta <- cbind(obj@meta.data, obj@reductions$umap@cell.embeddings)

colnams.cluster = 'L2_label'
centers <- aggregate(meta[, c("umap_1", "umap_2")],list(cluster = meta[, colnams.cluster]),mean)
names(centers) = c(colnams.cluster,"umap_1", "umap_2")
centers[,colnams.cluster] = factor(centers[,colnams.cluster],levels = sort(centers[,colnams.cluster]))

p = ggplot(centers, aes_string("umap_1", "umap_2", label = colnams.cluster)) +
  #geom_point(aes_string(color = colnams.cluster), size = 0.2, show.legend = FALSE) +#
  scale_color_manual(values = L2_color) +
  guides(color = guide_legend(override.aes = list(size = 4),ncol = 3))+
  geom_text(size = 3) +
  theme_classic()

ggsave('dat_umap_level2_label.pdf',p,height = 5.5,width = 6)


#legend

colnams.cluster = 'Level1_label'
centers <- aggregate(meta[, c("umap_1", "umap_2")],list(cluster = meta[, colnams.cluster]),mean)
names(centers) = c(colnams.cluster,"umap_1", "umap_2")
centers[,colnams.cluster] = factor(centers[,colnams.cluster],levels = sort(centers[,colnams.cluster]))

p = ggplot(centers, aes_string("umap_1", "umap_2", label = colnams.cluster)) +
  geom_point(aes_string(color = colnams.cluster), size = 0.2, show.legend = T) +#
  scale_color_manual(values = L1_color) +
  guides(color = guide_legend(override.aes = list(size = 4),ncol = 1))+
  #geom_text(size = 3) +
  theme_classic()

ggsave('dat_umap_level1_legend.pdf',p,height = 8,width = 8)



colnams.cluster = 'Level2_label'
centers <- aggregate(meta[, c("umap_1", "umap_2")],list(cluster = meta[, colnams.cluster]),mean)
names(centers) = c(colnams.cluster,"umap_1", "umap_2")
centers[,colnams.cluster] = factor(centers[,colnams.cluster],levels = sort(centers[,colnams.cluster]))

p = ggplot(centers, aes_string("umap_1", "umap_2", label = colnams.cluster)) +
  geom_point(aes_string(color = colnams.cluster), size = 0.2, show.legend = T) +#
  scale_color_manual(values = L2_color) +
  guides(color = guide_legend(override.aes = list(size = 4),ncol = 3))+
  #geom_text(size = 3) +
  theme_classic()

ggsave('dat_umap_level2_legend.pdf',p,height = 7,width = 15)





# 绘制轮廓线
library(ggunchull)
library(ggrepel)

meta$Major = factor(meta$Major,levels = c("Endothelial","VSMC","Fibroblast","ASPC","Adipocyte","Neuron","Epithelial",
                                        "Myeloid","Neutrophil","B","Platelet","T","NK","Doublets"))

colnams.cluster = 'Major'
main_type_med <- meta %>% 
  dplyr::group_by(across(all_of(colnams.cluster)))%>% 
  summarise(umap_1 = median(umap_1)-1, umap_2 = median(umap_2)-1)

col = c("#A6CEE3","#FF62BC","#66C2A5","#EDB869","#84D7F7","#00B0F6",
        "#D8B6F7","#E78AC3","#FEF6AE","#FB9A99","#7FC28E","#FFD92F",
        "#91CDC1")


meta_tmp = sample_n(meta,0.05*nrow(meta))
p = ggplot(meta_tmp, aes(x = umap_1, y = umap_2)) +
  stat_unchull(aes_string(fill = colnams.cluster, color =colnams.cluster), 
               alpha = 0.05, size = 0.5, lty = 2, delta=0.25, show.legend=F) +
  #geom_point(aes_string(color = colnams.cluster), size = 0.2) +#, show.legend = FALSE
  scale_x_continuous(breaks = NULL) + 
  scale_y_continuous(breaks = NULL) +
  scale_color_manual(values = col) +
  guides(color = guide_legend(override.aes = list(size = 4)))+
  geom_text_repel(aes_string(label=colnams.cluster, color=colnams.cluster),#
                  fontface="bold", data=main_type_med, show.legend=F, size=4) +
  theme(
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(arrow = arrow(type = "closed")),
        axis.title = element_text(hjust = 0.05, size=12))

ggsave('data_umap_Major_outline.pdf',p,height = 5,width = 5.8)



