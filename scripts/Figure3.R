setwd('D:/OneDrive/Work/Thymus/Manuscript/Figure3/Plot')
#setwd('C:/Users/qxq/OneDrive/Work/Thymus/ScSn_integration/Figure/Figure3')

rm(list=ls())
library(dplyr)
library(reshape2)
library(stringr)
library(Seurat)
library(openxlsx)
library(ggplot2)

cell = 'snRNA_c10_Macro'
data = readRDS(str_glue('../../../ScSn_integration/Analysis/Macro/data_{cell}.rds'))

pdf(str_glue('data_umap_level2.pdf'),height = 3,width = 7)
DimPlot(data,cols = c("#D4EEF6", "#EEBED6", "#CD7560"),raster = F,group.by = 'Level2_label')+
  labs(title = '')+NoAxes()+
  guides(color = guide_legend(override.aes = list(size = 4)))
dev.off()


gene = c( #M1
          "IL1B", "TNF", "IL6", "CXCL10", "CXCL9","NOS2", "CD80","FCGR1A", 
          "S100A8", "S100A9", "GBP1", "GBP2", "GBP5",
          "CD86", "HLA-DRA", "HLA-DRB1","STAT1", "IRF5", "NFKBIA", "NLRP3", "PTGS2",
          #M2
          "CD163", "MRC1",  "MSR1","STAB1", "MERTK","CCL18", "TREM2","IL10", "ARG1",
          #亚群
          "MRC1", "CCL17", "CCL22", "IL1RN", "CLEC10A", "MSR1",#M2a（Wound Healing）
          "IL10", "CD86", "TNFSF14", "C1QA", "C1QB", "C1QC",#M2b（Immune-Regulatory）
          "CD163", "MERTK", "STAB1", "LRRC32",#M2c（Anti-inflammatory / Clearance）
          "VEGFA", "HIF1A", "PTX3",#M2d（Angiogenic）
          "SPP1", "CHI3L1", "ITGAX", "TREM2","LGALS3", "FABP5", "MMP9",#SPP1+ Macrophage（肿瘤/纤维化热点）
          "TREM2", "APOE", "APOC1", "LPL","GPNMB", "SLC40A1", "CD9", "CST7", #TREM2+ / Lipid-associated Macrophage (LAM/DAM)
          "IL1B", "TNF", "CXCL8", "CXCL10", "CXCL2","SOD2", "NLRP3",#Inflammatory Macrophage
          "C1QA", "C1QB", "C1QC", "LYVE1","MERTK", "SEPP1",#C1Q+ Tissue-resident Macrophage
          "LYVE1", "F13A1", "MRC1", "STAB1",#LYVE1+ Perivascular Macrophage
          "MKI67", "TOP2A", "BIRC5", "PCNA", "STMN1",#Poliferating Macrophage
          "ISG15", "IFIT1", "IFIT3", "MX1", "MX2","CXCL10", "GBP1", "GBP2", "GBP5"#IFN-response Macrophage
          
           ) %>% unique()

gene_set = 'All'
pdf(str_glue('data_{cell}_{gene_set}_dotplot.pdf'),height = 12,width = 4.5)
DotPlot(data,features = gene,scale = T,assay = 'RNA')+
  theme(panel.border = element_rect(colour = "black", fill=NA,size = 1))+
  coord_flip()+
  labs(x = 'Gene',y = 'Cell cluster')+
  theme(axis.text.x = element_text(angle = 45,vjust =1,hjust = 1))
dev.off()



gene = c( #M1
  "CD86", "HLA-DRA", "HLA-DRB1","STAT1", "IRF5", "NFKBIA", "NLRP3", "PTGS2",
  #M2
  "CD163", "MRC1",  "MSR1","STAB1", "MERTK","CCL18",
  #亚群
  "C1QA", "C1QB", "C1QC",#M2b（Immune-Regulatory）
  "SPP1", "MMP9", "CHI3L1", "ITGAX", "TREM2","LGALS3", "FABP5",#SPP1+ Macrophage（肿瘤/纤维化热点）
  "TREM2", "APOE", "APOC1", "LPL","GPNMB", "SLC40A1",  #TREM2+ / Lipid-associated Macrophage (LAM/DAM)
  "IL1B", "TNF", "CXCL8", "SOD2", "NLRP3",#Inflammatory Macrophage
  "LYVE1", "F13A1", "MRC1", "STAB1"#LYVE1+ Perivascular Macrophage

) %>% unique()



gene = c( "CD68", "HLA-DRA", "HLA-DRB1","STAT1", 
          "CD163", "MRC1",  "MSR1","STAB1", "CCL18",
          "C1QA", "C1QB", "C1QC",#M2b（Immune-Regulatory）
          "SELENOP","APOE","IL18","CCL18",
          
          "CD36","LYVE1","CCL2","CCL3","CCL4","CCL8",'FGF13',
          "PDGFRA","COL1A1","COL15A1","ADIPOQ","PLIN1",
          "SPP1", "MMP9", 
          "TREM2" #TREM2+ / Lipid-associated Macrophage (LAM/DAM)
) %>% unique()


data$L2_label = factor(data$L2_label,levels = rev(sort(unique(data$L2_label))))
pdf(str_glue('data_{cell}_dotplot.pdf'),height = 3,width = 7)
DotPlot(data,features = gene,scale = T,assay = 'RNA',group.by = 'L2_label')+
  theme(panel.border = element_rect(colour = "black", fill=NA,size = 1))+
  #coord_flip()+
  labs(x = 'Gene',y = 'Cell cluster')+
  theme(axis.text.x = element_text(angle = 90,vjust =0.5,hjust = 1,face = "italic"),
        axis.line = element_blank(),
        legend.direction = "horizontal",
        legend.position = "bottom",
        panel.background = element_rect(fill = NA,colour = 'black',size = 0.5))
dev.off()


ncol = 6
png(str_glue('data_{cell}_feature.png'),res = 300,
    height = ceiling(length(gene)/ncol)*2.5,width = ncol*3,units = 'in')
FeaturePlot(data,features = gene,ncol = ncol)
dev.off()



#绘图Feature
gene = c( "CD68","CD163","APOE","IL18","CD36",
          "LYVE1","PDGFRA","COL15A1","ADIPOQ","PLIN1"
) %>% unique()

ncol = 5
# png(str_glue('data_{cell}_feature_selected.png'),res = 300,
#     height = ceiling(length(gene)/ncol)*2.5,width = ncol*3,units = 'in')
# FeaturePlot(data,features = gene,ncol = ncol)+
#   theme(plot.title = element_text(face = "italic"))#只斜体了最后一个
# dev.off()



# 基因斜体
plots=list()
for (i in 1:length(gene)){
  plots[[i]]=FeaturePlot(data,features = gene[i])+NoAxes()+
  theme(plot.title = element_text(face = "italic"),
        #panel.border = element_rect(fill = NA,color = "black",size=1.5,linetype = "solid"),
        #legend.position = 'none'
        )}
  
p <- wrap_plots(plots, ncol = 5)

png(str_glue('data_{cell}_feature_selected.png'),res = 300,
    height = ceiling(length(gene)/ncol)*2.5,width = ncol*3,units = 'in')
print(p)
dev.off()


pdf(str_glue('data_{cell}_feature_selected.pdf'),
    height = ceiling(length(gene)/ncol)*2.5,width = ncol*3)
print(p)
dev.off()
