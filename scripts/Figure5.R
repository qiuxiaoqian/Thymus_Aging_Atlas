setwd('D:/OneDrive/Work/Thymus/Manuscript/Figure5')
#setwd('C:/Users/qxq/OneDrive/Work/Thymus/Manuscript/Figure5')

rm(list=ls())
library(dplyr)
library(reshape2)
library(stringr)
library(Seurat)
library(openxlsx)
library(ggplot2)

res_path = 'D:/OneDrive/Work/Thymus/ScSn_integration/Analysis/Macro/DEG'
compare_group = c('c26',"c25")

for (compare_group in list(c('c26',"c25"),c('c27',"c25"),c('c27',"c26"))){
  
  
  marker_res= read.xlsx(str_glue('{res_path}/{compare_group[1]}vs.{compare_group[2]}.xlsx'))
  
  logfc_col = str_glue("{compare_group[1]}/{compare_group[2]}.avg_log2FC")
  padj_col  = str_glue("{compare_group[1]}/{compare_group[2]}.p_val_adj")
  
  
  df <- marker_res %>%
    mutate(
      logfc = .data[[logfc_col]],
      padj  = .data[[padj_col]],
      log10_p = -log10(padj),
      sig = case_when(
        padj < 0.05 & logfc > 0.25 ~ "Up",
        padj < 0.05 & logfc < -0.25 ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  
  # ---- 抽样 NS 中 logFC 在 [-5,5] 的一半 ----
  df_NS_sampled <- df %>%
    filter(sig == "NS", logfc > -5, logfc < 5) %>%
    slice_sample(prop = 0.2)
  
  df_other <- df %>% filter(!(sig == "NS" & logfc > -5 & logfc < 5))
  df_plot <- bind_rows(df_NS_sampled, df_other)
  
  
  
  # ---- 分别取上/下调 top10 ----
  top_labels_up <- df %>%
    filter(sig == "Up") %>%
    arrange(padj) %>%
    slice(1:10)
  
  top_labels_down <- df %>%
    filter(sig == "Down") %>%
    arrange(padj) %>%
    slice(1:10)
  
  top_labels_real <- bind_rows(top_labels_up, top_labels_down) %>%
    distinct(gene, .keep_all = TRUE)
  
  
  #一定标记APOE,CD36,PDGFRA
  force_genes <- c("APOE", "CD36", "PDGFRA",'FGF13',"CCL2",'CCL3','CCL4','CCL8')
  
  
  #加入df_plot
  df_plot <- bind_rows(df_plot, df %>% filter(gene %in% force_genes)) %>% distinct(gene, .keep_all = TRUE)
  df_plot$sig = factor(df_plot$sig,levels = c("Up", "Down", "NS"))
  
  
  # 加入label
  # force_genes 中，但不在 top10 的
  force_only_labels <- df %>%
    filter(gene %in% force_genes) %>%
    filter(!(gene %in% top_labels_real$gene))
  
  # 黑圈点（只针对 force_genes）
  highlight <- df %>% filter(gene %in% force_genes)
  

  
  library(ggrepel)
  p = ggplot(df_plot, aes(x = logfc, y = log10_p, color = sig)) +
    geom_point(size = 1.2, alpha = 0.7,shape = 16) +
    geom_point(data = highlight,aes(x = logfc, y = log10_p),size = 1.2,
               color = "black",stroke = 1,shape = 1,show.legend = FALSE)+
    # 🔴🔵 Top genes：按 sig 着色
    geom_text_repel(data = top_labels_real, aes(label = gene, color = sig),
                    size = 3, max.overlaps = 50,show.legend = FALSE) +
    # ⚫ force_genes 但非 top：黑色字体
    geom_text_repel(data = force_only_labels,aes(label = gene),
                    color = "black",size = 3,max.overlaps = 50,show.legend = FALSE) +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey")) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed",size = 0.2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",size = 0.2) +
    guides(color = guide_legend(override.aes = list(size = 3)))+
    labs(
      x = str_glue("log2 Fold Change ({compare_group[1]} vs {compare_group[2]})"),
      y = "-log10(Adjusted P-value)",
      color = "Significance",
      title = str_glue("Volcano Plot: {compare_group[1]} vs {compare_group[2]}")
    ) +
    theme_classic(base_size = 12)+
    theme(axis.line = element_line(size = 0.5))
  
  ggsave(str_glue('{compare_group[1]}_vs_{compare_group[2]}.pdf'), p, width = 5.5, height = 3.5)
  
  
  
}



#FGF13的表达
dat = readRDS('D:/OneDrive/Work/Thymus/ScSn_integration/anno2/Thymus_scsn_anno2_finalname.rds')
dat = subset(dat,Level1 != 'Doublets')
dat = subset(dat,Platform =='snRNA')
#data = readRDS('D:/OneDrive/Work/Thymus/ScSn_integration/Analysis/Macro/data_snRNA_c10_Macro.rds')

L1_color_used = c("#A6CEE3", "#1F78B4", "#EDB869", "#FB9A99", "#B2DF8A", "#84D7F7", 
                  "#D8B6F7", "#E3EDC0", "#E78AC3", "#D4EEF6", "#FC8D62", "#8DA0CB", 
                  "#66C2A5", "#6F89D3", "#5050FF", "#F9C9BD", "#F0E685", "#33A02C", 
                  "#F07590", "#4DBBD5", "#FFD92F", "#E31A1C", "#E16CB7", "#B3E2CD", 
                  "#EEBED6", "#7FC28E", "#5DB4F3", "#FFF2AE", "#B589BC", "#CD7560", 
                  "#FDBF6F", "#DCDCDC") %>% unique()


celltype_metadata = read.xlsx('D:/OneDrive/Work/Thymus/ScSn_integration/data/celltype_metadata.xlsx')
dat$Level1 = factor(dat$Level1,levels = unique(celltype_metadata$Level1))

gene_list = c('CCL2','CCL3','CCL4','CCL8','FGF13','ANGPTL4')
for (gene in gene_list){
  png(str_glue('vlnplot/dat_Level1_{gene}_vlnplot.png'),res = 300,height = 5,width = 10,units = 'in')
  p = VlnPlot(dat,features = gene,group.by  = "Level1",cols = L1_color_used)+
    labs(x = 'Cell type')+
    theme(legend.position = 'none')
  print(p)
  dev.off()
  
  pdf(str_glue('vlnplot/dat_Level1_{gene}_vlnplot.pdf'),height = 5,width = 10)
  p = VlnPlot(dat,features = gene,group.by  = "Level1",cols = L1_color_used)+
    labs(x = 'Cell type')+
    theme(legend.position = 'none')
  print(p)
  dev.off()
  
}


#dotplot
celltype_metadata = read.xlsx('D:/OneDrive/Work/Thymus/ScSn_integration/data/celltype_metadata.xlsx')
dat$Level1 = factor(dat$Level1,levels = unique(celltype_metadata$Level1))

gene_list = c('CCL2','CCL3','CCL4','CCL8')#,'FGF13','ANGPTL4'

pdf(str_glue('dat_CCL_dotplot_Level1.pdf'),height = 3,width = 10)
DotPlot(dat,features = gene_list,scale = T,assay = 'RNA',group.by  = "Level1")+
  theme(panel.border = element_rect(colour = "black", fill=NA,size = 1),
        axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
        axis.text.y = element_text(face = "italic"))+
  coord_flip()+labs(x = 'Gene',y = '')
dev.off()


pdf(str_glue('dat_CCL_dotplot_Level2.pdf'),height = 5,width = 16)
DotPlot(dat,features = gene_list,scale = T,assay = 'RNA',group.by  = "Level2_label")+
  theme(panel.border = element_rect(colour = "black", fill=NA,size = 1),
        axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),
        axis.text.y = element_text(,face = "italic"))+
  coord_flip()+labs(x = 'Gene',y = '')
dev.off()


