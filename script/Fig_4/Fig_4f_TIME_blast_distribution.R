rm(list = ls())

# libPaths --------------------------------------------------------------------
# .libPaths(c(
#   "/share/home/bioinfo/anaconda3/envs/R-4.2.3/lib/R/library"
# ))


.libPaths()
options(warn = 1)

# Load Package ------------------------------------------------------------

# devtools::install_github("satijalab/seurat", ref = "develop")


library(readr) #R 4.2.2 #R 4.2.3
library(ggplot2)
#R 4.2.2 #R 4.2.3
library(tidyverse)
#R 4.2.2
library(openxlsx)
#R 4.2.3
library(readxl)
#R 4.2.3
library(dplyr)
#R 4.2.2 #R 4.2.3

library(pheatmap) #R 4.2.2 #R 4.2.3
library(ggrastr) #R 4.2.3
library(ggpubr) #R 4.2.2 #R 4.2.3

{
  library(survminer)
  library(survival)
}


# input path --------------------------------------------------------------

if (TRUE) {
  ## output
  output_dir = "/share/home/wyz/Project/result/2025/202505/20250527_1_cibersortx_immune_cell_MDS_sAML_AML_merge/20250528_2_blast_barplot_MDS_AML_blast_map_kmean_group"
  
}


# input path --------------------------------------------------------------
{
  cell_f = "/share/home/wyz/Project/result/2025/202505/20250527_1_cibersortx_immune_cell_MDS_sAML_AML_merge/20250528_1_blast_barplot_MDS_sAML_AML_blast_map_kmean_group/20250528_1_blast_barplot_MDS_sAML_AML_like_map_kmean_group_2239_sample_add_BM_Blass_416_keep.xls"
}

{
  ## read file
  cell_df = read.table(cell_f,
                       sep = "\t",
                       header = T,
                       check.names = FALSE) # 阻止列名中非法符号-号变为.号 # row.names = 1,
  print(head(cell_df))
  dim(cell_df) #
  
}

{
  colnames(cell_df) <- gsub("-", "_", colnames(cell_df))
  print(colnames(cell_df))
}



# -------------------------------------------------------------------------


# -------------------------------------------------------------------------

{
  cell_df$Group = cell_df$kmean_k4_mark
  cell_df$Final_Blast_Percent = (cell_df$BM_Blast_Percent) / 100
  
  table(cell_df$Group)
  cell_df$Group
  
  cell_df$Group = factor(
    cell_df$Group ,
    levels = c(
      "k4_3_CD8T_Low",
      "k4_4_CD8_Int_CD4_High",
      "k4_2_CD8T_High",
      "k4_1_CD8T_Int_CD4T_Low"
    )
  )
  table(cell_df$Group)
}




# -------------------------------------------------------------------------

{
  # Kruskal-Wallis
  kruskal_test_result <- kruskal.test(Final_Blast_Percent ~ Group, data = cell_df)
  
  kruskal_test_result$p.value
}



# -------------------------------------------------------------------------
{
  cns.plot.format = theme(
    plot.background = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    
    axis.text = element_text(color = "black", size = 7),
    axis.title = element_text(color = "black", size = 7),
    plot.title = element_text(color = "black", size = 7),
    
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(color = "black", size = 7),
    legend.title = element_text(color = "black", size = 7)
  )
}


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
{
  #
  plot_object =
    ggplot(data = cell_df, aes(x = Final_Blast_Percent, y = Group)) +
    geom_boxplot(
      width = 0.7,
      fill = "#29abe2",
      colour = "black",
      outlier.shape = NA,
      fatten = 2
    ) +
    
    stat_boxplot(geom = 'errorbar', width = 0.2) +
    
    labs(x = "BM Blasts(%)", y = NULL, title = "MDS,sAML,AML") +
    theme(
      legend.position = "right",
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    ) +
    
    theme(axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7)) +
    scale_x_continuous(
      breaks = c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
      labels = scales::percent
    ) +
    geom_vline(
      xintercept = c(0, 0.05, 0.1, 0.2),
      linetype = "dashed",
      #
      color = "black",
      #
      linewidth = 0.2          #
    ) +
    
    cns.plot.format
  
  output_f_path = file.path(output_dir,
                            "20250528_2_blast_barplot_MDS_sAML_AML_blast_map_kmean_group.pdf")
  ggsave(
    filename = output_f_path,
    plot = plot_object,
    device = pdf,
    width = 180,
    height = 100,
    units = "mm"
  )
}