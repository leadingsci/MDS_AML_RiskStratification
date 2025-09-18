rm(list = ls())

.libPaths()
options(warn = 1)

# Load Package ------------------------------------------------------------

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
  library(survminer) # 加载包
  library(survival) # 加载包
}


# input path --------------------------------------------------------------

if (TRUE) {
  ## output
  output_dir = "/share/home/wyz/Project/result/2025/202506/20250605_20250520_1_cibersortx_tumor_cell_kmean_5448_MDS_sAML_AML_merge/20250605_2_blast_barplot_MDS_sAML_AML_like_map_kmean_group_split_blast"
}


# input path --------------------------------------------------------------


{
  ## experimental group
  cell_f = "/share/home/wyz/Project/result/2025/202506/20250605_20250520_1_cibersortx_tumor_cell_kmean_5448_MDS_sAML_AML_merge/20250605_2_blast_barplot_MDS_sAML_AML_like_map_kmean_group_split_blast/20250522_1_blast_5448_sample_add_BM_Blass_540_keep_mark.xlsx"
}


{
  cell_df = read_excel(cell_f)
  cell_df
}

{
  colnames(cell_df) <- gsub("-", "_", colnames(cell_df))
  print(colnames(cell_df))
}



# -------------------------------------------------------------------------


# -------------------------------------------------------------------------

{
  cell_df$Group = cell_df$kmean_k5_mark
  cell_df$Final_Blast_Percent = as.numeric(cell_df$BM_Blast_Percent) / 100
  
  table(cell_df$Group)
  table(cell_df$Final_Blast_Percent)
  cell_df$Group
  table(cell_df$kmean_k5_mark)
  
  cell_df$Group = factor(
    cell_df$Group ,
    levels = c(
      "k5_4_GMP",
      "k5_1_Mature",
      "k5_3_Prog",
      "k5_2_Primitive",
      "k5_5_intermediate"
    )
  )
  table(cell_df$Group)
}

{
  table(cell_df$BM_Blast_Percent_mark)
  cell_df$BM_Blast_Percent_mark  = factor(cell_df$BM_Blast_Percent_mark, levels = c("BM<20", "BM>=20"))
}




# -------------------------------------------------------------------------

{
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
  plot_object =
    ggplot(data = cell_df, aes(x = Final_Blast_Percent, y = Group)) +
    geom_boxplot(
      width = 0.7,
      fill = "#29abe2",
      #
      colour = "black",
      #
      outlier.shape = NA,
      #
      fatten = 2               #
    ) +
    
    stat_boxplot(geom = 'errorbar', width = 0.2) +
    labs(x = "BM Blasts(%)", y = NULL, title = "MDS,sAML,AML Kruskal-Wallis (P=3.569215e-30)") +
    theme(
      legend.position = "right",
      axis.line = element_line(color = "black"),
      #
      axis.ticks = element_line(color = "black")          #
    ) +
    
    theme(axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7)) +
    scale_x_continuous(
      breaks = c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
      labels = scales::percent,
      limits = c(0.1, 0.2)              # Limit the range of the X-axis
      # limits = c(0.2, 1)              # Limit the range of the X-axis
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
                            "20250704_1_blast_barplot_MDS_sAML_AML_like_map_kmean_group.pdf")
  ggsave(
    filename = output_f_path,
    plot = plot_object,
    device = pdf,
    width = 180,
    height = 100,
    units = "mm"
  )
  
}