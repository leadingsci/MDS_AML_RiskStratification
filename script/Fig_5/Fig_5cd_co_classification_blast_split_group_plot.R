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
  output_dir = "/share/home/wyz/Project/result/2025/202506/20250606_3_cibersortx_immune_cell/20250607_3_MDS_AML_tumor_immune_cell_blast_split_group_plot"
  
}

{
  cell_f = "/share/home/wyz/Project/result/2025/202506/20250606_3_cibersortx_immune_cell/20250607_1_MDS_AML_tumor_immune_cell_alluvial_plot_data/20250607_1_MDS_AML_tumor_immune_cell_alluvial_plot_data_inter_4795_sample.xls"
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


{
  ## remove NA
  table(cell_df$BM_Blast_Percent)
  cell_df_keep <- cell_df %>% filter(!is.na(BM_Blast_Percent))
  
  cell_df_keep = cell_df_keep %>% filter(BM_Blast_Percent != ".")
  cell_df_keep = cell_df_keep %>% filter(BM_Blast_Percent != "https://mirrors.e-ducation.cn/CRAN/")
  
  cell_df_keep = cell_df_keep %>% filter(BM_Blast_Percent != "0 (<5%)")
  cell_df_keep = cell_df_keep %>% filter(BM_Blast_Percent != "0.5 (5% - 10%)")
  cell_df_keep = cell_df_keep %>% filter(BM_Blast_Percent != "1.5 (11% -20%)")
  cell_df_keep = cell_df_keep %>% filter(BM_Blast_Percent != "")
  
  table(cell_df_keep$BM_Blast_Percent)
  cell_df_keep$BM_Blast_Percent_1 = as.numeric(cell_df_keep$BM_Blast_Percent) /
    100
  dim(cell_df_keep)
  
}

{
  cell_df_keep
  
  output_f_path = file.path(
    output_dir,
    paste0(
      "20250607_3_MDS_AML_tumor_to_immune_cell_blast_split_group_plot.xls"
    )
  )
  write.table(
    cell_df_keep,
    file = output_f_path,
    sep = "\t",
    row.names = TRUE,
    col.names = TRUE,
    quote = FALSE
  )
}


{
  unique(cell_df_keep$tumor_to_immune_mark)
  
  # >   unique(cell_df_keep$tumor_to_immune_mark)
  # [1] "Primitive_CD8T_Int_CD4T_Low"    "Mature_CD8T_Int_CD4T_Low"       "GMP_CD8_Int_CD4_High"           "GMP_CD8T_Int_CD4T_Low"          "intermediate_CD8T_Int_CD4T_Low" "Mature_CD8T_Low"
  # [7] "Primitive_CD8T_Low"             "Primitive_CD8T_High"            "GMP_CD8T_Low"                   "GMP_CD8T_High"                  "Mature_CD8_Int_CD4_High"        "intermediate_CD8T_Low"
  # [13] "Primitive_CD8_Int_CD4_High"     "intermediate_CD8_Int_CD4_High"  "Prog_CD8T_Low"                  "Mature_CD8T_High"               "Prog_CD8T_Int_CD4T_Low"         "Prog_CD8_Int_CD4_High"
  # [19] "Prog_CD8T_High"                 "intermediate_CD8T_High"
  
  ## 排序
  order_list = c(
    "Prog_CD8T_High",
    "GMP_CD8T_High",
    "Prog_CD8_Int_CD4_High",
    "intermediate_CD8T_High",
    "Mature_CD8T_High",
    "Mature_CD8_Int_CD4_High",
    "Mature_CD8T_Int_CD4T_Low",
    "GMP_CD8T_Int_CD4T_Low" ,
    "intermediate_CD8T_Int_CD4T_Low",
    "GMP_CD8_Int_CD4_High",
    "GMP_CD8T_Low",
    "intermediate_CD8T_Low",
    "intermediate_CD8_Int_CD4_High",
    "Primitive_CD8T_Low",
    "Prog_CD8T_Low",
    "Primitive_CD8T_Int_CD4T_Low",
    "Mature_CD8T_Low",
    "Primitive_CD8_Int_CD4_High",
    "Prog_CD8T_Int_CD4T_Low",
    "Primitive_CD8T_High"
  )
  
  length(order_list)
  
  cell_df_keep$tumor_to_immune_mark = factor(cell_df_keep$tumor_to_immune_mark, levels = order_list)
}

{
  unique(cell_df_keep$immune_to_tumor_mark)
  
  # [1] "CD8T_Int_CD4T_Low_Primitive"    "CD8T_Int_CD4T_Low_Mature"       "CD8_Int_CD4_High_GMP"           "CD8T_Int_CD4T_Low_GMP"          "CD8T_Int_CD4T_Low_intermediate" "CD8T_Low_Mature"
  # [7] "CD8T_Low_Primitive"             "CD8T_High_Primitive"            "CD8T_Low_GMP"                   "CD8T_High_GMP"                  "CD8_Int_CD4_High_Mature"        "CD8T_Low_intermediate"
  # [13] "CD8_Int_CD4_High_Primitive"     "CD8_Int_CD4_High_intermediate"  "CD8T_Low_Prog"                  "CD8T_High_Mature"               "CD8T_Int_CD4T_Low_Prog"         "CD8_Int_CD4_High_Prog"
  # [19] "CD8T_High_Prog"                 "CD8T_High_intermediate"
  
  
  
  ## 排序
  order_list2 = c(
    "CD8T_High_Prog",
    "CD8T_High_GMP",
    "CD8_Int_CD4_High_Prog",
    "CD8T_High_intermediate",
    "CD8T_High_Mature",
    "CD8_Int_CD4_High_Mature",
    "CD8T_Int_CD4T_Low_Mature" ,
    "CD8T_Int_CD4T_Low_GMP" ,
    "CD8T_Int_CD4T_Low_intermediate",
    "CD8_Int_CD4_High_GMP",
    "CD8T_Low_GMP",
    "CD8T_Low_intermediate",
    "CD8_Int_CD4_High_intermediate",
    "CD8T_Low_Primitive",
    "CD8T_Low_Prog",
    "CD8T_Int_CD4T_Low_Primitive",
    "CD8T_Low_Mature",
    "CD8_Int_CD4_High_Primitive",
    "CD8T_Int_CD4T_Low_Prog",
    "CD8T_High_Primitive"
    
  )
  
  length(order_list2)
  
  cell_df_keep$immune_to_tumor_mark = factor(cell_df_keep$immune_to_tumor_mark, levels = order_list2)
}


# -------------------------------------------------------------------------
{
  # CNS ggplot
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


# -------------------------------------------------------------------------
{
  plot_object =
    ggplot(data = cell_df_keep, aes(x = BM_Blast_Percent_1, y = immune_to_tumor_mark)) +
    geom_boxplot(
      width = 0.7,
      fill = "#29abe2",
      colour = "black",
      outlier.shape = NA,
      fatten = 2
    ) +
    
    stat_boxplot(geom = 'errorbar', width = 0.2) +
    labs(x = "BM Blasts(%)", y = NULL, title = "immune to tumor group(blast greater than or equal to 20)") +
    theme(
      legend.position = "right",
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    ) +
    
    theme(axis.text.x = element_text(size = 7),
          axis.text.y = element_text(size = 7)) +
    scale_x_continuous(
      breaks = c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
      labels = scales::percent,
      limits = c(0.2, 1)
      # limits = c(0,0.2)
    ) +
    geom_vline(
      xintercept = c(0, 0.05, 0.1, 0.2),
      linetype = "dashed",
      color = "black",
      linewidth = 0.2
    ) +
    cns.plot.format
  output_f_path = file.path(
    output_dir,
    "20250607_8_MDS_AML_immune_to_tumor_cell_blast_split_group_data_20_100.pdf"
  )
  ggsave(
    filename = output_f_path,
    plot = plot_object,
    device = pdf,
    width = 180,
    height = 120,
    units = "mm"
  )
}