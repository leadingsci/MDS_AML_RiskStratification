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


# input path --------------------------------------------------------------

if (TRUE) {
  ## output
  output_dir = "/share/home/wyz/Project/result/2025/202506/20250606_3_cibersortx_immune_cell/20250606_1_scRNA_subtype_MDS_sAML_AML_immune_percent_survival_barplot_k4"
}

# input path --------------------------------------------------------------

{
  cell_f = "/share/home/wyz/Project/result/2025/202505/20250527_1_cibersortx_immune_cell_MDS_sAML_AML_merge/20250527_4_immune_cell_percent_kmean_split_group_del_GSE146173/20250527_3_immune_cell_percent_kmean_split_group_kmean_4798_sample.xls"
}

{
  ## read file
  cell_df = read.table(cell_f,
                       sep = "\t",
                       header = T,
                       check.names = FALSE) # 阻止列名中非法符号-号变为.号# row.names = 0,
  print(head(cell_df))
  dim(cell_df)
  
}

{
  colnames(cell_df) <- gsub("-", "_", colnames(cell_df))
  print(colnames(cell_df))
}



# -------------------------------------------------------------------------
if (TRUE) {
  library(dplyr)
  cell_df_T_sum = cell_df %>% mutate(T_sum = CD4_T_Percent  + CD8_T_Percent)
  head(cell_df_T_sum)
}



{
  cell_list = c("CD4_T_Percent", "CD8_T_Percent", "NK_Percent", "B_Percent")
  
  library(tidyr)
  
  
  cell_df_long = cell_df_T_sum %>%
    pivot_longer(cols = cell_list,
                 names_to = "Cell_Type",
                 values_to = "Cell_Value")
  
  cell_df_long
  print(colnames(cell_df_long))
  cell_df_long[, c("Sample_ID", "kmean_k4", "Cell_Type", "Cell_Value", "T_sum")]
}


{
  cell_df_long = cell_df_long %>%
    mutate(Cell_Type = str_replace(Cell_Type, "_Percent$", ""))
  
  cell_df_long[, c("Sample_ID", "kmean_k4", "Cell_Type", "Cell_Value")]
}


# -------------------------------------------------------------------------
{
  cell_list = c("CD4_T", "CD8_T", "NK", "B")
  
  cell_df_long$Cell_Type = factor(cell_df_long$Cell_Type, levels = cell_list)
  print(cell_df_long$Cell_Type)
}

# -------------------------------------------------------------------------

{
  # save
  output_f_path = file.path(
    output_dir,
    "20250527_1_MDS_AML_immune_percent_survival_barplot_4798_sample.xls"
  )
  write.table(
    cell_df_long,
    file = output_f_path,
    sep = "\t",
    row.names = F,
    quote = FALSE
  )
  
}

table(cell_df_long$kmean_k4)



# -------------------------------------------------------------------------
# k4_1
# -------------------------------------------------------------------------
{
  cell_select_df = cell_df_long[cell_df_long$kmean_k4 == "1", ]
}


# -------------------------------------------------------------------------

{
  library(dplyr)
  cell_select_df_sorted <- cell_select_df %>%
    arrange(desc(T_sum))
  dim(cell_select_df_sorted)
}

{
  ## factor
  cell_select_df_sorted$Sample_ID = factor(cell_select_df_sorted$Sample_ID,
                                           levels = unique(cell_select_df_sorted$Sample_ID))
}



# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
{
  color_list = c(
    "CD4_T" = "#ca8282",
    "CD8_T" = "#e19793",
    "NK" = "#66abda",
    "B" = "#f5d966"
  )
}



# -------------------------------------------------------------------------

{
  group_name = "k4_1"
  sample_num = 971
  
  plot_object =
    ggplot(cell_select_df_sorted,
           aes(x = Sample_ID, y = Cell_Value, fill = Cell_Type)) +
    geom_bar(
      stat = "identity",
      alpha = 1,
      width = 0.8,
      color = NA,
      size = 0
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), label = scales::percent) +
    scale_fill_manual(values = color_list) +
    labs(
      x = NULL,
      y = "Percentage",
      title = paste0(group_name, " (n=", sample_num, ")")
    ) +
    theme_bw() +
    theme(axis.title = element_text(size = 7)) +
    theme(
      plot.title = element_text(size = 7),
      axis.text.x = element_blank(),
      # axis.text.x = element_text(size = 3,angle=45,hjust=1,vjust=1),
      axis.text.y = element_text(size = 7)
    ) +
    theme(plot.margin = margin(1, 1, 1, 1.2, "cm")) +
    theme(panel.grid = element_blank()) +
    theme(
      legend.position = "top",
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 6),
      
      legend.key.size = unit(2, "mm"),
    )
  
  
  
  plot_object
  
  output_f_path = file.path(output_dir,
                            paste0("20250704_11_1_", group_name, "_barplot.pdf"))
  ggsave(
    filename = output_f_path,
    plot = plot_object,
    device = pdf,
    width = 0.1 * sample_num,
    height = 120,
    units = "mm"
  )
  
}



# -------------------------------------------------------------------------
# k4_2
# -------------------------------------------------------------------------
{
  cell_select_df = cell_df_long[cell_df_long$kmean_k4 == "2", ]
  
  
}


# -------------------------------------------------------------------------

{
  library(dplyr)
  cell_select_df_sorted <- cell_select_df %>%
    arrange(desc(T_sum))
  
  
  
  # print(cell_select_df_sorted$HSC_Prog_sum)
  dim(cell_select_df_sorted)
}

{
  ## 固定X轴的顺序
  cell_select_df_sorted$Sample_ID = factor(cell_select_df_sorted$Sample_ID,
                                           levels = unique(cell_select_df_sorted$Sample_ID))
  
}



# -------------------------------------------------------------------------

# 画图
# -------------------------------------------------------------------------
{
  # 配色
  color_list = c(
    "HSC_like" = "#d5211d",
    "Prog_like" = "#ed7a1c",
    "GMP_like" = "#3976ab",
    "Mono_like" = "#4da54a",
    "cDC_like" = "#8e4b97"
  )
  
}


# -------------------------------------------------------------------------

{
  group_name = "k4_2"
  sample_num = 1218
  plot_object =
    ggplot(cell_select_df_sorted,
           aes(x = Sample_ID, y = Cell_Value, fill = Cell_Type)) +
    geom_bar(
      stat = "identity",
      alpha = 1,
      width = 0.8,
      color = NA,
      size = 0
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), label = scales::percent) +
    scale_fill_manual(values = color_list) +
    labs(
      x = NULL,
      y = "Percentage",
      title = paste0(group_name, " (n=", sample_num, ")")
    ) +
    theme_bw() +
    theme(axis.title = element_text(size = 7)) +
    theme(
      plot.title = element_text(size = 7),
      axis.text.x = element_blank(),
      # axis.text.x = element_text(size = 3,angle=45,hjust=1,vjust=1),
      axis.text.y = element_text(size = 7)
    ) +
    theme(plot.margin = margin(1, 1, 1, 1.2, "cm")) +
    theme(panel.grid = element_blank()) +
    theme(
      legend.position = "top",
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 6),
      legend.key.size = unit(2, "mm"),
    )
  
  plot_object
  
  output_f_path = file.path(output_dir,
                            paste0("20250704_12_2_", group_name, "_barplot.pdf"))
  ggsave(
    filename = output_f_path,
    plot = plot_object,
    device = pdf,
    width = 0.1 * sample_num,
    height = 120,
    units = "mm"
  )
}



# -------------------------------------------------------------------------
# k4_3
# -------------------------------------------------------------------------
{
  cell_select_df = cell_df_long[cell_df_long$kmean_k4 == "3", ]
  
}


# -------------------------------------------------------------------------

{
  library(dplyr)
  cell_select_df_sorted <- cell_select_df %>%
    arrange(desc(T_sum))
  print(cell_select_df_sorted$HSC_Prog_sum)
  dim(cell_select_df_sorted)
}

{
  ## factor
  cell_select_df_sorted$Sample_ID = factor(cell_select_df_sorted$Sample_ID,
                                           levels = unique(cell_select_df_sorted$Sample_ID))
  
}



# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
{
  color_list = c(
    "HSC_like" = "#d5211d",
    "Prog_like" = "#ed7a1c",
    "GMP_like" = "#3976ab",
    "Mono_like" = "#4da54a",
    "cDC_like" = "#8e4b97"
  )
  
}


# -------------------------------------------------------------------------

{
  group_name = "k4_3"
  sample_num = 1895
  plot_object =
    ggplot(cell_select_df_sorted,
           aes(x = Sample_ID, y = Cell_Value, fill = Cell_Type)) +
    geom_bar(
      stat = "identity",
      alpha = 1,
      width = 0.8,
      color = NA,
      size = 0
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), label = scales::percent) +
    scale_fill_manual(values = color_list) +
    labs(
      x = NULL,
      y = "Percentage",
      title = paste0(group_name, " (n=", sample_num, ")")
    ) +
    theme_bw() +
    theme(axis.title = element_text(size = 7)) +
    theme(
      plot.title = element_text(size = 7),
      axis.text.x = element_blank(),
      # axis.text.x = element_text(size = 3,angle=45,hjust=1,vjust=1),
      axis.text.y = element_text(size = 7)
    ) +
    theme(plot.margin = margin(1, 1, 1, 1.2, "cm")) +
    theme(panel.grid = element_blank()) +
    theme(
      legend.position = "top",
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 6),
      
      legend.key.size = unit(2, "mm"),
    )
  
  plot_object
  
  output_f_path = file.path(output_dir,
                            paste0("20250704_13_2_", group_name, "_barplot.pdf"))
  ggsave(
    filename = output_f_path,
    plot = plot_object,
    device = pdf,
    width = 0.1 * sample_num,
    height = 120,
    units = "mm"
  )
  
}



# -------------------------------------------------------------------------
# k4_4
# -------------------------------------------------------------------------
{
  cell_select_df = cell_df_long[cell_df_long$kmean_k4 == "4", ]
}


# -------------------------------------------------------------------------

{
  library(dplyr)
  cell_select_df_sorted <- cell_select_df %>%
    arrange(desc(T_sum))
  
  print(cell_select_df_sorted$T_sum)
  dim(cell_select_df_sorted)
}

{
  ## 固定X轴的顺序
  cell_select_df_sorted$Sample_ID = factor(cell_select_df_sorted$Sample_ID,
                                           levels = unique(cell_select_df_sorted$Sample_ID))
  
}



# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
{
  color_list = c(
    "HSC_like" = "#d5211d",
    "Prog_like" = "#ed7a1c",
    "GMP_like" = "#3976ab",
    "Mono_like" = "#4da54a",
    "cDC_like" = "#8e4b97"
  )
  
}


# -------------------------------------------------------------------------
{
  group_name = "k4_4"
  sample_num = 714
  plot_object =
    ggplot(cell_select_df_sorted,
           aes(x = Sample_ID, y = Cell_Value, fill = Cell_Type)) +
    geom_bar(
      stat = "identity",
      alpha = 1,
      width = 0.8,
      color = NA,
      size = 0
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), label = scales::percent) +
    scale_fill_manual(values = color_list) +
    labs(
      x = NULL,
      y = "Percentage",
      title = paste0(group_name, " (n=", sample_num, ")")
    ) +
    theme_bw() +
    theme(axis.title = element_text(size = 7)) +
    theme(
      plot.title = element_text(size = 7),
      axis.text.x = element_blank(),
      # axis.text.x = element_text(size = 3,angle=45,hjust=1,vjust=1),
      axis.text.y = element_text(size = 7)
    ) +
    theme(plot.margin = margin(1, 1, 1, 1.2, "cm")) +
    theme(panel.grid = element_blank()) +
    theme(
      legend.position = "top",
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 6),
      
      legend.key.size = unit(2, "mm"),
    )
  
  plot_object
  
  output_f_path = file.path(output_dir,
                            paste0("20250704_14_2_", group_name, "_barplot.pdf"))
  ggsave(
    filename = output_f_path,
    plot = plot_object,
    device = pdf,
    width = 0.1 * sample_num,
    height = 120,
    units = "mm"
  )
}