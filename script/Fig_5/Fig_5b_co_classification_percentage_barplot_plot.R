rm(list = ls())



.libPaths()
options(warn = 1)

# Load Package ------------------------------------------------------------

# devtools::install_github("satijalab/seurat", ref = "develop")
{
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
  
}



# input path --------------------------------------------------------------

if (TRUE) {
  ## output
  output_dir = "/share/home/wyz/Project/result/2025/202507/20250708_1_kmean_group_sign_barplot/20250708_8_kmean_group_sign_barplot_tumor_immune"
}

# input path --------------------------------------------------------------

{
  ## experimental group
  #cell_f = "/share/home/wyz/Project/result/2025/202505/20250520_1_cibersortx_tumor_cell_MDS_sAML_AML_merge/20250520_3_tumor_cell_percent_kmean_split_group/20250520_3_MDS_sAML_sAML_tumor_cell_percent_kmean_split_group.xls"
  cell_f = "/share/home/wyz/Project/result/2025/202506/20250606_3_cibersortx_immune_cell/20250607_1_MDS_AML_tumor_immune_cell_alluvial_plot_data/20250607_1_MDS_AML_tumor_immune_cell_alluvial_plot_data_inter_4795_sample.xls"
}

{
  ## read file
  cell_df = read.table(cell_f,
                       sep = "\t",
                       header = T,
                       check.names = FALSE) # 阻止列名中非法符号-号变为.号# row.names = 0,
  print(head(cell_df))
  #print(count_df)
  dim(cell_df) # [1] 5448   59
  
}


{
  colnames(cell_df) <- gsub("-", "_", colnames(cell_df))
  print(colnames(cell_df))
}


# -------------------------------------------------------------------------

{
  ## 选择列
  select_list = c("Sample_ID", "DiseaseClass", "tumor_to_immune_mark")
  
  cell_keep_df = cell_df %>% select(select_list)
  dim(cell_keep_df)
  head(cell_keep_df)
}


{
  table(cell_keep_df$DiseaseClass)
  
}

# -------------------------------------------------------------------------
{
  cell_keep_df = cell_keep_df %>% filter(DiseaseClass %in% c("MDS", "AML"))
  
  ## 固定顺序
  cell_keep_df$DiseaseClass = factor(cell_keep_df$DiseaseClass, levels = rev(c("MDS", "AML")))
  table(cell_keep_df$DiseaseClass)
  
}
{
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
  
  cell_keep_df$tumor_to_immune_mark = factor(cell_keep_df$tumor_to_immune_mark, levels = order_list)
}


{
  # 计算每个 kmean_k5_mark 的比例
  data_summary <- cell_keep_df %>%
    group_by(DiseaseClass, tumor_to_immune_mark) %>%
    summarise(count = n()) %>%
    group_by(DiseaseClass) %>%
    mutate(percentage = count / sum(count) * 100)
  data_summary
  
}


{
  output_f_path = file.path(
    output_dir,
    "20250708_1_MDS_sAML_AML_tumor_immune_cell_clinical_data_simple.xls"
  )
  write.table(
    data_summary,
    file = output_f_path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
  
}

if (FALSE) {
  plot_object = ggbarplot(
    data_summary,
    x = "tumor_to_immune_mark",
    y = "percentage",
    color = "DiseaseClass",
    fill = "DiseaseClass",
    # palette = "npg",
    position = position_dodge(0.8)
  ) +
    scale_fill_manual(values = c("MDS" = "#a3ba87", "AML" = "#6ca7b2")) +  # 手动设置MDS为蓝色，AML为红色
    scale_color_manual(values = c("MDS" = "#a3ba87", "AML" = "#6ca7b2")) + # 手动设置color的颜色
    rotate_x_text(angle = 45) +
    labs(x = NULL, y = "Percentage(%)", title = "Tumor to Immune Group") +  # 添加标签和标题
    # scale_y_continuous(breaks = seq(0, 100, by = 25),limits = c(0, 100)) +# 设置Y轴的刻度为 0, 25, 50, 75, 100
    scale_y_continuous(breaks = seq(0, 30, by = 10), limits = c(0, 30)) +
    # 设置Y轴的刻度为 0, 25, 50, 75, 100
    theme(legend.position = "right") +  # 调整legend位置
    # geom_text(aes(label = round(percentage, 1)), position = position_dodge(0.8), vjust = -0.5,size.unit = "pt",size=5)  # 添加柱子上的Y轴数值标签
    geom_text(
      aes(label = round(percentage, 1)),
      position = position_dodge(0.8),
      vjust = 0.5,
      size.unit = "pt",
      size = 5
    )  # 添加柱子上的Y轴数值标签
  
  
  output_f_path = file.path(
    output_dir,
    "20250708_2_MDS_sAML_AML_tumor_immune_cell_clinical_data_simple_w220mm_h140.pdf"
  )
  ggsave(
    filename = output_f_path,
    plot = plot_object,
    device = pdf,
    width = 220,
    height = 140,
    units = "mm"
  )
}


{
  plot_object = ggbarplot(
    data_summary,
    x = "tumor_to_immune_mark",
    y = "percentage",
    color = "DiseaseClass",
    fill = "DiseaseClass",
    # palette = "npg",
    position = position_dodge(0.8)
  ) +
    scale_fill_manual(values = c("MDS" = "#a3ba87", "AML" = "#6ca7b2")) +  # 手动设置MDS为蓝色，AML为红色
    scale_color_manual(values = c("MDS" = "#a3ba87", "AML" = "#6ca7b2")) + # 手动设置color的颜色
    # rotate_x_text(angle = 45) +
    coord_flip() +  # 旋转坐标轴，使柱子沿Y轴显示
    labs(x = NULL, y = "Percentage(%)", title = "Tumor to Immune Group") +  # 添加标签和标题
    # scale_y_continuous(breaks = seq(0, 100, by = 25),limits = c(0, 100)) +# 设置Y轴的刻度为 0, 25, 50, 75, 100
    scale_y_continuous(breaks = seq(0, 30, by = 10), limits = c(0, 30)) +
    # 设置Y轴的刻度为 0, 25, 50, 75, 100
    theme(legend.position = "right") +  # 调整legend位置
    # geom_text(aes(label = round(percentage, 1)), position = position_dodge(0.8), vjust = -0.5,size.unit = "pt",size=5)  # 添加柱子上的Y轴数值标签
    geom_text(
      aes(label = round(percentage, 1)),
      position = position_dodge(0.8),
      vjust = 0.5,
      size.unit = "pt",
      size = 5
    )  # 添加柱子上的Y轴数值标签
  
  
  output_f_path = file.path(
    output_dir,
    "20250708_1_2_MDS_sAML_AML_tumor_immune_cell_clinical_data_simple_w230mm_h140.pdf"
  )
  ggsave(
    filename = output_f_path,
    plot = plot_object,
    device = pdf,
    width = 230,
    height = 140,
    units = "mm"
  )
}