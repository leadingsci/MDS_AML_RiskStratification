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
  library(survminer) # 加载包
  library(survival) # 加载包
}


# input path --------------------------------------------------------------

if (TRUE) {
  ## sample info
  #sample_sheet_f = "/share/home/wyz/Lab/Code/ProjectCode/2024/202403/20240305_scRNA_AML/sample_info/scRNA_Concat_Colname_Ref_20240316.xlsx"
  
  ## input_dir
  # input_dir = "/share/home/wyz/Project/result/20240320_scRNA_AML_Annotation/20240321_scRNA_GSE_GSE116256_annotation/scRNA_SMU_GSE116256_PRJNA720840_GSE205490_HRA000084_AML_MDS_83_origin_rds_combined_qc_STD_Harmony_res_anno_GES116256.rds"
  
  ## output
  
  output_dir = "/share/home/wyz/Project/result/2025/202505/20250527_1_cibersortx_immune_cell_MDS_sAML_AML_merge/20250528_4_MDS_AML_immune_cell_os_time"
  
}


# input path --------------------------------------------------------------

{
  ## experimental group
  # cell_f = "/share/home/wyz/Project/result/2025/202503/20250301_2_merge_MDS_sAML_AML_ref_immune_cell_re_cibersortx/20250303_15_GSE58831_Chip_STD_cibersortx_to_percent_k7/20250225_1_GSE58831_Clinical_SampleInfo_status_with_survival_sorted_cell_percent_k3_k7.xls"
}

{
  cell_f = "/share/home/wyz/Project/result/2025/202505/20250527_1_cibersortx_immune_cell_MDS_sAML_AML_merge/20250528_4_MDS_AML_immune_cell_os_time/20250529_3_MDS_AML_immune_cell_HR_data.xls"
}

{
  ## read file
  cell_df = read.table(cell_f,
                       sep = "\t",
                       header = T,
                       check.names = FALSE) # 阻止列名中非法符号-号变为.号 # row.names = 1,
  print(head(cell_df))
  #print(count_df)
  dim(cell_df) # [1] 342  56
  
}

{
  ## 对列名，存在-号的，俊改为_
  # 替换列名中的 - 为 _
  colnames(cell_df) <- gsub("-", "_", colnames(cell_df))
  print(colnames(cell_df))
}

{
  head(cell_df)
  unique(cell_df$Blast_Percent_group)
  cell_df$Blast_Percent_group = factor(
    cell_df$Blast_Percent_group,
    levels = c(
      "blast0_2",
      "blast2_5",
      "blast5_10",
      "blast10_20",
      "blast20_50",
      "blast50_100"
    )
  )
  
  cell_df["Group"] = cell_df["kmean_k4_mark"]
  unique(cell_df$Group)
  # cell_df$Group = factor(cell_df$Group ,levels = c("k5_3_Prog","k5_2_Primitive","k5_4_GMP","k5_5_intermediate","k5_1_Mature"))
  # cell_df$Group = factor(cell_df$Group ,levels = c("k4_2_CD8T_High","k4_1_CD8T_Int_CD4T_Low","k4_3_CD8_Int_CD4_High","k4_4_CD8T_Low"))
  # cell_df$Group = factor(cell_df$Group ,levels = c("k4_4_CD8T_Low","k4_2_CD8T_High","k4_1_CD8T_Int_CD4T_Low","k4_3_CD8_Int_CD4_High"))
  cell_df$Group = factor(
    cell_df$Group ,
    levels = c(
      "k4_2_CD8T_High",
      "k4_1_CD8T_Int_CD4T_Low",
      "k4_3_CD8_Int_CD4_High",
      "k4_4_CD8T_Low"
    )
  )
}



# -------------------------------------------------------------------------

{
  ## 将Group和Blast_Percent_group列的文本相连，生成新列Group_k_blast
  cell_df$Group_Kmean_Blast = paste(cell_df$Group, cell_df$Blast_Percent_group, sep = "_")
  
  ## cell_df 按Group和Blast_Percent_group这两列进行排序
  cell_df_sorted = cell_df %>% arrange(Group, Blast_Percent_group)
  head(cell_df_sorted)
}


{
  ## 剔除Group_Kmean_Blast计数小于3的行
  # 按Group_Kmean_Blast分组，筛选行数≥3的组
  cell_df_sorted_keep <- cell_df_sorted %>%
    group_by(Group_Kmean_Blast) %>%
    filter(n() >= 3) %>%  # n()统计每组行数
    ungroup()             # 取消分组
  
  table(cell_df_sorted_keep$Group_Kmean_Blast)
}
# -------------------------------------------------------------------------

{
  # 计算 Kruskal-Wallis 检验的 p 值
  kruskal_test_result <- kruskal.test(Survival_years ~ Blast_Percent_group, data = cell_df_sorted_keep)
  
  # 输出 p 值
  kruskal_test_result$p.value
  #[1] 4.876502e-08
}



# -------------------------------------------------------------------------
{
  # CNS ggplot格式
  cns.plot.format = theme(
    # Step1: 设置图片的背景
    # 1)设置整个图片的背景为空，这样即可获得透明背景
    plot.background = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    #2)设置不显示panel的线条，不显示panel的背景，但给panel添加边界，黑色，宽度0.5mm
    #panel.border=element_rect(color="black",linewidth=0.5,fill=NA),
    
    # Step2: 设置坐标轴和坐标轴刻度线的颜色及粗细
    #1)设置横纵坐标轴为空，因为前面已经指定了panel的边框，这里就不需要再加坐标轴了，否则会出现叠加的线条。
    #axis.line=element_blank(),   # 注释，代表保留X和Y轴的刻度轴边框
    #2)指定坐标轴刻度线为黑色，宽度0.5mm
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    
    # Step3: 设置文字的字体、字号和颜色
    # 1)指定坐标轴的文字为黑色，字号7pt
    axis.text = element_text(color = "black", size = 7),
    # 2)指定坐标轴标题文字为黑色，字号7pt
    axis.title = element_text(color = "black", size = 7),
    # 3)指定图片标题为黑色，字号7pt
    plot.title = element_text(color = "black", size = 7),
    
    # Step4: 对legend进行设置
    # 1) 指定legend的背景设为空
    legend.background = element_blank(),
    # 2）指定legend的字号和颜色
    legend.key = element_blank(),
    legend.text = element_text(color = "black", size = 7),
    legend.title = element_text(color = "black", size = 7)
  )
}


# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
{
  head(cell_df_sorted_keep)
  table(cell_df_sorted_keep$Survival_days)
  table(cell_df_sorted_keep$Group)
  table(cell_df_sorted_keep$Blast_Percent_group)
}



# -------------------------------------------------------------------------

{
  # 生成箱线图（基于原始数据cell_df，无需预处理）
  plot_object = ggplot(cell_df_sorted_keep,
                       aes(x = Survival_years, y = Group, fill = Blast_Percent_group)) +  # 按Blast百分比填充颜色
    
    # 核心修改：使用箱线图几何对象
    geom_boxplot(
      position = position_dodge(preserve = "single", width = 1),
      # 并排显示
      width = 0.7,
      # 箱体宽度
      color = "black",
      # 边框颜色
      outlier.size = 1.5,
      # 离群点大小
      outlier.shape = NA,
      # 空心圆点 # 21
      alpha = 1,
      # 透明度
      coef = 0              # 修改须线长度 (默认1.5，改为2时须线会变长)
    ) +
    # 颜色配置
    scale_fill_manual(
      name = "Blast Percentage",
      values = c(
        "blast0_2" = "#4da54a",
        "blast2_5" = "#3976ab",
        "blast5_10" = "#9d552a",
        "blast10_20" = "#8e4b97",
        "blast20_50" = "#e37fae",
        "blast50_100" = "#d5211d"
      )
    ) +
    # 坐标轴设置
    scale_x_continuous(limits = c(0, 7), breaks = seq(0, 10, 1)) +
    labs(x = "25-75% OS range (years)", y = "Taxonomy", title = "Survival Distribution by Group and Blast Percentage Kruskal-Wallis (P=1.637907e-08)") +
    # 主题优化
    theme_classic(base_size = 7) +
    theme(
      axis.text = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.5),
      legend.position = "right",
      legend.key.size = unit(0.6, "cm"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  output_f_path = file.path(output_dir, "20250529_6_MDS_AML_immune_cell_os_time.pdf")
  ggsave(
    filename = output_f_path,
    plot = plot_object,
    device = pdf,
    width = 180,
    height = 100,
    units = "mm"
  )
}

# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
{
  head(cell_df_sorted_keep)
  cell_df_sorted_keep_select = cell_df_sorted_keep %>% filter(Group == "k5_sub_k3_Prog_high_Primitive")
  table(cell_df_sorted_keep_select$Blast_Percent_group)
}
