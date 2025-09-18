rm(list = ls())

.libPaths()
options(warn = 1)

# Load Package ------------------------------------------------------------
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
  output_dir = "/share/home/wyz/Project/result/2025/202507/20250708_1_kmean_group_sign_barplot/20250708_6_kmean_group_sign_barplot_tumor"
}

# input path --------------------------------------------------------------

{
  cell_f = "/share/home/wyz/Project/result/2025/202507/20250708_1_kmean_group_sign_barplot/20250708_1_kmean_group_sign_barplot_tumor/20250708_1_blast_barplot_MDS_sAML_AML_like_map_kmean_group_5448_sample_info_and_kmean.xls"
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

{
  select_list = c("Sample_ID", "DiseaseClass", "kmean_k5_mark")
  
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
  
  cell_keep_df$kmean_k5_mark = factor(
    cell_keep_df$kmean_k5_mark ,
    levels = c(
      "k5_4_GMP",
      "k5_1_Mature",
      "k5_3_Prog",
      "k5_2_Primitive",
      "k5_5_intermediate"
    )
  )
  
}


{
  # Calculate the proportion of each kmean_k5_mark
  data_summary <- cell_keep_df %>%
    group_by(DiseaseClass, kmean_k5_mark) %>%
    summarise(count = n()) %>%
    group_by(DiseaseClass) %>%
    mutate(percentage = count / sum(count) * 100)
  data_summary
}


{
  output_f_path = file.path(output_dir,
                            "20250708_1_MDS_sAML_AML_tumor_cell_clinical_data_simple.xls")
  write.table(
    data_summary,
    file = output_f_path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
}

{
  data_summary$DiseaseClass
  data_summary$kmean_k5_mark
}

{
  plot_object = ggbarplot(
    data_summary,
    x = "kmean_k5_mark",
    y = "percentage",
    color = "DiseaseClass",
    fill = "DiseaseClass",
    position = position_dodge(0.8)
  ) +
    scale_fill_manual(values = c("MDS" = "#a3ba87", "AML" = "#6ca7b2")) +
    scale_color_manual(values = c("MDS" = "#a3ba87", "AML" = "#6ca7b2")) +
    coord_flip() +
    labs(x = NULL, y = "Percentage(%)", title = "Tumor Kmean Group") +  #
    scale_y_continuous(breaks = seq(0, 50, by = 10), limits = c(0, 50)) +
    #
    theme(legend.position = "right") +  #
    geom_text(
      aes(label = round(percentage, 1)),
      position = position_dodge(0.8),
      vjust = 0.5,
      size.unit = "pt",
      size = 5
    )
  
  output_f_path = file.path(
    output_dir,
    "20250708_1_1_MDS_sAML_AML_tumor_cell_clinical_data_simple_w150mm_h120mm.pdf"
  )
  ggsave(
    filename = output_f_path,
    plot = plot_object,
    device = pdf,
    width = 150,
    height = 120,
    units = "mm"
  )
}