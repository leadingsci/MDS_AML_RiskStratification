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
  library(ggpubr) #R 4.2.2
  
}



# input path --------------------------------------------------------------

if (TRUE) {
  ## output
  output_dir = "/share/home/wyz/Project/result/2025/202507/20250708_1_kmean_group_sign_barplot/20250708_7_kmean_group_sign_barplot_immune"
}

# input path --------------------------------------------------------------

{
  cell_f = "/share/home/wyz/Project/result/2025/202507/20250708_1_kmean_group_sign_barplot/20250708_4_kmean_group_sign_barplot_immune/20250708_1_blast_barplot_MDS_sAML_AML_like_map_immune_kmean_group_4787_sample_info_and_kmean.xls"
}

{
  ## read file
  cell_df = read.table(cell_f,
                       sep = "\t",
                       header = T,
                       check.names = FALSE)
  print(head(cell_df))
  dim(cell_df)
  
}


{
  colnames(cell_df) <- gsub("-", "_", colnames(cell_df))
  print(colnames(cell_df))
}


# -------------------------------------------------------------------------

{
  ##
  select_list = c("Sample_ID", "DiseaseClass", "kmean_k4_mark")
  
  cell_keep_df = cell_df %>% select(select_list)
  dim(cell_keep_df)
  head(cell_keep_df)
}


{
  table(cell_keep_df$DiseaseClass)
  # AML  MDS sAML
  # 3370 1415   13
}

# -------------------------------------------------------------------------
{
  cell_keep_df = cell_keep_df %>% filter(DiseaseClass %in% c("MDS", "AML"))
  
  ## 固定顺序
  cell_keep_df$DiseaseClass = factor(cell_keep_df$DiseaseClass, levels = rev(c("MDS", "AML")))
  table(cell_keep_df$DiseaseClass)
  
  cell_keep_df$kmean_k4_mark = factor(
    cell_keep_df$kmean_k4_mark ,
    levels = c(
      "k4_2_CD8T_High",
      "k4_1_CD8T_Int_CD4T_Low",
      "k4_3_CD8_Int_CD4_High",
      "k4_4_CD8T_Low"
    )
  )
  
}


{
  data_summary <- cell_keep_df %>%
    group_by(DiseaseClass, kmean_k4_mark) %>%
    summarise(count = n()) %>%
    group_by(DiseaseClass) %>%
    mutate(percentage = count / sum(count) * 100)
  data_summary
  
}



{
  output_f_path = file.path(output_dir,
                            "20250708_4_MDS_sAML_AML_immune_cell_clinical_data_simple.xls")
  write.table(
    data_summary,
    file = output_f_path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
}


if (TRUE) {
  plot_object = ggbarplot(
    data_summary,
    x = "kmean_k4_mark",
    y = "percentage",
    color = "DiseaseClass",
    fill = "DiseaseClass",
    # palette = "npg",
    position = position_dodge(0.8)
  ) +
    scale_fill_manual(values = c("MDS" = "#a3ba87", "AML" = "#6ca7b2")) +  #
    scale_color_manual(values = c("MDS" = "#a3ba87", "AML" = "#6ca7b2")) + #
    # rotate_x_text(angle = 45) +
    coord_flip() +  #
    labs(x = NULL, y = "Percentage(%)", title = "Immune Kmean Group") +  #
    scale_y_continuous(breaks = seq(0, 50, by = 10), limits = c(0, 50)) +
    #
    theme(legend.position = "right") +  #
    geom_text(
      aes(label = round(percentage, 1)),
      position = position_dodge(0.8),
      vjust = 0.5,
      size.unit = "pt",
      size = 5
    )  #
  
  output_f_path = file.path(
    output_dir,
    "20250708_1_3_MDS_sAML_AML_immune_cell_clinical_data_simple_w170mm_h120mm.pdf"
  )
  ggsave(
    filename = output_f_path,
    plot = plot_object,
    device = pdf,
    width = 170,
    height = 120,
    units = "mm"
  )
}
