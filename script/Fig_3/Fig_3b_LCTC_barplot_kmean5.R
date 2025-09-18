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
  output_dir = "/share/home/wyz/Project/result/2025/202506/20250605_20250520_1_cibersortx_tumor_cell_kmean_5448_MDS_sAML_AML_merge/20250605_1_scRNA_subtype_MDS_sAML_AML_tumor_percent_survival_barplot_k5"
}

# input path --------------------------------------------------------------

{
  ## experimental group
  cell_f = "/share/home/wyz/Project/result/2025/202505/20250520_1_cibersortx_tumor_cell_MDS_sAML_AML_merge/20250520_3_tumor_cell_percent_kmean_split_group/20250520_3_MDS_sAML_sAML_tumor_cell_percent_kmean_split_group.xls"
}

{
  ## read file
  cell_df = read.table(cell_f,
                       sep = "\t",
                       header = T,
                       check.names = FALSE)
  print(head(cell_df))
  #print(count_df)
  dim(cell_df) # [1] 5448   41
  
}

{
  colnames(cell_df) <- gsub("-", "_", colnames(cell_df))
  print(colnames(cell_df))
}



# -------------------------------------------------------------------------
{
  # Calculate the sum of the "HSC" and "Prog" columns for each sample
  library(dplyr)
  cell_df_HSC_Prog_sum = cell_df %>% mutate(HSC_Prog_sum = HSC_like_Percent  + Prog_like_Percent)
  head(cell_df_HSC_Prog_sum)
}

{
  # Calculate the total sum of HSC_Prog_sum for each of the k2 to k8 groups
  # Sorting for each group
  
  cell_df_HSC_Prog_sum_k2 <- cell_df_HSC_Prog_sum %>%
    group_by(kmean_k2) %>%
    mutate(HSC_Prog_sum_kmean_k2 = sum(HSC_Prog_sum)) %>%
    ungroup()
  
  cell_df_HSC_Prog_sum_k3 <- cell_df_HSC_Prog_sum_k2 %>%
    group_by(kmean_k3) %>%
    mutate(HSC_Prog_sum_kmean_k3 = sum(HSC_Prog_sum)) %>%
    ungroup()
  
  cell_df_HSC_Prog_sum_k4 <- cell_df_HSC_Prog_sum_k3 %>%
    group_by(kmean_k4) %>%
    mutate(HSC_Prog_sum_kmean_k4 = sum(HSC_Prog_sum)) %>%
    ungroup()
  
  cell_df_HSC_Prog_sum_k5 <- cell_df_HSC_Prog_sum_k4 %>%
    group_by(kmean_k5) %>%
    mutate(HSC_Prog_sum_kmean_k5 = sum(HSC_Prog_sum)) %>%
    ungroup()
  
  cell_df_HSC_Prog_sum_k6 <- cell_df_HSC_Prog_sum_k5 %>%
    group_by(kmean_k6) %>%
    mutate(HSC_Prog_sum_kmean_k5 = sum(HSC_Prog_sum)) %>%
    ungroup()
  
  cell_df_HSC_Prog_sum_k7 <- cell_df_HSC_Prog_sum_k6 %>%
    group_by(kmean_k7) %>%
    mutate(HSC_Prog_sum_kmean_k7 = sum(HSC_Prog_sum)) %>%
    ungroup()
  
  cell_df_HSC_Prog_sum_k8 <- cell_df_HSC_Prog_sum_k7 %>%
    group_by(kmean_k8) %>%
    mutate(HSC_Prog_sum_kmean_k8 = sum(HSC_Prog_sum)) %>%
    ungroup()
  
  cell_df_HSC_Prog_sum_k9 <- cell_df_HSC_Prog_sum_k8 %>%
    group_by(kmean_k9) %>%
    mutate(HSC_Prog_sum_kmean_k9 = sum(HSC_Prog_sum)) %>%
    ungroup()
  
  
  cell_df_HSC_Prog_sum_k9[, c("Sample_ID",
                              "kmean_k8",
                              "HSC_Prog_sum",
                              "HSC_Prog_sum_kmean_k8")]
  
}



# -------------------------------------------------------------------------
{
  ## Change the long table to a wide table
  cols_to_normalize = c("HSC_like",
                        "Prog_like",
                        "GMP_like",
                        "ProMono_like",
                        "Mono_like",
                        "cDC_like")
  cell_list = paste0(cols_to_normalize, "_Percent")
  library(tidyr)
  
  
  cell_df_long = cell_df_HSC_Prog_sum_k9 %>%
    pivot_longer(cols = cell_list,
                 names_to = "Cell_Type",
                 values_to = "Cell_Value")
  
  cell_df_long
  print(colnames(cell_df_long))
  cell_df_long[, c("Sample_ID",
                   "kmean_k8",
                   "Cell_Type",
                   "Cell_Value",
                   "HSC_Prog_sum")]
  
  
}

{
  ## rename
  cell_df_long = cell_df_long %>%
    mutate(Cell_Type = str_replace(Cell_Type, "_Percent$", ""))
  
  cell_df_long[, c("Sample_ID",
                   "kmean_k8",
                   "Cell_Type",
                   "Cell_Value",
                   "HSC_Prog_sum")]
}


# -------------------------------------------------------------------------
{
  ## factor
  cell_list = c("HSC_like",
                "Prog_like",
                "GMP_like",
                "ProMono_like",
                "Mono_like",
                "cDC_like")
  
  cell_df_long$Cell_Type = factor(cell_df_long$Cell_Type, levels = cell_list)
  print(cell_df_long$Cell_Type)
}


# -------------------------------------------------------------------------

{
  # save
  output_f_path = file.path(
    output_dir,
    "20250605_1_scRNA_subtype_MDS_AML_tumor_percent_survival_barplot.xls"
  )
  write.table(
    cell_df_long,
    file = output_f_path,
    sep = "\t",
    row.names = F,
    quote = FALSE
  )
  
}


# -------------------------------------------------------------------------

# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
# k5_1
# -------------------------------------------------------------------------
{
  cell_select_df = cell_df_long[cell_df_long$kmean_k5 == "1", ]
  
  cell_select_df[, c("Sample_ID",
                     "kmean_k5",
                     "Cell_Type",
                     "Cell_Value",
                     "HSC_Prog_sum")]
}


# -------------------------------------------------------------------------

{
  ## sort
  library(dplyr)
  cell_select_df_sorted <- cell_select_df %>%
    arrange(desc(HSC_Prog_sum))
  
  cell_select_df_sorted[, c("Sample_ID",
                            "kmean_k5",
                            "Cell_Type",
                            "Cell_Value",
                            "HSC_Prog_sum")]
  
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
    "HSC_like" = "#ca8282",
    "Prog_like" = "#e19793",
    "GMP_like" = "#898989",
    "ProMono_like" = "#66abda",
    "Mono_like" = "#668aa3",
    "cDC_like" = "#f5d966"
  )
  
}



# -------------------------------------------------------------------------

{
  group_name = "k5_1"
  sample_num = 1652
  
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
      axis.text.y = element_text(size = 7)
    ) +
    theme(plot.margin = margin(1, 1, 1, 1.2, "cm")) + #
    theme(panel.grid = element_blank()) +  #
    theme(
      legend.position = "top",
      #
      legend.title = element_text(size = 6),
      #
      legend.text = element_text(size = 6),
      #
      legend.key.size = unit(2, "mm"),
      #
    )
  
  plot_object
  
  output_f_path = file.path(output_dir,
                            paste0("20250605_1_2_", group_name, "_barplot_0.1.pdf"))
  ggsave(
    filename = output_f_path,
    plot = plot_object,
    device = pdf,
    width = 0.1 * sample_num,
    height = 120,
    units = "mm",
    limitsize = FALSE
  )
}



# -------------------------------------------------------------------------
# k5_2
# -------------------------------------------------------------------------
{
  cell_select_df = cell_df_long[cell_df_long$kmean_k5 == "2", ]
  cell_select_df[, c("Sample_ID",
                     "kmean_k5",
                     "Cell_Type",
                     "Cell_Value",
                     "HSC_Prog_sum")]
}


# -------------------------------------------------------------------------

{
  library(dplyr)
  cell_select_df_sorted <- cell_select_df %>%
    arrange(desc(HSC_Prog_sum))
  
  cell_select_df_sorted[, c("Sample_ID",
                            "kmean_k5",
                            "Cell_Type",
                            "Cell_Value",
                            "HSC_Prog_sum")]
  
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
  group_name = "k5_2"
  sample_num = 976
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
      axis.text.y = element_text(size = 7)
    ) +
    theme(plot.margin = margin(1, 1, 1, 1.2, "cm")) + #
    theme(panel.grid = element_blank()) +  #
    theme(
      legend.position = "top",
      #
      legend.title = element_text(size = 6),
      #
      legend.text = element_text(size = 6),
      #
      legend.key.size = unit(2, "mm"),
      #
    )
  
  plot_object
  
  output_f_path = file.path(output_dir,
                            paste0("20250605_2_2_", group_name, "_barplot_0.1.pdf"))
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
# k5_3
# -------------------------------------------------------------------------
{
  cell_select_df = cell_df_long[cell_df_long$kmean_k5 == "3", ]
  
  cell_select_df[, c("Sample_ID",
                     "kmean_k5",
                     "Cell_Type",
                     "Cell_Value",
                     "HSC_Prog_sum")]
}


# -------------------------------------------------------------------------

{
  library(dplyr)
  cell_select_df_sorted <- cell_select_df %>%
    arrange(desc(HSC_Prog_sum))
  
  cell_select_df_sorted[, c("Sample_ID",
                            "kmean_k5",
                            "Cell_Type",
                            "Cell_Value",
                            "HSC_Prog_sum")]
  
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
  group_name = "k5_3"
  sample_num = 1765
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
      axis.text.y = element_text(size = 7)
    ) +
    theme(plot.margin = margin(1, 1, 1, 1.2, "cm")) + #
    theme(panel.grid = element_blank()) +  #
    theme(
      legend.position = "top",
      #
      legend.title = element_text(size = 6),
      #
      legend.text = element_text(size = 6),
      #
      legend.key.size = unit(2, "mm"),
      #
    )
  
  plot_object
  
  output_f_path = file.path(output_dir,
                            paste0("20250520_3_1_", group_name, "_barplot_0.1.pdf"))
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
# k5_4
# -------------------------------------------------------------------------
{
  cell_select_df = cell_df_long[cell_df_long$kmean_k5 == "4", ]
  
  cell_select_df[, c("Sample_ID",
                     "kmean_k5",
                     "Cell_Type",
                     "Cell_Value",
                     "HSC_Prog_sum")]
}


# -------------------------------------------------------------------------

{
  library(dplyr)
  cell_select_df_sorted <- cell_select_df %>%
    arrange(desc(HSC_Prog_sum))               #
  
  cell_select_df_sorted[, c("Sample_ID",
                            "kmean_k5",
                            "Cell_Type",
                            "Cell_Value",
                            "HSC_Prog_sum")]
  
  print(cell_select_df_sorted$HSC_Prog_sum)
  dim(cell_select_df_sorted)
}

{
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
  group_name = "k5_4"
  sample_num = 605
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
    scale_y_continuous(expand = c(0, 0), label = scales::percent) + #
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
      axis.text.y = element_text(size = 7)
    ) +
    theme(plot.margin = margin(1, 1, 1, 1.2, "cm")) + #
    theme(panel.grid = element_blank()) +  #
    theme(
      legend.position = "top",
      #
      legend.title = element_text(size = 6),
      #
      legend.text = element_text(size = 6),
      #
      legend.key.size = unit(2, "mm"),
      #
    )
  
  plot_object
  
  output_f_path = file.path(output_dir,
                            paste0("20250520_4_1_", group_name, "_barplot_0.1.pdf"))
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
# k5_5
# -------------------------------------------------------------------------
{
  cell_select_df = cell_df_long[cell_df_long$kmean_k5 == "5", ]
  
  cell_select_df[, c("Sample_ID",
                     "kmean_k5",
                     "Cell_Type",
                     "Cell_Value",
                     "HSC_Prog_sum")]
}


# -------------------------------------------------------------------------

{
  library(dplyr)
  cell_select_df_sorted <- cell_select_df %>%
    arrange(desc(HSC_Prog_sum))               #
  
  cell_select_df_sorted[, c("Sample_ID",
                            "kmean_k5",
                            "Cell_Type",
                            "Cell_Value",
                            "HSC_Prog_sum")]
  
  print(cell_select_df_sorted$HSC_Prog_sum)
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
  group_name = "k5_5"
  sample_num = 450
  plot_object =
    ggplot(cell_select_df_sorted,
           aes(x = Sample_ID, y = Cell_Value, fill = Cell_Type)) +
    geom_bar(
      stat = "identity",
      alpha = 1,
      width = 0.8,
      color = NA,
      size = 0
    ) +  #
    scale_x_discrete(expand = c(0, 0)) +  #
    scale_y_continuous(expand = c(0, 0), label = scales::percent) + #
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
      axis.text.y = element_text(size = 7)
    ) +
    theme(plot.margin = margin(1, 1, 1, 1.2, "cm")) + #
    theme(panel.grid = element_blank()) +  #
    theme(
      legend.position = "top",
      #
      legend.title = element_text(size = 6),
      #
      legend.text = element_text(size = 6),
      #
      #
      legend.key.size = unit(2, "mm"),
      #
    )
  
  plot_object
  
  output_f_path = file.path(output_dir,
                            paste0("20250520_5_1_", group_name, "_barplot_0.1.pdf"))
  ggsave(
    filename = output_f_path,
    plot = plot_object,
    device = pdf,
    width = 0.3 * sample_num,
    height = 120,
    units = "mm"
  )
}