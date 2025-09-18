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
  library(survminer) #
  library(survival) #
}


# input path --------------------------------------------------------------
if (TRUE) {
  output_dir = "/share/home/wyz/Project/result/2025/202505/20250521_1_cibersortx_tumor_cell_MDS_sAML_AML_merge/20250522_3_MDS_AML_tumor_cell_HR"
}


# input path --------------------------------------------------------------
{
  cell_f = "/share/home/wyz/Project/result/2025/202505/20250521_1_cibersortx_tumor_cell_MDS_sAML_AML_merge/20250522_3_MDS_AML_tumor_cell_HR/20250522_3_MDS_AML_tumor_cell_HR_data.xls"
}

{
  ## read file
  cell_df = read.table(cell_f,
                       sep = "\t",
                       header = T,
                       check.names = FALSE)
  print(head(cell_df))
  dim(cell_df) #
  
}

{
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
  
  cell_df["Group"] = cell_df["kmean_k5_mark"]
  unique(cell_df$Group)
  cell_df$Group = factor(
    cell_df$Group ,
    levels = c(
      "k5_3_Prog",
      "k5_2_Primitive",
      "k5_4_GMP",
      "k5_5_intermediate",
      "k5_1_Mature"
    )
  )
}



# -------------------------------------------------------------------------
{
  cell_df$Group_Kmean_Blast = paste(cell_df$Group, cell_df$Blast_Percent_group, sep = "_")
  
  cell_df_sorted = cell_df %>% arrange(Group, Blast_Percent_group)
  head(cell_df_sorted)
}


{
  cell_df_sorted_keep <- cell_df_sorted %>%
    group_by(Group_Kmean_Blast) %>%
    filter(n() >= 3) %>%  #
    ungroup()             #
  
  table(cell_df_sorted_keep$Group_Kmean_Blast)
}
# -------------------------------------------------------------------------

{
  # Kruskal-Wallis test
  kruskal_test_result <- kruskal.test(Survival_years ~ Blast_Percent_group, data = cell_df_sorted_keep)
  
  # p
  kruskal_test_result$p.value
  #[1] 1.637907e-08
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
{
  head(cell_df_sorted_keep)
  table(cell_df_sorted_keep$Survival_days)
  table(cell_df_sorted_keep$Group)
  table(cell_df_sorted_keep$Blast_Percent_group)
}

# -------------------------------------------------------------------------
{
  plot_object = ggplot(cell_df_sorted_keep,
                       aes(x = Survival_years, y = Group, fill = Blast_Percent_group)) +
    geom_boxplot(
      position = position_dodge(preserve = "single", width = 1),
      width = 0.7,
      color = "black",
      #
      outlier.size = 1.5,
      #
      outlier.shape = NA,
      #
      alpha = 1,
      #
      coef = 0              #
    ) +
    #
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
    #
    scale_x_continuous(limits = c(0, 7), breaks = seq(0, 10, 1)) +
    labs(x = "25-75% OS range (years)", y = "Taxonomy", title = "Survival Distribution by Group and Blast Percentage Kruskal-Wallis (P=1.637907e-08)") +
    #
    theme_classic(base_size = 7) +
    theme(
      axis.text = element_text(color = "black"),
      axis.line = element_line(linewidth = 0.5),
      legend.position = "right",
      legend.key.size = unit(0.6, "cm"),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  output_f_path = file.path(output_dir, "20250522_4_MDS_AML_tumor_cell_os_time.pdf")
  ggsave(
    filename = output_f_path,
    plot = plot_object,
    device = pdf,
    width = 180,
    height = 100,
    units = "mm"
  )
}