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
  output_dir = "/share/home/wyz/Project/result/2025/202507/20250710_1_bluk_rnaseq_split_group_barplot/20250804_2_bluk_rnaseq_split_group_barplot_data_tumor_blast"
  
}

# input path --------------------------------------------------------------

{
  blast_f = "/share/home/wyz/Project/result/2025/202507/20250710_1_bluk_rnaseq_split_group_barplot/20250730_1_bluk_rnaseq_split_group_barplot_data/20250710_1_bluk_rnaseq_split_group_barplot_blast_data.xls"
  day_f = "/share/home/wyz/Project/result/2025/202507/20250710_1_bluk_rnaseq_split_group_barplot/20250730_1_bluk_rnaseq_split_group_barplot_data/20250730_1_bluk_rnaseq_split_group_barplot_Survival_days_data.xls"
}

{
  ## read file
  blast_df = read.table(blast_f,
                        sep = "\t",
                        header = T,
                        check.names = FALSE) #
  print(head(blast_df))
  dim(blast_df) #
  
}


{
  ## read file
  day_df = read.table(day_f,
                      sep = "\t",
                      header = T,
                      check.names = FALSE) #
  print(head(day_df))
  dim(day_df) #
  
}

{
  ## day_df
  library(dplyr)
  head(day_df)
  
  day_df = day_df %>% filter(Group %in% c("OS_Year_0_1", "OS_Year_1_5", "OS_Year_5_up"))
  table(day_df$Group)
}


{
  table(blast_df$Group)
  table(day_df$Group)
  
  unique(blast_df$Group)
  unique(day_df$Group)
  
  head(blast_df)
  head(day_df)
}




# -------------------------------------------------------------------------
{
  merge_df = blast_df
}

{
  ## save
  output_f_path = file.path(output_dir,
                            "20250710_1_bluk_rnaseq_split_group_barplot_data.xls")
  write.table(
    merge_df,
    file = output_f_path,
    sep = "\t",
    row.names = F,
    quote = FALSE
  )
}

# -------------------------------------------------------------------------
{
  head(merge_df)
  
  cell_list = c(
    "HSC_like_Percent",
    "Prog_like_Percent",
    "GMP_like_Percent",
    "ProMono_like_Percent",
    "Mono_like_Percent",
    "cDC_like_Percent"
  )
}


{
  ## long df to wide df
  library(tidyr)
  
  df_long <- merge_df %>%
    pivot_longer(cols = cell_list,
                 names_to = "Cell_Type",
                 values_to = "Value")
  head(df_long)
  print(dim(df_long))
  
  table(df_long$Group)
  unique(df_long$Group)
  df_long[is.na(df_long$Group), ]
}

# -------------------------------------------------------------------------

{
  table(df_long$Group)
  df_long$Group = factor(
    df_long$Group,
    levels = c(
      "Blast_0_2",
      "Blast_2_5",
      "Blast_5_10",
      "Blast_10_20",
      "Blast_20_40",
      "Blast_40_80",
      "Blast_80_100"
    )
  )
  
  table(df_long$Group)
  unique(df_long$Group)
}

{
  df_long$Cell_Type_Mark <- gsub("_Percent", "", df_long$Cell_Type)
  cell_list2 = c("HSC_like",
                 "Prog_like",
                 "GMP_like",
                 "ProMono_like",
                 "Mono_like",
                 "cDC_like")
  df_long$Cell_Type_Mark = factor(df_long$Cell_Type_Mark , levels = cell_list2)
  head(df_long)
  head(df_long)
  
}


{
  ## save
  output_f_path = file.path(output_dir,
                            "20250710_2_bluk_rnaseq_split_group_barplot_data_long_df.xls")
  write.table(
    df_long,
    file = output_f_path,
    sep = "\t",
    row.names = F,
    quote = FALSE
  )
  
}

# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
{
  ## line plot
}

{
  head(df_long)
  df_avg  = df_long %>% group_by(Group, Cell_Type_Mark) %>% summarize(Value = mean(Value, na.rm = TRUE), groups = "drop")
  df_avg
}

{
  library(dplyr)
  df_summary <- df_long %>%
    group_by(Cell_Type_Mark, Group) %>%
    summarize(
      Mean_Value = mean(Value, na.rm = TRUE),
      SD_Value = sd(Value, na.rm = TRUE),
      n = n(),
      #
      SEM_Value = SD_Value / sqrt(n)
    )  #
  
  df_summary
  
}

{
  plot_object = ggplot(
    df_summary,
    aes(
      x = Group,
      y = Mean_Value,
      group = Cell_Type_Mark,
      color = Cell_Type_Mark
    )
  ) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    geom_errorbar(
      aes(ymin = Mean_Value - SEM_Value, ymax = Mean_Value + SEM_Value),
      width = .2,
      position = position_dodge(0.00)
    ) +
    geom_hline(
      data = subset(df_summary, Group == "Blast_0_2"),
      aes(yintercept = Mean_Value, color = Cell_Type_Mark),
      linetype = "dashed",
      size = 0.5
    ) +  #
    labs(x = NULL, y = "Percent", title = "Tumor Cells Percent"  # ) +
         
         scale_y_continuous(
           labels = scales::percent,
           breaks = seq(0, 1, by = 0.05),
           limits = c(0, 0.5),
         ) +
           theme(
             legend.position = "right",
             #
             panel.background = element_rect(fill = "white", color = NA),
             #
             plot.background = element_rect(fill = "white", color = NA),
             #
             strip.text = element_text(face = "bold"),
             axis.text.x = element_text(angle = 45, hjust = 1),
             axis.line.x = element_line(color = "black", size = 0.5),
             axis.line.y = element_line(color = "black", size = 0.5)
           ) +
           scale_color_manual(
             values = c(
               "HSC_like" = "#1163a6",
               "Prog_like" = '#b156c1',
               "GMP_like" = "#25af31",
               "ProMono_like" = '#7f9c21',
               "Mono_like" = "#65cdc3",
               "cDC_like" = "#f6b9a4"
             )
           )
         
         plot_object
         output_f_path = file.path(output_dir,
                                   "20250730_3_9_bluk_rnaseq_split_group_barplot_data_long_df.pdf")
         ggsave(
           filename = output_f_path,
           plot = plot_object,
           device = pdf,
           width = 210,
           height = 160,
           units = "mm"
         )
}