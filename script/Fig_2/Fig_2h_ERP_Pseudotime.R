.libPaths()

rm(list = ls())

# init waring
options(warn = 1)


# Load Package ------------------------------------------------------------

library(tidyverse)

library(readr)
library(dplyr)

library(ggplot2)

# library(openxlsx);
library(readxl)

library(RColorBrewer)

library(Seurat)

# library(future);


library(SummarizedExperiment)
library(scCustomize)
library(viridis)
library(gridExtra)
library(ComplexHeatmap)

library(ggplotify)



if (TRUE) {
  ## output
  output_dir = "/share/home/wyz/Project/result/2025/202508/20240801_scRNA_GSE_Rice_summary"
}




# Load Package ------------------------------------------------------------
{
  input_f_path = "/share/home/wyz/Project/result/20240401_scRNA_AML_Annotation/20240410_scRNA_GSE_Rice_HSC_monocle/seurat_object_STD_select_HSC_MEP_ERP_re/scRNA_GSE_Rice_monocle_Pseudotime_re_20240412_1_1.xlsx"
  
  ptime_df = readxl::read_excel(
    path = input_f_path,
    col_name = T,
    col_types = "text",
    na = ""
  )
  dim(ptime_df)
  print(ptime_df)
  head(ptime_df)
  class(ptime_df)
  colnames(ptime_df)
  ptime_df[, c("Size_Factor", "num_genes_expressed", "Pseudotime", "State")]
}



{
  ptime_df$Pseudotime = as.numeric(ptime_df$Pseudotime)
  
  
  ptime_df <- ptime_df %>%
    dplyr::mutate(Pseudotime_streched = (Pseudotime - min(Pseudotime)) / (max(Pseudotime) - min(Pseudotime)) * 100)
  
  ptime_df[, c(
    "Size_Factor",
    "num_genes_expressed",
    "Pseudotime",
    "State",
    "Pseudotime_streched"
  )]
  
}

{
  # factor
  unique(ptime_df$DiseaseClass)
  ptime_df$DiseaseClass = factor(ptime_df$DiseaseClass,
                                 levels = c("Normal", "MDS", "tAML", "AML"))
  
  unique(ptime_df$CellType_Rice_20240401_1)
  ptime_df$CellType_Rice_20240401_1 = factor(ptime_df$CellType_Rice_20240401_1,
                                             levels = c("HSC/MPP", "MEP", "ERP"))
  
  unique(ptime_df$Sample_ID)
  ptime_df$Sample_ID = factor(ptime_df$Sample_ID, levels = unique(ptime_df$Sample_ID))
}


{
  ## 选取细胞类型
  ptime_df_ERP = subset(ptime_df, subset = CellType_Rice_20240401_1 == "ERP")
  unique(ptime_df_ERP$CellType_Rice_20240401_1)
  
  ## 取样本均值
  ptime_df_ERP_sample_mean = ptime_df_ERP %>% group_by(DiseaseClass, Sample_ID, celltype) %>%
    mutate(Mean_Pseudotime = mean(Pseudotime_streched, na.rm = TRUE))
  
  ## 去重
  ptime_df_ERP_sample_mean = select(ptime_df_ERP_sample_mean,
                                    DiseaseClass,
                                    Sample_ID,
                                    celltype,
                                    Mean_Pseudotime) %>% distinct()
  
  ptime_df_ERP_sample_mean
  table(ptime_df_ERP_sample_mean$Mean_Pseudotime)
  
  
  ## 排序
  ptime_df_ERP_sample_mean$DiseaseClass = factor(ptime_df_ERP_sample_mean$DiseaseClass,
                                                 levels = c("Normal", "MDS", "tAML", "AML"))
  unique(ptime_df_ERP_sample_mean$Sample_ID)
  
  output_f_path = file.path(
    output_dir,
    "STAT-scRNA_GSE_Rice_monocle_Pseudotime_re_ERP_mean_20240412_1_1.xlsx"
  )
  write.table(
    ptime_df_ERP_sample_mean,
    file = output_f_path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
}


# -------------------------------------------------------------------------
# ERP -------------------------------------------------------------------------
{
  seurat_object_plot = ggplot(data = ptime_df_ERP,
                              aes(x = Pseudotime_streched, y = DiseaseClass, color = DiseaseClass)) +
    geom_boxplot(outlier.shape = NA, fill = NA) +
    
    stat_boxplot(geom = 'errorbar', width = 0.3) +
    
    # geom_jitter(aes(color=DiseaseClass), size=3, alpha=0.9, width = 0.1) +
    # scale_color_npg() +
    scale_color_manual(values = c('#00A087FF', '#4DBBD5FF', '#F39B7FFF', '#E64B35FF')) +
    scale_x_continuous(limits = c(60, 120),
                       breaks = seq(10, 100, by = 10)) +
    geom_signif(
      comparisons = list(
        c("Normal", "MDS"),
        c("Normal", "tAML"),
        c("Normal", "AML"),
        c("MDS", "tAML"),
        c("tAML", "AML"),
        c("MDS", "AML")
      ),
      map_signif_level = T,
      test = wilcox.test,
      #
      # test = t.test, #
      step_increase = 0.02,
      size = 0.2,
      color = "black",
      textsize = 4
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, vjust = 0.5),
      # axis.ticks.x = element_blank(),
      # axis.text.x = element_blank(),
      axis.text.x = element_text(angle = 0, hjust = 1)
    ) +
    
    labs(title = "ERP - wilcox.test", x = "Pseudotime(streched) of ERP", y =
           "DiseaseClass")
  
  seurat_object_plot
  
  output_f_path = file.path(output_dir, "20250801_2_1_scRNA_cell_ERP.pdf")
  ggsave(
    filename = output_f_path,
    plot = seurat_object_plot,
    device = pdf,
    width = 250,
    height = 148,
    units = "mm"
  )
}