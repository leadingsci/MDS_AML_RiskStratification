rm(list = ls())


.libPaths(c("/share/home/bioinfo/anaconda3/envs/R-4.2.3/lib/R/library"))
.libPaths()

# init waring
options(warn = 1)


# Load Package ------------------------------------------------------------
{
  library(tidyverse)
  
  library(readr)
  library(dplyr)
  
  library(ggplot2)
  
  # library(openxlsx);
  library(readxl)
  
  library(RColorBrewer)
  library(ggpubr)
  library(Seurat)
  
  # library(future);
  library(SummarizedExperiment)
  library(scCustomize)
  library(viridis)
  library(gridExtra)
  library(ComplexHeatmap)
  
  library(ggplotify)
}






if (TRUE) {
  ## output
  output_dir = "/share/home/wyz/Project/result/2025/202508/20250829_1_scRNA_GSE_Rice_GMP_Neu_Annotation/20250829_1_scRNA_GSE_Rice_GMP_Neu_Module_GSF3_Pathway_Neu"
  
}


if (TRUE) {
  input_f_path_GSE = "/share/home/wyz/Project/result/20240415_scRNA_GSE_Rice_GMP_Neu_Annotation/scRNA_SMU_Rice_14_GSE116256_Rice_1_HRA000084_Rice_4_origin_rds_combined_qc_STD_Harmony_res_select_GMP_Neu_harmony_res.rds"
  seurat_object_STD_select_Neu_Neu <- readRDS(input_f_path_GSE)
  levels(seurat_object_STD_select_Neu_Neu)
  table(seurat_object_STD_select_Neu_Neu@meta.data$Sample_ID)
}



{
  ## Neu
  seurat_object_STD_select_Neu = subset(x = seurat_object_STD_select_Neu_Neu, idents = c("Neu1", "Neu2"))
  levels(seurat_object_STD_select_Neu)
  table(seurat_object_STD_select_Neu@meta.data$Sample_ID)
}

{
  head(seurat_object_STD_select_Neu@meta.data)
  table(seurat_object_STD_select_Neu@meta.data$CellType_Rice_GMP_Neu_20240415_1)
}


# gene list ------------------------------------------------------------
go_f_1 = "/share/home/wyz/Lab/Code/ProjectCode/2024/202404/sample_info/GSF3_Pathway.txt"
go_1_df = read.table(file = go_f_1, sep = "\t", header = FALSE)

head(go_1_df)
dim(go_1_df)

#  --------------------------------------------------------
go_f_1_df_keep = go_1_df[go_1_df$V1 %in% rownames(seurat_object_STD_select_Neu@assays$RNA@counts), ]
head(go_f_1_df_keep)
class(go_f_1_df_keep)
dim(go_f_1_df_keep) #


# add Module ------------------------------------------------------------
# # AddModuleScore  ---------------------------------------------------
seurat_object_STD_select_Neu =
  AddModuleScore(
    object = seurat_object_STD_select_Neu,
    features = list(go_f_1_df_keep),
    #
    ctrl = 100,
    # default
    name = 'GSF3_Pathway_Score',
    #
    seed = 1,
  )


# Module ------------------------------------------------------------

head(seurat_object_STD_select_Neu@meta.data)

#  ------------------------------------------------------------

#  DataFrame
Neu_score_result_cbind <- cbind(
  Sample_ID = seurat_object_STD_select_Neu@meta.data$Sample_ID,
  DiseaseClass = seurat_object_STD_select_Neu@meta.data$DiseaseClass,
  RNA_snn_res.0.2 = seurat_object_STD_select_Neu@meta.data$RNA_snn_res.0.2,
  CellType_Rice_GMP_Neu_20240415_1 = as.vector(
    seurat_object_STD_select_Neu@meta.data$CellType_Rice_GMP_Neu_20240415_1
  ),
  GSF3_Pathway_Score = seurat_object_STD_select_Neu@meta.data$GSF3_Pathway_Score1
) #

head(Neu_score_result_cbind)

#
Neu_score_result_cbind <- as.data.frame(Neu_score_result_cbind)
head(Neu_score_result_cbind)

output_f_path = file.path(output_dir,
                          "STAT-GSE_Rice_Neu_GSF3_Pathway_Score1_20240416_1.xls") #
write.table(
  Neu_score_result_cbind,
  file = output_f_path,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)


# data ------------------------------------------------------------
{
  # factor ------------------------------------------------------------
  table(Neu_score_result_cbind$DiseaseClass)
  Neu_score_result_cbind$DiseaseClass = factor(Neu_score_result_cbind$DiseaseClass,
                                               levels = c("Normal", "MDS", "tAML", "AML"))
  levels(Neu_score_result_cbind$DiseaseClass)
  
  #
  head(Neu_score_result_cbind)
  Neu_score_result_cbind$GSF3_Pathway_Score
  Neu_score_result_cbind$GSF3_Pathway_Score = as.numeric(Neu_score_result_cbind$GSF3_Pathway_Score)
  
  Neu_score_result_diseaseclass_mean_df <- Neu_score_result_cbind %>%
    group_by(DiseaseClass) %>%
    mutate(Average_GSF3_Pathway_Score = mean(GSF3_Pathway_Score))
  
  head(Neu_score_result_diseaseclass_mean_df)
  print(Neu_score_result_diseaseclass_mean_df)
  dim(Neu_score_result_diseaseclass_mean_df)
  
  output_f_path = file.path(output_dir,
                            "STAT-GSE_Rice_Neu_GSF3_Pathway_Score_mean_20240416_1.xls")  #
  write.table(
    Neu_score_result_diseaseclass_mean_df,
    file = output_f_path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
  
  Neu_score_result_mean_df = Neu_score_result_diseaseclass_mean_df
  head(Neu_score_result_mean_df)
  table(Neu_score_result_diseaseclass_mean_df$Average_GSF3_Pathway_Score)
  
  
  
}
# pathway ------------------------------------------------------------
{
  range(Neu_score_result_mean_df$Average_GSF3_Pathway_Score)
  range(Neu_score_result_mean_df$GSF3_Pathway_Score)
  
  Neu_score_result_mean_df$Average_GSF3_Pathway_Score
  class(Neu_score_result_mean_df$Average_GSF3_Pathway_Score)
  summary(Neu_score_result_mean_df$Average_GSF3_Pathway_Score)
  
  min_value = min(Neu_score_result_mean_df$Average_GSF3_Pathway_Score)
  max_value = max(Neu_score_result_mean_df$Average_GSF3_Pathway_Score)
  colors_list =  c("#3952A1",
                   "#57A0E4",
                   "#B5D3F7",
                   "#FBFCFE",
                   "#FF946D",
                   "#D62E1E",
                   "#D02823")
  breaks_list = seq(min_value, max_value, length.out = 7)
  labels_list = round(seq(min_value, max_value, length.out = 7), 3)
  limits_list = c(min_value - 0.002 , max_value + 0.002)
  
  seurat_object_plot =
    ggviolin(
      data = Neu_score_result_mean_df,
      x = "DiseaseClass",
      y = "GSF3_Pathway_Score",
      fill = "Average_GSF3_Pathway_Score",
      show.legend = TRUE
    ) +
    #scale_fill_gradientn(colors = viridis(7))
    scale_fill_gradientn(
      colors = colors_list,
      breaks = breaks_list,
      labels = labels_list,
      limits = limits_list,
      # 设置刻度范
      guide = guide_colorbar(
        barwidth = 15,
        barheight = 1,
        title.position = "top"
      )
    ) +
    guides(
      fill = guide_colorbar(
        title = "Average",
        #
        direction = "vertical",
        #
        title.position = "top",
        #
        title.hjust = 0.5
      )
    ) +           #
    theme(
      legend.position = "right",
      #
      legend.box = "vertical",
      #
      legend.justification = "center",
      #
      legend.margin = margin(
        t = 0,
        r = 10,
        b = 0,
        l = 0
      )
    ) +  #
    
    
    geom_boxplot(
      width = 0.15,
      outlier.color = NA,
      fill = NA,
      alpha = 1,
      linewidth = 0.5
    ) +
    
    geom_signif(
      comparisons = list(
        c("Normal", "MDS"),
        c("Normal", "tAML"),
        c("Normal", "AML"),
        c("MDS", "tAML")
      ),
      map_signif_level = T,
      test = wilcox.test,
      #
      step_increase = 0.08,
      size = 0.2,
      textsize = 5
    ) +
    labs(
      x = "Disease Class",
      y = "Module score",
      title = "Neu - Signaling by CSF3 (G-CSF) Pathway(30 genes)",
      fill = "Average"
    ) +  # 添加标签和标题
    theme(legend.position = "right")  #
  
  seurat_object_plot
  
  output_f_path = file.path(output_dir,
                            "QC-GSE_Rice_Neu_GSF3_Pathway_Score_20240416_1_3_1.pdf")
  ggsave(
    filename = output_f_path,
    plot = seurat_object_plot,
    device = pdf,
    width = 210,
    height = 148,
    units = "mm"
  )
  
}

# -------------------------------------------------------------------------
# Sample Average Score -------------------------------------------------------------------------
{
  #
  library(dplyr)
  score_result_diseaseclass_mean_df = Neu_score_result_mean_df
  
  score_result_sample_id_mean_df <- score_result_diseaseclass_mean_df %>%
    group_by(DiseaseClass, Sample_ID) %>%
    mutate(CellType_Sample_ID_Average_Pathway_Score = mean(GSF3_Pathway_Score))
  
  head(score_result_sample_id_mean_df)
  dim(score_result_sample_id_mean_df)
  table(score_result_sample_id_mean_df$Sample_ID)
  unique(score_result_sample_id_mean_df$Sample_ID)
  table(score_result_sample_id_mean_df$CellType_Sample_ID_Average_Pathway_Score)
  colnames(score_result_sample_id_mean_df)
  
  score_result_sample_id_mean_keep_df = score_result_sample_id_mean_df
  
  
  
  
}



# Statistical chart of inter-group differences in sample means  ------------------------------------------------------------
{
  head(score_result_sample_id_mean_keep_df)
  score_result_sample_id_mean_keep_df_select =  score_result_sample_id_mean_keep_df %>% select(Sample_ID,
                                                                                               DiseaseClass,
                                                                                               CellType_Sample_ID_Average_Pathway_Score)
  score_result_sample_id_mean_keep_df_select
  
  score_result_sample_id_mean_keep_df_deldup = score_result_sample_id_mean_keep_df_select %>% distinct()
  score_result_sample_id_mean_keep_df_deldup
  
  output_f_path = file.path(output_dir,
                            "STAT-GSE_Rice_Neu_Score_mean_Sample_ID_20250618_1.xls")
  write.table(
    score_result_sample_id_mean_keep_df_deldup,
    file = output_f_path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
}


{
  library(ggpubr)
  library(RColorBrewer)
  
  disease_colors = c(
    "CD4_T_01_FYB1" = "#e37fae",
    "CD4_T_02_CCL5" = "#19977f",
    "CD4_T_03_MT_ND4L" = "#808bac",
    "CD4_T_04_ISG15" = "#4da54a",
    "CD4_T_05_APOOL" = "#ed7a1c",
    "CD4_Treg_01_ITGB1" = "#379ec1",
    "CD8_T_01_NBEAL1" = "#d74b33",
    "CD8_T_02_LMNA" = "#3b5181",
    "CD8_T_03_HBA2" = "#8bc7ba",
    "CD8_T_04_CD8B" = "#9d552a",
    "CD8_T_05_CCL4" = "#8e4b97",
    "CD8_T_06_CCR7" = "#e7957b",
    "CD8_T_07_S100A11" = "#3976ab",
    "CD8_T_08_ATP5MC2" = "#d5211d",
    "CD8_T_09_MT_ND2" = "#00ffff",
    "CD8_T_10_EEF1G_PCNP" = "#ff00ff",
    "NK_01_IL7R" = "#d74b33",
    "NK_02_MT_CYB" = "#3b5181",
    "NK_03_GZMK" = "#8bc7ba",
    "NK_04_GZMH" = "#9d552a",
    "NK_05_CD3D" = "#8e4b97",
    "NK_06_SPON2" = "#e7957b"
  )
  
  {
    table(score_result_sample_id_mean_keep_df$DiseaseClass)
    disease_colors = c(
      "Normal" = "#4da54a" ,
      "MDS" = "#3b5181" ,
      "tAML" = "#8e4b97" ,
      "AML" = "#d5211d"
    )
    
  }
}


# -------------------------------------------------------------------------


# -------------------------------------------------------------------------

{
  custom_theme <- theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5),
          axis.ticks.x = element_blank(),
          # axis.text.x = element_blank(),)
          library(ggsci)
          
          # step
          seurat_object_plot = ggplot(
            data = score_result_sample_id_mean_keep_df_deldup,
            aes(x = DiseaseClass, y = CellType_Sample_ID_Average_Pathway_Score, color =
                  DiseaseClass)
          ) +
            geom_boxplot(outlier.size = 0.2) +
            stat_boxplot(geom = 'errorbar', width = 0.2) +
            scale_color_manual(values = disease_colors, name = "Group") +   # Apply custom colors to ellipses and other elements
            scale_fill_manual(values = disease_colors, name = "Group") +   # Apply custom colors to ellipses and other elements
            custom_theme +
            theme(axis.text.x = element_text(
              angle = 45,
              hjust = 1,
              size = 7
            )) + #
            labs(
              x = NULL,
              y = "Module score",
              title = "GSF3 Score in Neu cells",
              fill = "Average"
            )   #
          seurat_object_plot
          
          output_f_path = file.path(output_dir,
                                    "20250829_1_3_scRNA_GSE_Neu_GSF3_Score_diseaseclass_mean.pdf")
          ggsave(
            filename = output_f_path,
            plot = seurat_object_plot,
            device = pdf,
            width = 120,
            height = 90,
            units = "mm"
          )
}