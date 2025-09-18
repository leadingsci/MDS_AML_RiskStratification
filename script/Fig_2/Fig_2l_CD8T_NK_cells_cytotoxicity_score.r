.libPaths()

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
library(ggpubr)
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
  output_dir = "/share/home/wyz/Project/result/20240509_scRNA_GSE_TNK/20240509_scRNA_GSE_TNK_module_score"
  
}


if (TRUE) {
  input_f_path_GSE = "/share/home/wyz/Project/result/20240429_scRNA_GSE_TNK/20240506_scRNA_GSE_CellPhoneDB/scRNA_SMU_GSE116256_PRJNA720840_GSE205490_HRA000084_AML_MDS_84_origin_rds_combined_qc_STD_Harmony_res1_select_20240507.rds"
  
  seurat_object_STD <- readRDS(input_f_path_GSE)
  levels(seurat_object_STD)
  table(seurat_object_STD@meta.data$Sample_ID)
  head(seurat_object_STD@meta.data)
  table(seurat_object_STD@meta.data$GSE116256_TransferLabel_merge_new)
}



{
  ## selecet
  seurat_object_keep = subset(x = seurat_object_STD, subset = GSE116256_TransferLabel_merge_new == c("CD8 T", "NK"))
  levels(seurat_object_keep)
  table(seurat_object_keep@meta.data$GSE116256_TransferLabel_merge_new)
}



# gene lit ------------------------------------------------------------
go_f_path = "/share/home/wyz/Lab/Code/ProjectCode/2024/202405/sample_info/cytotoxicity-associated genes.xlsx"
go_1_df = readxl::read_excel(go_f_path, col_names = TRUE)

head(go_1_df)
dim(go_1_df)

# filter gene list --------------------------------------------------------
go_f_1_df_keep = go_1_df[go_1_df$GeneName %in% rownames(seurat_object_keep@assays$RNA@counts), ]
head(go_f_1_df_keep)
dim(go_f_1_df_keep)
class(go_f_1_df_keep)



# add Module------------------------------------------------------------

# # AddModuleScore  ---------------------------------------------------
seurat_object_keep_module =
  AddModuleScore(
    object = seurat_object_keep,
    features = list(go_f_1_df_keep$GeneName),
    #
    ctrl = 100,
    # 默认
    name = 'Cytotoxicity_Scores',
    #
    seed = 1,
  )

#  Module ------------------------------------------------------------

head(seurat_object_keep_module@meta.data)

# df ------------------------------------------------------------
score_result_cbind <- cbind(
  Sample_ID = seurat_object_keep_module@meta.data$Sample_ID,
  DiseaseClass = seurat_object_keep_module@meta.data$DiseaseClass,
  RNA_snn_res.1 = seurat_object_keep_module@meta.data$RNA_snn_res.1,
  CellType  =  as.character(
    seurat_object_keep_module@meta.data$GSE116256_TransferLabel_merge_new
  ),
  Pathway_Score = seurat_object_keep_module@meta.data$Cytotoxicity_Scores1
) #

head(score_result_cbind)

score_result_cbind_df <- as.data.frame(score_result_cbind)
head(score_result_cbind_df)

output_f_path = file.path(output_dir,
                          "scRNA_GSE_CD8_NK_Cytotoxicity_Scores_20240509_1.xls") #
write.table(
  score_result_cbind_df,
  file = output_f_path,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)


#  --------------------------------------------------------------------

# data ------------------------------------------------------------
{
  # factor ------------------------------------------------------------
  table(score_result_cbind_df$DiseaseClass)
  score_result_cbind_df$DiseaseClass = factor(score_result_cbind_df$DiseaseClass,
                                              levels = c("Normal", "MDS", "tAML", "AML"))
  levels(score_result_cbind_df$DiseaseClass)
  
  head(score_result_cbind_df)
  score_result_cbind_df$Pathway_Score = as.numeric(score_result_cbind_df$Pathway_Score)
  
  score_result_diseaseclass_mean_df <- score_result_cbind_df %>%
    group_by(DiseaseClass) %>%
    mutate(Average_Pathway_Score = mean(Pathway_Score))
  
  head(score_result_diseaseclass_mean_df)
  dim(score_result_diseaseclass_mean_df)
  
  output_f_path = file.path(
    output_dir,
    "scRNA_GSE_CD8_NK_Cytotoxicity_Score_diseaseclass_mean_20240509_1.xls"
  )  # 需要修改
  write.table(
    score_result_diseaseclass_mean_df,
    file = output_f_path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
  
  table(score_result_diseaseclass_mean_df$Average_Pathway_Score)
  
}
------------------------------------------------------------{
  range(score_result_diseaseclass_mean_df$Average_Pathway_Score)
  range(score_result_diseaseclass_mean_df$Pathway_Score)
  
  summary(score_result_diseaseclass_mean_df$Average_Pathway_Score)
  
  min_value = min(score_result_diseaseclass_mean_df$Average_Pathway_Score)
  max_value = max(score_result_diseaseclass_mean_df$Average_Pathway_Score)
  
  colors_list =  c("#3952A1",
                   "#57A0E4",
                   "#B5D3F7",
                   "#FBFCFE",
                   "#FF946D",
                   "#D62E1E",
                   "#D02823")
  breaks_list = seq(min_value, max_value, length.out = 5)
  labels_list = round(seq(min_value, max_value, length.out = 5), 3)
  limits_list = c(min_value - 0.0001 , max_value + 0.0001)
  
  
  seurat_object_plot =
    ggviolin(
      data = score_result_diseaseclass_mean_df,
      x = "DiseaseClass",
      y = "Pathway_Score",
      fill = "Average_Pathway_Score",
      show.legend = TRUE
    ) +
    scale_fill_gradientn(
      colors = colors_list,
      breaks = breaks_list,
      labels = labels_list,
      limits = limits_list,
      #
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
      width = 0.05,
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
      title = "Cytotoxicity Score(4 genes) in CD8 T/NK cells",
      fill = "Average"
    ) +  #
    theme(legend.position = "right")  #
  
  seurat_object_plot
  
  output_f_path = file.path(
    output_dir,
    "scRNA_GSE_CD8_NK_Cytotoxicity_Score_diseaseclass_mean_20240509_1.pdf"
  )
  ggsave(
    filename = output_f_path,
    plot = seurat_object_plot,
    device = pdf,
    width = 210,
    height = 148,
    units = "mm"
  )
  
}