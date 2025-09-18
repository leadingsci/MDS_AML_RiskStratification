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
  output_dir = "/share/home/wyz/Project/result/2025/202502/20250206_1_unrice_scRNA_TNK_exp/20250207_2_scRNA_GSE_TNK_module_score_TNK"
  
}


if (TRUE) {
  input_f_path_GSE = "/share/home/wyz/Project/result/20240429_scRNA_GSE_TNK/20240506_scRNA_GSE_CellPhoneDB/scRNA_SMU_GSE116256_PRJNA720840_GSE205490_HRA000084_AML_MDS_84_origin_rds_combined_qc_STD_Harmony_res1_select_20240507.rds"
  
  # 读入
  seurat_object_STD <- readRDS(input_f_path_GSE)
  levels(seurat_object_STD)
  table(seurat_object_STD@meta.data$Sample_ID)
  head(seurat_object_STD@meta.data)
  table(seurat_object_STD@meta.data$GSE116256_TransferLabel_merge_new)
}



{
  ## 选择Neu子集
  seurat_object_keep = subset(x = seurat_object_STD, subset = GSE116256_TransferLabel_merge_new %in% c("CD4 T", "CD8 T", "NK"))
  levels(seurat_object_keep)
  table(seurat_object_keep@meta.data$GSE116256_TransferLabel_merge_new)
}



# 获得基因列表 ------------------------------------------------------------
## 读取表格
if (FALSE) {
  # go_f_path = "/share/home/wyz/Lab/Code/ProjectCode/2024/202405/sample_info/Immunosuppression-related_genes.xlsx"
  # go_1_df = readxl::read_excel(go_f_path,col_names = TRUE)
  #
  # # 近保留在seurat中存在的基因 --------------------------------------------------------
  # go_f_1_df_keep = go_1_df[go_1_df$GeneName %in% rownames(seurat_object_keep@assays$RNA@counts),]
  
}

## 读取文本
if (TRUE) {
  go_f_1 = "/share/home/wyz/Lab/Code/ProjectCode/2024/202405/sample_info/GO0043065.txt"
  go_1_df = read.table(file = go_f_1,
                       sep = "\t",
                       header = FALSE)
  
  # keep gene  --------------------------------------------------------
  go_f_1_df_keep = go_1_df[go_1_df$V1 %in% rownames(seurat_object_keep@assays$RNA@counts), ]
  
}
head(go_1_df)
dim(go_1_df) # [1] 725  16
head(go_f_1_df_keep)
class(go_f_1_df_keep)
dim(go_f_1_df_keep) # [1] 703  16


# add Module分数 ------------------------------------------------------------

# # AddModuleScore 评分 ---------------------------------------------------
# 分开每个细胞类型来做
seurat_object_keep_module =
  AddModuleScore(
    object = seurat_object_keep,
    features = list(go_f_1_df_keep$V1),
    #
    ctrl = 100,
    # default
    name = 'Apoptosis_Score',
    #
    seed = 1,
  )

# 查看 Module分数 ------------------------------------------------------------

head(seurat_object_keep_module@meta.data)

# save ------------------------------------------------------------

#  DataFrame
score_result_cbind <- cbind(
  Sample_ID = seurat_object_keep_module@meta.data$Sample_ID,
  DiseaseClass = seurat_object_keep_module@meta.data$DiseaseClass,
  RNA_snn_res.1 = seurat_object_keep_module@meta.data$RNA_snn_res.1,
  CellType  =  as.character(
    seurat_object_keep_module@meta.data$GSE116256_TransferLabel_merge_new
  ),
  Pathway_Score = seurat_object_keep_module@meta.data$Apoptosis_Score
) # 需要更改通路

head(score_result_cbind)

#
score_result_cbind_df <- as.data.frame(score_result_cbind)
head(score_result_cbind_df)

output_f_path = file.path(output_dir, "scRNA_GSE_TNK_Apoptosis_Score_20240509_1.xls") # 需要更改通路
write.table(
  score_result_cbind_df,
  file = output_f_path,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)


{
  # 20250205 rename -------------------------------------------------------------------------
  score_result_cbind_df$DiseaseClass[score_result_cbind_df$DiseaseClass == "Normal"] <- "HC"
  score_result_cbind_df$DiseaseClass[score_result_cbind_df$DiseaseClass == "tAML"] <- "sAML"
  head(score_result_cbind_df)
  table(score_result_cbind_df$DiseaseClass)
}



# data  ------------------------------------------------------------
{
  # factor ------------------------------------------------------------
  table(score_result_cbind_df$DiseaseClass)
  score_result_cbind_df$DiseaseClass = factor(score_result_cbind_df$DiseaseClass,
                                              levels = c("HC", "MDS", "sAML", "AML"))
  levels(score_result_cbind_df$DiseaseClass)
  
  # as.numeric
  head(score_result_cbind_df)
  score_result_cbind_df$Pathway_Score = as.numeric(score_result_cbind_df$Pathway_Score)
  
  score_result_diseaseclass_mean_df <- score_result_cbind_df %>%
    group_by(DiseaseClass) %>%
    mutate(Average_Pathway_Score = mean(Pathway_Score))
  
  head(score_result_diseaseclass_mean_df)
  dim(score_result_diseaseclass_mean_df)
  
  output_f_path = file.path(
    output_dir,
    "scRNA_GSE_TNK_CD4_Immunosuppression_Score_diseaseclass_mean_20240509_1.xls"
  )  #
  write.table(
    score_result_diseaseclass_mean_df,
    file = output_f_path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
  
  #
  table(score_result_diseaseclass_mean_df$Average_Pathway_Score)
  
  
}

# -------------------------------------------------------------------------



#
## 20250207
{
  # CNS ggplot
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


# single cell pathway  ------------------------------------------------------------
{
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
  limits_list = c(min_value - 0.0002 , max_value + 0.0002)
  
  
  seurat_object_plot =
    ggviolin(
      data = score_result_diseaseclass_mean_df,
      x = "DiseaseClass",
      y = "Pathway_Score",
      fill = "Average_Pathway_Score",
      show.legend = TRUE,
      size = 0.3
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
      width = 0.1,
      outlier.color = NA,
      fill = "white",
      alpha = 1,
      linewidth = 0.1
    ) +
    
    geom_signif(
      comparisons = list(
        c("HC", "MDS"),
        c("MDS", "sAML"),
        c("sAML", "AML"),
        c("HC", "sAML"),
        c("MDS", "AML"),
        c("HC", "AML")
      ),
      map_signif_level = T,
      test = wilcox.test,
      #
      step_increase = 0.03,
      size = 0.2,
      tip_length = 0.01,
      textsize = 1.5
    ) +
    labs(
      x = NULL,
      y = "Module score",
      title = "Apoptosis Score (703 genes) in T/NK cells",
      fill = "Average"
    ) +  # 添加标签和标题
    theme(legend.position = "right") +  #
    cns.plot.format +
    theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.justification = "center",
      legend.margin = margin(
        t = 0,
        r = 10,
        b = 0,
        l = 0
      ),
      legend.key.width = unit(4, "mm"),
      #
      legend.key.height = unit(5, "mm")  #
    )
  
  seurat_object_plot
  
  output_f_path = file.path(output_dir,
                            "scRNA_GSE_TNK_Apoptosis_Score_diseaseclass_mean_20240509_1.pdf")
  ggsave(
    filename = output_f_path,
    plot = seurat_object_plot,
    device = pdf,
    width = 80,
    height = 50,
    units = "mm"
  )
  
}

# -------------------------------------------------------------------------

# -------------------------------------------------------------------------


# -------------------------------------------------------------------------

# Sample Average Score -------------------------------------------------------------------------
{
  library(dplyr)
  score_result_sample_id_mean_df <- score_result_diseaseclass_mean_df %>%
    group_by(Sample_ID) %>%
    mutate(Sample_ID_Average_Pathway_Score = mean(Pathway_Score))
  
  head(score_result_sample_id_mean_df)
  dim(score_result_sample_id_mean_df)
  table(score_result_sample_id_mean_df$Sample_ID)
  unique(score_result_sample_id_mean_df$Sample_ID)
  table(score_result_sample_id_mean_df$Sample_ID_Average_Pathway_Score)
  colnames(score_result_sample_id_mean_df)
  
  score_result_sample_id_mean_keep_df = score_result_sample_id_mean_df
}

# dotplot  ------------------------------------------------------------
{
  head(score_result_sample_id_mean_keep_df)
  score_result_sample_id_mean_keep_df_select = select(
    score_result_sample_id_mean_keep_df,
    Sample_ID,
    DiseaseClass,
    Sample_ID_Average_Pathway_Score
  )
  score_result_sample_id_mean_keep_df_deldup = score_result_sample_id_mean_keep_df_select %>% distinct()
  score_result_sample_id_mean_keep_df_deldup
  
  output_f_path = file.path(output_dir,
                            "STAT-GSE_Rice_TNK_APOPTOSIS_Score_mean_Sample_ID_20240521_1.xls")
  write.table(
    score_result_sample_id_mean_keep_df_deldup,
    file = output_f_path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
}


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------

{
  library(ggpubr)
  library(RColorBrewer)
}


{
  # Define a custom color palette for the DiseaseClass levels
  disease_colors <- c(
    "HC" = "#4da54a",
    #
    "MDS" = "#3b5181",
    # Red
    "sAML" = "#8e4b97",
    # Green
    "AML" = "#d5211d"   # Blue
  )
}


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

{
  #
  seurat_object_plot = ggplot(
    data = score_result_sample_id_mean_keep_df_deldup,
    aes(x = DiseaseClass, y = Sample_ID_Average_Pathway_Score, color =
          DiseaseClass)
  ) +
    geom_boxplot(size = 0.2, outlier.shape = NA) +
    
    # 加入errorbar
    stat_boxplot(geom = 'errorbar', width = 0.2) +
    geom_jitter(
      aes(color = DiseaseClass),
      size = 0.2,
      alpha = 1,
      width = 0.1
    ) +
    scale_color_manual(values = disease_colors, name = "Group") +   # Apply custom colors to ellipses and other elements
    scale_fill_manual(values = disease_colors, name = "Group") +   # Apply custom colors to ellipses and other elements
    
    # custom_theme +
    labs(
      x = NULL,
      y = "Module score",
      title = "Apoptosis Score (703 genes) in T/NK cells",
      fill = "Group"
    ) +  #
    
    theme(# legend.position = "none",
      axis.line = element_line(color = "black"),
      #
      axis.ticks = element_line(color = "black")          # ) +
      
      
      cns.plot.format +
        
        geom_signif(
          comparisons = list(
            c("HC", "MDS"),
            c("MDS", "sAML"),
            c("sAML", "AML"),
            c("HC", "sAML"),
            c("MDS", "AML"),
            c("HC", "AML")
          ),
          map_signif_level = T,
          test = wilcox.test,
          #
          step_increase = 0.03,
          size = 0.2,
          tip_length = 0.01,
          textsize = 1.5,
          color = "black",
          #
          annotations = NULL  #
        )
      
      seurat_object_plot
      
      output_f_path = file.path(output_dir,
                                "scRNA_GSE_TNK_Apoptosis_Score_diseaseclass_mean_20240523_1.pdf")
      ggsave(
        filename = output_f_path,
        plot = seurat_object_plot,
        device = pdf,
        width = 75,
        height = 40,
        units = "mm"
      )
      
}

# -------------------------------------------------------------------------
{
  #
  custom_theme <- theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, vjust = 0.5),
          axis.ticks.x = element_blank(),)
  
  #
  seurat_object_plot = ggplot(
    data = score_result_sample_id_mean_keep_df_deldup,
    aes(x = DiseaseClass, y = Sample_ID_Average_Pathway_Score, color = DiseaseClass)
  ) +
    geom_boxplot() +
    #
    stat_boxplot(geom = 'errorbar', width = 0.3) +
    #
    geom_jitter(
      aes(color = DiseaseClass),
      size = 3,
      alpha = 0.9,
      width = 0.1
    ) +
    #
    scale_color_npg() +
    
    geom_signif(
      comparisons = list(
        c("Normal", "MDS"),
        c("Normal", "tAML"),
        c("Normal", "AML"),
        c("MDS", "tAML"),
        c("tAML", "AML"),
        c("MDS", "AML")
      ),
      map_signif_level = F,
      test = wilcox.test,
      #
      step_increase = 0.08,
      size = 0.2,
      color = "black",
      textsize = 5
    ) +
    custom_theme +
    labs(
      x = NULL,
      y = "Module score",
      title = "Apoptosis Score (703 genes) in T/NK cells",
      fill = "Average"
    )  #
  seurat_object_plot
  
  output_f_path = file.path(
    output_dir,
    "QC-scRNA_GSE_TNK_Apoptosis_Score_diseaseclass_mean_20240523_2.pdf"
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