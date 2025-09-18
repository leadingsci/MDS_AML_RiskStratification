.libPaths(
  c(
    "/share/home/bioinfo/anaconda3/envs/R-4.2.3/lib/R/library",
    "/share/apps/Anaconda3/2021.05/envs/R-4.2/lib/R/library",
    "/share/home/wyz/miniconda3/envs/r423base/lib/R/library",
    "/share/home/wyz/miniconda3/envs/r423base_ipynb/lib/R/library"
  )
)
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
  output_dir = "/share/home/wyz/Project/result/20240509_scRNA_GSE_TNK/20240516_scRNA_GSE_Rice_summary/20240516_scRNA_GSE_Rice_HSC_redar_res0.2"
  
}


if (TRUE) {
  # read
  input_f_path = "/share/home/wyz/Project/result/20240401_scRNA_AML_Annotation/20240409_scRNA_GSE_Rice_HSC_Annotation/scRNA_SMU_Rice_14_GSE116256_Rice_1_HRA000084_Rice_4_origin_rds_combined_qc_STD_Harmony_res_select_HSC_harmony_res.rds"
  
  seurat_object_STD <- readRDS(input_f_path)
  print(seurat_object_STD)
  levels(seurat_object_STD)
  
  
}

{
  ## select
  seurat_object_STD_HSC = subset(seurat_object_STD, subset = CellType_HSC_20240409_2 %in% c("HSC/MPP1", "HSC/MPP2"))
  levels(seurat_object_STD_HSC)
}

{
  # res 0.2
  Idents(object = seurat_object_STD_HSC) <- "CellType_HSC_20240409_2"
  levels(seurat_object_STD_HSC)
  
}


{
  # 保存
  ## 保存rds
  if (TRUE) {
    # output
    file_name = "scRNA_SMU_Rice_14_GSE116256_Rice_1_HRA000084_Rice_4_origin_rds_combined_qc_STD_Harmony_res_select_HSC_harmony_res_20240516.rds"
    
    saveRDS(seurat_object_STD_HSC, file.path(output_dir, file_name))
    
    # read
    #seurat_object_STD = readRDS(file.path(output_dir, file_name))
  }
  
}

# -------------------------------------------------------------------------

# 画图
## ummap图
if (TRUE) {
  col_list = c(
    "#d20000",
    "#e3b8a5",
    "#3b62a4",
    "#5caf34",
    "#907366",
    "#87c9c1",
    "#889b22",
    "#a855cb",
    "#b9617f"
  )
  
  seurat_object_plot =
    DimPlot(
      seurat_object_STD_HSC,
      reduction = "std.umap.harmony.hsc",
      group.by = c("CellType_HSC_20240409_2"),
      # split.by =c("Source"),
      combine = TRUE,
      label.size = 5,
      pt.size = 0.1,
      label = TRUE,
      ncol = 1,
      cols = col_list[1:2],
      raster = FALSE
    ) +
    ggtitle("UAMP")
  
  seurat_object_plot
  
  output_f_path = file.path(
    output_dir,
    "QC-GSE_Rice_CellType_Rice_CellType_HSC_20240409_2_umap_20240516_1.pdf"
  )
  ggsave(
    filename = output_f_path,
    plot = seurat_object_plot,
    device = pdf,
    width = 168,
    height = 148,
    units = "mm"
  )
  
}


# -------------------------------------------------------------------------


## read
gene_sheet_f =  "/share/home/wyz/Lab/Code/ProjectCode/2024/202404/sample_info/HSC_stemness_DatafileS4_GSEA_gene_sets.xlsx"

gene_sheet_df = readxl::read_excel(
  path = gene_sheet_f,
  col_name = T,
  col_types = "text",
  na = "",
  sheet = "Gene_symbol"
)
dim(gene_sheet_df)
print(gene_sheet_df)
class(gene_sheet_df)

## 5 gene set
Stemness_gene_df <- subset(gene_sheet_df, Type == "Stemness")
Erythroid_gene_df <- subset(gene_sheet_df, Type == "Erythroid")
Megakaryocytic_gene_df <- subset(gene_sheet_df, Type == "Megakaryocytic")
Myeloid_gene_df <- subset(gene_sheet_df, Type == "Myeloid")
Lymphoid_gene_df <- subset(gene_sheet_df, Type == "Lymphoid")

print(Stemness_gene_df)
print(Stemness_gene_df$Gene_symbol)


# # 确保基因列表存在于对象中 ---------------------------------------------------
Stemness_gene_df_keep = Stemness_gene_df[Stemness_gene_df$Gene_symbol %in% rownames(seurat_object_STD_HSC@assays$RNA@counts), ]
Erythroid_gene_df_keep = Erythroid_gene_df[Erythroid_gene_df$Gene_symbol %in% rownames(seurat_object_STD_HSC@assays$RNA@counts), ]
Megakaryocytic_gene_df_keep = Megakaryocytic_gene_df[Megakaryocytic_gene_df$Gene_symbol %in% rownames(seurat_object_STD_HSC@assays$RNA@counts), ]
Myeloid_gene_df_keep = Myeloid_gene_df[Myeloid_gene_df$Gene_symbol %in% rownames(seurat_object_STD_HSC@assays$RNA@counts), ]
Lymphoid_gene_df_keep = Lymphoid_gene_df[Lymphoid_gene_df$Gene_symbol %in% rownames(seurat_object_STD_HSC@assays$RNA@counts), ]

dim(Stemness_gene_df_keep)



# # AddModuleScore 干性评分 ---------------------------------------------------
# 分开每个细胞类型来做
seurat_object_STD_HSC =
  AddModuleScore(
    object = seurat_object_STD_HSC,
    features = list(Stemness_gene_df_keep$Gene_symbol),
    ctrl = 100,
    # 默认
    name = 'Stemness_Score',
    seed = 1,
  )

seurat_object_STD_HSC =
  AddModuleScore(
    object = seurat_object_STD_HSC,
    features = list(Erythroid_gene_df_keep$Gene_symbol),
    ctrl = 100,
    # 默认
    name = 'Erythroid_Score',
    seed = 1,
  )


seurat_object_STD_HSC =
  AddModuleScore(
    object = seurat_object_STD_HSC,
    features = list(Megakaryocytic_gene_df_keep$Gene_symbol),
    ctrl = 100,
    # 默认
    name = 'Megakaryocytic_Score',
    seed = 1,
  )

seurat_object_STD_HSC =
  AddModuleScore(
    object = seurat_object_STD_HSC,
    features = list(Myeloid_gene_df_keep$Gene_symbol),
    ctrl = 100,
    # 默认
    name = 'Myeloid_Score',
    seed = 1,
  )

seurat_object_STD_HSC =
  AddModuleScore(
    object = seurat_object_STD_HSC,
    features = list(Lymphoid_gene_df_keep$Gene_symbol),
    ctrl = 100,
    # 默认
    name = 'Lymphoid_Score',
    seed = 1,
  )


seurat_object_STD_HSC@meta.data$Stemness_Score1


# save df
# row : cell name
# col : Stemness type
# value : score

# 创建 DataFrame
hsc_score_result_cbind <- cbind(
  Stemness = seurat_object_STD_HSC@meta.data$Stemness_Score1,
  Erythroid = seurat_object_STD_HSC@meta.data$Erythroid_Score1,
  Megakaryocytic = seurat_object_STD_HSC@meta.data$Megakaryocytic_Score1,
  Myeloid = seurat_object_STD_HSC@meta.data$Myeloid_Score1,
  Lymphoid = seurat_object_STD_HSC@meta.data$Lymphoid_Score1
)

head(hsc_score_result_cbind)

# 转为数据框
hsc_score_result_df <- as.data.frame(hsc_score_result_cbind)
head(hsc_score_result_df[0:5, 0:5])

# 添加列名
# 自定义行名
dim(hsc_score_result_df)
length(as.vector(colnames(seurat_object_STD_HSC)))
# 长度相等

rownames(hsc_score_result_df) <- as.vector(colnames(seurat_object_STD_HSC))
head(hsc_score_result_df[0:5, 0:5])


hsc_score_result_df <- hsc_score_result_df %>%
  rownames_to_column(var = "Cell_ID")
hsc_score_result_df


output_f_path = file.path(output_dir, "STAT-GSE_Rice_HSC_stem_Sorce_20240516.xls")
write.table(
  hsc_score_result_df,
  file = output_f_path,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)



## Calculate the average score of the differentiation preference of each cell for each subgroup
head(seurat_object_STD_HSC@meta.data)
seurat_object_STD_HSC_score_df = seurat_object_STD_HSC@meta.data[, c(
  "CellType_HSC_20240409_2",
  "Stemness_Score1",
  "Erythroid_Score1",
  "Megakaryocytic_Score1",
  "Myeloid_Score1",
  "Lymphoid_Score1"
)]
head(seurat_object_STD_HSC_score_df)

# Specify the classification of a certain column and calculate the mean value for each classification.
seurat_object_STD_HSC_mean_df <- seurat_object_STD_HSC_score_df %>%
  group_by(CellType_HSC_20240409_2) %>%
  summarise_all(mean)

print(seurat_object_STD_HSC_mean_df)

# rename
colnames(seurat_object_STD_HSC_mean_df) = c("subClass",
                                            "Stemness",
                                            "Erythroid",
                                            "Megakaryocytic",
                                            "Myeloid",
                                            "Lymphoid")

output_f_path = file.path(output_dir,
                          "STAT-GSE_Rice_HSC_stem_Sorce_mean_20240516.xls")
write.table(
  seurat_object_STD_HSC_mean_df,
  file = output_f_path,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)

# Build radar chart data tables for each disease group ---------------------------------------------------------------------
seurat_object_STD_HSC_score_group_df = seurat_object_STD_HSC@meta.data[, c(
  "CellType_HSC_20240409_2",
  "DiseaseClass",
  "Stemness_Score1",
  "Erythroid_Score1",
  "Megakaryocytic_Score1",
  "Myeloid_Score1",
  "Lymphoid_Score1"
)]
head(seurat_object_STD_HSC_score_group_df)
dim(seurat_object_STD_HSC_score_group_df)

# Calculate the mean of the groups
seurat_object_STD_HSC_score_group_mean_df <-
  seurat_object_STD_HSC_score_group_df %>%
  group_by(DiseaseClass, CellType_HSC_20240409_2) %>%
  summarise_all(mean)

head(seurat_object_STD_HSC_score_group_mean_df)
print(seurat_object_STD_HSC_score_group_mean_df$CellType_HSC_20240409_2)
dim(seurat_object_STD_HSC_score_group_mean_df)

# 重命名列名
colnames(seurat_object_STD_HSC_score_group_mean_df) = c(
  "DiseaseClass",
  "SubClass",
  "Stemness",
  "Erythroid",
  "Megakaryocytic",
  "Myeloid",
  "Lymphoid"
)

output_f_path = file.path(output_dir,
                          "STAT-GSE_Rice_HSC_stem_Sorce_group_mean_20240516.xls")
write.table(
  seurat_object_STD_HSC_score_group_mean_df,
  file = output_f_path,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)


# dotplot ---------------------------------------------------------------------
{
  #  ---------------------------------------------------------------------
  # Group the scores according to the diseases
  head(seurat_object_STD_HSC@meta.data)
  seurat_object_STD_HSC_desease_score_df = seurat_object_STD_HSC@meta.data[, c(
    "DiseaseClass",
    "CellType_HSC_20240409_2",
    "Stemness_Score1",
    "Erythroid_Score1",
    "Megakaryocytic_Score1",
    "Myeloid_Score1",
    "Lymphoid_Score1"
  )]
  head(seurat_object_STD_HSC_desease_score_df)
  
  # 按疾病分组
  table(seurat_object_STD_HSC_desease_score_df$DiseaseClass)
  
  seurat_object_STD_HSC_Normal_score_df <- subset(x = seurat_object_STD_HSC_desease_score_df, DiseaseClass == "Normal")
  seurat_object_STD_HSC_MDS_score_df <- subset(x = seurat_object_STD_HSC_desease_score_df, DiseaseClass == "MDS")
  seurat_object_STD_HSC_tAML_score_df <- subset(x = seurat_object_STD_HSC_desease_score_df, DiseaseClass == "tAML")
  seurat_object_STD_HSC_AML_score_df <- subset(x = seurat_object_STD_HSC_desease_score_df, DiseaseClass == "AML")
  
  head(seurat_object_STD_HSC_Normal_score_df)
  dim(seurat_object_STD_HSC_Normal_score_df)
  
}

{
  # Specify the classification of a certain column and calculate the enrichment score and mean value for each classification.
  seurat_object_STD_HSC_Normal_score_mean_df <- seurat_object_STD_HSC_Normal_score_df[, -1] %>%
    group_by(CellType_HSC_20240409_2) %>%
    summarise_all(mean)
  
  seurat_object_STD_HSC_MDS_score_mean_df <- seurat_object_STD_HSC_MDS_score_df[, -1] %>%
    group_by(CellType_HSC_20240409_2) %>%
    summarise_all(mean)
  
  seurat_object_STD_HSC_tAML_score_mean_df <- seurat_object_STD_HSC_tAML_score_df[, -1] %>%
    group_by(CellType_HSC_20240409_2) %>%
    summarise_all(mean)
  
  seurat_object_STD_HSC_AML_score_mean_df <- seurat_object_STD_HSC_AML_score_df[, -1] %>%
    group_by(CellType_HSC_20240409_2) %>%
    summarise_all(mean)
  
  print(seurat_object_STD_HSC_Normal_score_mean_df)
  
  # rename
  colnames(seurat_object_STD_HSC_Normal_score_mean_df) = c("SubClass",
                                                           "Stemness",
                                                           "Erythroid",
                                                           "Megakaryocytic",
                                                           "Myeloid",
                                                           "Lymphoid")
  colnames(seurat_object_STD_HSC_MDS_score_mean_df) = c("SubClass",
                                                        "Stemness",
                                                        "Erythroid",
                                                        "Megakaryocytic",
                                                        "Myeloid",
                                                        "Lymphoid")
  colnames(seurat_object_STD_HSC_tAML_score_mean_df) = c("SubClass",
                                                         "Stemness",
                                                         "Erythroid",
                                                         "Megakaryocytic",
                                                         "Myeloid",
                                                         "Lymphoid")
  colnames(seurat_object_STD_HSC_AML_score_mean_df) = c("SubClass",
                                                        "Stemness",
                                                        "Erythroid",
                                                        "Megakaryocytic",
                                                        "Myeloid",
                                                        "Lymphoid")
  
  seurat_object_STD_HSC_Normal_score_mean_df
}



{
  ## Calculate the enrichment scores for the data frames of 4 diseases
  seurat_object_STD_HSC_Normal_score_mean_df$DiseaseClass =  rep("Normal", 2)
  seurat_object_STD_HSC_MDS_score_mean_df$DiseaseClass = rep("MDS", 2)
  seurat_object_STD_HSC_tAML_score_mean_df$DiseaseClass = rep("tAML", 2)
  seurat_object_STD_HSC_AML_score_mean_df$DiseaseClass = rep("AML", 2)
  
  seurat_object_STD_HSC_Normal_score_mean_df
  
  seurat_object_STD_HSC_Normal_score_mean_df$SubGroup = paste0(
    seurat_object_STD_HSC_Normal_score_mean_df$DiseaseClass,
    "_",
    seurat_object_STD_HSC_Normal_score_mean_df$SubClass
  )
  seurat_object_STD_HSC_MDS_score_mean_df$SubGroup = paste0(
    seurat_object_STD_HSC_MDS_score_mean_df$DiseaseClass,
    "_",
    seurat_object_STD_HSC_MDS_score_mean_df$SubClass
  )
  seurat_object_STD_HSC_tAML_score_mean_df$SubGroup = paste0(
    seurat_object_STD_HSC_tAML_score_mean_df$DiseaseClass,
    "_",
    seurat_object_STD_HSC_tAML_score_mean_df$SubClass
  )
  seurat_object_STD_HSC_AML_score_mean_df$SubGroup = paste0(
    seurat_object_STD_HSC_AML_score_mean_df$DiseaseClass,
    "_",
    seurat_object_STD_HSC_AML_score_mean_df$SubClass
  )
  
  seurat_object_STD_HSC_Normal_score_mean_df
  
  HSC_merge_score_mean_class_df <- rbind(
    seurat_object_STD_HSC_Normal_score_mean_df,
    seurat_object_STD_HSC_MDS_score_mean_df,
    seurat_object_STD_HSC_tAML_score_mean_df,
    seurat_object_STD_HSC_AML_score_mean_df
  )
  HSC_merge_score_mean_class_df
  
  
  output_f_path = file.path(output_dir,
                            "STAT-GSE_Rice_HSC_stem_Sorce_mean_merge_20240516.xls")
  write.table(
    HSC_merge_score_mean_class_df,
    file = output_f_path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
}


{
  # set
  set.seed(123)
  # Correct by column
  HSC_merge_score_mean_class_df_scale = apply(HSC_merge_score_mean_class_df[2:6], 2, scale)
  HSC_merge_score_mean_class_df_scale
  
  HSC_merge_score_mean_class_df_scale = cbind(HSC_merge_score_mean_class_df[, c(7, 1, 8)],
                                              HSC_merge_score_mean_class_df_scale)
  HSC_merge_score_mean_class_df_scale
  
  output_f_path = file.path(output_dir,
                            "STAT-GSE_Rice_HSC_stem_Sorce_mean_scale_merge_20240516.xls")
  write.table(
    HSC_merge_score_mean_class_df_scale,
    file = output_f_path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
}




{
  # Obtain the enrichment scores of each subgroup in the four groups
  Normal_HSC1_Stemness_Score = subset(seurat_object_STD_HSC_Normal_score_df,
                                      CellType_HSC_20240409_2 == "HSC/MPP1")$Stemness_Score1
  Normal_HSC2_Stemness_Score = subset(seurat_object_STD_HSC_Normal_score_df,
                                      CellType_HSC_20240409_2 == "HSC/MPP2")$Stemness_Score1
  
  MDS_HSC1_Stemness_Score = subset(seurat_object_STD_HSC_MDS_score_df,
                                   CellType_HSC_20240409_2 == "HSC/MPP1")$Stemness_Score1
  MDS_HSC2_Stemness_Score = subset(seurat_object_STD_HSC_MDS_score_df,
                                   CellType_HSC_20240409_2 == "HSC/MPP2")$Stemness_Score1
  
  tAML_HSC1_Stemness_Score = subset(seurat_object_STD_HSC_tAML_score_df,
                                    CellType_HSC_20240409_2 == "HSC/MPP1")$Stemness_Score1
  tAML_HSC2_Stemness_Score = subset(seurat_object_STD_HSC_tAML_score_df,
                                    CellType_HSC_20240409_2 == "HSC/MPP2")$Stemness_Score1
  
  AML_HSC1_Stemness_Score = subset(seurat_object_STD_HSC_AML_score_df,
                                   CellType_HSC_20240409_2 == "HSC/MPP1")$Stemness_Score1
  AML_HSC2_Stemness_Score = subset(seurat_object_STD_HSC_AML_score_df,
                                   CellType_HSC_20240409_2 == "HSC/MPP2")$Stemness_Score1
}

{
  Normal_HSC1_Erythroid_Score = subset(seurat_object_STD_HSC_Normal_score_df,
                                       CellType_HSC_20240409_2 == "HSC/MPP1")$Erythroid_Score1
  Normal_HSC2_Erythroid_Score = subset(seurat_object_STD_HSC_Normal_score_df,
                                       CellType_HSC_20240409_2 == "HSC/MPP2")$Erythroid_Score1
  MDS_HSC1_Erythroid_Score = subset(seurat_object_STD_HSC_MDS_score_df,
                                    CellType_HSC_20240409_2 == "HSC/MPP1")$Erythroid_Score1
  MDS_HSC2_Erythroid_Score = subset(seurat_object_STD_HSC_MDS_score_df,
                                    CellType_HSC_20240409_2 == "HSC/MPP2")$Erythroid_Score1
  tAML_HSC1_Erythroid_Score = subset(seurat_object_STD_HSC_tAML_score_df,
                                     CellType_HSC_20240409_2 == "HSC/MPP1")$Erythroid_Score1
  tAML_HSC2_Erythroid_Score = subset(seurat_object_STD_HSC_tAML_score_df,
                                     CellType_HSC_20240409_2 == "HSC/MPP2")$Erythroid_Score1
  AML_HSC1_Erythroid_Score = subset(seurat_object_STD_HSC_AML_score_df,
                                    CellType_HSC_20240409_2 == "HSC/MPP1")$Erythroid_Score1
  AML_HSC2_Erythroid_Score = subset(seurat_object_STD_HSC_AML_score_df,
                                    CellType_HSC_20240409_2 == "HSC/MPP2")$Erythroid_Score1
  
}

{
  Normal_HSC1_Megakaryocytic_Score = subset(seurat_object_STD_HSC_Normal_score_df,
                                            CellType_HSC_20240409_2 == "HSC/MPP1")$Megakaryocytic_Score1
  Normal_HSC2_Megakaryocytic_Score = subset(seurat_object_STD_HSC_Normal_score_df,
                                            CellType_HSC_20240409_2 == "HSC/MPP2")$Megakaryocytic_Score1
  MDS_HSC1_Megakaryocytic_Score = subset(seurat_object_STD_HSC_MDS_score_df,
                                         CellType_HSC_20240409_2 == "HSC/MPP1")$Megakaryocytic_Score1
  MDS_HSC2_Megakaryocytic_Score = subset(seurat_object_STD_HSC_MDS_score_df,
                                         CellType_HSC_20240409_2 == "HSC/MPP2")$Megakaryocytic_Score1
  tAML_HSC1_Megakaryocytic_Score = subset(seurat_object_STD_HSC_tAML_score_df,
                                          CellType_HSC_20240409_2 == "HSC/MPP1")$Megakaryocytic_Score1
  tAML_HSC2_Megakaryocytic_Score = subset(seurat_object_STD_HSC_tAML_score_df,
                                          CellType_HSC_20240409_2 == "HSC/MPP2")$Megakaryocytic_Score1
  AML_HSC1_Megakaryocytic_Score = subset(seurat_object_STD_HSC_AML_score_df,
                                         CellType_HSC_20240409_2 == "HSC/MPP1")$Megakaryocytic_Score1
  AML_HSC2_Megakaryocytic_Score = subset(seurat_object_STD_HSC_AML_score_df,
                                         CellType_HSC_20240409_2 == "HSC/MPP2")$Megakaryocytic_Score1
}

{
  Normal_HSC1_Myeloid_Score = subset(seurat_object_STD_HSC_Normal_score_df,
                                     CellType_HSC_20240409_2 == "HSC/MPP1")$Myeloid_Score1
  Normal_HSC2_Myeloid_Score = subset(seurat_object_STD_HSC_Normal_score_df,
                                     CellType_HSC_20240409_2 == "HSC/MPP2")$Myeloid_Score1
  MDS_HSC1_Myeloid_Score = subset(seurat_object_STD_HSC_MDS_score_df,
                                  CellType_HSC_20240409_2 == "HSC/MPP1")$Myeloid_Score1
  MDS_HSC2_Myeloid_Score = subset(seurat_object_STD_HSC_MDS_score_df,
                                  CellType_HSC_20240409_2 == "HSC/MPP2")$Myeloid_Score1
  tAML_HSC1_Myeloid_Score = subset(seurat_object_STD_HSC_tAML_score_df,
                                   CellType_HSC_20240409_2 == "HSC/MPP1")$Myeloid_Score1
  tAML_HSC2_Myeloid_Score = subset(seurat_object_STD_HSC_tAML_score_df,
                                   CellType_HSC_20240409_2 == "HSC/MPP2")$Myeloid_Score1
  AML_HSC1_Myeloid_Score = subset(seurat_object_STD_HSC_AML_score_df,
                                  CellType_HSC_20240409_2 == "HSC/MPP1")$Myeloid_Score1
  AML_HSC2_Myeloid_Score = subset(seurat_object_STD_HSC_AML_score_df,
                                  CellType_HSC_20240409_2 == "HSC/MPP2")$Myeloid_Score1
  
}

{
  Normal_HSC1_Lymphoid_Score = subset(seurat_object_STD_HSC_Normal_score_df,
                                      CellType_HSC_20240409_2 == "HSC/MPP1")$Lymphoid_Score1
  Normal_HSC2_Lymphoid_Score = subset(seurat_object_STD_HSC_Normal_score_df,
                                      CellType_HSC_20240409_2 == "HSC/MPP2")$Lymphoid_Score1
  MDS_HSC1_Lymphoid_Score = subset(seurat_object_STD_HSC_MDS_score_df,
                                   CellType_HSC_20240409_2 == "HSC/MPP1")$Lymphoid_Score1
  MDS_HSC2_Lymphoid_Score = subset(seurat_object_STD_HSC_MDS_score_df,
                                   CellType_HSC_20240409_2 == "HSC/MPP2")$Lymphoid_Score1
  tAML_HSC1_Lymphoid_Score = subset(seurat_object_STD_HSC_tAML_score_df,
                                    CellType_HSC_20240409_2 == "HSC/MPP1")$Lymphoid_Score1
  tAML_HSC2_Lymphoid_Score = subset(seurat_object_STD_HSC_tAML_score_df,
                                    CellType_HSC_20240409_2 == "HSC/MPP2")$Lymphoid_Score1
  AML_HSC1_Lymphoid_Score = subset(seurat_object_STD_HSC_AML_score_df,
                                   CellType_HSC_20240409_2 == "HSC/MPP1")$Lymphoid_Score1
  AML_HSC2_Lymphoid_Score = subset(seurat_object_STD_HSC_AML_score_df,
                                   CellType_HSC_20240409_2 == "HSC/MPP2")$Lymphoid_Score1
  
}

{
  ## Group-wise P-value statistics
  ## Each corresponding subgroup of each group was compared with the corresponding subgroup of the normal group
  # MDS vs Normal
  Stemness_HSC1_MDS_vs_Normal = wilcox.test(Normal_HSC1_Stemness_Score, MDS_HSC1_Stemness_Score)
  Stemness_HSC2_MDS_vs_Normal = wilcox.test(Normal_HSC2_Stemness_Score, MDS_HSC2_Stemness_Score)
  
  Erythroid_HSC1_MDS_vs_Normal = wilcox.test(Normal_HSC1_Erythroid_Score, MDS_HSC1_Erythroid_Score)
  Erythroid_HSC2_MDS_vs_Normal = wilcox.test(Normal_HSC2_Erythroid_Score, MDS_HSC2_Erythroid_Score)
  
  Megakaryocytic_HSC1_MDS_vs_Normal = wilcox.test(Normal_HSC1_Megakaryocytic_Score,
                                                  MDS_HSC1_Megakaryocytic_Score)
  Megakaryocytic_HSC2_MDS_vs_Normal = wilcox.test(Normal_HSC2_Megakaryocytic_Score,
                                                  MDS_HSC2_Megakaryocytic_Score)
  
  Myeloid_HSC1_MDS_vs_Normal = wilcox.test(Normal_HSC1_Myeloid_Score, MDS_HSC1_Myeloid_Score)
  Myeloid_HSC2_MDS_vs_Normal = wilcox.test(Normal_HSC2_Myeloid_Score, MDS_HSC2_Myeloid_Score)
  
  Lymphoid_HSC1_MDS_vs_Normal = wilcox.test(Normal_HSC1_Lymphoid_Score, MDS_HSC1_Lymphoid_Score)
  Lymphoid_HSC2_MDS_vs_Normal = wilcox.test(Normal_HSC2_Lymphoid_Score, MDS_HSC2_Lymphoid_Score)
  
}

{
  # tAML vs Normal
  Stemness_HSC1_tAML_vs_Normal = wilcox.test(Normal_HSC1_Stemness_Score, tAML_HSC1_Stemness_Score)
  Stemness_HSC2_tAML_vs_Normal = wilcox.test(Normal_HSC2_Stemness_Score, tAML_HSC2_Stemness_Score)
  
  Erythroid_HSC1_tAML_vs_Normal = wilcox.test(Normal_HSC1_Erythroid_Score, tAML_HSC1_Erythroid_Score)
  Erythroid_HSC2_tAML_vs_Normal = wilcox.test(Normal_HSC2_Erythroid_Score, tAML_HSC2_Erythroid_Score)
  
  Megakaryocytic_HSC1_tAML_vs_Normal = wilcox.test(Normal_HSC1_Megakaryocytic_Score,
                                                   tAML_HSC1_Megakaryocytic_Score)
  Megakaryocytic_HSC2_tAML_vs_Normal = wilcox.test(Normal_HSC2_Megakaryocytic_Score,
                                                   tAML_HSC2_Megakaryocytic_Score)
  
  Myeloid_HSC1_tAML_vs_Normal = wilcox.test(Normal_HSC1_Myeloid_Score, tAML_HSC1_Myeloid_Score)
  Myeloid_HSC2_tAML_vs_Normal = wilcox.test(Normal_HSC2_Myeloid_Score, tAML_HSC2_Myeloid_Score)
  
  Lymphoid_HSC1_tAML_vs_Normal = wilcox.test(Normal_HSC1_Lymphoid_Score, tAML_HSC1_Lymphoid_Score)
  Lymphoid_HSC2_tAML_vs_Normal = wilcox.test(Normal_HSC2_Lymphoid_Score, tAML_HSC2_Lymphoid_Score)
}


{
  # AML vs Normal
  Stemness_HSC1_AML_vs_Normal = wilcox.test(Normal_HSC1_Stemness_Score, AML_HSC1_Stemness_Score)
  Stemness_HSC2_AML_vs_Normal = wilcox.test(Normal_HSC2_Stemness_Score, AML_HSC2_Stemness_Score)
  
  Erythroid_HSC1_AML_vs_Normal = wilcox.test(Normal_HSC1_Erythroid_Score, AML_HSC1_Erythroid_Score)
  Erythroid_HSC2_AML_vs_Normal = wilcox.test(Normal_HSC2_Erythroid_Score, AML_HSC2_Erythroid_Score)
  
  Megakaryocytic_HSC1_AML_vs_Normal = wilcox.test(Normal_HSC1_Megakaryocytic_Score,
                                                  AML_HSC1_Megakaryocytic_Score)
  Megakaryocytic_HSC2_AML_vs_Normal = wilcox.test(Normal_HSC2_Megakaryocytic_Score,
                                                  AML_HSC2_Megakaryocytic_Score)
  
  Myeloid_HSC1_AML_vs_Normal = wilcox.test(Normal_HSC1_Myeloid_Score, AML_HSC1_Myeloid_Score)
  Myeloid_HSC2_AML_vs_Normal = wilcox.test(Normal_HSC2_Myeloid_Score, AML_HSC2_Myeloid_Score)
  
  Lymphoid_HSC1_AML_vs_Normal = wilcox.test(Normal_HSC1_Lymphoid_Score, AML_HSC1_Lymphoid_Score)
  Lymphoid_HSC2_AML_vs_Normal = wilcox.test(Normal_HSC2_Lymphoid_Score, AML_HSC2_Lymphoid_Score)
  
}

## Construct a list of P-value vectors
{
  Stemness_list = c(
    0.05,
    0.05,
    Stemness_HSC1_MDS_vs_Normal[["p.value"]],
    Stemness_HSC2_MDS_vs_Normal[["p.value"]],
    
    Stemness_HSC1_tAML_vs_Normal[["p.value"]],
    Stemness_HSC2_tAML_vs_Normal[["p.value"]],
    
    Stemness_HSC1_AML_vs_Normal[["p.value"]],
    Stemness_HSC2_AML_vs_Normal[["p.value"]]
  )
  
  Erythroid_list = c(
    0.05,
    0.05,
    Erythroid_HSC1_MDS_vs_Normal[["p.value"]],
    Erythroid_HSC2_MDS_vs_Normal[["p.value"]],
    Erythroid_HSC1_tAML_vs_Normal[["p.value"]],
    Erythroid_HSC2_tAML_vs_Normal[["p.value"]],
    
    Erythroid_HSC1_AML_vs_Normal[["p.value"]],
    Erythroid_HSC2_AML_vs_Normal[["p.value"]]
  )
  
  Megakaryocytic_list = c(
    0.05,
    0.05,
    Megakaryocytic_HSC1_MDS_vs_Normal[["p.value"]],
    Megakaryocytic_HSC2_MDS_vs_Normal[["p.value"]],
    Megakaryocytic_HSC1_tAML_vs_Normal[["p.value"]],
    Megakaryocytic_HSC2_tAML_vs_Normal[["p.value"]],
    Megakaryocytic_HSC1_AML_vs_Normal[["p.value"]],
    Megakaryocytic_HSC2_AML_vs_Normal[["p.value"]]
  )
  
  Myeloid_list = c(
    0.05,
    0.05,
    Myeloid_HSC1_MDS_vs_Normal[["p.value"]],
    Myeloid_HSC2_MDS_vs_Normal[["p.value"]],
    Myeloid_HSC1_tAML_vs_Normal[["p.value"]],
    Myeloid_HSC2_tAML_vs_Normal[["p.value"]],
    Myeloid_HSC1_AML_vs_Normal[["p.value"]],
    Myeloid_HSC2_AML_vs_Normal[["p.value"]]
  )
  
  Lymphoid_list = c(
    0.05,
    0.05,
    Lymphoid_HSC1_MDS_vs_Normal[["p.value"]],
    Lymphoid_HSC2_MDS_vs_Normal[["p.value"]],
    Lymphoid_HSC1_tAML_vs_Normal[["p.value"]],
    Lymphoid_HSC2_tAML_vs_Normal[["p.value"]],
    Lymphoid_HSC1_AML_vs_Normal[["p.value"]],
    Lymphoid_HSC2_AML_vs_Normal[["p.value"]]
  )
  
}

{
  ## Create a data frame
  HSC_score_pavlue_df = data.frame(
    SubGroup = HSC_merge_score_mean_class_df_scale$SubGroup,
    Stemness_socre = HSC_merge_score_mean_class_df_scale$Stemness,
    Erythroid_socre = HSC_merge_score_mean_class_df_scale$Erythroid,
    Megakaryocytic_socre = HSC_merge_score_mean_class_df_scale$Megakaryocytic,
    Myeloid_socre = HSC_merge_score_mean_class_df_scale$Myeloid,
    Lymphoid_socre = HSC_merge_score_mean_class_df_scale$Lymphoid,
    
    Stemness_pvalue = Stemness_list,
    Erythroid_pvalue = Erythroid_list,
    Megakaryocytic_pvalue = Megakaryocytic_list,
    Myeloid_pvalue = Myeloid_list,
    Lymphoid_pvalue = Lymphoid_list,
    
    Stemness_qvalue = -log10(Stemness_list),
    Erythroid_qvalue = -log10(Erythroid_list),
    Megakaryocytic_qvalue = -log10(Megakaryocytic_list),
    Myeloid_qvalue = -log10(Myeloid_list),
    Lymphoid_qvalue = -log10(Lymphoid_list)
  )
  
  print(HSC_score_pavlue_df)
  
  output_f_path = file.path(output_dir, "STAT-GSE_Rice_HSC_sorce_pavlue.xls")
  write.table(
    HSC_score_pavlue_df,
    file = output_f_path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
  
}

{
  ## Convert the wide table into a long table
  HSC_score_pavlue_df_long <-
    HSC_score_pavlue_df %>%
    pivot_longer(cols = starts_with(
      c(
        "Stemness_socre",
        "Erythroid_socre",
        "Megakaryocytic_socre",
        "Myeloid_socre",
        "Lymphoid_socre"
      )
    ),
    names_to = "direction_type",
    values_to = "direction_mean_scale_score") %>%
    pivot_longer(cols = starts_with(
      c(
        "Stemness_pvalue",
        "Erythroid_pvalue",
        "Megakaryocytic_pvalue",
        "Myeloid_pvalue",
        "Lymphoid_pvalue"
      )
    ),
    names_to = "direction_pvalue_type",
    values_to = "direction_pvalue") %>%
    pivot_longer(cols = starts_with(
      c(
        "Stemness_qvalue",
        "Erythroid_qvalue",
        "Megakaryocytic_qvalue",
        "Myeloid_qvalue",
        "Lymphoid_qvalue"
      )
    ),
    names_to = "direction_qvalue_type",
    values_to = "direction_qvalue")
  
  HSC_score_pavlue_df_long
  
  output_f_path = file.path(output_dir, "STAT-GSE_Rice_HSC_sorce_pavlue_df_long.xls")
  write.table(
    HSC_score_pavlue_df_long,
    file = output_f_path,
    sep = "\t",
    row.names = TRUE,
    col.names = TRUE,
    quote = FALSE
  )
  
}



# -------------------------------------------------------------------------

{
  # dotplot --------------------------------------------------------------------
  # direction_type
  HSC_score_pavlue_df_long$direction_type = factor(HSC_score_pavlue_df_long$direction_type , levels = rev(unique(
    HSC_score_pavlue_df_long$direction_type
  )))
  HSC_score_pavlue_df_long$SubGroup = factor(HSC_score_pavlue_df_long$SubGroup,
                                             levels = unique(HSC_score_pavlue_df_long$SubGroup))
  
  HSC_score_pavlue_df_long
  HSC_score_pavlue_df_long$direction_type
  HSC_score_pavlue_df_long$SubGroup
  
  range(HSC_score_pavlue_df_long$direction_qvalue)  #[1] 0.02181361        Inf
  range(HSC_score_pavlue_df_long$direction_mean_scale_score) # [1] -1.854837  1.651418
  
  # as.numeric
  HSC_score_pavlue_df_long <-
    HSC_score_pavlue_df_long %>%
    mutate(across(
      c(
        direction_mean_scale_score,
        direction_pvalue,
        direction_qvalue
      ),
      as.numeric
    ))
  
}

{
  p1 = ggplot(HSC_score_pavlue_df_long,
              aes(x = SubGroup, y = direction_type)) +
    geom_point(aes(size = direction_qvalue, color = direction_mean_scale_score))
  p1
  
}

{
  p2 <- p1 +
    scale_size_continuous(range = c(0, 8)) +
    
    # color
    scale_color_gradientn(
      colors = c("#462CF7", "#A88FEF", "#FEF8FB", "#E2A890", "#CE2725"),
      limits = c(min((HSC_score_pavlue_df_long$direction_mean_scale_score)
      ), max((HSC_score_pavlue_df_long$direction_mean_scale_score)
      )),
      name = "Normalized\nEnrichment\nScore",
      breaks = c(-2, -1, 0, 1, 2)
    ) +
    theme_classic() +
    theme(
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      x = 'Group',
      y = 'Direction',
      color =  "Normalized\nEnrichment\nScore",
      size = '-log10(Pvalue)'
    ) +
    
    # rename
    scale_y_discrete(
      labels = c(
        "Lymphoid_socre" = "Lymphoid",
        "Myeloid_socre" = "Myeloid",
        "Megakaryocytic_socre" = "Megakaryocytic",
        "Erythroid_socre" = "Erythroid",
        "Stemness_socre" = "Stemness"
      )
    )
  
  
  p2
  
  output_f_path = file.path(output_dir,
                            "QC-GSE_Rice_HSC_sorce_pavlue_dotplot_20240409_2.pdf")
  ggsave(
    filename = output_f_path,
    plot = p2,
    device = pdf,
    width = 210,
    height = 148,
    units = "mm"
  )
}



# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

{
  ## no add Normal
  HSC_score_pavlue_df_long_filtered = HSC_score_pavlue_df_long %>% filter(!str_detect(SubGroup , "Normal"))
  HSC_score_pavlue_df_long_filtered
}

{
  p1 = ggplot(HSC_score_pavlue_df_long_filtered,
              aes(x = SubGroup, y = direction_type)) + #建立映射
    geom_point(aes(size = direction_qvalue, color = direction_mean_scale_score)) #绘制散点
  p1
  
}

{
  p2 <- p1 +
    scale_size_continuous(range = c(0, 8)) + #气泡大小范围调整
    
    # color
    scale_color_gradientn(
      colors = c("#462CF7", "#A88FEF", "#FEF8FB", "#E2A890", "#CE2725"),
      limits = c(min((HSC_score_pavlue_df_long$direction_mean_scale_score)
      ), max((HSC_score_pavlue_df_long$direction_mean_scale_score)
      )),
      name = "Normalized\nEnrichment\nScore",
      breaks = c(-2, -1, 0, 1, 2)
    ) +
    theme_classic() +
    theme(
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 12),
      #
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      x = 'Group',
      y = 'Direction',
      color =  "Normalized\nEnrichment\nScore",
      size = '-log10(Pvalue)'
    ) +
    
    #
    scale_y_discrete(
      labels = c(
        "Lymphoid_socre" = "Lymphoid",
        "Myeloid_socre" = "Myeloid",
        "Megakaryocytic_socre" = "Megakaryocytic",
        "Erythroid_socre" = "Erythroid",
        "Stemness_socre" = "Stemness"
      )
    )
  
  
  p2
  
  output_f_path = file.path(output_dir,
                            "QC-GSE_Rice_HSC_sorce_pavlue_filtered_dotplot_20240409_2.pdf")
  ggsave(
    filename = output_f_path,
    plot = p2,
    device = pdf,
    width = 210,
    height = 148,
    units = "mm"
  )
}