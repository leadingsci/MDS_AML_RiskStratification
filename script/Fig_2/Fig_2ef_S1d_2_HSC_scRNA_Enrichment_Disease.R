rm(list = ls())

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
  output_dir = "/share/home/wyz/Project/result/2025/202508/20250805_1_scRNA_GSE_Rice_HSC_redar/20250805_2_scRNA_GSE_Rice_HSC_redar_res0.2_2"
  
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
  # save
  if (TRUE) {
    # output
    # file_name = "scRNA_SMU_Rice_14_GSE116256_Rice_1_HRA000084_Rice_4_origin_rds_combined_qc_STD_Harmony_res_select_HSC_harmony_res_20240516.rds"
    file_name = "scRNA_Harmony_res_select_HSC_harmony_res_20240516.rds"
    saveRDS(seurat_object_STD_HSC, file.path(output_dir, file_name))
    
    # read
    #seurat_object_STD = readRDS(file.path(output_dir, file_name))
  }
  
}

# -------------------------------------------------------------------------

## umap
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

## 5 gene sets of differentiated types
Stemness_gene_df <- subset(gene_sheet_df, Type == "Stemness")
Erythroid_gene_df <- subset(gene_sheet_df, Type == "Erythroid")
Megakaryocytic_gene_df <- subset(gene_sheet_df, Type == "Megakaryocytic")
Myeloid_gene_df <- subset(gene_sheet_df, Type == "Myeloid")
Lymphoid_gene_df <- subset(gene_sheet_df, Type == "Lymphoid")

print(Stemness_gene_df)
print(Stemness_gene_df$Gene_symbol)


#  ---------------------------------------------------
Stemness_gene_df_keep = Stemness_gene_df[Stemness_gene_df$Gene_symbol %in% rownames(seurat_object_STD_HSC@assays$RNA@counts), ]
Erythroid_gene_df_keep = Erythroid_gene_df[Erythroid_gene_df$Gene_symbol %in% rownames(seurat_object_STD_HSC@assays$RNA@counts), ]
Megakaryocytic_gene_df_keep = Megakaryocytic_gene_df[Megakaryocytic_gene_df$Gene_symbol %in% rownames(seurat_object_STD_HSC@assays$RNA@counts), ]
Myeloid_gene_df_keep = Myeloid_gene_df[Myeloid_gene_df$Gene_symbol %in% rownames(seurat_object_STD_HSC@assays$RNA@counts), ]
Lymphoid_gene_df_keep = Lymphoid_gene_df[Lymphoid_gene_df$Gene_symbol %in% rownames(seurat_object_STD_HSC@assays$RNA@counts), ]

dim(Stemness_gene_df_keep)



# # AddModuleScore  ---------------------------------------------------
seurat_object_STD_HSC =
  AddModuleScore(
    object = seurat_object_STD_HSC,
    features = list(Stemness_gene_df_keep$Gene_symbol),
    # 官方输入基因名为list
    ctrl = 100,
    # default
    name = 'Stemness_Score',
    seed = 1,
  )

seurat_object_STD_HSC =
  AddModuleScore(
    object = seurat_object_STD_HSC,
    features = list(Erythroid_gene_df_keep$Gene_symbol),
    # 官方输入基因名为list
    ctrl = 100,
    # default
    name = 'Erythroid_Score',
    seed = 1,
  )


seurat_object_STD_HSC =
  AddModuleScore(
    object = seurat_object_STD_HSC,
    features = list(Megakaryocytic_gene_df_keep$Gene_symbol),
    # 官方输入基因名为list
    ctrl = 100,
    # default
    name = 'Megakaryocytic_Score',
    seed = 1,
  )

seurat_object_STD_HSC =
  AddModuleScore(
    object = seurat_object_STD_HSC,
    features = list(Myeloid_gene_df_keep$Gene_symbol),
    # 官方输入基因名为list
    ctrl = 100,
    # default
    name = 'Myeloid_Score',
    seed = 1,
  )

seurat_object_STD_HSC =
  AddModuleScore(
    object = seurat_object_STD_HSC,
    features = list(Lymphoid_gene_df_keep$Gene_symbol),
    # 官方输入基因名为list
    ctrl = 100,
    # default
    name = 'Lymphoid_Score',
    seed = 1,
  )


seurat_object_STD_HSC@meta.data$Stemness_Score1


# save df
hsc_score_result_cbind <- cbind(
  Stemness = seurat_object_STD_HSC@meta.data$Stemness_Score1,
  Erythroid = seurat_object_STD_HSC@meta.data$Erythroid_Score1,
  Megakaryocytic = seurat_object_STD_HSC@meta.data$Megakaryocytic_Score1,
  Myeloid = seurat_object_STD_HSC@meta.data$Myeloid_Score1,
  Lymphoid = seurat_object_STD_HSC@meta.data$Lymphoid_Score1
)

head(hsc_score_result_cbind)

# to df
hsc_score_result_df <- as.data.frame(hsc_score_result_cbind)
head(hsc_score_result_df[0:5, 0:5])

dim(hsc_score_result_df)
length(as.vector(colnames(seurat_object_STD_HSC)))

rownames(hsc_score_result_df) <- as.vector(colnames(seurat_object_STD_HSC))
head(hsc_score_result_df[0:5, 0:5])



hsc_score_result_df <- hsc_score_result_df %>%
  rownames_to_column(var = "Cell_ID")
hsc_score_result_df

output_f_path = file.path(output_dir, "20250801_1_HSC_cell_pathway_score.xls")
write.table(
  hsc_score_result_df,
  file = output_f_path,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE
)



## Calculate the average score of the cell differentiation preference for each subgroup.
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

#  Rename column names
colnames(seurat_object_STD_HSC_mean_df) = c("subClass",
                                            "Stemness",
                                            "Erythroid",
                                            "Megakaryocytic",
                                            "Myeloid",
                                            "Lymphoid")


output_f_path = file.path(output_dir,
                          "20250801_2_HSC_cell_pathway_score_HSC_mean_score.xls")
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

colnames(seurat_object_STD_HSC_score_group_mean_df) = c(
  "DiseaseClass",
  "SubClass",
  "Stemness",
  "Erythroid",
  "Megakaryocytic",
  "Myeloid",
  "Lymphoid"
)

output_f_path = file.path(
  output_dir,
  "20250801_3_HSC_cell_pathway_score_HSC_Desease_mean_score_redar_data.xls"
)
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
  # Group the scores according to the diseases ---------------------------------------------------------------------
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
  
  table(seurat_object_STD_HSC_desease_score_df$DiseaseClass)
  
  seurat_object_STD_HSC_Normal_score_df <- subset(x = seurat_object_STD_HSC_desease_score_df, DiseaseClass == "Normal")
  seurat_object_STD_HSC_MDS_score_df <- subset(x = seurat_object_STD_HSC_desease_score_df, DiseaseClass == "MDS")
  seurat_object_STD_HSC_tAML_score_df <- subset(x = seurat_object_STD_HSC_desease_score_df, DiseaseClass == "tAML")
  seurat_object_STD_HSC_AML_score_df <- subset(x = seurat_object_STD_HSC_desease_score_df, DiseaseClass == "AML")
  
  head(seurat_object_STD_HSC_Normal_score_df)
  dim(seurat_object_STD_HSC_Normal_score_df)
  
}

{
  # Specify the classification of a certain column and calculate the enrichment score and mean value for each classification
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
  
  
  output_f_path = file.path(
    output_dir,
    "20250801_4_HSC_cell_pathway_score_HSC_Desease_mean_score_dotplot_data.xls"
  )
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
  # 按列进行矫正
  HSC_merge_score_mean_class_df_scale = apply(HSC_merge_score_mean_class_df[2:6], 2, scale)
  HSC_merge_score_mean_class_df_scale
  
  
  HSC_merge_score_mean_class_df_scale = cbind(HSC_merge_score_mean_class_df[, c(7, 1, 8)],
                                              HSC_merge_score_mean_class_df_scale)
  HSC_merge_score_mean_class_df_scale
  
  output_f_path = file.path(
    output_dir,
    "20250801_5_HSC_cell_pathway_score_HSC_Desease_mean_score_dotplot_data_scale.xls"
  )
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
  # Obtain the enrichment scores of each subgroup for the four groups
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
  
  output_f_path = file.path(
    output_dir,
    "20250801_6_HSC_cell_pathway_score_HSC_Desease_mean_score_dotplot_data_scale_pvalue.xls"
  )
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
  ## long df
  
  head(HSC_score_pavlue_df)
  # 1)
  HSC_score_pavlue_df_long_NES = HSC_score_pavlue_df %>%
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
    values_to = "direction_mean_scale_score")
  head(HSC_score_pavlue_df_long_NES)
  colnames(HSC_score_pavlue_df_long_NES)
  
  # 2）
  HSC_score_pavlue_df_long_pvalue = HSC_score_pavlue_df_long_NES %>%
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
    values_to = "direction_pvalue")
  head(HSC_score_pavlue_df_long_pvalue)
  colnames(HSC_score_pavlue_df_long_pvalue)
  
  
  # 3）-log10(Pvalue)
  HSC_score_pavlue_df_long_qvalue = HSC_score_pavlue_df_long_pvalue %>%
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
  head(HSC_score_pavlue_df_long_qvalue)
  colnames(HSC_score_pavlue_df_long_qvalue)
  
  # output_f_path = file.path(output_dir,"STAT-GSE_Rice_HSC_sorce_pavlue_df_long.xls")
  output_f_path = file.path(
    output_dir,
    "20250801_8_long3_HSC_cell_pathway_score_HSC_Desease_mean_score_dotplot_data_scale_pvalue_long.xls"
  )
  write.table(
    HSC_score_pavlue_df_long_qvalue,
    file = output_f_path,
    sep = "\t",
    row.names = TRUE,
    col.names = TRUE,
    quote = FALSE
  )
}


{
  library(dplyr)
  # Filtering condition: The values in all three columns are exactly the same (after removing the suffix "_qvalue" and "_pvalue", they remain the same)
  filtered_df <- HSC_score_pavlue_df_long_qvalue[sapply(1:nrow(HSC_score_pavlue_df_long_qvalue), function(i) {
    # Extract values from three columns
    type <- HSC_score_pavlue_df_long_qvalue$direction_type[i]
    ptype <- HSC_score_pavlue_df_long_qvalue$direction_pvalue_type[i]
    qtype <- HSC_score_pavlue_df_long_qvalue$direction_qvalue_type[i]
    
    # Remove suffix for comparison
    base_type <- gsub("_socre$", "", type)
    base_ptype <- gsub("_pvalue$", "", ptype)
    base_qtype <- gsub("_qvalue$", "", qtype)
    
    # Check if the three base types are the same
    base_type == base_ptype && base_ptype == base_qtype
  }), ]
  
  filtered_df
  dim(filtered_df)
  
  output_f_path = file.path(
    output_dir,
    "20250801_9_long3_deldup_HSC_cell_pathway_score_HSC_Desease_mean_score_dotplot_data_scale_pvalue_long.xls"
  )
  write.table(
    filtered_df,
    file = output_f_path,
    sep = "\t",
    row.names = TRUE,
    col.names = TRUE,
    quote = FALSE
  )
}

{
  HSC_score_pavlue_df_long = filtered_df
  HSC_score_pavlue_df_long
}

{
  HSC_score_pavlue_df_long_adj = HSC_score_pavlue_df_long
}

{
  # If the value in the "direction_pvalue" column is less than 0.0001, then replace it with 0.00001.
  library(dplyr)
  
  HSC_score_pavlue_df_long_adj <- HSC_score_pavlue_df_long_adj %>%
    mutate(direction_pvalue = ifelse(direction_pvalue < 0.0001, 0.000099, direction_pvalue))
  HSC_score_pavlue_df_long_adj$direction_pvalue
  
  # 赋值给qvalue
  HSC_score_pavlue_df_long_adj$direction_qvalue = -log10(HSC_score_pavlue_df_long_adj$direction_pvalue)
  HSC_score_pavlue_df_long_adj$direction_qvalue
  
  
}

{
  library(dplyr)
  
  # 假设你的数据框名称为 df
  HSC_score_pavlue_df_long_adj <- HSC_score_pavlue_df_long_adj %>%
    mutate(
      Pvalue_Mark = case_when(
        direction_pvalue < 0.0001 ~ "****",
        # P<0.0001
        direction_pvalue < 0.001  ~ "***",
        # P<0.001
        direction_pvalue < 0.01   ~ "**",
        # P<0.01
        direction_pvalue < 0.05   ~ "*",
        # P<0.05
        direction_pvalue >= 0.05  ~ "ns",
        # P<=0.05
        TRUE ~ ""  # else
      )
    )
  
  HSC_score_pavlue_df_long_adj$Pvalue_Mark
}

{
  output_f_path = file.path(
    output_dir,
    "20250801_10_long3_deldup_limit_max_HSC_cell_pathway_score_HSC_Desease_mean_score_dotplot_data_scale_pvalue_long.xls"
  )
  write.table(
    HSC_score_pavlue_df_long_adj,
    file = output_f_path,
    sep = "\t",
    row.names = TRUE,
    col.names = TRUE,
    quote = FALSE
  )
}

{
  HSC_score_pavlue_df_long = HSC_score_pavlue_df_long_adj
}

# data ok -------------------------------------------------------------------------

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
  
  # Convert the columns "Score1", "Score2", and "Score3" to numeric type
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
  HSC_score_pavlue_df_long
}

{
  p1 = ggplot(HSC_score_pavlue_df_long,
              aes(x = SubGroup, y = direction_type)) +
    geom_point(aes(size = direction_qvalue, color = direction_mean_scale_score))
  p1
  
}

{
  p2 <- p1 +
    scale_color_gradientn(
      colors = c("#462CF7", "#A88FEF", "#FEF8FB", "#E2A890", "#CE2725"),
      # limits = c(min((HSC_score_pavlue_df_long$direction_mean_scale_score)),max((HSC_score_pavlue_df_long$direction_mean_scale_score))),
      limits = c(-2, 2),
      name = "Normalized\nEnrichment\nScore",
      breaks = c(-2, -1, 0, 1, 2)
    ) +
    scale_size_continuous(
      range = c(0, 16),
      breaks = c(0, 1, 1.3, 2, 3, 4),
      #
      name = '-log10(Pvalue)'  #
    ) +
    theme_classic() +
    theme(
      legend.text = element_text(size = 12),
      #
      legend.title = element_text(size = 12),
      #
      axis.text = element_text(size = 12),
      #
      axis.title = element_text(size = 12),
      #
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) + #
    labs(
      x = 'Group',
      y = 'Direction',
      #xy轴标题修改
      color =  "Normalized\nEnrichment\nScore",
      size = '-log10(Pvalue)'
    ) + #
    
    # 重命名纵坐标的刻度标签
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
  
  output_f_path = file.path(
    output_dir,
    "20250801_8_7_HSC_cell_pathway_score_HSC_Desease_mean_score_dotplot_data_scale_pvalue_long.pdf"
  )
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
  ## no show Normal group
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
    scale_color_gradientn(
      colors = c("#462CF7", "#A88FEF", "#FEF8FB", "#E2A890", "#CE2725"),
      # limits = c(min((HSC_score_pavlue_df_long$direction_mean_scale_score)),max((HSC_score_pavlue_df_long$direction_mean_scale_score))),
      limits = c(-2, 2),
      name = "Normalized\nEnrichment\nScore",
      breaks = c(-2, -1, 0, 1, 2)
    ) +
    scale_size_continuous(
      range = c(0, 16),
      #
      breaks = c(0, 1, 1.3, 2, 3, 4),
      #
      name = '-log10(Pvalue)'  #
    ) +
    theme_classic() +
    theme(
      legend.text = element_text(size = 12),
      #
      legend.title = element_text(size = 12),
      #
      axis.text = element_text(size = 12),
      #
      axis.title = element_text(size = 12),
      #
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) + #
    labs(
      x = 'Group',
      y = 'Direction',
      #
      color =  "Normalized\nEnrichment\nScore",
      size = '-log10(Pvalue)'
    ) + #
    
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
  output_f_path = file.path(
    output_dir,
    "20250801_10_1_HSC_cell_pathway_score_HSC_Desease_mean_score_dotplot_data_scale_pvalue_long.pdf"
  )
  ggsave(
    filename = output_f_path,
    plot = p2,
    device = pdf,
    width = 210,
    height = 148,
    units = "mm"
  )
}
