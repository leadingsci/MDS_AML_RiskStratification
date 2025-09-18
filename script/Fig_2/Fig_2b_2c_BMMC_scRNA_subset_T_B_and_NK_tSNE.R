rm(list = ls())

# libPaths --------------------------------------------------------------------


.libPaths()
options(warn = 1)

# Load Package ------------------------------------------------------------

# devtools::install_github("satijalab/seurat", ref = "develop")
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
  
  library(Seurat)
  #R 4.2.2 #R 4.2.3
  
  # library(harmony); #R 4.2.3
  # library(clustree) #R 4.2.3
  # library(future); #R 4.2.2 #R 4.2.3
  
  # library(Matrix)
  # packageVersion("Matrix")
  
  library(ComplexHeatmap)
  
}




# input path --------------------------------------------------------------

if (TRUE) {
  ## sample info
  # sample_sheet_f = "/share/home/wyz/Lab/Code/ProjectCode/2024/202403/20240305_scRNA_AML/sample_info/scRNA_Concat_Colname_Ref_20240316.xlsx"
  
  ## input_dir
  # input_dir = "/share/home/wyz/Project/result/20240320_scRNA_AML_Annotation/20240321_scRNA_GSE_GSE116256_annotation/scRNA_SMU_GSE116256_PRJNA720840_GSE205490_HRA000084_AML_MDS_83_origin_rds_combined_qc_STD_Harmony_res_anno_GES116256.rds"
  
  ## output
  output_dir = "/share/home/wyz/Project/result/2025/202503/20250307_1_UMAP_immune/20250801_1_UMAP_immune_color"
  
}

# read file ---------------------------------------------------------------


# -------------------------------------------------------------------------
if (TRUE) {
  input_f_path = "/share/home/wyz/Project/result/2025/202502/20250221_1_MDS_GSE58831_immune_cell_re_cibersortx/20250221_1_T_sub_cell_qc/scRNA_res1_select_TNK_res1_select_cd4cd8_delna_20250221.rds"
  seurat_object_T <- readRDS(input_f_path)
  print(seurat_object_T)
  head(seurat_object_T@meta.data)
  table(seurat_object_T@meta.data$CellType_cd4cd8_1)
  table(seurat_object_T@meta.data$GSE116256_TransferLabel)
  table(seurat_object_T@meta.data$cluster_mark_unique)
  levels(seurat_object_T@meta.data$cluster_mark_unique)
}



{
  seurat_object_T@meta.data$cluster_mark_unique_delna = droplevels(seurat_object_T@meta.data$cluster_mark_unique)
  
  length(table(seurat_object_T@meta.data$cluster_mark_unique_delna))
}

if (TRUE) {
  input_f_path = "/share/home/wyz/Project/result/2025/202502/20250221_1_MDS_GSE58831_immune_cell_re_cibersortx/20250221_2_B_sub_cell_qc/scRNA_select_B_NK_anno_20241204_20250221.rds"
  seurat_object_B_NK <- readRDS(input_f_path)
  print(seurat_object_B_NK@meta.data)
  table(seurat_object_B_NK@meta.data$GSE116256_TransferLabel)
  table(seurat_object_B_NK@meta.data$cluster_mark_unique)
}


{
  seurat_object_B_NK@meta.data$cluster_mark_unique_delna = droplevels(seurat_object_B_NK@meta.data$cluster_mark_unique)
  table(seurat_object_B_NK@meta.data$cluster_mark_unique_delna)
  length(table(seurat_object_B_NK@meta.data$cluster_mark_unique_delna))
}



# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

{
  ## rename
  
  ## 1. remove factor
  table(seurat_object_T$cluster_mark_unique)
  levels(seurat_object_T$cluster_mark_unique)
  
  ## 2. transfer to character
  seurat_object_T$cluster_mark_unique_20250320 = as.character(seurat_object_T$cluster_mark_unique)
  
  ## 3.rename
  seurat_object_T$cluster_mark_unique_20250320[seurat_object_T$cluster_mark_unique_20250320 == "CD8_T_10_EEF1G_PCNP"] <- "CD8_T_10_EEF1G"
  seurat_object_T$cluster_mark_unique_20250320[seurat_object_T$cluster_mark_unique_20250320 == "CD4_T_03_MT-ND4L"] <- "CD4_T_03_MT_ND4L"
  seurat_object_T$cluster_mark_unique_20250320[seurat_object_T$cluster_mark_unique_20250320 == "CD8_T_09_MT-ND2"] <- "CD8_T_09_MT_ND2"
  
  ## 4. levels
  seurat_object_T$cluster_mark_unique_20250320 = factor(
    seurat_object_T$cluster_mark_unique_20250320,
    levels = c(
      "CD4_T_01_FYB1",
      "CD4_T_02_CCL5",
      "CD4_T_03_MT_ND4L",
      "CD4_T_04_ISG15",
      "CD4_T_05_APOOL",
      "CD4_Treg_01_ITGB1",
      "CD8_T_01_NBEAL1",
      "CD8_T_02_LMNA",
      "CD8_T_03_HBA2",
      "CD8_T_04_CD8B",
      "CD8_T_05_CCL4",
      "CD8_T_06_CCR7",
      "CD8_T_07_S100A11",
      "CD8_T_08_ATP5MC2",
      "CD8_T_09_MT_ND2",
      "CD8_T_10_EEF1G"
      
    )
  )
  
  ## 5. levels
  levels(seurat_object_T$cluster_mark_unique_20250320)
  table(seurat_object_T$cluster_mark_unique_20250320)
  
}

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

{
  ## rename
  
  ## 1. remove factor
  table(seurat_object_B_NK$cluster_mark_unique)
  levels(seurat_object_B_NK$cluster_mark_unique)
  
  ## 2. transfer to character
  seurat_object_B_NK$cluster_mark_unique_20250320 = as.character(seurat_object_B_NK$cluster_mark_unique)
  
  ## 3.rename
  seurat_object_B_NK$cluster_mark_unique_20250320[seurat_object_B_NK$cluster_mark_unique_20250320 == "B_03_HLA-DRA"] <- "B_02_HLA-DRA"
  seurat_object_B_NK$cluster_mark_unique_20250320[seurat_object_B_NK$cluster_mark_unique_20250320 == "B_05_STMN1"] <- "B_03_STMN1"
  seurat_object_B_NK$cluster_mark_unique_20250320[seurat_object_B_NK$cluster_mark_unique_20250320 == "B_06_DNTT"] <- "B_04_DNTT"
  seurat_object_B_NK$cluster_mark_unique_20250320[seurat_object_B_NK$cluster_mark_unique_20250320 == "B_07_MS4A1_CD79A"] <- "B_05_MS4A1_CD79A"
  
  ## 4. levels
  seurat_object_B_NK$cluster_mark_unique_20250320 = factor(
    seurat_object_B_NK$cluster_mark_unique_20250320,
    levels = c(
      "NK_01_IL7R",
      "NK_02_MT-CYB",
      "NK_03_GZMK",
      "NK_04_GZMH",
      "NK_05_CD3D",
      "NK_06_SPON2",
      "B_01_MS4A1_HLA-DRA",
      "B_02_HLA-DRA",
      "B_03_STMN1",
      "B_04_DNTT",
      "B_05_MS4A1_CD79A"
    )
  )
  
  ## 5. levels
  levels(seurat_object_B_NK$cluster_mark_unique_20250320)
  table(seurat_object_B_NK$cluster_mark_unique_20250320)
  
}


{
  seurat_object_T
  seurat_object_B_NK
}

# -------------------------------------------------------------------------


## save
if (TRUE) {
  # output
  file_name = paste0("scRNA_STD_Harmony_res1_select_TNK_res1_select_cd4cd8_delna_20250321.rds")
  print(file_name)
  saveRDS(seurat_object_T, file.path(output_dir, file_name))
  
  # read
  #seurat_object_STD = readRDS(file.path(output_dir, file_name))
}


if (TRUE) {
  # output
  file_name = paste0(
    "scRNA_STD_Harmony_res1_anno_GES116256_select_B_NK_anno_20241204_20250321.rds"
  )
  print(file_name)
  saveRDS(seurat_object_B_NK, file.path(output_dir, file_name))
}
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------

{
  ## color
  library(ggsci)
  library(RColorBrewer)
  library(scales)
  
  col_list <- colorRampPalette((pal_npg("nrc")(9)))(16)
  show_col(col_list)
  
}

{
  color_mapping16_format = c(
    "CD8_T_01_NBEAL1" = "#fccdac",
    "CD8_T_02_LMNA" = "#7670b2",
    "CD8_T_03_HBA2" = "#6b3e97",
    "CD8_T_04_CD8B" = "#189d77",
    "CD8_T_05_CCL4" = "#d76127",
    "CD8_T_06_CCR7" = "#c9b2d5",
    "CD8_T_07_S100A11" = "#e62a89",
    "CD8_T_08_ATP5MC2" = "#65a644",
    "CD8_T_09_MT_ND2" = "#e4ab23",
    "CD8_T_10_EEF1G" = "#a4772b",
    
    "CD4_T_01_FYB1" = "#a6cde1",
    "CD4_T_02_CCL5" = "#1f78b4",
    "CD4_T_03_MT_ND4L" = "#b4d789",
    "CD4_T_04_ISG15" = "#329f48",
    "CD4_T_05_APOOL" = "#f59898",
    "CD4_Treg_01_ITGB1" = "#f47e20"
  )
  
  color_mapping16_format
  
  table(seurat_object_T$cluster_mark_unique_20250320)
  levels(seurat_object_T$cluster_mark_unique_20250320)
}


# T umap-------------------------------------------------------------------------

seurat_object_plot = NULL

if (TRUE) {
  seurat_object_plot =
    DimPlot(
      seurat_object_T,
      reduction = "std.umap.harmony.tnk.cd48",
      group.by = c("cluster_mark_unique_20250320"),
      combine = TRUE,
      label.size = 2,
      pt.size = 0.001,
      label = TRUE,
      ncol = 1,
      cols = color_mapping16_format,
      raster = FALSE
    ) +
    ggtitle("UMAP_CD4_T_CD8_T") +
    theme(
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 4),
      legend.key.size = unit(0.1, "cm"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      #axis.text.x = element_blank(),
      #axis.text.y = element_blank(),
      title = element_text(size = 8),
      plot.title = element_text(hjust = 0.5)
    )
  seurat_object_plot
  
  output_f_path = file.path(output_dir, "20250801_1_3_UMAP_CD4_T_CD8_T.pdf")
  ggsave(
    filename = output_f_path,
    plot = seurat_object_plot,
    device = pdf,
    width = 120,
    height = 100,
    units = "mm"
  )
}

# T tsne-------------------------------------------------------------------------
if (TRUE) {
  seurat_object_plot =
    DimPlot(
      seurat_object_T,
      reduction = "std.tsne.harmony.tnk.cd48",
      group.by = c("cluster_mark_unique_20250320"),
      combine = TRUE,
      label.size = 2,
      pt.size = 0.001,
      label = TRUE,
      ncol = 1,
      cols = color_mapping16_format,
      raster = FALSE
    ) +
    ggtitle("tSNE_CD4_T_CD8_T") +
    theme(
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 4),
      legend.key.size = unit(0.1, "cm"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      #axis.text.x = element_blank(),
      #axis.text.y = element_blank(),
      title = element_text(size = 8),
      plot.title = element_text(hjust = 0.5)
    )
  seurat_object_plot
  
  output_f_path = file.path(output_dir, "20250801_1_5_tsne_CD4_T_CD8_T.pdf")
  ggsave(
    filename = output_f_path,
    plot = seurat_object_plot,
    device = pdf,
    width = 120,
    height = 100,
    units = "mm"
  )
}


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
{
  seurat_object_B_NK
  table(seurat_object_B_NK@meta.data$cluster_mark_unique_delna)
  # NK_01_IL7R       NK_02_MT-CYB         NK_03_GZMK         NK_04_GZMH         NK_05_CD3D        NK_06_SPON2 B_01_MS4A1_HLA-DRA       B_03_HLA-DRA
  # 6417               3667               4733               8688               8239               7152               6842               6681
  # B_05_STMN1          B_06_DNTT   B_07_MS4A1_CD79A
  # 1967               1462                726
  length(table(seurat_object_B_NK@meta.data$cluster_mark_unique_delna))
}
{
  color_mapping11_format = c(
    "NK_01_IL7R" = "#cfb3d6",
    "NK_02_MT-CYB" = "#f8cce0",
    "NK_03_GZMK" = "#905da5",
    "NK_04_GZMH" = "#898b35",
    "NK_05_CD3D" = "#b28072",
    "NK_06_SPON2" = "#a5292a",
    "B_01_MS4A1_HLA-DRA" = "#fab263",
    "B_02_HLA-DRA" = "#bebdbd",
    "B_03_STMN1" = "#1f78b4",
    "B_04_DNTT" = "#a6cde1",
    "B_05_MS4A1_CD79A" = "#68c2a3"
  )
  
}


# B umap-------------------------------------------------------------------------


if (TRUE) {
  seurat_object_plot =
    DimPlot(
      seurat_object_B_NK,
      reduction = "std.umap.harmony.b.nk",
      group.by = c("cluster_mark_unique_20250320"),
      combine = TRUE,
      label.size = 2,
      pt.size = 0.001,
      label = TRUE,
      ncol = 1,
      cols = color_mapping11_format,
      raster = FALSE
    ) +
    ggtitle("UMAP_B_NK") +
    theme(
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 4),
      legend.key.size = unit(0.1, "cm"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      #axis.text.x = element_blank(),
      #axis.text.y = element_blank(),
      title = element_text(size = 8),
      plot.title = element_text(hjust = 0.5)
    )
  seurat_object_plot
  
  output_f_path = file.path(output_dir, "20250801_1_2_UMAP_B_NK.pdf")
  ggsave(
    filename = output_f_path,
    plot = seurat_object_plot,
    device = pdf,
    width = 120,
    height = 100,
    units = "mm"
  )
}

# B tsne-------------------------------------------------------------------------


if (TRUE) {
  seurat_object_plot =
    DimPlot(
      seurat_object_B_NK,
      reduction = "std.tsne.harmony.b.nk",
      group.by = c("cluster_mark_unique_20250320"),
      combine = TRUE,
      label.size = 2,
      pt.size = 0.001,
      label = TRUE,
      ncol = 1,
      cols = color_mapping11_format,
      raster = FALSE
    ) +
    ggtitle("tSNE_B_NK") +
    theme(
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 4),
      legend.key.size = unit(0.1, "cm"),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      #axis.text.x = element_blank(),
      #axis.text.y = element_blank(),
      title = element_text(size = 8),
      plot.title = element_text(hjust = 0.5)
    )
  seurat_object_plot
  
  output_f_path = file.path(output_dir, "20250801_1_2_tsne_B_NK.pdf")
  ggsave(
    filename = output_f_path,
    plot = seurat_object_plot,
    device = pdf,
    width = 120,
    height = 100,
    units = "mm"
  )
}