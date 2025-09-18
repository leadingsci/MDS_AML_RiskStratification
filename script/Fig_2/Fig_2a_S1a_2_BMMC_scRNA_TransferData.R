# -----------------------------------------------------------------------------
# Author: [Yongzhang Yu]
# Description: [Fig_2a_S1a_2_BMMC_scRNA_TransferData]
# -----------------------------------------------------------------------------


# init --------------------------------------------------------------------

rm(list = ls())

# Load Package ------------------------------------------------------------
library(readr) #R 4.2.2 #R 4.2.3
library(tidyverse)
#R 4.2.2
library(openxlsx)
#R 4.2.3
library(readxl)
#R 4.2.3
library(dplyr)
#R 4.2.2 #R 4.2.3

library(ggplot2)
#R 4.2.2 #R 4.2.3
library(ggrastr) #R 4.2.3
library(ggpubr) #R 4.2.2 #R 4.2.3

library(Seurat)
#R 4.2.2 #R 4.2.3
options(warn = 1)

# input path --------------------------------------------------------------
if (TRUE) {
  ## sample info
  sample_sheet_f = "/share/home/wyz/Lab/Code/ProjectCode/2024/202403/20240305_scRNA_AML/sample_info/scRNA_Concat_Colname_Ref_20240316.xlsx"
  
  ## output
  output_dir = "/share/home/wyz/Project/result/20240320_scRNA_AML_Annotation/20240321_scRNA_GSE_GSE116256_annotation/"
  
}

# read file ---------------------------------------------------------------
## object
if (TRUE) {
  # 读入文件
  input_f_path_GSE = "/share/home/wyz/Project/result/20240316_scRNA_AML/20240320_scRNA_GSE_AML_MDS/scRNA_SMU_GSE116256_PRJNA720840_GSE205490_HRA000084_AML_MDS_83_origin_rds_combined_qc_STD_Harmony_res.rds"
  
  seurat_object_GSE <- readRDS(input_f_path_GSE)
  print(seurat_object_GSE)
  print(unique(seurat_object_GSE@meta.data[["Sample_ID"]]))
  print(unique(seurat_object_GSE@meta.data[["Source"]]))
}


## 获取参考的对象
if (TRUE) {
  seurat_ref =  subset(seurat_object_GSE, Source %in% c("GSE116256"))
  print(seurat_ref)
  
  seurat_query = seurat_object_GSE
  print(seurat_query)
  
  # find anchors
  anchors <- FindTransferAnchors(
    reference = seurat_ref,
    query = seurat_query,
    dims = 1:30,
    normalization.method = "LogNormalize"
  )
  print(anchors)
  
  # transfer labels
  # Specify the column name of the reference cell type that uses the reference dataset
  predictions <- TransferData(anchorset = anchors, refdata = seurat_ref$CellTypeRef)
  print(predictions)
  
  # Insert the annotation results into the meta
  seurat_query = AddMetaData(
    object = seurat_query,
    metadata = predictions$predicted.id,
    col.name = "GSE116256_TransferLabel"
  )
  # table(seurat_query@meta.data$Sample_ID)
}


## save
if (TRUE) {
  ## 保存数据框
  output_f_path = file.path(output_dir, "Annotation-sample_transfer_label.xls")
  print(output_f_path)
  write.table(
    predictions,
    file = output_f_path,
    sep = "\t",
    row.names = TRUE,
    col.names = TRUE,
    quote = FALSE
  )
}

## save rds
if (TRUE) {
  # output
  file_name = paste0(
    "scRNA_SMU_GSE116256_PRJNA720840_GSE205490_HRA000084_AML_MDS_83_origin_rds_combined_qc_STD_Harmony_res_anno_GES116256",
    ".rds"
  )
  print(file_name)
  file.path(output_dir, file_name)
  saveRDS(seurat_query, file.path(output_dir, file_name))
  
  ## read
  # seurat_object_STD_harmony = readRDS(file.path(output_dir, file_name))
}



# plot --------------------------------------------------------------------

seurat_query@meta.data$RNA_snn_res.1 = factor(x = seurat_query@meta.data$RNA_snn_res.1,
                                              levels = seq(0, 42))
print(seurat_query@meta.data$RNA_snn_res.1)


## umap split source
if (TRUE) {
  seurat_object_plot =
    DimPlot(
      seurat_query,
      reduction = "std.umap.harmony",
      group.by = c("RNA_snn_res.1", "GSE116256_TransferLabel"),
      split.by = c("Source"),
      combine = TRUE,
      label.size = 2,
      pt.size = 0.5,
      label = TRUE,
      ncol = 3,
      raster = FALSE
    ) +
    ggtitle("QC-STD_harmony_umap_integrated_res1")
  seurat_object_plot
  
  output_f_path = file.path(output_dir,
                            "QC-STD_umap_integrated_res1_split_Source.pdf")
  ggsave(
    filename = output_f_path,
    plot = seurat_object_plot,
    device = pdf,
    width = 297,
    height = 210,
    units = "mm"
  )
}


## umap split source
if (TRUE) {
  seurat_object_plot =
    DimPlot(
      seurat_query,
      reduction = "std.umap.harmony",
      group.by = c("GSE116256_TransferLabel"),
      split.by = c("Source"),
      combine = TRUE,
      label.size = 4,
      pt.size = 0.5,
      label = TRUE,
      ncol = 5,
      raster = FALSE
    ) +
    ggtitle("QC-STD_harmony_umap_integrated_res1")
  seurat_object_plot
  
  output_f_path = file.path(output_dir,
                            "QC-STD_umap_integrated_res1_split_Source_2.pdf")
  ggsave(
    filename = output_f_path,
    plot = seurat_object_plot,
    device = pdf,
    width = 297,
    height = 210,
    units = "mm"
  )
}