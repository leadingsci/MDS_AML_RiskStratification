# -----------------------------------------------------------------------------
# Author: [Yongzhang Yu]
# Description: [Fig_2a_S1a_1_BMMC_scRNA_Data_Integration]
# -----------------------------------------------------------------------------


# init --------------------------------------------------------------------

rm(list = ls())


# R-4.2.3
.libPaths(c("/share/home/bioinfo/anaconda3/envs/R-4.2.3/lib/R/library"))

.libPaths()
packageVersion("Seurat")
options(warn = 1)


# Load Package ------------------------------------------------------------
if (TRUE) {
  library(tidyverse)
  
  library(readr)
  library(dplyr)
  
  library(ggplot2)
  
  library(openxlsx)
  
  library(readxl)
  
  
  library(Seurat)
  
  library(harmony)
  
  library(clustree)
  library(future)
  
  
  library(pheatmap)
}



# input path --------------------------------------------------------------

if (TRUE) {
  ## sample info
  sample_sheet_f = "/share/home/wyz/Lab/Code/ProjectCode/2024/202403/20240305_scRNA_AML/sample_info/scRNA_Concat_Colname_Ref_20240316.xlsx"
  
  ## input_dir
  input_dir = "/share/home/wyz/Project/result/20240316_scRNA_AML/20240320_scRNA_GSE_AML_MDS/"
  
  ## output
  output_dir =  "/share/home/wyz/Project/result/20240316_scRNA_AML/20240320_scRNA_GSE_AML_MDS/"
  
}

# Part 1 ---------------------------------------------------------------------

# read file ---------------------------------------------------------------
# Normal Batch ---------------------------------------------------------------
## object_1
if (TRUE) {
  input_f_path_SMU_Normal = "/share/home/wyz/Project/result/20240305_scRNA_AML/20240306_scRNA_SMU_normal/scRNA_SMU_Normal_1_origin_rds_combined_qc.rds"
  
  seurat_object_SMU_Normal <- readRDS(input_f_path_SMU_Normal)
  print(seurat_object_SMU_Normal)
  
}

## object_2
if (TRUE) {
  # 读入文件
  input_f_path_GSE116256_Normal = "/share/home/wyz/Project/result/20240305_scRNA_AML/20240305_scRNA_GSE116256_normal_2/scRNA_GSE116256_Normal_1_origin_rds_combined_qc.rds"
  seurat_object_GSE116256_Normal <- readRDS(input_f_path_GSE116256_Normal)
  print(seurat_object_GSE116256_Normal)
  
}

## object_3
if (TRUE) {
  # 读入文件
  seurat_object_GSE120221_Normal = "/share/home/wyz/Project/result/20240305_scRNA_AML/20240311_scRNA_GSE120221_normal_13/scRNA_GSE120221_Normal_13_origin_rds_combined_qc.rds"
  
  seurat_object_GSE120221_Normal <- readRDS(seurat_object_GSE120221_Normal)
  print(seurat_object_GSE120221_Normal)
  
}


## object_4
if (TRUE) {
  input_f_path_HRA000084_Normal = "/share/home/wyz/Project/result/20240305_scRNA_AML/20240307_scRNA_SMU_HRA000084_normal/scRNA_HRA000084_Normal_2_origin_rds_combined_qc.rds"
  seurat_object_HRA000084_Normal <- readRDS(input_f_path_HRA000084_Normal)
  print(seurat_object_HRA000084_Normal)
}

## object_5
if (TRUE) {
  input_f_path_GSE223060_Normal = "/share/home/wyz/Project/result/20240305_scRNA_AML/20240311_scRNA_GSE223060_normal_3/scRNA_GSE223060_Normal_3_origin_rds_combined_qc.rds"
  seurat_object_GSE223060_Normal <- readRDS(input_f_path_GSE223060_Normal)
  print(seurat_object_GSE223060_Normal)
  
}

# AML Batch ---------------------------------------------------------------

## object_1
if (TRUE) {
  # 读入文件
  input_f_path_SMU = "/share/home/wyz/Project/result/20240316_scRNA_AML/20240316_scRNA_SMU_AML_MDS/scRNA_SMU_AML_MDS_24_origin_rds_combined_qc.rds"
  
  seurat_object_SMU <- readRDS(input_f_path_SMU)
  print(seurat_object_SMU)
}

## object_2
if (TRUE) {
  # 读入文件
  input_f_path_GSE116256 = "/share/home/wyz/Project/result/20240316_scRNA_AML/20240316_scRNA_GSE116256_AML/scRNA_GSE116256_AML_33_origin_rds_combined_qc.rds"
  seurat_object_GSE116256 <- readRDS(input_f_path_GSE116256)
  print(seurat_object_GSE116256)
}


## object_3
if (TRUE) {
  # 读入文件
  input_f_path_PRJNA720840_GSE205490_HRA000084 = "/share/home/wyz/Project/result/20240316_scRNA_AML/20240316_scRNA_PRJNA720840_GSE205490_HRA000084_AML_MDS/scRNA_PRJNA720840_GSE205490_HRA000084_AML_MDS_7_origin_rds_combined_qc.rds"
  seurat_object_PRJNA720840_GSE205490_HRA000084 <- readRDS(input_f_path_PRJNA720840_GSE205490_HRA000084)
  print(seurat_object_PRJNA720840_GSE205490_HRA000084)
}

## object_to_list
if (TRUE) {
  sample_merge_rds_list <- list(
    seurat_object_SMU_Normal,
    seurat_object_GSE116256_Normal,
    seurat_object_GSE120221_Normal,
    seurat_object_HRA000084_Normal,
    seurat_object_GSE223060_Normal,
    seurat_object_SMU,
    seurat_object_GSE116256,
    seurat_object_PRJNA720840_GSE205490_HRA000084
  )
  print(sample_merge_rds_list)
}



## Merge
if (TRUE) {
  seurat_object <- merge(
    x = sample_merge_rds_list[[1]],
    y = c(sample_merge_rds_list[-1]),
    add.cell.ids = c(
      "SMU_Normal",
      "GSE116256_Normal",
      "GSE120221_Normal",
      "HRA000084_Normal",
      "GSE223060_Normal",
      "SMU",
      "GSE116256",
      "PRJNA720840_GSE205490_HRA000084"
    ),
    merge.data = TRUE,
    project = "SMU_seurat_object_PRJNA720840_GSE205490_HRA000084"
  )
  print(seurat_object)
  # 39899 features across 243568 samples within 1 assay
  print(table(seurat_object@meta.data$Source))
}


# select code protein gene -------------------------------------------------------------

# Extract counts
counts_dgCMatrix <- GetAssayData(object = seurat_object, slot = "counts")
print(counts_dgCMatrix)
# 39899 x 243568 sparse Matrix of class "dgCMatrix"

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero_dgCMatrix <- counts_dgCMatrix > 0
print(nonzero_dgCMatrix)
# 39899 x 243568 sparse Matrix of class "lgCMatrix"

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes_logivctor <- Matrix::rowSums(nonzero_dgCMatrix) >= 10
print(keep_genes_logivctor)
head(keep_genes_logivctor)
print(class(keep_genes_logivctor))
print(length(keep_genes_logivctor))
# [1] 39899


# Only keeping those genes expressed in more than 10 cells
filtered_counts_dgCMatrix <- counts_dgCMatrix[keep_genes_logivctor, ]
print(filtered_counts_dgCMatrix)
# 30510 x 243568 sparse Matrix of class "dgCMatrix"



## get code gene list
gene_sheet_f  = "/share/home/wyz/Lab/Database/2023/Gencode/hg38/20240320_genecode_protein_coding_gene/gencode.v44.basic.annotation.split.gene.protein_coding.xlsx"

## read
gene_sheet_df = readxl::read_excel(
  path = gene_sheet_f,
  col_name = T,
  col_types = "text",
  na = ""
)
dim(gene_sheet_df)
print(gene_sheet_df)
class(gene_sheet_df)

# get column names
colname_list = colnames(gene_sheet_df)
class(colname_list)
length(colname_list)
print(colname_list)

## gene_name length
pc_gene_list = as.vector(gene_sheet_df$gene_name)
print(length(pc_gene_list))
# [1] 20046
print(pc_gene_list)


# conduction vector
filtered_counts_dgCMatrix_genename = filtered_counts_dgCMatrix@Dimnames[[1]]
print(filtered_counts_dgCMatrix_genename)
print(length(filtered_counts_dgCMatrix_genename))
# [1] 30510

keep_pc_genes_logi = filtered_counts_dgCMatrix_genename %in% pc_gene_list
head(keep_pc_genes_logi)
print(class(keep_pc_genes_logi))
table(keep_pc_genes_logi)
# FALSE  TRUE
# 13891 16619

names(keep_pc_genes_logi) = filtered_counts_dgCMatrix_genename
print(keep_pc_genes_logi)
print(class(keep_pc_genes_logi))
# [1] "logical"
table(keep_pc_genes_logi)



filtered_counts_keep_pcgene_dgCMatrix <- filtered_counts_dgCMatrix[keep_pc_genes_logi, ]
print(filtered_counts_keep_pcgene_dgCMatrix)
# 16619 x 243568 sparse Matrix of class "dgCMatrix"

# Reassign to filtered Seurat object
seurat_object_filtered <- CreateSeuratObject(filtered_counts_keep_pcgene_dgCMatrix, meta.data = seurat_object@meta.data)

print(seurat_object_filtered)
# 16619 features across 243568 samples within 1 assay


# Part 2 ---------------------------------------------------------------------

# normalization ---------------------------------------------------------------------
if (TRUE) {
  seurat_object_STD = NormalizeData(
    seurat_object_filtered,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )
  seurat_object_STD = FindVariableFeatures(
    seurat_object_STD,
    assay = "RNA",
    selection.method = "vst",
    nfeatures = 3000
  )
  
  # clean up
  if (TRUR) {
    rm(seurat_object_SMU_Normal)
    rm(seurat_object_GSE116256_Normal)
    rm(seurat_object_GSE120221_Normal)
    rm(seurat_object_HRA000084_Normal)
    rm(seurat_object_GSE223060_Normal)
    rm(seurat_object_SMU)
    
    rm(seurat_object_GSE116256)
    rm(seurat_object_PRJNA720840_GSE205490_HRA000084)
    
    
    rm(sample_merge_rds_list)
    rm(seurat_object)
  }
  
  
  if (TRUE) {
    library(future)
    nbrOfWorkers()
    # plan(multisession,workers=100)
    plan(multisession)
    nbrOfWorkers()
    # options(future.globals.maxSize = "1000GB")
    options(future.globals.maxSize = 1000 * 1024^9)
    getOption("future.globals.maxSize")
  }
  
  # long time
  seurat_object_STD = ScaleData(
    seurat_object_STD,
    assay = "RNA",
    features = rownames(seurat_object_STD),
    vars.to.regress = "percent.mt",
    model.use = "linear"
  )
  
  #
  seurat_object_STD = RunPCA(
    seurat_object_STD,
    assay = "RNA",
    reduction.name = "pca",
    npcs = 50,
    verbose = TRUE
  )
  seurat_object_STD = RunUMAP(
    seurat_object_STD,
    dims = 1:30,
    reduction = "pca",
    reduction.name = "std.umap.unintegrated",
    assay = "RNA",
    verbose = TRUE
  )
}


if (TRUE) {
  # output
  file_name = paste0(
    "scRNA_SMU_GSE116256_PRJNA720840_GSE205490_HRA000084_AML_MDS_83_origin_rds_combined_qc_STD",
    ".rds"
  )
  print(file_name)
  saveRDS(seurat_object_STD, file.path(output_dir, file_name))
  
  # read
  # seurat_object_STD = readRDS(file.path(output_dir, file_name))
}
print(seurat_object_STD)

# Integration before
## umap1
if (TRUE) {
  seurat_object_plot =
    DimPlot(
      seurat_object_STD,
      reduction = "std.umap.unintegrated",
      group.by = c("Source"),
      combine = TRUE,
      label.size = 4,
      pt.size = 0.2,
      label = TRUE,
      ncol = 1,
      raster = FALSE
    ) +
    ggtitle("QC-STD_map_unintegrated_Sample_ID")
  seurat_object_plot
  
  output_f_path = file.path(output_dir, "QC-STD_harmony_umap_unintegrated_Source.pdf")
  ggsave(
    filename = output_f_path,
    plot = seurat_object_plot,
    device = pdf,
    width = 297,
    height = 210,
    units = "mm"
  )
}

## umap2
if (TRUE) {
  seurat_object_plot =
    DimPlot(
      seurat_object_STD,
      reduction = "std.unintegrated.umap",
      group.by = c("Source"),
      split.by = c("Sample_ID"),
      combine = TRUE,
      label.size = 4,
      pt.size = 0.2,
      label = TRUE,
      ncol = 5,
      raster = FALSE
    ) +
    ggtitle("QC-STD_umap_unintegrated_Sample_ID")
  seurat_object_plot
  
  output_f_path = file.path(output_dir,
                            "QC-STD_harmony_umap_unintegrated_Sample_ID.pdf")
  ggsave(
    filename = output_f_path,
    plot = seurat_object_plot,
    device = pdf,
    width = 297,
    height = 210,
    units = "mm"
  )
}


# Part 3 ---------------------------------------------------------------------
# Data Integration ---------------------------------------------------------------------
if (TRUE) {
  ## harmony
  seurat_object_STD_harmony <- RunHarmony(
    seurat_object_STD,
    group.by.vars = "Sample_ID",
    assay.use = "RNA",
    max_iter = 20
  )
}


## save rds
if (TRUE) {
  # output
  file_name = paste0(
    "scRNA_SMU_GSE116256_PRJNA720840_GSE205490_HRA000084_AML_MDS_83_origin_rds_combined_qc_STD_Harmony",
    ".rds"
  )
  print(file_name)
  file.path(output_dir, file_name)
  saveRDS(seurat_object_STD_harmony,
          file.path(output_dir, file_name))
  
  ## read
  seurat_object_STD_harmony = readRDS(file.path(output_dir, file_name))
  print(seurat_object_STD_harmony)
}


## dimensionality reduction
if (TRUE) {
  seurat_object_STD_harmony <- FindNeighbors(
    seurat_object_STD_harmony,
    dims = 1:30,
    reduction = "harmony",
    assay = "RNA",
    verbose = TRUE
  )
  
  resolution_seq = append(seq(0.4, 2, 0.1), c(3, 4, 5))
  length(resolution_seq)
  
  if (TRUE) {
    library(future)
    availableCores()
    nbrOfWorkers()
    plan(multisession)
    nbrOfWorkers()
  }
  
  seurat_object_STD_harmony <- FindClusters(
    seurat_object_STD_harmony,
    resolution = resolution_seq,
    assay = "RNA",
    verbose = TRUE
  )
  
  
  ## Find the best resolution
  p_cutree_res = clustree(seurat_object_STD_harmony@meta.data,
                          prefix = "RNA_snn_res.",
                          node_label = "RNA_snn_res.")
  p_cutree_res
  
  
  output_f_path = file.path(output_dir,
                            "QC-seurat_object_STD_harmony_p_cutree_res.pdf")
  ggsave(
    filename = output_f_path,
    plot = p_cutree_res,
    device = pdf,
    width = 297,
    height = 210,
    units = "mm"
  )
  
  
  p_cutree_res_nolabel = clustree(seurat_object_STD_harmony@meta.data, prefix = "RNA_snn_res.")
  p_cutree_res_nolabel
  
  
  output_f_path = file.path(output_dir,
                            "QC-seurat_object_STD_harmony_p_cutree_res_nolabel.pdf")
  ggsave(
    filename = output_f_path,
    plot = p_cutree_res,
    device = pdf,
    width = 297,
    height = 210,
    units = "mm"
  )
  
  seurat_object_STD_harmony <- RunUMAP(
    seurat_object_STD_harmony,
    dims = 1:30,
    reduction = "harmony",
    reduction.name = "std.umap.harmony",
    assay = "RNA",
    verbose = TRUE
  )
  seurat_object_STD_harmony <- RunTSNE(
    seurat_object_STD_harmony,
    dims = 1:30,
    reduction = "harmony",
    reduction.name = "std.tsne.harmony",
    assay = "RNA",
    verbose = TRUE
  )
}


## save rds
if (TRUE) {
  # output
  file_name = paste0(
    "scRNA_SMU_GSE116256_PRJNA720840_GSE205490_HRA000084_AML_MDS_83_origin_rds_combined_qc_STD_Harmony_res",
    ".rds"
  )
  print(file_name)
  saveRDS(seurat_object_STD_harmony,
          file.path(output_dir, file_name))
}

## umap res1
if (TRUE) {
  seurat_object_plot =
    DimPlot(
      seurat_object_STD_harmony,
      reduction = "std.umap.harmony",
      group.by = c("RNA_snn_res.1"),
      split.by = c("Source"),
      combine = TRUE,
      label.size = 4,
      pt.size = 0.5,
      label = TRUE,
      ncol = 4,
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


## umap split sample

table(seurat_object_STD_harmony@meta.data$Sample_ID)

if (TRUE) {
  seurat_object_plot =
    DimPlot(
      seurat_object_STD_harmony,
      reduction = "std.umap.harmony",
      group.by = c("RNA_snn_res.1"),
      split.by = c("Sample_ID"),
      combine = TRUE,
      label.size = 2,
      pt.size = 0.5,
      label = TRUE,
      ncol = 12,
      raster = FALSE
    ) + NoLegend() +
    ggtitle("QC-STD_harmony_umap_integrated_res1")
  seurat_object_plot
  
  output_f_path = file.path(output_dir,
                            "QC-STD_umap_integrated_res1_split_Sample_ID_ALL.pdf")
  ggsave(
    filename = output_f_path,
    plot = seurat_object_plot,
    device = pdf,
    width = 297,
    height = 210,
    units = "mm"
  )
}


## umap split DiseaseClass
if (TRUE) {
  # 调整顺序
  seurat_object_STD_harmony@meta.data$DiseaseClass = factor(
    x = seurat_object_STD_harmony@meta.data$DiseaseClass,
    levels = c("Normal", "MDS", "tAML", "AML")
  )
  print(seurat_object_STD_harmony@meta.data$DiseaseClass)
  print(table(seurat_object_STD_harmony@meta.data$DiseaseClass))
  
  seurat_object_plot =
    DimPlot(
      seurat_object_STD_harmony,
      reduction = "std.umap.harmony",
      group.by = c("RNA_snn_res.1"),
      split.by = c("DiseaseClass"),
      combine = TRUE,
      label.size = 4,
      pt.size = 0.5,
      label = TRUE,
      ncol = 4,
      raster = FALSE
    ) +
    ggtitle("QC-STD_harmony_umap_integrated_res1")
  seurat_object_plot
  
  output_f_path = file.path(output_dir,
                            "QC-STD_umap_integrated_res1_split_DiseaseClass.pdf")
  ggsave(
    filename = output_f_path,
    plot = seurat_object_plot,
    device = pdf,
    width = 297,
    height = 210,
    units = "mm"
  )
}


## tsne
if (TRUE) {
  seurat_object_plot =
    DimPlot(
      seurat_object_STD_harmony,
      reduction = "std.tsne.harmony",
      group.by = c("RNA_snn_res.1"),
      split.by = c("Source"),
      combine = TRUE,
      label.size = 4,
      pt.size = 0.5,
      label = TRUE,
      ncol = 4,
      raster = FALSE
    ) +
    ggtitle("QC-STD_harmony_tsne_integrated_res1")
  seurat_object_plot
  
  output_f_path = file.path(output_dir, "QC-STD_tsne_integrated_res1.pdf")
  ggsave(
    filename = output_f_path,
    plot = seurat_object_plot,
    device = pdf,
    width = 297,
    height = 210,
    units = "mm"
  )
}


## tsne split sample
if (TRUE) {
  seurat_object_plot =
    DimPlot(
      seurat_object_STD_harmony,
      reduction = "std.tsne.harmony",
      group.by = c("RNA_snn_res.1"),
      split.by = c("Sample_ID"),
      combine = TRUE,
      label.size = 2,
      pt.size = 0.5,
      label = TRUE,
      ncol = 12,
      raster = FALSE
    ) + NoLegend() +
    ggtitle("QC-STD_harmony_tsne_integrated_res1")
  seurat_object_plot
  
  output_f_path = file.path(output_dir,
                            "QC-STD_tsne_integrated_res1_split_Sample_ID.pdf")
  ggsave(
    filename = output_f_path,
    plot = seurat_object_plot,
    device = pdf,
    width = 297,
    height = 210,
    units = "mm"
  )
}

## tsne split Source
if (TRUE) {
  seurat_object_plot =
    DimPlot(
      seurat_object_STD_harmony,
      reduction = "std.tsne.harmony",
      group.by = c("RNA_snn_res.1"),
      split.by = c("Source"),
      combine = TRUE,
      label.size = 4,
      pt.size = 0.2,
      label = FALSE,
      ncol = 1,
      raster = FALSE
    ) +
    ggtitle("QC-STD_harmony_tsne_integrated_res1.6")
  seurat_object_plot
  
  output_f_path = file.path(output_dir,
                            "QC-STD_tsne_integrated_res1_split_Source.pdf")
  ggsave(
    filename = output_f_path,
    plot = seurat_object_plot,
    device = pdf,
    width = 297,
    height = 210,
    units = "mm"
  )
}


# Part 4 ------------------------------------------------------------

# Split by cell type ------------------------------------------------------------
# Draw a clustering heatmap -------------------------------------------------------------------

if (TRUE) {
  ## 1）Read in the gene table
  sample_sheet_f = "/share/home/wyz/Lab/Code/ProjectCode/2024/202403/20240305_scRNA_AML/sample_info/2019cell_marker.xlsx"
  
  ## Read in the sample information table
  sample_sheet_df = readxl::read_excel(
    path = sample_sheet_f,
    col_name = T,
    col_types = "text",
    na = ""
  )
  dim(sample_sheet_df)
  print(sample_sheet_df)
  class(sample_sheet_df)
  
  ## 2）factor sequence
  seurat_object_STD_harmony_temp = seurat_object_STD_harmony
  seurat_object_STD_harmony_temp@meta.data[["CellTypeRef"]] = factor(
    x = seurat_object_STD_harmony_temp@meta.data[["CellTypeRef"]],
    levels = c(
      "HSC",
      "Prog",
      "GMP",
      "ProMono",
      "Mono",
      "cDC",
      "pDC",
      "earlyEry",
      "lateEry",
      "ProB",
      "B",
      "Plasma",
      "T",
      "CTL",
      "NK"
    )
  )
  
  ## 3）plot
  seurat_object_plot =
    DoHeatmap(
      seurat_object_STD_harmony_temp,
      features = as.character(unique(sample_sheet_df$Marker基因)),
      group.by = "CellTypeRef",
      assay = 'RNA',
      size = 5.5,
      draw.lines = TRUE,
      lines.width = 2,
      group.bar.height = 0.02,
      combine = TRUE,
      group.colors = c("#00BFC4", "#AB82FF", "#00CD00", "#C77CFF")
    ) +
    scale_fill_gradientn(colors = c("white", "LightSkyBlue1", "firebrick3")) +
    ggtitle("QC-SCT_harmony_CellTypeRef_DoHeatmap")
  
  
  output_f_path = file.path(output_dir, "QC-SCT_harmony_CellTypeRef_DoHeatmap.pdf")
  ggsave(
    filename = output_f_path,
    plot = seurat_object_plot,
    device = pdf,
    width = 297,
    height = 210,
    units = "mm"
  )
}