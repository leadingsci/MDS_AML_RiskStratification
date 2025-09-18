

# init --------------------------------------------------------------------

rm(list = ls())


# libPaths template
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


# max core ------------------------------------------------------------------
library(future)
if (nbrOfWorkers() < 40) {
  plan(multisession, workers = nbrOfWorkers())
} else{
  plan(multisession, workers = 40)
}


# input path --------------------------------------------------------------

if (TRUE) {
  ## sample info
  sample_sheet_f = "/share/home/wyz/Lab/Code/ProjectCode/2024/202403/20240305_scRNA_AML/sample_info/scRNA_Concat_Colname_Ref_20240316.xlsx"
  
  ## input_dir
  input_dir = "/share/home/wyz/Project/result/20240316_scRNA_AML/20240320_scRNA_GSE_Rice/"
  
  ## output
  output_dir =  "/share/home/wyz/Project/result/20240316_scRNA_AML/20240320_scRNA_GSE_Rice/"
  
}

# Part 1 ---------------------------------------------------------------------

# read file ---------------------------------------------------------------
## object_1
if (TRUE) {
  input_f_path_SMU = "/share/home/wyz/Project/result/20240316_scRNA_AML/20240316_scRNA_SMU_Rice/scRNA_SMU_Rice_13_origin_rds_combined_qc.rds"
  
  seurat_object_SMU <- readRDS(input_f_path_SMU)
  print(seurat_object_SMU)
  
}

## object_2
if (TRUE) {
  # 读入文件
  input_f_path_GSE116256 = "/share/home/wyz/Project/result/20240316_scRNA_AML/20240318_scRNA_GSE116256_Rice/scRNA_GSE116256_Rice_1_origin_rds_combined_qc.rds"
  seurat_object_GSE116256 <- readRDS(input_f_path_GSE116256)
  print(seurat_object_GSE116256)
}


## object_3
if (TRUE) {
  # 读入文件
  input_f_path_HRA000084 = "/share/home/wyz/Project/result/20240316_scRNA_AML/20240316_scRNA_HRA000084_Rice/scRNA_HRA000084_Rice_4_origin_rds_combined_qc.rds"
  seurat_object_HRA000084 <- readRDS(input_f_path_HRA000084)
  print(seurat_object_HRA000084)
}




## object_to_list
if (TRUE) {
  sample_merge_rds_list <- list(seurat_object_SMU,
                                seurat_object_GSE116256,
                                seurat_object_HRA000084)
  print(sample_merge_rds_list)
}

## Merge
if (TRUE) {
  seurat_object <- merge(
    x = sample_merge_rds_list[[1]],
    y = c(sample_merge_rds_list[-1]),
    add.cell.ids = c("SMU", "GSE116256", "HRA000084"),
    merge.data = TRUE,
    project = "SMU_GSE116256_HRA000084"
  )
}

# select code protein gene -------------------------------------------------------------
# select gene -------------------------------------------------------------

# Extract counts
counts_dgCMatrix <- GetAssayData(object = seurat_object, slot = "counts")
print(counts_dgCMatrix)
# 39899 x 243568 sparse Matrix of class "dgCMatrix"

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero_dgCMatrix <- counts_dgCMatrix > 0
print(nonzero_dgCMatrix)
# 38714 x 67452 sparse Matrix of class "lgCMatrix"

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
# 至少10存在10个转录本的基因
keep_genes_logivctor <- Matrix::rowSums(nonzero_dgCMatrix) >= 10
print(keep_genes_logivctor)
head(keep_genes_logivctor)
print(class(keep_genes_logivctor))
print(length(keep_genes_logivctor))
# [1] 38714
table(keep_genes_logivctor)
# FALSE  TRUE
# 11679 27035

# Only keeping those genes expressed in more than 10 cells
filtered_counts_dgCMatrix <- counts_dgCMatrix[keep_genes_logivctor, ]
print(filtered_counts_dgCMatrix)
# 27035 x 67452 sparse Matrix of class "dgCMatrix"



## get code protein gene
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

## gene_name length
colname_list = colnames(gene_sheet_df)
class(colname_list)
length(colname_list)
print(colname_list)

# conduction vector
pc_gene_list = as.vector(gene_sheet_df$gene_name)
print(length(pc_gene_list))
# [1] 20046
print(pc_gene_list)

#
filtered_counts_dgCMatrix_genename = filtered_counts_dgCMatrix@Dimnames[[1]]
print(filtered_counts_dgCMatrix_genename)
print(length(filtered_counts_dgCMatrix_genename))
# [1] 27035

keep_pc_genes_logi = filtered_counts_dgCMatrix_genename %in% pc_gene_list
head(keep_pc_genes_logi)
print(class(keep_pc_genes_logi))
print(length(keep_pc_genes_logi))
# [1] 27035


names(keep_pc_genes_logi) = filtered_counts_dgCMatrix_genename
print(keep_pc_genes_logi)
print(class(keep_pc_genes_logi))
# [1] "logical"
table(keep_pc_genes_logi)
# FALSE  TRUE
# 10376 16659


filtered_counts_keep_pcgene_dgCMatrix <- filtered_counts_dgCMatrix[keep_pc_genes_logi, ]
print(filtered_counts_keep_pcgene_dgCMatrix)
# 16659 x 67452 sparse Matrix of class "dgCMatrix"

# Reassign to filtered Seurat object
seurat_object_filtered <- CreateSeuratObject(filtered_counts_keep_pcgene_dgCMatrix, meta.data = seurat_object@meta.data)

print(seurat_object_filtered)
# 16659 features across 67452 samples within 1 assay


# Part 2 ---------------------------------------------------------------------
# normalization ---------------------------------------------------------------------
if (TRUE) {
  # normalization
  
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
  
  if (TRUE) {
    library(future)
    nbrOfWorkers()
    plan(multisession)
    nbrOfWorkers()
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
    "scRNA_SMU_Rice_14_GSE116256_Rice_1_HRA000084_Rice_4_origin_rds_combined_qc_STD",
    ".rds"
  )
  print(file_name)
  saveRDS(seurat_object_STD, file.path(output_dir, file_name))
  
  # read
  # seurat_object_STD = readRDS(file.path(output_dir, file_name))
}


# Integration before
## uamp
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

## uamp
if (TRUE) {
  seurat_object_plot =
    DimPlot(
      seurat_object_STD,
      reduction = "std.umap.unintegrated",
      group.by = c("Source"),
      split.by = c("Sample_ID"),
      combine = TRUE,
      label.size = 4,
      pt.size = 0.5,
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
  table(seurat_object_STD@meta.data[["Sample_ID"]])
  
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
    "scRNA_SMU_Rice_14_GSE116256_Rice_1_HRA000084_Rice_4_origin_rds_combined_qc_STD_Harmony",
    ".rds"
  )
  print(file_name)
  file.path(output_dir, file_name)
  saveRDS(seurat_object_STD_harmony,
          file.path(output_dir, file_name))
  
  ## read
  # seurat_object_STD_harmony = readRDS(file.path(output_dir, file_name))
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
    cluster.name = "harmony.integrated.cluster",
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
    plot = p_cutree_res_nolabel,
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
    "scRNA_SMU_Rice_14_GSE116256_Rice_1_HRA000084_Rice_4_origin_rds_combined_qc_STD_Harmony_res",
    ".rds"
  )
  print(file_name)
  saveRDS(seurat_object_STD_harmony,
          file.path(output_dir, file_name))
}

## uamp res1
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
      ncol = 3,
      raster = FALSE
    ) +
    ggtitle("QC-STD_harmony_umap_integrated_res1")
  seurat_object_plot
  
  output_f_path = file.path(output_dir, "QC-STD_umap_integrated_res1.pdf")
  ggsave(
    filename = output_f_path,
    plot = seurat_object_plot,
    device = pdf,
    width = 297,
    height = 210,
    units = "mm"
  )
}


## uamp split celltype
if (TRUE) {
  seurat_object_plot =
    DimPlot(
      seurat_object_STD_harmony,
      reduction = "std.umap.harmony",
      group.by = c("CellTypeRef"),
      split.by = c("Source"),
      combine = TRUE,
      label.size = 4,
      pt.size = 0.5,
      label = TRUE,
      ncol = 3,
      raster = FALSE
    ) +
    ggtitle("QC-STD_harmony_umap_integrated_res1.5")
  seurat_object_plot
  
  output_f_path = file.path(output_dir,
                            "QC-STD_umap_integrated_res1_split_CellTypeRef.pdf")
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
      ncol = 3,
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


## tsne split CellTypeRef
if (TRUE) {
  seurat_object_plot =
    DimPlot(
      seurat_object_STD_harmony,
      reduction = "std.tsne.harmony",
      group.by = c("CellTypeRef"),
      split.by = c("Source"),
      combine = TRUE,
      label.size = 4,
      pt.size = 0.5,
      label = TRUE,
      ncol = 3,
      raster = FALSE
    ) +
    ggtitle("QC-STD_harmony_tsne_integrated_res1.3")
  seurat_object_plot
  
  output_f_path = file.path(output_dir,
                            "QC-STD_tsne_integrated_res1.3_split_CellTypeRef.pdf")
  ggsave(
    filename = output_f_path,
    plot = seurat_object_plot,
    device = pdf,
    width = 297,
    height = 210,
    units = "mm"
  )
}


# 按样本进行split -----------------------------------------------------------
## umap
if (TRUE) {
  seurat_object_plot =
    DimPlot(
      seurat_object_STD_harmony,
      reduction = "STD_integrated_umap_harmony",
      group.by = c("CellTypeRef"),
      split.by = "Sample_ID",
      combine = TRUE,
      label.size = 4,
      pt.size = 0.5,
      label = TRUE,
      ncol = 2,
      raster = FALSE
    ) +
    
    ggtitle("QC-STD_harmony_umap_integrated_split_Sample_ID") + theme(# legend.text = element_text(size = 4),  # 调整图例文本字体大小
      # legend.title = element_text(size = 4)   # 调整图例标题字体大小)
      
      output_f_path = file.path(output_dir,
                                "QC-STD_harmony_umap_integrated_split_Sample_ID.pdf")
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
      reduction = "STD_integrated_tsne_harmony",
      group.by = c("CellTypeRef"),
      split.by = "Sample_ID",
      combine = TRUE,
      label.size = 4,
      pt.size = 0.5,
      label = TRUE,
      ncol = 2,
      raster = FALSE
    ) +
    
    ggtitle("QC-STD_harmony_tsne_integrated_split_Sample_ID") + theme(# legend.text = element_text(size = 4),  # 调整图例文本字体大小
      # legend.title = element_text(size = 4)   # 调整图例标题字体大小)
      
      output_f_path = file.path(output_dir,
                                "QC-STD_harmony_tsne_integrated_split_Sample_ID.pdf")
      ggsave(
        filename = output_f_path,
        plot = seurat_object_plot,
        device = pdf,
        width = 297,
        height = 210,
        units = "mm"
      )
      
}