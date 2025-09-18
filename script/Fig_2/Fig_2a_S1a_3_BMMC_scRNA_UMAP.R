# -----------------------------------------------------------------------------
# Author: [Yongzhang Yu]
# Description: [Fig_2a_S1a_3_BMMC_scRNA_UMAP]
# -----------------------------------------------------------------------------

{
  rm(list = ls())
  
  .libPaths()
  options(warn = 1)
  
  # Load Package ------------------------------------------------------------
  # devtools::install_github("satijalab/seurat", ref = "develop")
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
  #library(ggpubr) #R 4.2.2 #R 4.2.3
  
  library(Seurat)
  #R 4.2.2 #R 4.2.3
  
  library(harmony)
  #R 4.2.3
  library(clustree) #R 4.2.3
  # library(future); #R 4.2.2 #R 4.2.3
  
  library(Matrix)
  packageVersion("Matrix")
  
  library(ComplexHeatmap)
  
}


# input path --------------------------------------------------------------

if (TRUE) {
  ## output
  output_dir = "/share/home/wyz/Project/result/2025/202501/20250118_2_unrice_scRNA_subtype/20250801_1_unrice_scRNA_subtype_color"
}


# Load seurat object ------------------------------------------------------------

if (TRUE) {
  # read file
  input_f_path = "/share/home/wyz/Project/result/20240320_scRNA_AML_Annotation/20240321_scRNA_GSE_GSE116256_annotation/scRNA_SMU_GSE116256_PRJNA720840_GSE205490_HRA000084_AML_MDS_83_origin_rds_combined_qc_STD_Harmony_res_anno_GES116256.rds"
  
  seurat_object <- readRDS(input_f_path)
  print(seurat_object)
  levels(seurat_object)
  head(seurat_object@meta.data)
  table(seurat_object@meta.data$GSE116256_TransferLabel)
  # >   table(seurat_object@meta.data$GSE116256_TransferLabel)
  #
  # B          cDC     cDC-like          CTL     earlyEry          GMP     GMP-like          HSC     HSC-like      lateEry         Mono    Mono-like           NK          pDC       Plasma         ProB
  # 15331         2959         9365        11325         4037         2734        12573         1600         8108        13254        18595        14110        40727         1205         3498         1849
  # Prog    Prog-like      ProMono ProMono-like            T
  # 7592        17234         2561         2109        52802
  
  table(seurat_object@meta.data$Sample_ID)
  table(seurat_object@meta.data$RNA_snn_res.1)
  # 0     1    10    11    12    13    14    15    16    17    18    19     2    20    21    22    23    24    25    26    27    28    29     3    30    31    32    33    34    35    36    37    38    39
  # 31129 28793  6391  6140  6078  4207  4125  3535  3419  3080  3018  2295 28516  2162  1787  1696  1151   826   817   763   654   565   490 26345   292   270   248   198   136   131    86    85    81    74
  # 4    40    41    42     5     6     7     8     9
  # 20218    38    13     2 17565 13692  9285  6660  6512
}


{
  ## 统计样本信息
  table(seurat_object@meta.data$Source)
  
  seurat_object_SMU = subset(seurat_object, subset = Source == "SMU")
  
  meta_SMU = seurat_object_SMU@meta.data
  
  table(meta_SMU$Sample_ID)
  unique(meta_SMU$Sample_ID)
  length(unique(meta_SMU$Sample_ID))
  
  table(meta_SMU$DiseaseClass)
  
  
  table(seurat_object_SMU@meta.data$DiseaseClass)
  table(seurat_object_SMU@meta.data$Disease,
        seurat_object_SMU@meta.data$Sample_ID)
}

{
  seurat_object_public = subset(seurat_object, subset = Source != "SMU")
  meta_public = seurat_object_public@meta.data
  
  table(meta_public$Source)
  unique(meta_public$Sample_ID)
  
  unique(meta_public %>% filter(DiseaseClass == "Normal") %>% select(c("Sample_ID")))
  dim(unique(
    meta_public %>% filter(DiseaseClass == "Normal") %>% select(c("Sample_ID"))
  ))
  
  unique(meta_public %>% filter(DiseaseClass == "MDS") %>% select(c("Sample_ID")))
  dim(unique(meta_public %>% filter(DiseaseClass == "MDS") %>% select(c("Sample_ID"))))
  
  unique(meta_public %>% filter(DiseaseClass == "tAML") %>% select(c("Sample_ID")))
  dim(unique(meta_public %>% filter(DiseaseClass == "tAML") %>% select(c("Sample_ID"))))
  
  unique(meta_public %>% filter(DiseaseClass == "AML") %>% select(c("Sample_ID")))
  dim(unique(meta_public %>% filter(DiseaseClass == "AML") %>% select(c("Sample_ID"))))
}


# -------------------------------------------------------------------------


{
  table(seurat_object_SMU@meta.data$DiseaseClass)
  seurat_object_SMU_MDS_sAML_AML = subset(seurat_object_SMU, subset = DiseaseClass %in% c("MDS", "tAML", "AML"))
  table(seurat_object_SMU_MDS_sAML_AML@meta.data$DiseaseClass)
  seurat_object_SMU_MDS_sAML_AML
}

{
  seurat_object_SMU_HC = subset(seurat_object_SMU, subset = DiseaseClass %in% c("Normal"))
  table(seurat_object_SMU_HC@meta.data$DiseaseClass)
  seurat_object_SMU_HC
}



{
  table(seurat_object_public@meta.data$DiseaseClass)
  seurat_object_public_MDS_sAML_AML = subset(seurat_object_public, subset =
                                               DiseaseClass %in% c("MDS", "tAML", "AML"))
  table(seurat_object_public_MDS_sAML_AML@meta.data$DiseaseClass)
  seurat_object_public_MDS_sAML_AML
}

{
  seurat_object_public_HC = subset(seurat_object_public, subset = DiseaseClass %in% c("Normal"))
  table(seurat_object_public_HC@meta.data$DiseaseClass)
  seurat_object_public_HC
}


# -------------------------------------------------------------------------
{
  #
  seurat_object@meta.data$GSE116256_TransferLabel_20250118 = seurat_object@meta.data$GSE116256_TransferLabel
  
  ## Rename the "like" cells
  seurat_object@meta.data$GSE116256_TransferLabel_20250118 = gsub("HSC-like",
                                                                  "HSC",
                                                                  seurat_object@meta.data$GSE116256_TransferLabel_20250118)
  seurat_object@meta.data$GSE116256_TransferLabel_20250118 = gsub("Prog-like",
                                                                  "Prog",
                                                                  seurat_object@meta.data$GSE116256_TransferLabel_20250118)
  seurat_object@meta.data$GSE116256_TransferLabel_20250118 = gsub("GMP-like",
                                                                  "GMP",
                                                                  seurat_object@meta.data$GSE116256_TransferLabel_20250118)
  seurat_object@meta.data$GSE116256_TransferLabel_20250118 = gsub(
    "ProMono-like",
    "ProMono",
    seurat_object@meta.data$GSE116256_TransferLabel_20250118
  )
  seurat_object@meta.data$GSE116256_TransferLabel_20250118 = gsub("Mono-like",
                                                                  "Mono",
                                                                  seurat_object@meta.data$GSE116256_TransferLabel_20250118)
  seurat_object@meta.data$GSE116256_TransferLabel_20250118 = gsub("cDC-like",
                                                                  "cDC",
                                                                  seurat_object@meta.data$GSE116256_TransferLabel_20250118)
  
  table(seurat_object@meta.data$GSE116256_TransferLabel_20250118)
  # HSC     Prog      GMP  ProMono     Mono      cDC      pDC earlyEry  lateEry     ProB        B   Plasma        T      CTL       NK
  # 9708    24826    15307     4670    32705    12324     1205     4037    13254     1849    15331     3498    52802    11325    40727
}



# -------------------------------------------------------------------------

{
  ## factor
  seurat_object@meta.data$GSE116256_TransferLabel_20250118 = factor(
    seurat_object@meta.data$GSE116256_TransferLabel_20250118,
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
}



# -------------------------------------------------------------------------
{
  seurat_object_HSC = subset(seurat_object, subset = GSE116256_TransferLabel_20250118 ==
                               "HSC")
  seurat_object_Prog = subset(seurat_object, subset = GSE116256_TransferLabel_20250118 ==
                                "Prog")
  seurat_object_GMP = subset(seurat_object, subset = GSE116256_TransferLabel_20250118 ==
                               "GMP")
  
  table(seurat_object_HSC@meta.data$RNA_snn_res.1)
  table(seurat_object_Prog@meta.data$RNA_snn_res.1)
  table(seurat_object_GMP@meta.data$RNA_snn_res.1)
}

# -------------------------------------------------------------------------
{
  ## res1 to new celltype
  celltype_mapping = c(
    "0" = "T",
    "1" = "NK",
    "2" = "Mono",
    "3" = "Prog",
    "4" = "T",
    "5" = "NK",
    "6" = "B",
    "7" = "GMP",
    "8" = "lateEry",
    "9" = "ProB",
    "10" = "cDC",
    "11" = "cDC",
    "12" = "earlyEry",
    "13" = "lateEry",
    "14" = "Mono",
    "15" = "Plasma",
    "16" = "HSC",
    "17" = "lateEry",
    "18" = "T",
    "19" = "ProMono",
    "20" = "B",
    "21" = "Mono",
    "22" = "pDC",
    "23" = "earlyEry",
    "24" = "T",
    "25" = "T",
    "26" = "B",
    "27" = "Prog",
    "28" = "Mono",
    "29" = "lateEry",
    "30" = "lateEry",
    "31" = "lateEry",
    "32" = "NK",
    "33" = "Mono",
    "34" = "Mono",
    "35" = "GMP",
    "36" = "Plasma",
    "37" = "GMP",
    "38" = "T",
    "39" = "Prog",
    "40" = "B",
    "41" = "Prog",
    "42" = "T"
  )
  
  
  # Add a new 'celltype' column based on the 'res1' column
  seurat_object$GSE116256_TransferLabel_20250120 <- celltype_mapping[as.character(seurat_object@meta.data$RNA_snn_res.1)]
  
  table(seurat_object$GSE116256_TransferLabel_20250120)
  
  # Convert numeric 'res' column to the corresponding celltype using the mapping
  # seurat_object$celltype <- cluster_to_celltype[as.character(seurat_object@meta.data$res)]
  
}


{
  ## factor
  seurat_object@meta.data$GSE116256_TransferLabel_20250120 = factor(
    seurat_object@meta.data$GSE116256_TransferLabel_20250120,
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
}




# -------------------------------------------------------------------------



{
  # celltype color
  celltype_colors = c(
    "HSC" = "#e11f26",
    "Prog" = "#fccdac",
    "GMP" = "#e62a89",
    "ProMono" = "#a6cde1",
    "Mono" = "#1f78b4",
    "cDC" = "#b4d789",
    "pDC" = "#329f48",
    "earlyEry" = "#f47e20",
    "lateEry" = "#b05a28",
    "ProB" = "#c9b2d5",
    "B" = "#6b3e97",
    "Plasma" = "#f1c9df",
    "T" = "#b4decb",
    "CTL" = "#b4decb",
    "NK" = "#e4ab23"
  )
}


# -------------------------------------------------------------------------

# -------------------------------------------------------------------------

if (TRUE) {
  ## UMAP
  seurat_object_plot =
    DimPlot(
      seurat_object,
      reduction = "std.umap.harmony",
      group.by = c("GSE116256_TransferLabel_20250120"),
      combine = TRUE,
      label.size = 2,
      pt.size = 0.01,
      label = TRUE,
      ncol = 1,
      cols = celltype_colors,
      raster = FALSE
    ) +
    ggtitle("unrice_scRNA_UMAP") +
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
  
  output_f_path = file.path(output_dir, "20250118_9_unrice_scRNA_subtype_UMAP.pdf")
  ggsave(
    filename = output_f_path,
    plot = seurat_object_plot,
    device = pdf,
    width = 110,
    height = 90,
    units = "mm"
  )
}



if (TRUE) {
  ## UMAP
  seurat_object_plot =
    DimPlot(
      seurat_object,
      reduction = "std.tsne.harmony",
      group.by = c("GSE116256_TransferLabel_20250120"),
      combine = TRUE,
      label.size = 2,
      pt.size = 0.01,
      label = TRUE,
      ncol = 1,
      cols = celltype_colors,
      # 使用自定义颜色
      raster = FALSE
    ) +
    ggtitle("unrice_scRNA_UMAP") +
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
  
  output_f_path = file.path(output_dir, "20250118_10_unrice_scRNA_subtype_tSNE.pdf")
  ggsave(
    filename = output_f_path,
    plot = seurat_object_plot,
    device = pdf,
    width = 110,
    height = 90,
    units = "mm"
  )
}