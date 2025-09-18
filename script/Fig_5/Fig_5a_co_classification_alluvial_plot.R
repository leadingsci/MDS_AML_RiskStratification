
# -----------------------------------------------------------------------------
# Author: [Yongzhang Yu]
# Description: [Fig_5a_co_classification_alluvial_plot]
# -----------------------------------------------------------------------------


  
  {
    rm(list = ls())
    
    ## libpath
    .libPaths()
    options(warn = 1)
    
  }
  
  
  # Load Package ------------------------------------------------------------
  {
    # library(tidyverse);
    library(readr)
    library(dplyr)
    
    library(ggplot2)
    
    # library(openxlsx);
    library(readxl)
    
    library(RColorBrewer)
    # library(ggpubr)
    library(dplyr)
    
  }
  
  library(ggalluvial)
  
  # input -------------------------------------------------------------------
  if (TRUE) {
    ## input_dir
    cell_f = "/share/home/wyz/Project/result/2025/202506/20250606_3_cibersortx_immune_cell/20250607_1_MDS_AML_tumor_immune_cell_alluvial_plot_data/20250607_1_MDS_AML_tumor_immune_cell_alluvial_plot_data_inter_4795_sample.xls"
    
    ## output
    output_dir = "/share/home/wyz/Project/result/2025/202506/20250606_3_cibersortx_immune_cell/20250607_1_MDS_AML_tumor_immune_cell_alluvial_plot_data"
  }
  
  
  {
    ## read file
    cell_df = read.table(cell_f,
                         sep = "\t",
                         header = T,
                         check.names = FALSE)
    print(head(cell_df))
    dim(cell_df) # [1]  16619 189307
  }
  
  {
    cell_df$Sample_Count = 1
    head(cell_df)
  }
  
  {
    ## 固定顺序
    table(cell_df$kmean_k5_mark_short)
    cell_df$kmean_k5_mark_short
    cell_df$kmean_k5_mark_short = factor(
      cell_df$kmean_k5_mark_short,
      levels = c("GMP", "Mature", "Prog", "Primitive", "intermediate")
    )
    
    
    table(cell_df$kmean_k4_mark_short)
    cell_df$kmean_k4_mark_short
    cell_df$kmean_k4_mark_short = factor(
      cell_df$kmean_k4_mark_short,
      levels =  c(
        "CD8T_High",
        "CD8T_Int_CD4T_Low",
        "CD8_Int_CD4_High",
        "CD8T_Low"
      )
    )
    
    table(cell_df$DiseaseClass)
    unique(cell_df$DiseaseClass)
    # cell_df$DiseaseClass = factor(cell_df$DiseaseClass,levels = c("AML","sAML","MDS"))
    cell_df$DiseaseClass = factor(cell_df$DiseaseClass, levels = c("MDS", "sAML", "AML"))
  }
  
  
  # -------------------------------------------------------------------------
  
  {
    brewer.pal(9, "Set1")
    # [1] "#E41A1C" "#377EB8" "#4DAF4A" "#984EA3" "#FF7F00" "#FFFF33" "#A65628" "#F781BF" "#999999"
    color_list_2 = c("#4DAF4A", "#E41A1C", "#984EA3")
    color_list_3 = c("#E41A1C", "#984EA3", "#377EB8", "#4DAF4A", "#FF7F00")
    
    
    color_list_4 = c("#e19793", "#66abda", "#ca8282", "#f5d966")
    
    color_list_5 = c("#009e73", "#0072b2", "#cc79a7", "#f0e442")
  }
  
  # -------------------------------------------------------------------------
  
  
  
  # -------------------------------------------------------------------------
  # 2025-06-07
  
  # -------------------------------------------------------------------------
  
  
  # -------------------------------------------------------------------------
  
  {
    plot_object = ggplot(
      data = cell_df, aes(
        axis1 = DiseaseClass,
        # First variable on the X-axis
        axis2 = kmean_k5_mark_short,
        # Second variable on the X-axis
        axis3 = kmean_k4_mark_short,
        # Third variable on the X-axis
        y = Sample_Count
      )
    ) +
      geom_alluvium(
        aes(fill = kmean_k4_mark_short),
        curve_type = "cubic",
        width = 0,
        alpha = 0.9,
        color = NA
      ) + # fill kmean_k4_mark_short # DiseaseClass #0.01#,color = NA
      geom_stratum(alpha = 0.4, width = 0.2) + # aes(fill = kmean_k4_mark_short) #
      geom_text(
        stat = "stratum", aes(label = after_stat(stratum)),
        size = 2,
        alpha = 1
      ) + #
      scale_fill_manual(values = color_list_5) +
      # scale_fill_manual(values = colors) +
      theme_void() +
      # theme_minimal() +
      guides(fill = guide_legend(title = "Group"))  +
      ggtitle("tumor to immune")
    
    output_f_path = file.path(
      output_dir,
      "20250704_8_MDS_AML_tumor_immune_cell_alluvial_plot_alpha_0.4.pdf"
    )
    ggsave(
      filename  =  output_f_path,
      plot  =  plot_object,
      device  =  pdf,
      width  =  180,
      height  =  100,
      units  =  "mm"
    )
    
  }
  
  # -------------------------------------------------------------------------
  
  
  # -------------------------------------------------------------------------
  
  {
    plot_object = ggplot(
      data = cell_df, aes(
        axis1 = DiseaseClass,
        # First variable on the X-axis
        axis2 = kmean_k4_mark_short,
        # Second variable on the X-axis
        axis3 = kmean_k5_mark_short,
        # Third variable on the X-axis
        y = Sample_Count
      )
    ) +
      geom_alluvium(
        aes(fill = kmean_k5_mark_short),
        curve_type = "cubic",
        width = 0,
        alpha = 0.9,
        color = NA
      ) + # fill kmean_k4_mark_short # DiseaseClass #0.01#,color = NA
      geom_stratum(alpha = 1, width = 0.2) + # aes(fill = kmean_k4_mark_short) #
      geom_text(
        stat = "stratum", aes(label = after_stat(stratum)),
        size = 2,
        alpha = 1
      ) +
      
      scale_fill_manual(values = color_list_3) +
      # scale_fill_manual(values = colors) +
      theme_void() +
      # theme_minimal() +
      guides(fill = guide_legend(title = "Group"))  +
      ggtitle("immune to tumor")
    
    output_f_path = file.path(output_dir,
                              "20250607_4_MDS_AML_immune_to_tumor_cell_alluvial_plot.pdf")
    ggsave(
      filename  =  output_f_path,
      plot  =  plot_object,
      device  =  pdf,
      width  =  180,
      height  =  100,
      units  =  "mm"
    )
    
  }
  