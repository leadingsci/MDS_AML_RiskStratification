rm(list = ls())

.libPaths()
options(warn = 1)

# Load Package ------------------------------------------------------------
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
  
  library(survminer) # 加载包
  library(survival) # 加载包
}


# input path --------------------------------------------------------------
if (TRUE) {
  ## output
  output_dir = "/share/home/wyz/Project/result/2025/202505/20250527_1_cibersortx_immune_cell_MDS_sAML_AML_merge/20250527_6_immune_cell_kmean_MDS_sAML_AML_survival"
}


# input path --------------------------------------------------------------
{
  ## experimental group
  cell_f = "/share/home/wyz/Project/result/2025/202505/20250527_1_cibersortx_immune_cell_MDS_sAML_AML_merge/20250527_5_tumor_cell_kmean_MDS_sAML_AML_data/20250527_1_tumor_cell_kmean_MDS_sAML_AML_data_survival_1317_sample_del_unrice_MDS.xls"
}


{
  ## read file
  cell_df = read.table(
    cell_f,
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = FALSE
  )
  print(head(cell_df))
  dim(cell_df)
}

{
  colnames(cell_df) <- gsub("-", "_", colnames(cell_df))
  print(colnames(cell_df))
}

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

if (TRUE) {
  celltype = "kmean_k4"
  
  fit_Group <- survfit(Surv(Survival_days, Status_ggsurvfit) ~ kmean_k4, data = cell_df)
  
  #summary(fit_Group)
  
  surv_object = ggsurvplot(
    fit_Group,
    data = cell_df,
    surv.median.line = "hv",
    conf.int = FALSE,
    conf.int.alpha = 0.1,
    risk.table = TRUE,
    align = "hv",
    pval = TRUE,
    #palette = "hue",
    # palette = c("blue","red"),
    palette = c("#66abda", "#e19793", "#ca8282", "#f5d966"),
    title = paste0("MDS_AML_sAML : ", celltype),
    #
    xlab = "Follow up time(d)",
    #legend = c(0.8,0.75), #
    legend = "right",
    #
    legend.title = "celltype",
    #
    break.x.by = 365,
    #
    risk.table.y.text.col = TRUE,
    #
    risk.table.x.text.angle = 0,
    #
    risk.table.height = 0.2,
    #
    font.main = c(8, "bold", "black")
  )
  
  library(patchwork)
  plot_object = surv_object$plot / surv_object$table +
    plot_layout(heights = c(5, 1))
  #plot_object
  
  plot_f = paste0("20250704_3_1_MDS_AML_sAML_immune_percent_survival_",
                  celltype,
                  ".pdf")
  output_f_path = file.path(output_dir, plot_f)
  ggsave(
    filename = output_f_path,
    plot = plot_object,
    device = pdf,
    width = 170,
    height = 165,
    units = "mm"
  )
  
}

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------

{
  fig_Group3 <- pairwise_survdiff(Surv(Survival_days, Status_ggsurvfit) ~ kmean_k4, data = cell_df)
  pvalue_df = fig_Group3$p.value
  
  output_f_path = file.path(
    output_dir,
    "20250521_2_20250521_1_immune_cell_kmean_MDS_sAML_AML_survival_kmean4.xls"
  )
  write.table(
    pvalue_df,
    file = output_f_path,
    sep = "\t",
    row.names = T,
    quote = FALSE
  )
}