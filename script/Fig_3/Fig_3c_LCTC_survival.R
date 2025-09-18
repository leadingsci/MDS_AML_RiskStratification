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
  output_dir = "/share/home/wyz/Project/result/2025/202505/20250520_1_cibersortx_tumor_cell_MDS_sAML_AML_merge/20250520_5_tumor_cell_kmean_MDS_sAML_AML_survival"
}


# input path --------------------------------------------------------------
{
  ## experimental group
  cell_f = "/share/home/wyz/Project/result/2025/202505/20250520_1_cibersortx_tumor_cell_MDS_sAML_AML_merge/20250520_4_tumor_cell_kmean_MDS_sAML_AML_data/20250520_1_tumor_cell_kmean_MDS_sAML_AML_data_1528_sample.xls"
}


{
  ## read file
  cell_df = read.table(
    cell_f,
    sep = "\t",
    header = T,
    row.names = 1,
    check.names = FALSE
  ) # 阻止列名中非法符号-号变为.号
  print(head(cell_df))
  dim(cell_df) #
}

{
  colnames(cell_df) <- gsub("-", "_", colnames(cell_df))
  print(colnames(cell_df))
}


# -------------------------------------------------------------------------

if (TRUE) {
  celltype = "kmean_k5"
  fit_Group <- survfit(Surv(Survival_days, Status_ggsurvfit) ~ kmean_k5, data = cell_df) #
  summary(fit_Group)
  
  surv_object = ggsurvplot(
    fit_Group,
    data = cell_df,
    surv.median.line = "hv",
    conf.int = FALSE,
    conf.int.alpha = 0.1,
    risk.table = TRUE,
    align = "hv",
    pval = TRUE,
    palette =  c("#ca8282", "#66abda", "#e19793", "#898989", "#f5d966"),
    title = paste0("MDS_AML_sAML : ", celltype),
    #
    xlab = "Follow up time(d)",
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
    font.main = c(8, "bold", "black")  #
  )
  
  library(patchwork)
  plot_object = surv_object$plot / surv_object$table +
    plot_layout(heights = c(5, 1))
  
  plot_f = paste0("20250703_3_1_MDS_AML_sAML_tumor_percent_survival_",
                  celltype,
                  ".pdf")
  output_f_path = file.path(output_dir, plot_f)
  ggsave(
    filename = output_f_path,
    plot = plot_object,
    device = pdf,
    width = 170,
    height = 170,
    units = "mm"
  )
}

# -------------------------------------------------------------------------
{
  colnames(cell_df)
  "kmean_k2"          "kmean_k3"          "kmean_k4"          "kmean_k5"
  [17] "kmean_k6"          "kmean_k7"          "kmean_k8"          "kmean_k9"
}


{
  fig_Group3 <- pairwise_survdiff(Surv(Survival_days, Status_ggsurvfit) ~ kmean_k5, data = cell_df)
  pvalue_df = fig_Group3$p.value
  
  output_f_path = file.path(
    output_dir,
    "20250519_2_20250519_1_tumor_cell_kmean_MDS_sAML_AML_survival_kmean5.xls"
  )
  write.table(
    pvalue_df,
    file = output_f_path,
    sep = "\t",
    row.names = T,
    quote = FALSE
  )
}