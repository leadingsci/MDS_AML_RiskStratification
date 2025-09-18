
# -----------------------------------------------------------------------------
# Author: [Yongzhang Yu]
# Description: [Fig_6b_co_classification_risk_prediction_model_nomogram]
# -----------------------------------------------------------------------------


rm(list = ls())

.libPasts()
.libPaths(c(
  "/share/home/wyz/miniconda3/envs/r451base_rna/lib/R/library"
))
.libPaths()

library(rms)
library(survival)
setwd(
  "/share/home/wyz/Project/result/2025/202507/20250718_1_Risk_Model_Build/20250828_2_Risk_Model_Build_Union_420_nomogram"
) # Set working directory

output_dir = "/share/home/wyz/Project/result/2025/202507/20250718_1_Risk_Model_Build/20250828_2_Risk_Model_Build_Union_420_nomogram"


output_f = "/share/home/wyz/Project/result/2025/202507/20250718_1_Risk_Model_Build/20250723_6_Risk_Model_Build_Union_420/20250723_6_Risk_Model_Build_Union_420_all_workspace.Rdata"
# save.image(file=output_f)
load(output_f)

dd <- datadist(trainRiskOut)
options(datadist = "dd")


f <- cph(
  Surv(Survival_days, Status_Cox) ~ tumor_to_immune_mark,
  data = trainRiskOut,
  x = TRUE,
  y = TRUE,
  surv = TRUE,
  time.inc = 1
) # Using 365 days = 1 year

surv <- Survival(f)

nom <- nomogram(
  f,
  fun = list(function(x)
    surv(0.5, x), function(x)
      surv(1, x), function(x)
        surv(2, x), function(x)
          surv(3, x), function(x)
            surv(4, x), function(x)
              surv(5, x)),
  funlabel = c(
    "0.5-year survival",
    "1-year survival",
    "2-year survival",
    "3-year survival",
    "4-year survival",
    "5-year survival"
  ),
  lp = FALSE,
  maxscale = 100,
  fun.at = c(0.99, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05)
)

pdf(file = "20250720_1_15_Module_One_Riskscore_Group_Ref_Feature_Union_420.pdf",
    height = 8,
    width = 15)
plot(nom)
dev.off()


##########  C-index ###########
#install.packages("Hmisc")
library(Hmisc)

Cindex <- rcorrcens(trainRiskOut$Status_Cox ~ predict(mylog))
Cindex


library(survival)

multiCoxSum$concordance

> multiCoxSum$concordance
C     se(C)
# 0.7723346 0.0209568


################# C-index #######
# 0.5, the model has no predictive capability
# 0.5-0.7, poor accuracy
# 0.71-0.9, moderate accuracy
# &gt;0.9, high accuracy


### cutoff value
# Use "surv_cutpoint" to determine the optimal cut-off value
cutpoint_result <-
  surv_cutpoint(trainRiskOut,
                time = "Survival_days",
                event = "Status_Cox",
                variables = "tumor_to_immune_mark")