## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
# library(multiDEGGs)
# library(nestedcv)
# data("synthetic_metadata")
# data("synthetic_rnaseqData")
# 
# # Regularized linear model with interaction pairs only
# fit.glmnet <- nestcv.glmnet(
#   y = as.numeric(synthetic_metadata$response),
#   x = t(synthetic_rnaseqData),
#   modifyX = "multiDEGGs_filter",
#   modifyX_options = list(
#     keep_single_genes = FALSE,
#     nfilter = 20
#   ),
#   modifyX_useY = TRUE,
#   n_outer_folds = 5,
#   n_inner_folds = 6,
#   verbose = FALSE
# )
# 
# summary(fit.glmnet)

## ----echo=FALSE---------------------------------------------------------------
cat("Nested cross-validation with glmnet 
No filter 
Modifier:  multiDEGGs_filter 
Outer loop:  5-fold CV 
Inner loop:  6-fold CV 
100 observations, 14 predictors 
 
 
Final parameters: 
 lambda    alpha   
0.05894  0.10000   

Final coefficients: 
 (Intercept) TNF:TNFRSF1A    AKT2:MTOR   IL1B:IL1R2    FASLG:FAS TGFB3:TGFBR1  
    1.823874    -0.193020    -0.119887     0.052089    -0.035947    -0.033527  
MAP2K2:MAPK3  FANCD2:FAN1  
   -0.021308    -0.008862  
 
Result: 
       RMSE     R.squared   Pearson.r^2           MAE    
    0.47302       0.08148       0.09173       0.44153   " )

## ----eval=FALSE---------------------------------------------------------------
# # Random forest model including both pairs and individual genes
# fit.rf <- nestcv.train(
#   y = synthetic_metadata$response,
#   x = t(synthetic_rnaseqData),
#   method = "rf",
#   modifyX = "multiDEGGs_filter",
#   modifyX_options = list(
#     keep_single_genes = TRUE,
#     nfilter = 30
#   ),
#   modifyX_useY = TRUE,
#   n_outer_folds = 5,
#   n_inner_folds = 6,
#   verbose = FALSE
# )
# 
# fit.rf$summary
# 
# # Plot ROC on outer folds
# plot(fit.rf$roc)

## ----echo=FALSE---------------------------------------------------------------
cat("               Reference
Predicted       Non_responder Responder
  Non_responder            57         2
  Responder                 1        40

              AUC            Accuracy   Balanced accuracy   
           0.9979              0.9700              0.9676   " )

## ----eval=FALSE---------------------------------------------------------------
# # Dynamic selection with t-test for single genes
# fit.dynamic <- nestcv.glmnet(
#   y = as.numeric(synthetic_metadata$response),
#   x = t(synthetic_rnaseqData),
#   modifyX = "multiDEGGs_combined_filter",
#   modifyX_options = list(
#     filter_method = "ttest",
#     nfilter = 20,
#     dynamic_nfilter = TRUE,
#     keep_single_genes = FALSE
#   ),
#   modifyX_useY = TRUE,
#   n_outer_folds = 5,
#   n_inner_folds = 6,
#   verbose = FALSE
# )

## ----eval=FALSE---------------------------------------------------------------
# # Balanced selection with Wilcoxon-test importance
# fit.balanced <- nestcv.train(
#   y = synthetic_metadata$response,
#   x = t(synthetic_rnaseqData),
#   method = "rf",
#   modifyX = "multiDEGGs_combined_filter",
#   modifyX_options = list(
#     filter_method = "wilcoxon",
#     nfilter = 40,
#     dynamic_nfilter = FALSE
#   ),
#   modifyX_useY = TRUE,
#   n_outer_folds = 5,
#   n_inner_folds = 6,
#   verbose = FALSE
# )

## -----------------------------------------------------------------------------
citation("multiDEGGs")

