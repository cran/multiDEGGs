## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(multiDEGGs)
library(nestedcv)
data("synthetic_metadata")
data("synthetic_rnaseqData")

# Regularized linear model with interaction pairs only
fit.glmnet <- nestcv.glmnet(
  y = as.numeric(synthetic_metadata$response),
  x = t(synthetic_rnaseqData),
  modifyX = "multiDEGGs_filter",
  modifyX_options = list(
    keep_single_genes = FALSE,
    nfilter = 20
  ),
  modifyX_useY = TRUE,
  n_outer_folds = 5,
  n_inner_folds = 6,
  verbose = FALSE
)

summary(fit.glmnet)

## ----fig.width = 3, fig.height = 3--------------------------------------------
# Random forest model including both pairs and individual genes
fit.rf <- nestcv.train(
  y = synthetic_metadata$response,
  x = t(synthetic_rnaseqData),
  method = "rf",
  modifyX = "multiDEGGs_filter",
  modifyX_options = list(
    keep_single_genes = TRUE,
    nfilter = 30
  ),
  modifyX_useY = TRUE,
  n_outer_folds = 5,
  n_inner_folds = 6,
  verbose = FALSE
)

fit.rf$summary

# Plot ROC on outer folds
plot(fit.rf$roc)

## -----------------------------------------------------------------------------
# Dynamic selection with t-test for single genes
fit.dynamic <- nestcv.glmnet(
  y = as.numeric(synthetic_metadata$response),
  x = t(synthetic_rnaseqData),
  modifyX = "multiDEGGs_combined_filter",
  modifyX_options = list(
    filter_method = "ttest", 
    nfilter = 20,
    dynamic_nfilter = TRUE, 
    keep_single_genes = FALSE
  ),
  modifyX_useY = TRUE,
  n_outer_folds = 5,
  n_inner_folds = 6,
  verbose = FALSE
)

## -----------------------------------------------------------------------------
# Balanced selection with Wilcoxon-test importance
fit.balanced <- nestcv.train(
  y = synthetic_metadata$response,
  x = t(synthetic_rnaseqData),
  method = "rf",
  modifyX = "multiDEGGs_combined_filter",
  modifyX_options = list(
    filter_method = "wilcoxon", 
    nfilter = 40,
    dynamic_nfilter = FALSE
  ),
  modifyX_useY = TRUE,
  n_outer_folds = 5,
  n_inner_folds = 6,
  verbose = FALSE
)

## -----------------------------------------------------------------------------
citation("multiDEGGs")

