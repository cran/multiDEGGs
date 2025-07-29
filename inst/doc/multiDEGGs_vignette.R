## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load_data----------------------------------------------------------------
library(multiDEGGs)
data("synthetic_metadata")
data("synthetic_rnaseqData")
data("synthetic_proteomicData")
data("synthetic_OlinkData")

## -----------------------------------------------------------------------------
assayData_list <- list("RNAseq" = synthetic_rnaseqData,
                       "Proteomics" = synthetic_proteomicData,
                       "Olink" = synthetic_OlinkData)

deggs_object <- get_diffNetworks(assayData = assayData_list,
                                 metadata = synthetic_metadata,
                                 category_variable = "response",
                                 regression_method = "lm",
                                 padj_method = "bonferroni",
                                 verbose = FALSE,
                                 show_progressBar = FALSE,
                                 cores = 2)

## ----eval=FALSE---------------------------------------------------------------
# View_diffNetworks(deggs_object)

## ----warning=FALSE------------------------------------------------------------
get_multiOmics_diffNetworks(deggs_object, sig_threshold = 0.05)

## -----------------------------------------------------------------------------
deggs_object_oneOmic <- get_diffNetworks(assayData = synthetic_rnaseqData,
                                 metadata = synthetic_metadata,
                                 category_variable = "response",
                                 regression_method = "lm",
                                 padj_method = "bonferroni",
                                 verbose = FALSE,
                                 show_progressBar = FALSE,
                                 cores = 2)

get_sig_deggs(deggs_object_oneOmic, sig_threshold = 0.05)

## ----fig.width = 4.5, fig.height = 4, eval=FALSE------------------------------
# plot_regressions(deggs_object,
#                  assayDataName = "RNAseq",
#                  gene_A = "MTOR",
#                  gene_B = "AKT2",
#                  legend_position = "bottomright")

## -----------------------------------------------------------------------------
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
sessionInfo()

## -----------------------------------------------------------------------------
citation("multiDEGGs")

