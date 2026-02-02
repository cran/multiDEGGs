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
citation("multiDEGGs")

