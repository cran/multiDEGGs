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

## ----fig.width = 5, fig.height = 3.5------------------------------------------
plot_regressions(deggs_object,
                 assayDataName = "RNAseq",
                 gene_A = "MTOR", 
                 gene_B = "AKT2",
                 legend_position = "bottomright")

## -----------------------------------------------------------------------------
# Convert metadata into a named factor vector containing only the labels to 
# predict (to ensure compatibility with the nestedCV functions)
metadata_vector <- as.factor(synthetic_metadata[, "response"])

# Make sure the assay data you want to use for prediction is a matrix and
# transpose it, so features are in columns 
# (standard format for machine learning)
assayData <- as.matrix(t(synthetic_rnaseqData))

# NOTE: Make sure your vector is ALIGNED with assayData. 
# The order of the annotations in the metadata_vector must match the 
# samples in the rows of assayData

# Remove zero variance columns from data
assayData <- assayData[,apply(assayData, 2, var, na.rm=TRUE) != 0] 

# define a filtering function that extracts differential nodes and links
# using multiDEGGs:
DEGGs_modxy <- function(metadata, assayData, ...) {
  counts <- t(assayData)
  names(metadata_vector) <- rownames(assayData)
  
  deggs_object <- multiDEGGs::get_diffNetworks(
    assayData = counts,
    metadata = metadata_vector,
    percentile_vector = seq(0.25, 0.98, by = 0.05),
    use_qvalues = TRUE, 
    show_progressBar = FALSE,
    verbose = FALSE,
    cores = 1
  )
  
  pairs <- multiDEGGs::get_sig_deggs(deggs_object, 1, 0.05)
  
  # Take genes 2 by 2 from top pairs to lower ones
  keep_DEGGs <- unique(unlist(lapply(1:nrow(pairs), function(i) {
    row_n = c(pairs[i,1], pairs[i,2])
  })))
  
  # The following could be added if you want to set a maximum number of 
  # predictors to be selected: 
  # if (length(keep_DEGGs) > 50) {    # take only top 50 predictors 
  #   keep_DEGGs <- keep_DEGGs[1:50]
  #   pairs <- pairs[which(pairs$var1 %in% keep_DEGGs & 
  #                        pairs$var2 %in% keep_DEGGs), ]
  # }
  
  out <- list(keep_DEGGs = keep_DEGGs, pairs = pairs)
  class(out) <- "DEGGs_modxy"
  return(out)
}

## -----------------------------------------------------------------------------
# This custom predict function will add new columns to x (can be train or test)
predict.DEGGs_modxy <- function(DEGGs.object, newdata, filter = TRUE, 
                                interaction.type = "ratio",
                                sep = ":", ...) {
  if (length(DEGGs.object$keep) != 0) {
    pairs <- DEGGs.object$pairs
    x2a <- newdata[, pairs[, 1], drop = FALSE]
    x2b <- newdata[, pairs[, 2], drop = FALSE]
    
    if (interaction.type == "ratio") {
      x2 <- x2a/x2b
    } else {
      x2 <- x2a*x2b
    }
    
    colnames(x2) <- paste(colnames(x2a), colnames(x2b), sep = sep)
    
    if (filter) {
      keep <- DEGGs.object$keep[!is.na(DEGGs.object$keep)]
      x1 <- newdata[, keep]
      return(cbind(x1, x2))
    } else {
      return(cbind(newdata, x2))
    }
  }
  return(newdata)
}

## -----------------------------------------------------------------------------
sessionInfo()

## -----------------------------------------------------------------------------
citation("multiDEGGs")

