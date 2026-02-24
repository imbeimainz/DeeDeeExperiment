validDeeDeeExperiment <- function(object) {
  msg <- NULL

  ## check dea validity
  if (!is(getDEAInfo(object), "SimpleList")) {
    msg <- c(msg, "`dea` slot must be a SimpleList")
  }

  if (length(getDEAInfo(object)) > 0) {
    if (any(is.null(names(getDEAInfo(object))))) {
      msg <- c(msg, "`dea` must be a named list")
    }

    dea_names <- names(getDEAInfo(object))


    required_rowdata <- unlist(
      lapply(
        dea_names, function(arg) {
          c(
            paste0(arg, "_log2FoldChange"),
            paste0(arg, "_pvalue"),
            paste0(arg, "_padj")
          )
        }
      )
    )
    # this will catch almost every wrong thing
    if (!all(required_rowdata %in% colnames(rowData(object)))) {
      msg <- c(msg, "Some required columns were not found in the rowData")
    }
  }

  ## check fea validity
  if (!is(getFEAInfo(object), "SimpleList")) {
    msg <- c(msg, "`fea` must be a SimpleList")
  }

  if (length(getFEAInfo(object)) > 0) {
    if (any(is.null(names(getFEAInfo(object))))) {
      msg <- c(msg, "`fea` must be a named list")
    }
  }

  if (length(getFEAInfo(object)) > 0) {
    for (i in seq_len(length(getFEAInfo(object)))) {
      entry <- getFEAInfo(object)[[i]]
      if (!is(entry$original_object, "data.frame") &&
          !is(entry$original_object, "enrichResult") &&
          !is(entry$original_object, "gseaResult")) {
        msg <- c(msg, "FEA results should be either a data.frame,
                enrichResult, or a gseaResult object")
      }

      if (!(entry$fe_tool %in% c(
        "topGO", "clusterProfiler", "GeneTonic", "DAVID",
        "fgsea", "gsea", "enrichr", "gProfiler",
        "Not Specified"
      ))) {
        msg <- c(
          msg,
          "Some FEA entries have invalid or unrecognized `fea_tool` values")
      }
    }
  }


  if (is.null(msg)) {
    TRUE
  } else {
    msg
  }
}

#' @importFrom S4Vectors setValidity2
S4Vectors::setValidity2("DeeDeeExperiment", validDeeDeeExperiment)
