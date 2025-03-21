validDeeDeeExperiment <- function(object) {
  msg <- NULL

  if (!is(dea_info(object), "list")) {
    msg <- c(msg, "'dea' must be a list")
  }

  # if (length(dea_info(object)) == 0) {
  # msg <- c(msg, "'dea' must be a non-empty list")
  # }

  if (length(dea_info(object)) > 0) {
    if (any(is.null(names(dea_info(object))))) {
      msg <- c(msg, "'dea' must be a named list")
    }

    dea_names <- names(dea_info(object))

    required_rowdata <- unlist(
      lapply(
        dea_names, function(arg)
          c(
            paste0(arg, "_log2FoldChange"),
            paste0(arg, "_pvalue"),
            paste0(arg, "_padj")
          )
      )
    )
    if (!all(required_rowdata %in% colnames(rowData(object)))) {
      msg <- c(msg, "some required columns were not found in the rowData")
    }
  }

  if (is.null(msg)) {
    TRUE
  } else msg


}

#' @importFrom S4Vectors setValidity2
S4Vectors::setValidity2("DeeDeeExperiment", validDeeDeeExperiment)

