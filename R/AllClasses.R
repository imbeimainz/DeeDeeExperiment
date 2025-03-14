#' @rdname DeeDeeExperiment
#'
#' @exportClass DeeDeeExperiment
#'
#' @slot dea This slot is designed to hold the DE-related information. This is
#' internally being created upon importing from the list of DE results objects,
#' provided when instantiating the [DeeDeeExperiment].
#' @slot fea This slot is designed to hold Functional Enrichment related information.
#'
setClass("DeeDeeExperiment",
         contains = "RangedSummarizedExperiment",
         slots = representation(
           dea = "list",
           fea = "list"
         )
)
