#' DeeDeeExperiment
#'
#' `DeeDeeExperiment` is a an S4 class that allows integrating and managing transcriptomic
#' analysis results.
#'
#' `DeeDeeExpeirment` is an S4 class extending the `SummarizedExperiment` framework to
#' facilitate the integration and management of transcriptomic analysis results.
#' It introduces two dedicated slots to store Differential Expression (DE) analysis
#' results and Functional Enrichment analysis outcomes, providing a structured approach
#' for downstream analysis.
#'

#' @importFrom DESeq2 results DESeq vst
#' @importFrom edgeR topTags
#' @importFrom limma topTable
#' @importFrom stats na.omit
#' @importClassesFrom SummarizedExperiment RangedSummarizedExperiment
#' @importFrom S4Vectors metadata metadata<- DataFrame SimpleList
#' @importFrom SummarizedExperiment rowData rowData<- mcols assays
#' rowData rowData<- SummarizedExperiment
#' @importFrom utils packageVersion
#' @importFrom methods show as callNextMethod is new validObject
#'
#' @name DeeDeeExperiment-pkg
#' @docType package
"_PACKAGE"
NULL
