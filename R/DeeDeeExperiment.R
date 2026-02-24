#' @name DeeDeeExperiment
#'
#' @title The DeeDeeExperiment class
#'
#' @aliases
#' DeeDeeExperiment
#' DeeDeeExperiment-class
#'
#' @description
#' The `DeeDeeExperiment` class is integrate and manage omics analysis
#' results. It inherits from the `SingleCellExperiment` class, and additionally
#' stores DE-related/functional enrichment information via dedicated slots and
#' `rowData`.
#'
#' @param sce A `SingleCellExperiment` object, that will be used as a scaffold to
#' store the DE related information.
#' @param de_results A named list of DE results, in any of the formats supported
#' by the package (currently: results from `DESeq2`, `edgeR`, `limma`, `muscat`).
#' @param enrich_results A named list of functional enrichment results. Each
#' element can be either a data.frame (currently supports results from `topGO`,
#' `enrichR`, `gProfiler`, `fgsea`, `gsea`, `DAVID`, and output of `GeneTonic`
#' shakers), or an `enrichResult`/`gseaResult` objects (currently supports
#' `clusterProfiler`)
#'
#' @details
#' The `sce` parameter can be optionally left unspecified. If this is the case,
#' the resulting `DeeDeeExperiment` object will contain as features the ones
#' specified by the provided components of the object supplied via the
#' `de_results` parameter.
#'
#' The conversion of the components of the `de_results` list will be handled via
#' conversion functions to uniform the names and set of information which will
#' be stored in the returned `DeeDeeExperiment` object.
#' The names of the list will be used to define the `contrasts` for the
#' different DE analyses included, which will determine the way to access the
#' information stored in the `dea` slot of the `DeeDeeExperiment` object.
#'
#' The content of the `enrich_results` provided by the user will be validated to
#' ensure that it is properly formatted and correctly named. The FE tool can be
#' automatically detected, and based on that, the appropriate shaking method is
#' used to return a standardized format of the FEA results.
#' The names of the list will be used to attempt to associate each enrichment
#' result with a corresponding DE contrast stored in the `DeeDeeExperiment`
#' object, but it also can be defined by the user.
#'
#' Since a `DeeDeeExperiment` is also a `SummarizedExperiment` object, it can be
#' seamlessly provided downstream for visualization and in-depth exploration to
#' packages such as `iSEE` or similar.
#'
#'
#'
#' @return A `DeeDeeExperiment` object.
#' @export
#'
#' @author Najla Abassi, Lea Schwarz, and Federico Marini
#'
#' @examples
#' data("de_named_list", package = "DeeDeeExperiment")
#'
#' dde_onlyde <- DeeDeeExperiment(
#'   de_results = de_named_list
#' )
#'
#' # or, with a SE object as support - even without assay data available
#' library("SummarizedExperiment")
#'
#' rd_macrophage <- DataFrame(
#'   gene_id = rownames(de_named_list$ifng_vs_naive)
#' )
#' rownames(rd_macrophage) <- rownames(de_named_list$ifng_vs_naive)
#' se_macrophage_noassays <- SummarizedExperiment(
#'   assays = SimpleList(),
#'   rowData = rd_macrophage
#' )
#'
#' dde <- DeeDeeExperiment(
#'   se_macrophage_noassays,
#'   de_results = de_named_list
#' )
DeeDeeExperiment <- function(sce = SingleCellExperiment(),
                             de_results = NULL,
                             enrich_results = NULL) {

  if (!is.null(de_results)) {
    # capture variable name as a character
    entry_name <- deparse(substitute(de_results))
    de_results <- .check_de_results(de_results, entry_name)
  }

  if (.is_empty_sce(sce)) {
    # if empty sce and we pass de results, create a mock sce from it
    if (!is.null(de_results)) {
      cli::cli_alert_info(
        "creating a mock SCE from the rows of the DE result objects, if available"
      )

      if (any(vapply(de_results, function(x) is.null(rownames(x)), logical(1)))) {
        stop("Some elements in the de_results list do not have rownames!")
      }

      ## taking the union of all de_res elements
      ids <- unique(unlist(lapply(de_results, rownames)))

      rd_mock <- DataFrame(gene_id = ids, row.names = ids)

      sce <- SingleCellExperiment(assays = SimpleList(), rowData = rd_mock)

    } else if (is.null(de_results) &&
        is.null(enrich_results)) {
      # object <- new("DeeDeeExperiment",
      #               sce,
      #               dea = list(),
      #               fea = list())

      object <- .DeeDeeExperiment(sce,
                                  dea = SimpleList(),
                                  fea = SimpleList())

      # stash the package version
      metadata(object)[["version"]] <- packageVersion("DeeDeeExperiment")

      return(object)
    }
  }

  # if sce not empty
  if (!is(sce, "SingleCellExperiment")) {
    if (is(sce, "SummarizedExperiment")) {
      # check if it is SE and convert it into a RangedSE
      sce <- as(sce, "RangedSummarizedExperiment")

      # we'll build an sce from an se obj
      available_assays <- names(assays(sce))
      assay_list <- setNames(lapply(available_assays, function(x)
        assay(sce, x)),
        available_assays)

      se_to_sce <- SingleCellExperiment(
        assays = assay_list,
        colData = colData(sce),
        rowData = rowData(sce),
        metadata = metadata(sce)
      )
      sce <- se_to_sce
    }
  }

  sce_out <- sce

  dea_contrasts <- SimpleList()

  if (is.list(de_results) && identical(attr(de_results, "package"), "limma")) {
    cli::cli_alert_info(
      "Detected a limma result list for entry, importing accordingly."
    )

    limma_res <- .handle_limma_list(sce_out, de_results, entry_name)
    sce_out <- limma_res$sce
    dea_contrasts <- c(dea_contrasts, limma_res$dea_contrasts)

    # stop here the processing of dea?
    de_results <- list()
  }

  for (i in names(de_results)) {
    this_de <- de_results[[i]]

    # do different things according to what these objects are
    if (is(this_de, "DESeqResults")) {
      input_deseq2 <- .importDE_DESeq2(sce_out, this_de, i)
      sce_out <- input_deseq2$sce
      dea_contrasts[[i]] <- input_deseq2$dea_contrast

      # check for rowname mismatches
      rownames_x <- rownames(rowData(sce_out))
      rownames_y <- rownames(this_de)
      mismatched_rows <- sum(!rownames_x %in% rownames_y)

      affected_deas <- character()
      if (mismatched_rows > 0) {
        affected_deas <- c(affected_deas, i)
      }

      mismatch_percent <- (mismatched_rows / length(rownames_x)) * 100

      if (mismatch_percent > 50) {
        warning(
          "A Total number of ", mismatched_rows,
          " mismatched rows detected between `rownames(rowData(sce))` and  ",
          "rownames for the following dea element: ",
          i,
          " Unmatched genes will have NA values in rowData. ",
          ". Consider synchronizing your rownames in both se and",
          " de_results elements."
        )
      }
    } else if (is(this_de, "DGEExact") | is(this_de, "DGELRT")) {
      input_edgeR <- .importDE_edgeR(sce_out, this_de, i)
      sce_out <- input_edgeR$sce
      dea_contrasts[[i]] <- input_edgeR$dea_contrast

      # check for rowname mismatches
      rownames_x <- rownames(rowData(sce_out))
      rownames_y <- rownames(this_de)
      mismatched_rows <- sum(!rownames_x %in% rownames_y)

      affected_deas <- character()
      if (mismatched_rows > 0) {
        affected_deas <- c(affected_deas, i)
      }

      mismatch_percent <- (mismatched_rows / length(rownames_x)) * 100

      if (mismatch_percent > 50) {
        warning(
          "A Total number of ", mismatched_rows,
          " mismatched rows detected between `rownames(rowData(sce))`",
          " and rownames for the following dea element: ",
          i,
          " Unmatched genes will have NA values in rowData. ",
          " Consider synchronizing your rownames in both se and",
          " de_results elements."
        )
      }
    } else if (is(this_de, "MArrayLM")) {
      input_limma <- .importDE_limma(sce_out, this_de, i)
      sce_out <- input_limma$sce
      dea_contrasts[[i]] <- input_limma$dea_contrast

      # check for rowname mismatches
      rownames_x <- rownames(rowData(sce_out))
      rownames_y <- rownames(this_de)
      mismatched_rows <- sum(!rownames_x %in% rownames_y)

      affected_deas <- character()
      if (mismatched_rows > 0) {
        affected_deas <- c(affected_deas, i)
      }
      mismatch_percent <- (mismatched_rows / length(rownames_x)) * 100

      if (mismatch_percent > 50) {
        warning(
          "A Total number of ", mismatched_rows,
          " mismatched rows detected between `rownames(rowData(sce))`",
          " and rownames for the following dea element: ",
          i,
          " Unmatched genes will have NA values in rowData. ",
          " Consider synchronizing your rownames in both se and",
          " de_results elements."
        )
      }
    } else if (is(this_de, "data.frame")) {
      input_df <- .importDE_df(sce_out, this_de, i)
      sce_out <- input_df$sce
      dea_contrasts[[i]] <- input_df$dea_contrast

      # check for rowname mismatches
      rownames_x <- rownames(rowData(sce_out))
      rownames_y <- rownames(this_de)
      mismatched_rows <- sum(!rownames_x %in% rownames_y)

      affected_deas <- character()
      if (mismatched_rows > 0) {
        affected_deas <- c(affected_deas, i)
      }
      mismatch_percent <- (mismatched_rows / length(rownames_x)) * 100

      if (mismatch_percent > 50) {
        warning(
          "A Total number of ", mismatched_rows,
          " mismatched rows detected between `rownames(rowData(sce))`",
          " and rownames for the following dea element: ",
          i,
          " Unmatched genes will have NA values in rowData. ",
          " Consider synchronizing your rownames in both se and",
          " de_results elements."
        )
      }
    }
  }

  ## handle fea results

  fea_contrasts <- SimpleList()

  if (!is.null(enrich_results)) {
    # first check content
    # capture variable name as a char
    enrich_name <- deparse(substitute(enrich_results))
    enrich_results <- .check_enrich_results(enrich_results, enrich_name)

    # get de_name? try and link fea to dea by name
    for (fe in names(enrich_results)) {
      res_enrich <- enrich_results[[fe]]

      if (!is.null(de_results) && length(de_results) > 0) {
        # hoping here that the user names their results in a meaningful way
        matched_name <- .match_fe_to_de(fe, names(de_results))
        if (!is.na(matched_name) &&
            matched_name %in% names(de_results)) {
          de_res_name <- matched_name
          if (fe != matched_name) {

            cli::cli_alert_info(
              "FEA {.val {fe}} matched to DE contrast {.val {matched_name}}"
            )
          } else {
            cli::cli_alert_info(
              "FEA {.val {fe}} matched directly to DE contrast {.val {matched_name}}"
            )
          }
        } else {
          de_res_name <- NA_character_
          warning(
            "Could not match FEA '",
            fe,
            "' to any DE contrast.\n",
            "Available DE results: ",
            paste(names(de_results), collapse = ", "),
            "\n",
            " Consider naming your enrich_results starting with one of the",
            " following prefixes:",
            " 'topGO_', 'clusterProfiler_','GeneTonic_', 'DAVID_',",
            "'gsea_', 'fgsea_', 'enrichr_', 'gPro_',",
            "followed by the contrast name"
          )
        }
      } else {
        de_res_name <- NA_character_
        warning(
          "Could not match FEA '",
          fe,
          "' to a DE contrast because no DE results were provided.\n"
        )
      }

      fe_name <- fe # here goes the fea name

      # detect fea type (aka the packge used to generate the results)
      fe_tool <- .detect_fea_tool(res_enrich)


      res_enrich_shaken <- NULL # default

      res_enrich_shaken <- .shake_enrich_res(res_enrich, fe_tool)

      if (is.null(res_enrich_shaken)) {
        cli::cli_alert_info(
          "No shaking method available for this functional enrichment results.
                            Returning only the original object."
        )
      }

      fea_contrast <- list(
        de_name = de_res_name,
        fe_name = fe_name,
        shaken_results = res_enrich_shaken,
        original_object = res_enrich,
        fe_tool = fe_tool,
        fe_tool_version = if (
          fe_tool %in% loadedNamespaces()) packageVersion(fe_tool) else NA
      )

      fea_contrasts[[fe]] <- fea_contrast
    }
  }

  # object <- new("DeeDeeExperiment",
  #               sce_out,
  #               dea = dea_contrasts,
  #               fea = fea_contrasts
  # )

  object <- .DeeDeeExperiment(sce_out,
                              dea = dea_contrasts,
                              fea = fea_contrasts)

  # stash the package version
  metadata(object)[["version"]] <- packageVersion("DeeDeeExperiment")

  return(object)
}
