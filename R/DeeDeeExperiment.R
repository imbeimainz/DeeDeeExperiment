#' @name DeeDeeExperiment
#'
#' @title The DeeDeeExperiment class
#'
#' @aliases
#' DeeDeeExperiment
#' DeeDeeExperiment-class
#'
#' @description
#' The `DeeDeeExperiment` class is integrate and manage transcriptomic analysis results.
#' It inherits from the SummarizedExperiment class, and additionally stores
#' DE-related/functional enrichment information via dedicated slots and `colData`.
#'
#' @param se A `SummarizedExperiment` object, that will be used as a scaffold to
#' store the DE related information.
#' @param de_results A named list of DE results, in any of the formats supported by
#' the package (currently: results from DESeq2, edgeR, limma).
#' @param enrich_results A named list of functional enrichment results, in any of the
#' formats supported by the package (currently: results from topGO, clusterprofiler,
#' gsea, fgsea algorithms, or a data.frame generated with one of `GeneTonic` shakers)
#'
#' @details
#' The `se` parameter can be optionally left unspecified. If this is the case,
#' the resulting `DeeDeeExperiment` object will contain as features the ones
#' specified by the provided components of the object supplied via the
#' `de_results` parameter.
#'
#' The conversion of the components of the `de_results` list will be handled via
#' conversion functions to uniform the names and set of information which will
#' be stored in the returned `DeeDeeExperiment` object.
#' The names of the list will be used to define the `contrasts` for the different
#' DE analyses included, which will determine the way to access the information
#' stored in the `dea` slot of the `DeeDeeExperiment` object.
#'
#' The content of the `enrich_results` provided by the user will be validated to
#' ensure that it is properly formatted and correctly named. The FE tool can be automatically
#' detected, and based on that, the appropriate shaking method is used to return a standardized
#' format of the FEA results.
#' The names of the list will be used to attempt to associate each enrichment result
#' with a corresponding DE contrast stored in the `DeeDeeExperiment` object.
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
#' @author Najla Abassi, Lea Rothörl, and Federico Marini
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
#'   gene_id = rownames(de_named_list$ifng_vs_naive))
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
DeeDeeExperiment <- function(se = NULL,
                             de_results = NULL,
                             enrich_results = NULL) {
  # old <- S4Vectors:::disableValidity()
  # if (!isTRUE(old)) {
  #   S4Vectors:::disableValidity(TRUE)
  #   on.exit(S4Vectors:::disableValidity(old))
  # }


  # set up functional enrichment results list
  extracted_enrich_results <- list()

  if (!is.null(de_results)) {
    # capture variable name as a character
    entry_name <- deparse(substitute(de_results))
    de_results <- .check_de_results(de_results, entry_name)
  }

  if (!is.null(se)) {
    if (!is(se, "RangedSummarizedExperiment")) {
      # check if it is SE and convert it into a RangedSE
      if (is(se, "SummarizedExperiment")) { ### think again, should we use rather the SCE class?
        se <- as(se, "RangedSummarizedExperiment")
      } else {
        stop("'se' must be a RangedSummarizedExperiment object")
      }
    }

  }
  else {
    # if nothing is passed, return error
    if (length(de_results) == 0 & length(enrich_results) == 0) {
      stop("You have to provide at least an se object or a de_results object!")
    }
    # if no se passed but de_results is not empty, create a mock from it
    message("creating a mock SE from the rows of the DE result objects")
    # mock up the se from the de_results
    #first_de <- de_results[[1]]

    # independently of the class, the feature names are in the
    # rownames slot, TODO: check
    #ids <- rownames(first_de)

    ## check
    #stopifnot(!any(sapply(de_results, function(x) is.null(rownames(x)))))

    if(any(sapply(de_results, function(x) is.null(rownames(x))))) {
      stop("Some elements in the de_results list do not have rownames!")
    }

    ## taking rather the union of all de_res elements
    ids <- unique(unlist(lapply(de_results, rownames)))

    rd_mock <- DataFrame(gene_id = ids, row.names = ids)

    # way1
    se_mock <- SummarizedExperiment(assays = SimpleList(), rowData = rd_mock)
    # se_mock@NAMES <- NULL
    # rownames(se_mock) <- ids

    # no clue why this is strictly needed, but still it seems it is, if mocking up
    se <- as(se_mock, "RangedSummarizedExperiment")

  }


  # TODO: if no SE is really provided, instantiate some rownames, at least directly
  # from the rownames of the result objects
  # TODO: the row names are taken from the FIRST object in the de results then - or
  # from the union of all of them?


  if (is.null(de_results) &&
      is.null(enrich_results)) {
    object <- new("DeeDeeExperiment",
                  se,
                  dea = list(),
                  fea = list())

    # stash the package version
    metadata(object)[["version"]] <- packageVersion("DeeDeeExperiment")

    return(object)
  }


  # TODO: does not have to relate to an SE which has all the slots and all
  # ...


  # TODO: additional checks
  se_out <- se

  # here is where I add the names in the rowData to make all info matched
  # checks on the names
  #names(de_results)
  # if not there, "force add"
  # TODO

  dde_ids <- rownames(se_out)



  dea_contrasts <- list()

  for (i in names(de_results)) {
    this_de <- de_results[[i]]

    # do different things according to what these objects are
    if (is(this_de, "DESeqResults")) {
      input_deseq2 <- .importDE_DESeq2(se_out, this_de, i)
      se_out <- input_deseq2$se
      dea_contrasts[[i]] <- input_deseq2$dea_contrast

      # check for rowname mismatches
      rownames_x <- rownames(rowData(se_out))
      rownames_y <- rownames(this_de)
      mismatched_rows <- sum(!rownames_x %in% rownames_y)

      affected_deas <- character()
      if (mismatched_rows > 0) {
        affected_deas <- c(affected_deas, i)
      }

      mismatch_percent <- (mismatched_rows / length(rownames_x)) * 100

      if (mismatch_percent > 50) {
        warning(
          "A Total number of ", mismatched_rows," mistached rows detected between `rownames(rowData(se))` and rownames for the following dea element: ",
          i,
          ". Consider synchronizing your rownames in both se and de_results elements."
        )
      }


    } else if (is(this_de, "DGEExact") | is(this_de, "DGELRT")) {
      input_edgeR <- .importDE_edgeR(se_out, this_de, i)
      se_out <- input_edgeR$se
      dea_contrasts[[i]] <- input_edgeR$dea_contrast

      # check for rowname mismatches
      rownames_x <- rownames(rowData(se_out))
      rownames_y <- rownames(this_de)
      mismatched_rows <- sum(!rownames_x %in% rownames_y)

      affected_deas <- character()
      if (mismatched_rows > 0) {
        affected_deas <- c(affected_deas, i)
      }

      mismatch_percent <- (mismatched_rows / length(rownames_x)) * 100

      if (mismatch_percent > 50) {
        warning(
          "A Total number of ", mismatched_rows," mistached rows detected between `rownames(rowData(se))` and rownames for the following dea element: ",
          i,
          ". Consider synchronizing your rownames in both se and de_results elements."
        )
      }


    } else if (is(this_de, "MArrayLM")) {
      input_limma <- .importDE_limma(se_out, this_de, i)
      se_out <- input_limma$se
      dea_contrasts[[i]] <- input_limma$dea_contrast

      # check for rowname mismatches
      rownames_x <- rownames(rowData(se_out))
      rownames_y <- rownames(this_de)
      mismatched_rows <- sum(!rownames_x %in% rownames_y)

      affected_deas <- character()
      if (mismatched_rows > 0) {
        affected_deas <- c(affected_deas, i)
      }

      mismatch_percent <- (mismatched_rows / length(rownames_x)) * 100

      if (mismatch_percent > 50) {
        warning(
          "A Total number of ", mismatched_rows," mistached rows detected between `rownames(rowData(se))` and rownames for the following dea element: ",
          i,
          ". Consider synchronizing your rownames in both se and de_results elements."
        )
      }


    }
    # else if (is(this_de, "data.frame")) {
    #   input_custom <- .importDE_custom(se_out, this_de, i)
    #   se_out <- input_custom$se
    #   dea_contrasts[[i]] <- input_custom$dea_contrast
    #
    #   # check for rowname mismatches
    #   rownames_x <- rownames(rowData(se_out))
    #   rownames_y <- rownames(this_de)
    #   mismatched_rows <- sum(!rownames_x %in% rownames_y)
    #
    #   affected_deas <- character()
    #   if (mismatched_rows > 0) {
    #     affected_deas <- c(affected_deas, i)
    #   }
    #
    #   mismatch_percent <- (mismatched_rows / length(rownames_x)) * 100
    #
    #   if (mismatch_percent > 50) {
    #     warning(
    #       "A Total number of ", mismatched_rows," mistached rows detected between `rownames(rowData(se))` and rownames for the following dea element: ",
    #       i,
    #       ". Consider synchronizing your rownames in both se and de_results elements."
    #     )
    #   }
    # }
  }

  ## handle fea results

  fea_contrasts <- list()
  if (!is.null(enrich_results)) {
    # first check content
    enrich_name <- deparse(substitute(enrich_results)) # capture variable name as a char
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
            message("FEA '",
                    fe,
                    "' matched to DE contrast '", # in case of formatted name
                    matched_name,
                    "'")
          } else{
            message("FEA '",
                    fe,
                    "' matched **directly** to DE contrast '", # in case of the same name
                    matched_name, "'")
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
            "Consider naming your enrich_results starting with one of the following prefixes:",
            " 'topGO_', 'ClusterPro_','GeneTonic_', 'DAVID_','gsea_', 'fgsea_', 'enrichr_', 'gPro_',",
            "followed by the contrast name"
          )
        }
      } else {
        de_res_name <- NA_character_
        warning("Could not match FEA '",
                fe,
                "' to a DE contrast because no DE results were provided.\n")
      }

      fe_name <- fe # here goes the fea name

      # detect fea type (aka the packge used to generate the results)
      fe_tool <- .detect_fea_tool(res_enrich)


      # based on the fe_tool , we'll handle the shaking separately, so we always return res_enrich_shaken
      # unless that's an unknown tool used?
      # we can take any already shaken version

      res_enrich_shaken <- NULL # default

      if (fe_tool == "topGO") {
        # to be able to generate gtl objects we shouldn't convert enrich res into data.frame!!
        res_enrich_shaken <- GeneTonic::shake_topGOtableResult(res_enrich)

      } else if (fe_tool == "clusterPro") {
        if (is(res_enrich, "enrichResult")) {
          res_enrich_shaken <- GeneTonic::shake_enrichResult(res_enrich)
        }

      } else if (fe_tool == "GeneTonic") {
        res_enrich_shaken <- res_enrich

      }
      else if (fe_tool == "DAVID") {
        # we are not taking the output of the file!!  so we cannot
        # use genetonic shakers!!
        # create shakers for that
        res_enrich_shaken <- DeeDeeExperiment:::DeeDeefy_david(res_enrich)
      }

      else if (fe_tool == "fgsea") {
        res_enrich_shaken <- GeneTonic::shake_fgseaResult(res_enrich)
      }

      else if (fe_tool == "gsea") {
        if (is(res_enrich, "gseaResult")) {
          res_enrich_shaken <- GeneTonic::shake_gsenrichResult(res_enrich)
        }
      }

      else if (fe_tool == "enrichr") {
          res_enrich_shaken <- DeeDeeExperiment:::DeeDeefy_enrichr(res_enrich)
      }

      else if (fe_tool == "gProfiler") {
        res_enrich_shaken <- DeeDeeExperiment:::DeeDeefy_gprofiler(res_enrich)
      }


      if (is.null(res_enrich_shaken)) {
        message(
          "No shaking method available for this functional enrichment results.",
          " Returning only the original object."
        )
      }


      fea_contrast <- list(
        de_name = de_res_name,# links to de result
        fe_name = fe_name,
        shaken_results = res_enrich_shaken , # return shaken results for later use in genetonic
        original_object = res_enrich,
        fe_tool = fe_tool
      )


      fea_contrasts[[fe]] <- fea_contrast
    }

  }

  object <- new("DeeDeeExperiment",
                se_out,
                dea = dea_contrasts,
                fea = fea_contrasts)

  # stash the package version
  metadata(object)[["version"]] <- packageVersion("DeeDeeExperiment")

  return(object)


}
