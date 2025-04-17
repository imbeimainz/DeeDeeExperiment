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
#' @param enrich_results A named list of functional enrichment results
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
#' ### TODO: add description for fea
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
#' @author Lea Rothörl and Federico Marini
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
      if (is(se, "SummarizedExperiment")) {
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


    } else if (is(this_de, "data.frame")) {
      input_custom <- .importDE_custom(se_out, this_de, i)
      se_out <- input_custom$se
      dea_contrasts[[i]] <- input_custom$dea_contrast

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
  }

  ## handle fea results

  #### TODO: the list should be named

  fea_contrasts <- list()
  if (!is.null(enrich_results)) {
    # first check content
    enrich_name <- deparse(substitute(enrich_results)) # capture variable name as a char
    enrich_results <- .check_enrich_results(enrich_results, enrich_name)

    # get de_name? try and link fea to dea by name
    for (fe in names(enrich_results)) {
      res_enrich <- enrich_results[[fe]]

      if (!is.null(de_results) && length(de_results) > 0) {
        matched_name <- .match_fe_to_de(fe, names(de_results))
        if (!is.na(matched_name) && matched_name %in% names(de_results)) {
          de_res_name <- matched_name
          if (fe != matched_name) {
            message("FEA '", fe, "' matched to DE contrast '", matched_name, "'")
          }
          } else {
            de_res_name <- NA_character_
            warning(
              "Could not match FEA '", fe, "' to a DE contrast.\n", "Available DE results: ",
              paste(names(de_results), collapse = ", "), "\n",
              "Keep in mind that your enrich_results names should start with one of the following prefixes:",
              " 'GO_', 'ClusterPro_', 'KEGG_', 'Reactome_'"
            )
          }
    } else {
      de_res_name <- NA_character_
      warning("Could not match FEA '",
              fe,
              "' to a DE contrast because no DE results were provided.\n")
    }

      fe_name <- fe # here goes the fea name

      if ("GO.ID" %in% colnames(res_enrich)) {
        fe_type <- "TopGo"
      } else {
        fe_type <- "NULL" # placeholder for now
      }

      fea_contrast <- list(
        de_name = de_res_name, # links to de result
        fe_name = fe_name,
        original_object = res_enrich,
        GeneTonicList = NULL , # we'll put back a GT ready obj
        fe_type = fe_type
      )

      fea_contrasts[[fe]] <- fea_contrast
  }


    }



  # rowData(dde)[["new_rd"]] <- de_name

  object <- new("DeeDeeExperiment",
                se_out,
                dea = dea_contrasts,
                fea = fea_contrasts)

  # stash the package version
  metadata(object)[["version"]] <- packageVersion("DeeDeeExperiment")

  return(object)


}





# extends the rowData slot of the provided SE and returns also metadata

#' Import from `DESeq2` DE results
#'
#' @param se A `SummarizedExperiment` object
#' @param res_de A set of DE results, provided as `DESeqResults` as in the `DESeq2`
#' framework
#' @param de_name A character value, describing the contrast of interest. Will be
#' used to compose the column names in the `rowData` slot.
#'
#' @return A list, containing the updated `SummarizedExperiment` object, and the
#' standardized information on the DE analysis, as these are to be used in the
#' `DeeDee` framework.
#'
#' @noRd
#'
#' @examples
#' # TODO example
.importDE_DESeq2 <- function(se, res_de, de_name) {
  # checks TODO:
  # correct object format
  stopifnot(is(res_de, "DESeqResults"))
  # contain the right columns
  stopifnot(all(c("log2FoldChange", "pvalue", "padj") %in% colnames(res_de)))

  # contain the feature ids

  # p value different from NA respect the 0-1 interval
  stopifnot(all(na.omit(res_de$pvalue <= 1)) &
              all(na.omit(res_de$pvalue > 0)))

  #matched_ids <- match(rownames(res_de), rownames(se))

  matched_ids <- match(rownames(se), rownames(res_de)) # we align de res with se
  # only valid indices
  valid_matches <- !is.na(matched_ids)


  # Pre-fill rowData with NA
  rowData(se)[[paste0(de_name, "_log2FoldChange")]] <- NA
  rowData(se)[[paste0(de_name, "_pvalue")]]         <- NA
  rowData(se)[[paste0(de_name, "_padj")]]           <- NA


  # assign values only for matched indices, to have on both sides the
  # same length. we keep NA for unmatched genes
  rowData(se)[[paste0(de_name, "_log2FoldChange")]][valid_matches] <- res_de$log2FoldChange[matched_ids[valid_matches]]
  rowData(se)[[paste0(de_name, "_pvalue")]][valid_matches]         <- res_de$pvalue[matched_ids[valid_matches]]
  rowData(se)[[paste0(de_name, "_padj")]][valid_matches]           <- res_de$padj[matched_ids[valid_matches]]


  dea_contrast <- list(
    alpha = metadata(res_de)$alpha,
    lfcThreshold = metadata(res_de)$lfcThreshold,
    metainfo_logFC = mcols(res_de)$description[colnames(res_de) == "log2FoldChange"],
    metainfo_pvalue = mcols(res_de)$description[colnames(res_de) == "pvalue"],
    original_object = res_de,
    # object_name = deparse(substitute(res_de)),
    package = "DESeq2"
  )

  return(list(se = se, dea_contrast = dea_contrast))
}


#' Import from edgeR DE results
#'
#' @param se A SummarizedExperiment object
#' @param res_de A set of DE results, provided by the `edgeR` framework (either a
#' `DGEExact` or a `DGELRT` object).
#' @param de_name A character value, describing the contrast of interest. Will be
#' used to compose the column names in the rowData slot.
#'
#' @return A list, containing the updated SummarizedExperiment object, and the
#' standardized information on the DE analysis, as these are to be used in the
#' DeeDee framework.
#'
#' @noRd
#'
#' @examples
#' # TODO example
.importDE_edgeR <- function(se, res_de, de_name) {
  # checks object
  stopifnot(is(res_de, "DGEExact") || is(res_de, "DGELRT"))

  # extract columns
  res_tbl <- topTags(res_de, n = nrow(res_de), sort.by = "none")

  # p value different from NA respect the 0-1 interval
  stopifnot(all(na.omit(res_tbl$PValue <= 1)) &
              all(na.omit(res_tbl$PValue > 0)))

  # identify the logFC cols
  logFC_cols <- grep("^logFC", colnames(res_tbl), value = TRUE)


  matched_ids <- match(rownames(se), rownames(res_tbl)) # we align de res with se
  # only valid indices
  valid_matches <- !is.na(matched_ids)

  # pre-fill rowData with NA the assign the corresponding values only for matched
  # indices for logFC, accounting for the fact that the logFC column name in edgeR
  # depends on whether we have 1 or multiple contrasts
  for (i in logFC_cols) {
    rowData(se)[[paste0(de_name, "_log2FoldChange")]] <- NA
    # assign correspionding values
    rowData(se)[[paste0(de_name, "_log2FoldChange")]][valid_matches] <- res_tbl$table[[i]][matched_ids[valid_matches]]
  }

  # pre-fill rowData with NA the assign the corresponding values for matched indices for pval and padj
  rowData(se)[[paste0(de_name, "_pvalue")]]         <- NA
  rowData(se)[[paste0(de_name, "_padj")]]           <- NA


  # assign values only for matched indices, to have on both sides the
  # same length. we keep NA for unmatched genes
  rowData(se)[[paste0(de_name, "_pvalue")]][valid_matches]         <- res_tbl$table$PValue[matched_ids[valid_matches]]
  rowData(se)[[paste0(de_name, "_padj")]][valid_matches]           <- res_tbl$table$FDR[matched_ids[valid_matches]]

  dea_contrast <- list(
    alpha = NA,
    lfcThreshold = NA,
    metainfo_logFC = res_tbl$comparison,
    metainfo_pvalue = paste0("p-value adjusted using ", res_tbl$adjust.method),
    original_object = res_de,
    # object_name = deparse(substitute(res_tbl)),
    package = "edgeR"
  )

  return(list(se = se, dea_contrast = dea_contrast))
}



#' Import from `limma` DE results
#'
#' @param se A `SummarizedExperiment` object
#' @param res_de A set of DE results, provided in the `limma` framework (a `MArrayLM`
#' object).
#' @param de_name A character value, describing the contrast of interest. Will be
#' used to compose the column names in the `rowData` slot.
#'
#' @return A list, containing the updated `SummarizedExperiment` object, and the
#' standardized information on the DE analysis, as these are to be used in the
#' `DeeDee` framework.
#'
#' @noRd
#'
#' @examples
#' # TODO example
#' # ... limma_de <- lmFit
#' # will provide the outout of lmFit - see its examples
.importDE_limma <- function(se, res_de, de_name) {
  # checks object
  stopifnot(is(res_de, "MArrayLM"))

  # make sure there are at least 2 coefficients
  if (ncol(res_de$coefficients) < 2) {
    # we still need to manage the handling of 1 contrast
    warning(
      "The provided MArrayLM object has only ",
      ncol(res_de$coefficients),
      " coefficient(s). At least 2 are required."
    )
  }

  # extract columns
  res_tbl <- topTable(
    res_de,
    coef = 2,
    # this is forced internally, maybe offer more flexibility??
    number = nrow(res_de),
    sort.by = "none"
  )

  # p value different from NA respect the 0-1 interval
  stopifnot(all(na.omit(res_tbl$P.Value <= 1)) &
              all(na.omit(res_tbl$P.Value > 0)))

  #matched_ids <- match(rownames(res_tbl), rownames(se))

  matched_ids <- match(rownames(se), rownames(res_tbl)) # we align de res with se
  # only valid indices
  valid_matches <- !is.na(matched_ids)


  # Pre-fill rowData with NA
  rowData(se)[[paste0(de_name, "_log2FoldChange")]] <- NA
  rowData(se)[[paste0(de_name, "_pvalue")]]         <- NA
  rowData(se)[[paste0(de_name, "_padj")]]           <- NA


  # assign values only for matched indices, to have on both sides the
  # same length. we keep NA for unmatched genes
  rowData(se)[[paste0(de_name, "_log2FoldChange")]][valid_matches] <- res_tbl$logFC[matched_ids[valid_matches]]
  rowData(se)[[paste0(de_name, "_pvalue")]][valid_matches]         <- res_tbl$P.Value[matched_ids[valid_matches]]
  rowData(se)[[paste0(de_name, "_padj")]][valid_matches]           <- res_tbl$adj.P.Val[matched_ids[valid_matches]]

  dea_contrast <- list(
    alpha = NA,
    lfcThreshold = NA,
    metainfo_logFC = NA,
    metainfo_pvalue = NA,
    original_object = res_de,
    # object_name = deparse(substitute(res_tbl)),
    package = "limma"
  )

  return(list(se = se, dea_contrast = dea_contrast))


  # returns info (in the standardized manner)

}

# custom format can be a dataframe, can it be a list???

.importDE_custom <- function(se, res_de, de_name) {
  # checks object
  stopifnot(is(res_de, "data.frame"))

  # expected columns, mainly from DESeq, edgeR and limma?
  expected_columns <- list(
    logFC = c("log2FoldChange", "logFC"),
    pvalue = c("pvalue", "PValue", "P.Value"),
    padj = c("padj", "FDR", "adj.P.Val")
  )

  # find the matching column names in res_de
  matched_cols <- sapply(expected_columns, function(x) {
    match <- intersect(x, colnames(res_de))
    if (length(match) > 0) return(match[1])
    stop("Dataframe does not contain required columns: ", paste(x, collapse = ", "))
  })

  valid_matches <- rownames(se) %in% rownames(res_de)

  matched_ids <- match(rownames(se)[valid_matches], rownames(res_de))


  # Pre-fill rowData with NA
  rowData(se)[[paste0(de_name, "_log2FoldChange")]] <- NA
  rowData(se)[[paste0(de_name, "_pvalue")]]         <- NA
  rowData(se)[[paste0(de_name, "_padj")]]           <- NA

  # assign only matched indices. keep NA for unmatched genes
  rowData(se)[[paste0(de_name, "_log2FoldChange")]][valid_matches] <- res_de[[matched_cols["log2FoldChange"]]][matched_ids]
  rowData(se)[[paste0(de_name, "_pvalue")]][valid_matches]         <- res_de[[matched_cols["pvalue"]]][matched_ids]
  rowData(se)[[paste0(de_name, "_padj")]][valid_matches]           <- res_de[[matched_cols["padj"]]][matched_ids]


  dea_contrast <- list(
    alpha = NA,
    lfcThreshold = NA,
    metainfo_logFC = NA,
    metainfo_pvalue = NA,
    original_object = res_de,
    # object_name = deparse(substitute(res_tbl)),
    package = "custom input"
  )

  return(list(se = se, dea_contrast = dea_contrast))

}



.check_de_results <- function(x, entry_name = NULL) {
  ## checks the DE  input and processes it if it's 1 element
  ## if one single element is given, i.e not a list, it converts it into a list
  ## of length 1 and ensure it has a name
  if (is(x, "DGEExact") ||
      is(x, "DGELRT") || is(x, "MArrayLM") ||
      is(x, "DESeqResults")) {
    # convert into a named list
    x <- list(x)
    names(x) <- entry_name
  }

  ## if a list
  ok_types <- unlist(lapply(x, function(arg) {
    is(arg, "DESeqResults") || is(arg, "DGEExact") ||
      is(arg, "DGELRT") || is(arg, "MArrayLM") || is(arg, "data.frame")
  }))

  if (!all(ok_types)) {
    stop("All elements in the list must be of type DESeqResults, DGEExact, DGELRT, or MArrayLM. Alternatively, it can be a data.frame with at least a logFC, p-value and p-adjusted value columns.")
  }
  if (is.null(names(x)) || any(names(x) == "")) {
    stop("All elements in the provided de_results list must be named!")
  }
  return(x)
}


#' Checking the validity of the imported Enrichment results
#'
#' @param x de_results list
#' @param entry_name fea results name
#'
#' @returns a list of valid results elements
#' @noRd
#'
#' @examples
#' # will turn back a valid fea list
.check_enrich_results <- function(x, entry_name = NULL) {
  # check that:
  # you provided a name for your results
  if (is.null(entry_name)) {
    stop("You must provide a name for your enrichment results!")
  }

  # if results are not either a list or df throw an error
  # TODO later specific object types (like clusterpro results..)
  if (!is.list(x)) {
    # df are also lists :v
    stop("Enrichment results must be a list or a data frame!")
  }

  # check if results are either a df, or a list of dfs
  # if results is one df convert it into a named list
  if (is(x, "data.frame")) {
    x <- list(x)
    names(x) <- entry_name
  }

  # if a list
  ok_types <- unlist(lapply(x, function(arg) {
    is(arg, "data.frame")
  }))

  if (!all(ok_types)) {
    stop("All elements in the list must be of type data.frame!")
  }

  # check the columns for each df
  # for now i ll pretend we only work with res topGO
  required_enrich_cols <- c(
    "GO.ID",
    "Term",
    "Annotated",
    "Significant",
    "Expected",
    "Rank in p.value_classic",
    "p.value_elim",
    "p.value_classic",
    "genes"
  ) # do we have other possible names??? or simply select only
  # specific columns like GO.ID, term sometimes called description...

  for (i in seq_along(x)) {
    df <- x[[i]]
    missing_cols <- required_enrich_cols[!required_enrich_cols %in% colnames(df)]

    if (length(missing_cols) > 0) {
      stop("The following columns are missing: ",
           paste(missing_cols, collapse = ", "))
    }

  }

  return(x)
}




#' Find matching fea and dea results within a DeeDeeExperiment object
#'
#' @param fea_name name of fea to insert
#' @param dea_names names of available deas in DeeDeeExperiment
#' @param pattern acceptable prefixes for fea names, it is supposed to force the
#' user to call their result a specific way so that they can match their dea and
#' fea results
#'
#' @returns either the cleaned named, which is the corresponding dea name, or NA
#' if no match found
#' @noRd
#'
#' @examples
#' dea_names <- c("ctrl_vs_treat", "LPS", "IFNg")
#' .match_fe_to_de("GO_ctrl_vs_treat", dea_names)
.match_fe_to_de <- function(fea_name, dea_names,
                             pattern = "^(GO_|ClusterPro_|KEGG_|Reactome_)") {
  # what else can we put in the pattern??
  cleaned_name <- sub(pattern, "", fea_name, ignore.case = TRUE)
  if (cleaned_name %in% dea_names) {
    return(cleaned_name)
  } else {
    return(NA_character_) # Achtung this needs to be character
  }
}

deedee_import <- function(x) {
  # legacy code:

  # # ----------------------------- argument check ------------------------------
  # choices <- c("DESeq2", "edgeR", "limma")
  # checkmate::assertChoice(input_type, choices)
  #
  # if (input_type == "DESeq2") {
  #   checkmate::assertClass(data, "DESeqResults")
  #   logFC <- data$log2FoldChange
  #   pval <- data$padj
  #   input <- data.frame(logFC, pval)
  #   rownames(input) <- data@rownames
  # } else if (input_type == "edgeR") {
  #   checkmate::assertClass(data, "DGEExact")
  #   logFC <- data[["table"]][["logFC"]]
  #   pval <- data[["table"]][["PValue"]]
  #   input <- data.frame(logFC, pval)
  #   rownames(input) <- data[["genes"]][["genes"]]
  # } else if (input_type == "limma") {
  #   checkmate::assertDataFrame(data, types = "numeric")
  #   logFC <- data$logFC
  #   pval <- data$adj.P.Val
  #   input <- data.frame(logFC, pval)
  #   rownames(input) <- rownames(data)
  # }
  # return(input)
}
