#' Import from `DESeq2` DE results
#'
#' @param sce A `SingleCellExperiment` object
#' @param res_de A set of DE results, provided as `DESeqResults` as in the
#' `DESeq2` framework
#' @param de_name A character value, describing the contrast of interest. Will
#' be used to compose the column names in the `rowData` slot.
#'
#' @return A list, containing the updated `SingleCellExperiment` object, and the
#' standardized information on the DE analysis, as these are to be used in the
#' `DeeDee` framework.
#'
#' @noRd
.importDE_DESeq2 <- function(sce, res_de, de_name) {
  # correct object format
  stopifnot(is(res_de, "DESeqResults"))
  # contain the right columns
  stopifnot(all(c("log2FoldChange", "pvalue", "padj") %in% colnames(res_de)))

  # contain the feature ids

  # p value different from NA respect the 0-1 interval
  stopifnot(all(na.omit(res_de$pvalue <= 1)) &
              all(na.omit(res_de$pvalue >= 0)))

  matched_ids <- match(rownames(sce), rownames(res_de)) # we align de res with
  # se only valid indices
  valid_matches <- !is.na(matched_ids)

  sce <- .fill_rowdata_with_dea(sce = sce,
                                de_name = de_name,
                                de_res = res_de,
                                de_cols = c(logFC = "log2FoldChange",
                                            pval = "pvalue",
                                            padj = "padj"),
                                valid_matches = valid_matches,
                                matched_ids = matched_ids)

  if (is.null(metadata(sce)$singlecontrast)) {
    metadata(sce)$singlecontrast <- list()
  }
  metadata(sce)$singlecontrast[[de_name]] <- res_de

  dea_contrast <- list(
    alpha = metadata(res_de)$alpha,
    lfcThreshold = metadata(res_de)$lfcThreshold,
    metainfo_logFC = mcols(res_de)$description[colnames(res_de) == "log2FoldChange"],
    metainfo_pvalue = mcols(res_de)$description[colnames(res_de) == "pvalue"],
    original_object = list(
      metadata_storage = "singlecontrast",
      key   = de_name, #contrast name
      coef  = NULL),
    package = "DESeq2",
    package_version = packageVersion("DESeq2")
  )

  return(list(sce = sce, dea_contrast = dea_contrast))
}


#' Import from edgeR DE results
#'
#' @param sce A `SingleCellExperiment` object
#' @param res_de A set of DE results, provided by the `edgeR` framework (either
#' a `DGEExact` or a `DGELRT` object).
#' @param de_name A character value, describing the contrast of interest. Will
#' be used to compose the column names in the rowData slot.
#'
#' @return A list, containing the updated `SingleCellExperiment` object, and the
#' standardized information on the DE analysis, as these are to be used in the
#' `DeeDee` framework.
#'
#' @noRd
.importDE_edgeR <- function(sce, res_de, de_name) {
  # checks object
  stopifnot(is(res_de, "DGEExact") || is(res_de, "DGELRT"))

  # extract columns
  res_tbl <- topTags(res_de, n = nrow(res_de), sort.by = "none")

  # p value different from NA respect the 0-1 interval
  stopifnot(all(na.omit(res_tbl$PValue <= 1)) &
              all(na.omit(res_tbl$PValue >= 0)))

  # identify the logFC cols
  logFC_cols <- grep("^logFC", colnames(res_tbl), value = TRUE)


  matched_ids <- match(rownames(sce), rownames(res_tbl)) # we align de res with
  # se only valid indices
  valid_matches <- !is.na(matched_ids)

  # pre-fill rowData with NA the assign the corresponding values only for
  # matched indices for logFC, accounting for the fact that the logFC column
  # name in edgeR depends on whether we have 1 or multiple contrasts
  for (i in logFC_cols) {
    sce <- .fill_rowdata_with_dea(sce = sce,
                                  de_name = de_name,
                                  de_res = res_tbl,
                                  de_cols = c(logFC = i,
                                              pval = "PValue",
                                              padj = "FDR"),
                                  valid_matches = valid_matches,
                                  matched_ids = matched_ids)
  }

  if (is.null(metadata(sce)$singlecontrast)) {
    metadata(sce)$singlecontrast <- list()
  }
  metadata(sce)$singlecontrast[[de_name]] <- res_de

  dea_contrast <- list(
    alpha = NA,
    lfcThreshold = NA,
    metainfo_logFC = res_tbl$comparison,
    metainfo_pvalue = paste0("p-value adjusted using ", res_tbl$adjust.method),
    original_object = list(
      metadata_storage = "singlecontrast",
      key   = de_name, #contrast name
      coef  = NULL),
    package = "edgeR",
    package_version = packageVersion("edgeR")
  )

  return(list(sce = sce, dea_contrast = dea_contrast))
}



#' Import from `limma` DE results
#'
#' @param sce A `SingleCellExperiment` object
#' @param res_de A set of DE results, provided in the `limma` framework
#' (a `MArrayLM` object).
#' @param de_name A character value, describing the contrast of interest. Will
#' be used to compose the column names in the `rowData` slot.
#'
#' @return A list, containing the updated `SingleCellExperiment` object, and the
#' standardized information on the DE analysis, as these are to be used in the
#' `DeeDee` framework.
#'
#' @noRd
.importDE_limma <- function(sce, res_de, de_name) {
  # checks object
  stopifnot(is(res_de, "MArrayLM"))

  # make sure there are at least 2 coefficients
  if (ncol(res_de$coefficients) < 2) {
    # we still need to manage the handling of 1 contrast
    stop(
      "The provided MArrayLM object has only ",
      ncol(res_de$coefficients),
      " coefficient(s). At least 2 (intercept + 1 contrast) are required."
    )
  } else if (ncol(res_de$coefficients) > 2) {
    warning(
      "The provided MArrayLM object has ", ncol(res_de$coefficients) ,
      " coefficients. ",
      "Only coefficient 2 is being used for contrast '", de_name, "'.\n",
      "For proper handling of multiple contrasts, please use ",
      "`limma_list_for_dde()` and pass its output to a DeeDeeExperiment object."
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
              all(na.omit(res_tbl$P.Value >= 0)))

  matched_ids <- match(rownames(sce), rownames(res_tbl)) # we align de res with
  # se only valid indices
  valid_matches <- !is.na(matched_ids)

  sce <- .fill_rowdata_with_dea(sce = sce,
                                de_name = de_name,
                                de_res = res_tbl,
                                de_cols = c(logFC = "logFC",
                                            pval = "P.Value",
                                            padj = "adj.P.Val"),
                                valid_matches = valid_matches,
                                matched_ids = matched_ids)

  if (is.null(metadata(sce)$singlecontrast)) {
    metadata(sce)$singlecontrast <- list()
  }
  metadata(sce)$singlecontrast[[de_name]] <- res_de

  dea_contrast <- list(
    alpha = NA,
    lfcThreshold = NA,
    metainfo_logFC = NA,
    metainfo_pvalue = NA,
    original_object = list(
      metadata_storage = "singlecontrast",
      key   = de_name, #contrast name
      coef  = de_name),
    package = "limma",
    package_version = packageVersion("limma")
  )

  return(list(sce = sce, dea_contrast = dea_contrast))
}



#' Checking the validity of the imported DE results.
#' @param x de_results list
#' @param entry_name dea results name
#'
#' @returns a list of valid results elements
#' @noRd
.check_de_results <- function(x, entry_name = NULL) {
  ## checks the DE  input and processes it if it's 1 element
  ## if one single element is given, i.e not a list, it converts it into a list
  ## of length 1 and ensure it has a name
  if (is(x, "DGEExact") ||
      is(x, "DGELRT") || is(x, "MArrayLM") ||
      is(x, "DESeqResults") || is(x, "data.frame")) {
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
    stop("All elements in the list must be of type DESeqResults,",
         " DGEExact, DGELRT, or MArrayLM. Alternatively, it can be a data.frame",
         " with at least a 'log2FoldChange', 'pvalue' and 'padj' columns.")
  }
  if (is.null(names(x)) || any(names(x) == "")) {
    stop("All elements in the provided de_results list must be named!")
  }
  return(x)
}


#' Checking the validity of the imported Enrichment results.
#' This function will return a valid named fea list
#'
#' @param x fe_results list
#' @param entry_name fea results name
#'
#' @returns a list of valid results elements
#' @noRd
.check_enrich_results <- function(x, entry_name = NULL) {
  # check that:
  # you provided a name for your results

  # if (is.null(entry_name)) {
  #   stop("You must provide a name for your enrichment results!")
  # }

  # if results are not either a list or df or enrichResult or gseaResult obj
  # throw an error
  if (!(is(x, "data.frame") || is(x, "enrichResult") || is.list(x) ||
        is(x, "gseaResult"))) {
    stop(
      "Enrichment results must be a data frame, ",
      "an enrichResult object, a gseaResult object or a list of these elements!"
    )
  }

  # if results is not a list  (one df or enrichResult obj) put it into a
  # named list
  if (is(x, "data.frame") || is(x, "enrichResult") || is(x, "gseaResult")) {
    x <- list(x)
    names(x) <- entry_name
  }

  # check if the elements of the list are either data.frame or enrichResult obj
  # gost() returns a large list, so we can accept list
  x <- lapply(x, function(arg) {
    if (is(arg, "enrichResult") || is(arg, "data.frame") ||
        is(arg, "gseaResult")) {
      arg
    } else {
      stop("Elements in the list must be a data.frame or enrichResult",
           " or gseaResult object!")
    }
  })

  # check that all elements in the list have non_empty names
  if (is.null(names(x)) || any(names(x) == "")) {
    stop("All elements in the provided enrich_results list must be named!")
  }


  # check the columns for each df
  required_enrich_cols <- list(
    topGO = c("GO.ID", "Term", "Significant", "p.value_elim", "genes"),
    clusterProfiler = c("ID", "Description", "pvalue", "geneID", "Count"),
    GeneTonic = c("gs_id", "gs_description", "gs_pvalue", "gs_genes",
                  "gs_de_count"),
    DAVID = c(
      "Category", "Term", "Count", "X.", "PValue", "Genes", "List.Total",
      "Pop.Hits", "Pop.Total", "Fold.Enrichment", "Bonferroni", "Benjamini",
      "FDR"),
    fgsea = c("pathway", "pval", "padj", "ES", "NES", "size", "leadingEdge"),
    gsea = c("ID", "Description", "pvalue", "p.adjust", "core_enrichment"),
    enrichr = c(
      "Term", "Overlap", "P.value", "Adjusted.P.value", "Old.P.value",
      "Old.Adjusted.P.value", "Odds.Ratio", "Combined.Score", "Genes"
    ),
    gProfiler = c(
      "source", "term_name", "term_id", "term_size", "query_size",
      "intersection_size", "effective_domain_size"
    ) # taking only the intersection of both r and txt file outputs
  )


  for (i in names(x)) {
    df <- x[[i]]

    if (is(df, "enrichResult") || is(df, "gseaResult")) {
      cols <- colnames(df@result)
    }

    else {
      cols <- colnames(df)
    }

    matches <- vapply(required_enrich_cols, function(required_cols) {
      all(required_cols %in% cols)
    }, logical(1))

    if (!any(matches)) {
      stop(
        c(
          "Element `", i,
          "` does not contain the required columns for any known enrichment type! \n",
          "Please check that you re providing a valid enrichment result. \n",
          "Call `supported_fea_formats()` to see available formats"
        )
      )
    } ### long error msg?
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
.match_fe_to_de <- function(fea_name, dea_names,
                            pattern =
                              "^(topGO_|clusterProfiler_|GeneTonic_|DAVID_|gsea_|fgsea_|enrichr_|gPro_)") {

  cleaned_name <- sub(pattern, "", fea_name, ignore.case = TRUE)
  if (cleaned_name %in% dea_names) {
    return(cleaned_name)
  } else {
    return(NA_character_)
  }
}

#' detect the fe input type (e.g. topGO, clusterProfiler...)
#'
#' @param fe_res FE result table
#' @noRd
.detect_fea_tool <- function(fe_res) {
  stopifnot(is(fe_res, "data.frame") || is(fe_res, "enrichResult") ||
              is(fe_res, "gseaResult"))

  # detect fea type from what columns are found in fea
  required_enrich_cols <- list(
    topGO = c("GO.ID", "Term", "Significant", "p.value_elim", "genes"),
    clusterProfiler = c("ID", "Description", "pvalue", "geneID", "Count"),
    GeneTonic = c("gs_id", "gs_description", "gs_pvalue", "gs_genes"),
    DAVID = c(
      "Category", "Term", "Count", "X.", "PValue", "Genes", "List.Total",
      "Pop.Hits", "Pop.Total", "Fold.Enrichment", "Bonferroni", "Benjamini",
      "FDR"),
    fgsea = c("pathway", "pval", "padj", "ES", "NES", "size", "leadingEdge"),
    gsea = c("ID", "Description", "pvalue", "p.adjust", "core_enrichment"),
    enrichr = c(
      "Term", "Overlap", "P.value", "Adjusted.P.value", "Old.P.value",
      "Old.Adjusted.P.value", "Odds.Ratio", "Combined.Score", "Genes"
    ),
    gProfiler = c(
      "source", "term_name", "term_id", "term_size", "query_size",
      "intersection_size", "effective_domain_size"
    ) # taking only the intersection of both r
    # and txt file outputs
  )

  # extract result table if it is an enrichRes obj
  if (is(fe_res, "enrichResult")) {
    fe_res <- fe_res@result
  }

  # extract result table if it is a gseaResult obj
  if (is(fe_res, "gseaResult")) {
    fe_res <- fe_res@result
  }

  # get col names
  cols <- colnames(fe_res)

  # check for col name matches
  matches <- vapply(names(required_enrich_cols), function(tool) {
    all(required_enrich_cols[[tool]] %in% cols)
  }, logical(1))

  matched_tools <- names(matches)[matches]
  if (length(matched_tools) == 0) {
    return("Not Specified")
  }

  # just in case the user has a table with columns from 2 tools :v ?
  if (length(matched_tools) > 1) {
    warning(
      "Multiple FEA tool formats matched: ",
      paste(matched_tools, collapse = ", "),
      ". Returning all matches."
    )
  }

  return(matched_tools)
}



#' .DeeDeefy_david() , a slightly modified function based on the original shaker
#' for DAVID in GeneTonic. It takes a data.frame instead of the path to the
#' output
#'
#' @param david_output a data.frame of functional enrichment results exported
#' from DAVID
#'
#' @returns a data.frame in GeneTonic shaker standard format
#'
#' @noRd
.DeeDeefy_david <- function(david_output) {
  if (!is(david_output, "data.frame")) {
    stop("DAVID results should be a data.frame!")
  }

  exp_colnames <- c(
    "Category", "Term", "Count", "X.", "PValue", "Genes",
    "List.Total", "Pop.Hits", "Pop.Total", "Fold.Enrichment",
    "Bonferroni", "Benjamini", "FDR"
  )
  if (!all(exp_colnames %in% colnames(david_output))) {
    stop(
      "I could not find some of the usual column names from the DAVID output",
      " exported to file")
  }

  message("Found ", nrow(david_output),
          " gene sets in the file output from DAVID of which ",
          sum(david_output$PValue <= 0.05),
          " are significant (p-value <= 0.05).")
  message("Converting for usage within the DeeDeeExperiment framework...")

  mydf <- data.frame(
    gs_id = unlist(lapply(strsplit(david_output$Term, "~"),
                          function(arg) arg[[1]])),
    gs_description = unlist(lapply(strsplit(david_output$Term, "~"),
                                   function(arg) arg[[2]])),
    gs_pvalue = david_output$PValue,
    gs_genes = gsub(", ", ",", david_output$Genes),
    gs_de_count = david_output$Count,
    gs_bg_count = david_output$Pop.Hits,
    gs_ontology = david_output$Category,
    gs_generatio = david_output$Count / david_output$List.Total,
    gs_bgratio = david_output$Pop.Hits / david_output$Pop.Total,
    gs_foldenrich = david_output$Fold.Enrichment,
    gs_bonferroni = david_output$Bonferroni,
    gs_benjamini = david_output$Benjamini,
    gs_FDR = david_output$FDR,
    stringsAsFactors = FALSE
  )

  rownames(mydf) <- mydf$gs_id

  return(mydf)
}


#' .DeeDeefy_enrichr() , a slightly modified function based on the original
#' shaker for enrichR in GeneTonic. It takes a data.frame instead of the path to
#' the output
#'
#' @param enrichr_output a data.frame with the output of `enrichr`, related to a
#' specific set of genesets. Usually it is one of the members of the list
#' returned by the initial call to `enrichr`.
#'
#' @returns a data.frame in GeneTonic shaker standard format
#'
#' @noRd
.DeeDeefy_enrichr <- function(enrichr_output) {
  exp_colnames <- c(
    "Term", "Overlap", "P.value", "Adjusted.P.value",
    "Old.P.value", "Old.Adjusted.P.value", "Odds.Ratio",
    "Combined.Score", "Genes"
  )
  if (!all(exp_colnames %in% colnames(enrichr_output))) {
    stop(
      "I could not find some of the usual column names from the Enrichr output"
    )
  }

  message("Found ", nrow(enrichr_output),
          " gene sets in the file output from Enrichr of which ",
          sum(enrichr_output$P.value <= 0.05),
          " are significant (p-value <= 0.05).")
  message("Converting for usage within the DeeDeeExperiment framework...")

  mydf <- data.frame(
    gs_id = gsub("\\)", "", gsub("^.* \\(", "", enrichr_output$Term)),
    gs_description = gsub(" \\(GO.*$", "", enrichr_output$Term),
    gs_pvalue = enrichr_output$P.value,
    gs_genes = gsub(";", ",", enrichr_output$Genes),
    gs_de_count = as.numeric(
      unlist(lapply(strsplit(enrichr_output$Overlap, "/"),
                    function(arg) arg[[1]]))
    ),
    gs_bg_count = as.numeric(
      unlist(lapply(strsplit(enrichr_output$Overlap, "/"),
                    function(arg) arg[[2]]))
    ),
    gs_adj_pvalue = enrichr_output$Adjusted.P.value,
    stringsAsFactors = FALSE
  )

  rownames(mydf) <- mydf$gs_id

  return(mydf)
}

#' .DeeDeefy_gprofiler() , a slightly modified function based on the original
#' shaker for g:Profiler in GeneTonic. It takes a data.frame instead of the
#' path to the output
#'
#' @param gprofiler_output_df a data.frame of functional enrichment results
#' exported from g:Profiler
#' @param gprofiler_output a data.frame with the output of `gost()` in
#' `gprofiler2`.
#' Usually it is one of the members of the list returned by the initial call
#' to `gost()`
#'
#' @returns a data.frame in GeneTonic shaker standard format
#'
#' @noRd
.DeeDeefy_gprofiler <- function(gprofiler_output) {
  exp_colnames_textual <- c(
    "source", "term_name", "term_id", "adjusted_p_value",
    "negative_log10_of_adjusted_p_value", "term_size",
    "query_size", "intersection_size", "effective_domain_size",
    "intersections"
  )

  exp_colnames_rcall <- c(
    "query", "significant", "p_value", "term_size", "query_size",
    "intersection_size", "precision", "recall",
    "term_id", "source", "term_name", "effective_domain_size",
    "source_order", "parents", "evidence_codes", "intersection"
  )

  if (all(exp_colnames_textual %in% colnames(gprofiler_output))) {
    message("Found ", nrow(gprofiler_output),
            " gene sets in the file output from g:Profiler of which ",
            sum(gprofiler_output$adjusted_p_value <= 0.05),
            " are significant (p-value <= 0.05).")
    message("Converting for usage within the DeeDeeExperiment framework...")

    mydf <- data.frame(
      gs_id = gprofiler_output$term_id,
      gs_description = gprofiler_output$term_name,
      gs_pvalue = gprofiler_output$adjusted_p_value,
      gs_genes = gprofiler_output$intersections,
      gs_de_count = gprofiler_output$intersection_size,
      gs_bg_count = gprofiler_output$term_size,
      gs_adj_pvalue = gprofiler_output$adjusted_p_value,
      stringsAsFactors = FALSE
    )
  }

  else if (all(colnames(gprofiler_output) %in% exp_colnames_rcall)) {
    # using directly the output from the call from gprofiler2
    # if still a list, might need to select the appropriate element

    message("Found ", nrow(gprofiler_output),
            " gene sets in the file output from g:Profiler of which ",
            sum(gprofiler_output$p_value <= 0.05),
            " are significant (p-value <= 0.05).")
    message("Converting for usage within the DeeDeeExperiment framework...")

    mydf <- data.frame(
      gs_id = gprofiler_output$term_id,
      gs_description = gprofiler_output$term_name,
      gs_pvalue = gprofiler_output$p_value,
      gs_genes = gprofiler_output$intersection,
      gs_de_count = gprofiler_output$intersection_size,
      gs_bg_count = gprofiler_output$term_size,
      gs_adj_pvalue = gprofiler_output$p_value,
      gs_ontology = gprofiler_output$source,
      stringsAsFactors = FALSE
    )
  } else {
    stop(
      "I could not find some of the usual column names from the g:Profiler output.",
      " A possible reason could be that you did not specify `evcodes = TRUE`?",
      " This is required to fill in all the required fields of `res_enrich`"
    )
  }

  rownames(mydf) <- mydf$gs_id

  return(mydf)
}


#' .DeeDeefy_enrichResult() , a function based on the original shaker
#' for enrichResult objects in GeneTonic
#'
#' @param obj An `enrichResult` object, obtained via `clusterProfiler` (or also
#' via `reactomePA`)
#'
#' @returns a data.frame in GeneTonic shaker standard format
#'
#' @noRd
.DeeDeefy_enrichResult <- function(obj) {
  if (!is(obj, "enrichResult")) {
    stop("Provided object must be of class `enrichResult`")
  }

  if (is.null(obj@result$geneID)) {
    stop(
      "You are providing an object where the gene symbols are not specified, ",
      "this is required for running GeneTonic properly."
    )
  }

  message("Found ", nrow(obj@result),
          " gene sets in `enrichResult` object, of which ",
          nrow(as.data.frame(obj)), " are significant.")
  message("Converting for usage within the DeeDeeExperiment framework...")

  fullresults <- obj@result

  mydf <- data.frame(
    gs_id = fullresults$ID,
    gs_description = fullresults$Description,
    gs_pvalue = fullresults$pvalue,
    gs_genes = gsub("/", ",", fullresults$geneID),
    gs_de_count = fullresults$Count,
    gs_bg_count = unlist(lapply(strsplit(fullresults$BgRatio, "/"),
                                function(arg) arg[[1]])),
    gs_ontology = obj@ontology,
    gs_GeneRatio = fullresults$GeneRatio,
    gs_BgRatio = fullresults$BgRatio,
    gs_p.adjust = fullresults$p.adjust,
    gs_qvalue = fullresults$qvalue,
    stringsAsFactors = FALSE
  )

  rownames(mydf) <- mydf$gs_id

  return(mydf)
}



#' .DeeDeefy_topGOtableResult() , a function based on the original shaker
#' for topGOtableResult objects in GeneTonic
#'
#' @param obj An `topGOtableResult` object
#' @param p_value_column Character, specifying which column the p value for
#' enrichment has to be used. Example values are "p.value_elim" or
#' "p.value_classic"
#'
#' @returns a data.frame in GeneTonic shaker standard format
#'
#' @noRd
.DeeDeefy_topGOtableResult <- function(obj,
                                       p_value_column = "p.value_elim") {
  if (!all(
    c("GO.ID", "Term", "Annotated", "Significant", "Expected",
      "p.value_classic")
    %in% colnames(obj))) {
    stop(
      "The provided object must be of in the format specified by the",
      " `pcaExplorer::topGOtable` function or the `mosdef::run_topGO` function")
  }

  if (!p_value_column %in% colnames(obj)) {
    stop(
      "You specified a column for the p-value which is not contained in the",
      " provided object. \n",
      "Please check the colnames of your object in advance."
    )
  }

  if (!"genes" %in% colnames(obj)) {
    stop(
      "The column `genes` is not present in the provided object and",
      " is required for properly running GeneTonic.",
      "\nMaybe you did set `addGeneToTerms` to FALSE in the call to",
      " `pcaExplorer::topGOtable` or to `mosdef::run_topGO`?"
    )
  }

  # Thought: store somewhere the ontology if possible - in an extra column?
  message("Found ", nrow(obj), " gene sets in `topGOtableResult` object.")
  message("Converting for usage within the DeeDeeExperiment framework...")

  fullresults <- obj

  mydf <- data.frame(
    gs_id = fullresults$GO.ID,
    gs_description = fullresults$Term,
    gs_pvalue = fullresults[[p_value_column]],
    gs_genes = fullresults$genes,
    gs_de_count = fullresults$Significant,
    gs_bg_count = fullresults$Annotated,
    # gs_ontology = obj@ontology,
    gs_Expected = fullresults$Expected,
    stringsAsFactors = FALSE
  )

  rownames(mydf) <- mydf$gs_id

  return(mydf)
}


#' .DeeDeefy_gsenrichResult() , a function based on the original shaker
#' for gseaResult objects in GeneTonic
#'
#' @param obj An `gseaResult` object, obtained via `clusterProfiler`
#'
#' @returns a data.frame in GeneTonic shaker standard format
#'
#' @noRd
.DeeDeefy_gsenrichResult <- function(obj) {
  if (!is(obj, "gseaResult")) {
    stop("Provided object must be of class `gseaResult`")
  }

  if (is.null(obj@result$core_enrichment)) {
    stop(
      "You are providing an object where the `core_enrichment` is ",
      "not specified, this is required for running GeneTonic properly."
    )
  }

  message(
    "Using the content of the 'core_enrichment' column to generate the",
    " 'gs_genes' for GeneTonic...",
    " If you have that information available directly, please adjust the",
    " content accordingly.",
    "\n\nUsing the set of the 'core_enrichment' size to compute the 'gs_de_count'"
  )

  message("Found ", nrow(obj@result),
          " gene sets in `gseaResult` object, of which ",
          nrow(as.data.frame(obj)), " are significant.")
  message("Converting for usage within the DeeDeeExperiment framework...")

  fullresults <- obj@result

  mydf <- data.frame(
    gs_id = fullresults$ID,
    gs_description = fullresults$Description,
    gs_pvalue = fullresults$pvalue,
    gs_genes = gsub("/", ",", fullresults$core_enrichment),
    gs_de_count = lengths(strsplit(fullresults$core_enrichment, split = "/")),
    gs_bg_count = fullresults$setSize,
    gs_ontology = obj@setType,
    gs_NES = fullresults$NES,
    gs_p.adjust = fullresults$p.adjust,
    gs_qvalue = fullresults$qvalue,
    stringsAsFactors = FALSE
  )

  rownames(mydf) <- mydf$gs_id

  return(mydf)
}


#' .DeeDeefy_fgseaResult() , a function based on the original shaker
#' for fgsea output in GeneTonic
#'
#' @param fgsea_output a data.frame with the output of `fgsea()` in `fgsea`
#'
#' @returns a data.frame in GeneTonic shaker standard format
#'
#' @noRd
.DeeDeefy_fgseaResult <- function(fgsea_output) {
  if (!is(fgsea_output, "data.frame")) {
    stop("fgsea output should be a data.frame!")
  }
  exp_colnames <- c(
    "pathway", "pval", "padj", "ES", "NES",
    "size", "leadingEdge"
  )
  if (!all(exp_colnames %in% colnames(fgsea_output))) {
    stop(
      "I could not find some of the usual column names from the fgsea output.",
      " Maybe you performed additional processing/filtering steps?"
    )
  }

  if (!is(fgsea_output$leadingEdge, "list")) {
    stop("Expecting 'leadingEdge' column to be a list")
  }

  message("Found ", nrow(fgsea_output),
          " gene sets in the file output from fgsea of which ",
          sum(fgsea_output$padj <= 0.05), " are significant (p-value <= 0.05).")
  message("Converting for usage within the DeeDeeExperiment framework...")

  message(
    "Using the content of the 'leadingEdge' column to generate the 'gs_genes'",
    " for GeneTonic...",
    " If you have that information available directly, please adjust the",
    " content accordingly.",
    "\n\nUsing the set of the leadingEdge size to compute the 'gs_de_count'"
  )

  message(
    "\n\nfgsea is commonly returning no identifier for the gene sets used.",
    " Please consider editing the 'gs_id' field manually according to the gene ",
    " set you provided"
  )

  mydf <- data.frame(
    gs_id = fgsea_output$pathway,
    gs_description = fgsea_output$pathway,
    gs_pvalue = fgsea_output$pval,
    gs_genes = vapply(
      fgsea_output$leadingEdge,
      function(arg) paste(arg, collapse = ","), character(1)
    ),
    gs_de_count = lengths(fgsea_output$leadingEdge),
    gs_bg_count = fgsea_output$size,
    gs_NES = fgsea_output$NES,
    gs_adj_pvalue = fgsea_output$padj,
    stringsAsFactors = FALSE
  )

  rownames(mydf) <- mydf$gs_id

  # consider re-sorting by p-value?


  return(mydf)
}


.basic_str_wrap <- function(x, width = 80, ...) {
  paste(strwrap(x, width = width, ...), collapse = "\n")
}



#' Display available FEA formats
#' @returns a data.frame of possible FEA input formats
#'
#' @export
#' @examples
#' supported_fea_formats()
supported_fea_formats <- function() {
  data.frame(
    Format = c(
      "data.frame",
      "enrichResult",
      "gseaResult",
      "fgseaResult",
      "data.frame",
      "data.frame",
      "data.frame",
      "data.frame"
    ),
    Package = c(
      "topGO",
      "clusterProfiler",
      "clusterProfiler",
      "fgsea",
      "gprofiler2",
      "enrichR",
      "DAVID",
      "GeneTonic"
    )
  )
}


#' .is_empty_sce() , checks whether an sce oject is empty or not
#'
#' @param sce a `SingleCellExperiment` object
#'
#' @returns a logical, indicating whether the object is empty or not
#'
#' @noRd
.is_empty_sce <- function(sce) {
  !is.null(sce) &&
    is(sce, "SingleCellExperiment") &&
    nrow(sce) == 0 &&
    ncol(sce) == 0 &&
    length(assays(sce)) == 0
}


#' @noRd
.shake_enrich_res <- function(res_enrich, type) {
  switch(type,
  "topGO" = .DeeDeefy_topGOtableResult(res_enrich),
  "clusterProfiler" = if (is(res_enrich, "enrichResult")) {
    res_enrich_shaken <- .DeeDeefy_enrichResult(res_enrich)
  } else { NULL},
  "GeneTonic" = res_enrich,
  "DAVID" = .DeeDeefy_david(res_enrich),
  "fgsea" = .DeeDeefy_fgseaResult(res_enrich),
  "gsea" = if (is(res_enrich, "gseaResult")) {
    res_enrich_shaken <- .DeeDeefy_gsenrichResult(res_enrich)
  } else { NULL},
  "enrichr" = .DeeDeefy_enrichr(res_enrich),
  "gProfiler" = .DeeDeefy_gprofiler(res_enrich),
  NULL
)
}


#' .fill_rowdata_with_dea() fills the `rowData` of an `sce` with DEA results
#'
#' It pre-fills the columns (logFC, pval, padj) with NAs, and assigns the
#' corresponding DE statistics to the `rowData` of the `sce` object for a
#' given contrast
#'
#' @param sce a `SingleCellExperiment` object
#'
#' @param de_name a character vector indicating the name of the DE contrast
#' @param de_res a named list of DE results
#' @param de_cols a named vector specifying the name of the columns in de_res
#' for each DE statistics. Should have names: `"logFC"`, `"pval"`, `"padj"`
#' @param valid_matches a logical vector indicating which features/rows have
#' valid matches between `sce` and `de_res`
#' @param matched_ids integer vector specifying the positions of
#' matching features between `sce` and `de_res`
#' @return the updated `sce` object with new columns added to the `rowData`
#'
#' @noRd
.fill_rowdata_with_dea <- function(sce,
                                   de_name,
                                   de_res,
                                   de_cols = c(logFC = "log2FoldChange",
                                               pval = "pvalue",
                                               padj = "padj"),
                                   valid_matches,
                                   matched_ids) {

  suffixes <- c("_log2FoldChange", "_pvalue", "_padj")

  col_names <- c("logFC", "pval", "padj")

  for (i in seq_along(col_names)) {
    col_name <- paste0(de_name, suffixes[i])
    # pre-fill with NA
    rowData(sce)[[col_name]] <- NA

    # assign values only for matched indices
    rowData(sce)[[col_name]][valid_matches] <-
      de_res[[de_cols[col_names[i]]]][matched_ids[valid_matches]]
  }
  return(sce)
}



#' Import DE results as a `data.frame`, containing at least the following
#' statistics: "log2FoldChange", "pvalue", and "padj", with that exact name
#'
#' @param sce A `SingleCellExperiment` object
#' @param res_de A set of DE results, provided as `data.frame`
#' @param de_name A character value, describing the contrast of interest. Will
#' be used to compose the column names in the `rowData` slot.
#'
#' @return A list, containing the updated `SingleCellExperiment` object, and the
#' standardized information on the DE analysis, as these are to be used in the
#' `DeeDee` framework.
#'

#' @noRd
.importDE_df <- function(sce, res_de, de_name) {
  # correct object format
  stopifnot(is(res_de, "data.frame"))
  # contain the right columns
  stopifnot(all(c("log2FoldChange", "pvalue", "padj") %in% colnames(res_de)))

  # rownames exist
  if (is.null(rownames(res_de))){
    stop("`res_de` must have rownames (gene/feature IDs)")
  }

  # p value different from NA respect the 0-1 interval
  stopifnot(all(na.omit(res_de$pvalue <= 1)) &
              all(na.omit(res_de$pvalue >= 0)))

  matched_ids <- match(rownames(sce), rownames(res_de)) # we align de res with
  # se only valid indices
  valid_matches <- !is.na(matched_ids)

  sce <- .fill_rowdata_with_dea(sce = sce,
                                de_name = de_name,
                                de_res = res_de,
                                de_cols = c(logFC = "log2FoldChange",
                                            pval = "pvalue",
                                            padj = "padj"),
                                valid_matches = valid_matches,
                                matched_ids = matched_ids)

  pkg <- attr(res_de, "package")
  pkg_ver <- attr(res_de, "package_version")

  if (is.null(pkg)) pkg <- NA
  if (is.null(pkg_ver)) pkg_ver <- NA


  is_limma_multicontrast <- identical(pkg, "limma") &&
    !is.null(attr(res_de, "original_object"))

  if (!is_limma_multicontrast) {
    if (is.null(metadata(sce)$singlecontrast)) {
      metadata(sce)$singlecontrast <- list()
    }
    metadata(sce)$singlecontrast[[de_name]] <- res_de

    original_object <- list(
      metadata_storage = "singlecontrast",
      key   = de_name, #contrast name
      coef  = NULL
    )
  } else {
    original_object <- NULL # the other case is handled by constructor/method
  }


  dea_contrast <- list(
    alpha = NA,
    lfcThreshold = NA,
    metainfo_logFC = NA,
    metainfo_pvalue = NA,
    original_object = original_object,
    package = pkg,
    package_version = pkg_ver
  )

  return(list(sce = sce, dea_contrast = dea_contrast))
}



#' Convert a `MArrayLM` object with multiple contrasts into a list of DE results
#' tables compatible with DeeDeeExperiment
#'
#' This helper function extracts DE results for each contrast contained in a
#' `limma::MArrayLM` object and reformats them into a list of standardized data
#' frames suitable for integration in a `DeeDeeExperiment` object.
#' Each resulting data frame includes renamed columns:
#' `logFC` to `log2FoldChange` ,`P.Value` to `pvalue`, and `adj.P.Val` to `padj`
#'
#' @details
#' The function assumes that each column in `fit$coefficients` corresponds to
#' a contrast of interest. The names of the resulting list elements are taken
#' directly from `colnames(fit$coefficients)`. The names in the input object
#' must therefore accurately reflect the intended contrast names.
#'
#'
#' @param fit A `MArrayLM` object, as produced by the `limma` workflow
#' @param number An integer specifying the maximum number of genes to extract
#' per contrast
#' @param sort.by A character string specifying which statistic to rank the
#' genes by. It must be one of the values accepted by the `sort.by` argument in
#' the `limma::topTable()` function. Defaults to "none".
#'
#' @returns A named list of DE results tables, one per contrast. Each table
#' contains standardized columns and the list is annotated with metadata
#' indicating its package origin (`limma`).
#'
#' @seealso [limma::topTable()]
#'
#' @export
#'
#' @examples
#' data("de_limma", package = "DeeDeeExperiment")
#' new_limma_list <- limma_list_for_dde(de_limma)
limma_list_for_dde <- function(fit,
                               number = nrow(fit),
                               sort.by = "none") {
  # check type
  stopifnot(inherits(fit, "MArrayLM"))

  coef_names <- colnames(fit$coefficients)

  res_list <- setNames(lapply(coef_names, function(x) {
    res <- topTable(fit,
                    coef = x,
                    number = number,
                    sort.by = sort.by)

    rename_cols <- c("logFC" = "log2FoldChange",
                     "P.Value" = "pvalue",
                     "adj.P.Val" = "padj")
    intersect_cols <- intersect(names(rename_cols), colnames(res))
    colnames(res)[match(intersect_cols, colnames(res))] <-
      rename_cols[intersect_cols]

    # keep metadata
    attr(res, "package") <- "limma"
    attr(res, "package_version") <- as.character(packageVersion("limma"))
    attr(res, "original_object") <- fit

    res

  }), coef_names)

  attr(res_list, "package") <- "limma"
  attr(res_list, "package_version") <- as.character(packageVersion("limma"))

  cli::cli_alert_info(
    "Returning {length(res_list)} limma contrasts formatted for DeeDeeExperiment"
    )

  return(res_list)
}




#' Convert `muscat::pbDS()` results into a flat list of data frames compatible
#' with `DeeDeeExperiment`
#'
#' This helper function extracts and flattens the nested structure returned
#' by `muscat::pbDS()`, returning one table per contrast–cluster combination.
#' Each resulting data frame will have standardized column names
#' (`log2FoldChange`, `pvalue`, `padj`).
#'
#' The function is intended to simplify the integration of muscat  results into
#' a `DeeDeeExperiment` object. It automatically renames columns
#' (`logFC` to `log2FoldChange`, `p_val` to `pvalue`, and the selected adjusted
#' p-value column to `padj`) and annotates the resulting list with metadata
#' about the originating package.
#'
#' @details
#' The function checks that each contrast entry contains a valid `table`
#' component as expected from `pbDS()` output. Invalid or empty contrasts are
#' skipped with a warning message.
#' The names of the list elements in `res` must match contrast names found in
#' the `table` slot of each entry.
#'
#'
#' @param res A list, typically the output of `muscat::pbDS()` function,
#' containing one or more contrasts
#' @param padj_col A character string specifying which adjusted p-value column
#' to extract. It can be either "p_adj.loc" or "p_adj.glb".
#'
#' @returns A named list of data frames
#' @export
#'
#' @examples
#' data("muscat_pbDS_res", package = "DeeDeeExperiment")
#' new_muscat_list <- muscat_list_for_dde(list(`stim-ctrl` = muscat_res))
muscat_list_for_dde <- function(res, padj_col = c("p_adj.loc", "p_adj.glb")){

  padj_col <- match.arg(padj_col)

  # check type
  if (!is.list(res)) {
    stop("muscat output should be a list")
  }

  flat_list <- lapply(names(res), function(contrast_name){
    this_contrast <- res[[contrast_name]]

    # check for the correct structure for pbDS outpubt
    if (!("table" %in% names(this_contrast)) ||
        !(contrast_name %in% names(this_contrast$table)) ||
        !is.list(this_contrast$table[[contrast_name]])) {
      cli::cli_alert_warning(
        "Skipping contrast {.val {contrast_name}} . No valid  results tables found.")
      return(NULL)
    }

    tbls <- this_contrast$table[[contrast_name]]

    lapply(names(tbls), function(cell){
      df <- tbls[[cell]]

      stopifnot(is.data.frame(df))

      rename_cols <- c("logFC" = "log2FoldChange",
                       "p_val" = "pvalue")

      rename_cols[padj_col] <- "padj"

      intersect_cols <- intersect(names(rename_cols), colnames(df))
      colnames(df)[match(intersect_cols, colnames(df))] <-
        rename_cols[intersect_cols]

      # assign a unique name, combining the contrast and cluster name
      entry_name <- paste(contrast_name, cell, sep = "_")

      rownames(df) <- df$gene

      # keep metadata
      attr(df, "package") <- "muscat"
      attr(df, "package_version") <- as.character(packageVersion("muscat"))

      setNames(list(df), entry_name)

    })
  })

  flat_list <- unlist(unlist(flat_list, recursive = FALSE), recursive = FALSE)

  if (is.null(flat_list) || length(flat_list) == 0)
    stop("No valid muscat results found to convert")


  cli::cli_alert_info(
    "Returning {length(flat_list)} muscat contrast-cluster tables formatted for DeeDeeExperiment"
  )

  return(flat_list)

}



#' Handle limma multi-contrast results formatted for DeeDeeExperiment
#'
#'This internal helper takes the list of contrast-wise `data.frame`s produced by
#' `limma_list_for_dde()` and:
#'
#' * stores the original `MArrayLM` fit in `metadata(sce)$multicontrast`
#'   under the user-supplied `entry_name`, and
#' * imports each contrast table into `rowData(sce)` via `.importDE_df()`,
#'   assembling a named list of DEA contrast metadata to populate the `dea` slot
#'
#' The function assumes that `de_list` already has standardized column names
#' (`log2FoldChange`, `pvalue`, `padj`) and carries the attributes
#' `"package" = "limma"`, `"package_version"`, and `"original_object"` as set
#' by `limma_list_for_dde()`.
#'
#' @param sce A `singleCellExperiment` object in which the DE results will be
#' integrated
#'
#' @param de_list A named list of contrast-specific DE result tables, typically
#' the output of `limma_list_for_dde()`
#' @param entry_name A character string indicating the name under which the
#' original `MArrayLM` fit will be stored in `metadata(sce)$multicontrast`. This
#' is usually the object name as supplied by the user.
#'
#' @return A list containing the update `SingleCellExperiment` object, and
#' a named list of DEA contrasts ready to be merged into the dea slot of a dde
#' object
#'
#' @noRd
.handle_limma_list <- function(sce, de_list, entry_name) {
  #checks
  stopifnot(is.list(de_list))
  stopifnot(identical(attr(de_list, "package"), "limma"))

  first_de <- de_list[[1L]]
  fit <- attr(first_de, "original_object")

  if (is.null(metadata(sce)$multicontrast)) {
    metadata(sce)$multicontrast <- list()
  }

  if (!is.null(names(metadata(sce)$multicontrast)) &&
      entry_name %in% names(metadata(sce)$multicontrast)) {
    stop(
      "A limma multicontrast object named '", entry_name,
      "' already exists in metadata(sce)$multicontrast."
    )
  }

  metadata(sce)$multicontrast[[entry_name]] <- fit

  dea_contrasts <- list()

  for (i in names(de_list)) {
    sub_de <- de_list[[i]]

    imported <- .importDE_df(sce, sub_de, i)
    sce <- imported$sce

    dea_contrasts[[i]] <- imported$dea_contrast
    dea_contrasts[[i]]$package <- "limma"
    dea_contrasts[[i]]$package_version <- attr(sub_de, "package_version")

    dea_contrasts[[i]]$original_object <- list(
      metadata_storage = "multicontrast",
      key   = entry_name,  # name in metadata(sce)$multicontrast
      coef  = i          # which contrast inside the fit
    )
  }

  list(sce = sce, dea_contrasts = dea_contrasts)
}


#' Export DEA/FEA/ASSAY results from a `DeeDeeExperiment` to excel files
#'
#' Extracts DEA and/or FEA results stored in a `DeeDeeExperiment` and writes
#' them to excel files. Each contrast/result is written to a separate sheet.
#' DEA and FEA are written to separate files.
#'
#' @param x A `DeeDeeExperiment` object containing DEA/FEA results to be
#' extracted
#' @param res_type A character vector indicating the result to extract
#' (i.e. `"dea"`, `"fea"`, `"assay"`) from the corresponding slot
#' If set to `"all"`, DEAs, FEAs, as well as the specified assays are extracted
#' into separate excel files
#' @param res_format A character string specifying the DEA/FEAs output format
#' to be exported (i.e. `"minimal"` or `"original"`)
#' @param output_dir A character string specifying the directory where the excel
#' file will be written
#' @param file_name Optional character string specifying a file name. If `NULL`,
#' a default name is generated
#' @param assay_type Optional character vector specifying the assays to export.
#' It defaults to `"counts"`
#' @param force Logical. If `TRUE`, an existing file with the same name will be
#' overwritten
#'
#' @return An (invisible) character vector of file paths to the exported excel
#' files. The vector may contain one path `res_type = "dea"` or
#' `res_type = "fea" `or three paths `res_type = "all"`
#'
#' @export
#'
#' @examples
#' data("de_named_list", package = "DeeDeeExperiment")
#' data("topGO_results_list", package = "DeeDeeExperiment")
#' dde <- DeeDeeExperiment(de_results = de_named_list, enrich_results = topGO_results_list)
#'
#' export_result_for_dde(dde,
#' res_type = "dea",
#' res_format = "minimal", output_dir = tempdir(), force = TRUE)
#'
#' export_result_for_dde(dde,
#' res_type = c("dea", "fea"),
#' res_format = "original", output_dir = tempdir(), force = TRUE)
#'
#'
export_result_for_dde <- function(x,
                          res_type = c("assay","dea", "fea", "all"),
                          res_format = c("minimal", "original"),
                          output_dir = getwd(),
                          file_name = NULL,
                          assay_type = c("counts"),
                          force = FALSE) {
  # checks on the args
  if (!is(x, "DeeDeeExperiment")) {
    stop("`x` must be a `DeeDeeExperiment` object!")
  }

  if (!is.character(output_dir) || length(output_dir) != 1 ||
      output_dir == "") {
    stop("`output_dir` must be a non empty character string!")
  }

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  if (!is.null(file_name) && (!is.character(file_name) ||
                              length(file_name) != 1 || file_name == "")) {
    stop("`file_name` must be a non empty character string!")
  }

  if (!is.logical(force)) {
    stop("`force` must be logical!")
  }

  if(!is.character(assay_type)) {
    stop("`assay_type` must be a character vector!")
  }

  res_type <- match.arg(res_type, several.ok = TRUE)
  res_format <- match.arg(res_format)

  out_paths <- character()

  # export results
  if (any(res_type %in% c("dea", "all"))) {
    # get all deas
    deas <- getDEAList(x, format = res_format)
    out_paths <- c(out_paths, .write_res(list = deas,
                                         result_type = "DEA",
                                         res_format = res_format,
                                         output_dir = output_dir,
                                         file_name = file_name,
                                         force = force))
  }

  if (any(res_type %in% c("fea", "all"))) {
    # get all feas
    feas <- getFEAList(x, format = res_format)
    out_paths <- c(out_paths, .write_res(list = feas,
                                         result_type = "FEA",
                                         res_format = res_format,
                                         output_dir = output_dir,
                                         file_name = file_name,
                                         force = force))
  }


  if (any(res_type %in% c("assay", "all"))) {

    if (!all(assay_type %in% SummarizedExperiment::assayNames(x))) {
      missing <- setdiff(assay_type, SummarizedExperiment::assayNames(x))
      stop("assay not found in `x`: ", paste(missing, collapse = ", "))
    }

    assays <- setNames(
      lapply(assay_type, function(i){
      assay(x, i)}),
      assay_type)

    out_paths <- c(out_paths, .write_res(list = assays,
                                         result_type = "ASSAY",
                                         res_format = res_format,
                                         output_dir = output_dir,
                                         file_name = file_name,
                                         force = force))

  }


  if (length(out_paths) == 0) {
    stop("No results were exported!")
  }

  invisible(out_paths)

}


#' Write DEA/FEA/ASSAY results to an excel file
#'
#' @param list A named list of data frames (or coercible objects), where each
#' element represents one contrast/result/assay and will be written to a
#' separate sheet
#' @param result_type A character vector indicating the result to extract
#' (i.e. `"DEA"`, `"FEA"`, or `"ASSAY"`). It is used in the output file name
#' @param res_format A character string, specifying the DEA/FEAs output format
#' to export (i.e. `"minimal"` or `"original"`). It is used in the output file
#' name
#' @param output_dir A character string specifying the directory where the excel
#' file will be written
#' @param file_name Optional character string specifying a file name. If `NULL`,
#' a default name is generated
#' @param force Logical. If `TRUE`, an existing file with the same name will be
#' overwritten
#'
#' @return A character string giving the path to the written excel file,
#' or `NULL` if there were no results to export
#'
#' @noRd
.write_res <- function(list,
                       result_type,
                       res_format,
                       output_dir,
                       file_name,
                       force) {
  message(paste("Found", length(list), result_type, "results"))

  if (length(list) == 0) {
    warning(paste("No", result_type, "results to export."), call. = FALSE)
    return(NULL)
  }

  # add rownames as a column
  list <- lapply(list, function(df){
    rn <- rownames(df)
    if (!is.null(rn)) {
      col_name <- "id"
      col_name <- make.unique(c(names(df), col_name))[length(names(df)) + 1]
      df <- cbind(
        setNames(data.frame(rn, stringsAsFactors = FALSE), col_name),
        df)
    }
    df
  })

  sheet_name <- names(list)
  sheet_name <- make.unique(.clean_sheet_names(sheet_name))
  names(list) <- sheet_name

  prefix <- if (is.null(file_name)) {paste("dde", res_format, sep = "_")
  } else {
    file_name
  }

  base <- paste(prefix, result_type, sep = "_")
  out_file <- file.path(output_dir, paste0(base, ".xlsx"))

  if (file.exists(out_file) && !isTRUE(force)) {
    stop("File already exists: ", out_file,
         "\nSet `force = TRUE` to replace it.")
  }

  cli::cli_alert_success(
    "Writing results to:  {out_file} "
  )
  writexl::write_xlsx(list, path = out_file)

  out_file
}



#' Clean excel sheet names
#'
#' @param x A character vector of proposed Excel sheet names
#'
#' @return A character vector of excel-safe sheet names
#'
#' @noRd
.clean_sheet_names <- function(x) {
  # excel sheet name rules: max 31 chars; cannot contain : \ / ? * [ ]
  x <- gsub("[:\\\\/\\?\\*\\[\\]]", "_", x)
  substr(x, 1, 31)
}
