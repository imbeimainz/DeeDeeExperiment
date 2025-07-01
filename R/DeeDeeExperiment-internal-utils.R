#' Import from `DESeq2` DE results
#'
#' @param se A `SummarizedExperiment` object
#' @param res_de A set of DE results, provided as `DESeqResults` as in the
#' `DESeq2` framework
#' @param de_name A character value, describing the contrast of interest. Will
#' be used to compose the column names in the `rowData` slot.
#'
#' @return A list, containing the updated `SummarizedExperiment` object, and the
#' standardized information on the DE analysis, as these are to be used in the
#' `DeeDee` framework.
#'
#' @noRd
.importDE_DESeq2 <- function(se, res_de, de_name) {
  # correct object format
  stopifnot(is(res_de, "DESeqResults"))
  # contain the right columns
  stopifnot(all(c("log2FoldChange", "pvalue", "padj") %in% colnames(res_de)))

  # contain the feature ids

  # p value different from NA respect the 0-1 interval
  stopifnot(all(na.omit(res_de$pvalue <= 1)) &
              all(na.omit(res_de$pvalue > 0)))

  matched_ids <- match(rownames(se), rownames(res_de)) # we align de res with se
  # only valid indices
  valid_matches <- !is.na(matched_ids)


  # Pre-fill rowData with NA
  rowData(se)[[paste0(de_name, "_log2FoldChange")]] <- NA
  rowData(se)[[paste0(de_name, "_pvalue")]] <- NA
  rowData(se)[[paste0(de_name, "_padj")]] <- NA


  # assign values only for matched indices, to have on both sides the
  # same length. we keep NA for unmatched genes
  rowData(se)[[paste0(de_name, "_log2FoldChange")]][valid_matches] <-
    res_de$log2FoldChange[matched_ids[valid_matches]]
  rowData(se)[[paste0(de_name, "_pvalue")]][valid_matches] <-
    res_de$pvalue[matched_ids[valid_matches]]
  rowData(se)[[paste0(de_name, "_padj")]][valid_matches] <-
    res_de$padj[matched_ids[valid_matches]]


  dea_contrast <- list(
    alpha = metadata(res_de)$alpha,
    lfcThreshold = metadata(res_de)$lfcThreshold,
    metainfo_logFC = mcols(res_de)$description[colnames(res_de) == "log2FoldChange"],
    metainfo_pvalue = mcols(res_de)$description[colnames(res_de) == "pvalue"],
    original_object = res_de,
    package = "DESeq2",
    package_version = packageVersion("DESeq2")
  )

  return(list(se = se, dea_contrast = dea_contrast))
}


#' Import from edgeR DE results
#'
#' @param se A SummarizedExperiment object
#' @param res_de A set of DE results, provided by the `edgeR` framework (either
#' a `DGEExact` or a `DGELRT` object).
#' @param de_name A character value, describing the contrast of interest. Will
#' be used to compose the column names in the rowData slot.
#'
#' @return A list, containing the updated SummarizedExperiment object, and the
#' standardized information on the DE analysis, as these are to be used in the
#' DeeDee framework.
#'
#' @noRd
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


  matched_ids <- match(rownames(se), rownames(res_tbl)) # we align de res with
  # se only valid indices
  valid_matches <- !is.na(matched_ids)

  # pre-fill rowData with NA the assign the corresponding values only for
  # matched indices for logFC, accounting for the fact that the logFC column
  # name in edgeR depends on whether we have 1 or multiple contrasts
  for (i in logFC_cols) {
    rowData(se)[[paste0(de_name, "_log2FoldChange")]] <- NA
    # assign correspionding values
    rowData(se)[[paste0(de_name, "_log2FoldChange")]][valid_matches] <-
      res_tbl$table[[i]][matched_ids[valid_matches]]
  }

  # pre-fill rowData with NA the assign the corresponding values for matched
  # indices for pval and padj
  rowData(se)[[paste0(de_name, "_pvalue")]] <- NA
  rowData(se)[[paste0(de_name, "_padj")]] <- NA


  # assign values only for matched indices, to have on both sides the
  # same length. we keep NA for unmatched genes
  rowData(se)[[paste0(de_name, "_pvalue")]][valid_matches] <-
    res_tbl$table$PValue[matched_ids[valid_matches]]
  rowData(se)[[paste0(de_name, "_padj")]][valid_matches] <-
    res_tbl$table$FDR[matched_ids[valid_matches]]

  dea_contrast <- list(
    alpha = NA,
    lfcThreshold = NA,
    metainfo_logFC = res_tbl$comparison,
    metainfo_pvalue = paste0("p-value adjusted using ", res_tbl$adjust.method),
    original_object = res_de,
    package = "edgeR",
    package_version = packageVersion("edgeR")
  )

  return(list(se = se, dea_contrast = dea_contrast))
}



#' Import from `limma` DE results
#'
#' @param se A `SummarizedExperiment` object
#' @param res_de A set of DE results, provided in the `limma` framework
#' (a `MArrayLM` object).
#' @param de_name A character value, describing the contrast of interest. Will
#' be used to compose the column names in the `rowData` slot.
#'
#' @return A list, containing the updated `SummarizedExperiment` object, and the
#' standardized information on the DE analysis, as these are to be used in the
#' `DeeDee` framework.
#'
#' @noRd
.importDE_limma <- function(se, res_de, de_name) {
  # checks object
  stopifnot(is(res_de, "MArrayLM"))

  # make sure there are at least 2 coefficients
  if (ncol(res_de$coefficients) < 2) {
    # we still need to manage the handling of 1 contrast
    stop(
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

  matched_ids <- match(rownames(se), rownames(res_tbl)) # we align de res with
  # se only valid indices
  valid_matches <- !is.na(matched_ids)

  # Pre-fill rowData with NA
  rowData(se)[[paste0(de_name, "_log2FoldChange")]] <- NA
  rowData(se)[[paste0(de_name, "_pvalue")]] <- NA
  rowData(se)[[paste0(de_name, "_padj")]] <- NA


  # assign values only for matched indices, to have on both sides the
  # same length. we keep NA for unmatched genes
  rowData(se)[[paste0(de_name, "_log2FoldChange")]][valid_matches] <-
    res_tbl$logFC[matched_ids[valid_matches]]
  rowData(se)[[paste0(de_name, "_pvalue")]][valid_matches] <-
    res_tbl$P.Value[matched_ids[valid_matches]]
  rowData(se)[[paste0(de_name, "_padj")]][valid_matches] <-
    res_tbl$adj.P.Val[matched_ids[valid_matches]]

  dea_contrast <- list(
    alpha = NA,
    lfcThreshold = NA,
    metainfo_logFC = NA,
    metainfo_pvalue = NA,
    original_object = res_de,
    package = "limma",
    package_version = packageVersion("limma")
  )

  return(list(se = se, dea_contrast = dea_contrast))
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
    stop("All elements in the list must be of type DESeqResults,",
         " DGEExact, DGELRT, or MArrayLM. Alternatively, it can be a data.frame",
         " with at least a logFC, p-value and p-adjusted value columns.")
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

  if (is.null(entry_name)) {
    stop("You must provide a name for your enrichment results!")
  }

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
    return(NA_character_) # this needs to be character
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
    GeneRatio = fullresults$GeneRatio,
    BgRatio = fullresults$BgRatio,
    p.adjust = fullresults$p.adjust,
    qvalue = fullresults$qvalue,
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
    Expected = fullresults$Expected,
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
