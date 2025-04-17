#' @name DeeDeeExperiment-methods
#'
#' @title Methods for [DeeDeeExperiment] objects
#'
#' @aliases
#' dea_info
#' dea_info<-
#' dea_names
#' add_dea
#' remove_dea
#' dea
#' get_dea_list
#' fea_info
#' fea_info<-
#' fea_names
#' add_fea
#' remove_fea
#'
#' @description
#' The [DeeDeeExperiment()] class provides a family of methods to get
#' and set DE-related information and functional enrichment results in
#' [DeeDeeExperiment] objects.
#'
#' @param x A [DeeDeeExperiment()] object
#' @param value Replacement value for replacement methods.
#' @param dea A named list of DE results, in any of the formats supported by
#' the package (currently: results from DESeq2, edgeR, limma).
#' @param dea_name Character value, specifying the name of the DE analysis to
#' get or remove
#' @param verbose Logical, whether or not to display warnings. If TRUE, warnings/messages
#' will be displayed. If FALSE, the function runs silently
#' @param fea A named list of Functional Enrichment results
#' @param fea_name Character value, specifying the name of the functional enrichment
#' result to add or remove
#' @param fea_res A data frame or a named list of data frames containing fea results.
#' Each element should represent a table of enrichment terms (e.g., GO terms).
#' @param de_name A character string to explicitly specify the name of the de result this fea should be linked to.
#' If not provided, the function will attempt to match fea names to de results automatically.
#' @param fe_name A character string giving a name to the FE results.
#' It can only support `FDR` to adjust the threshold used for summarizing DE results.
#'
#' @return Return value varies depending on the individual methods, as described
#' below.
#'
#' @details
#' * `dea_info` and `dea_info<-` are the methods to get and set the `dea` information as a
#' whole. These methods return `DeeDeeExperiment` objects.
#' * `add_dea` and `remove_dea` are used to respectively add or remove DE-results
#' items. These methods also return `DeeDeeExperiment` objects, with updated
#' content in the `dea` slot.
#' * `dea_names` returns the names of the available DE contrasts in `DeeDeeExperiment` objects.
#' * `dea` and `get_dea_list` retrieve the `dea` information and provide
#' this as a `DataFrame` object (for a specific analysis) or as a list, with one
#' element for each reported analysis.
#' * `fea` and `fea<-` are the methods to get and set the `fea` information as a
#' whole. These methods return `DeeDeeExperiment` objects.
#' * `fea_names` returns the names of the available enrichment results in `DeeDeeExperiment` objects.
#' * `add_fea` and `remove_fea` are used to respectively add or remove functional
#' enrichment results items. These methods also return `DeeDeeExperiment` objects, with updated
#' content in the `fea` slot.
#' * `show` is the method to nicely print out the information of a `DeeDeeExperiment`
#' object.
#'
#' @examples
#' data("de_named_list", package = "DeeDeeExperiment")
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
#' # creating a `DeeDeeExperiment`
#' dde <- DeeDeeExperiment(
#'   se_macrophage_noassays,
#'   de_results = de_named_list
#' )
#' dde
#'
#' new_del <- list(
#'   ifng2 = de_named_list$ifng_vs_naive,
#'   ifngsalmo2 = de_named_list$ifngsalmo_vs_naive
#' )
#'
#' # add a new (set of) DE result(s)
#' dde_new <- add_dea(dde, new_del)
#' dde_new
#'
#' # removing DEAs
#' dde_removed <- remove_dea(dde, "ifng_vs_naive")
#' dde_removed
NULL


# dea slot - get & set ---------------------------------------------------------

#' @rdname DeeDeeExperiment-methods
#' @export
setMethod("dea_info",
          signature = "DeeDeeExperiment",
          definition = function(x) {
            x@dea
          })


#' @rdname DeeDeeExperiment-methods
#' @export
setReplaceMethod("dea_info",
                 signature = c("DeeDeeExperiment", "ANY"),
                 definition = function(x, value) {
                   x@dea <- value
                   validObject(x)
                   x
                 })


# dea info - add, remove, get --------------------------------------------------

#' @rdname DeeDeeExperiment-methods
#' @export
setMethod("dea_names",
          signature = "DeeDeeExperiment",
          definition = function(x){
            names(dea_info(x))
          }
          )

#' @rdname DeeDeeExperiment-methods
#' @export
setMethod("add_dea",
          signature = c("DeeDeeExperiment", "ANY"),
          definition = function(x, dea) {

            # capture name inside the env where the func is called
            entry_name <- deparse(substitute(dea, env = parent.frame()))

            # check and preocess dea
            dea <- .check_de_results(dea, entry_name)
            names(dea)
            names(dea_info(x))

            # dde must be a DeeDeeExp
            if (!is(x, "DeeDeeExperiment")) {
              stop("x must be DeeDeeExperiment object!")
            }
            # dea must be named list
            if (is.null(names(dea))) {
              stop("All elements in dea list must have names!")
            }

            # check that names are all unique, and do not overlap with the existing ones
            if (anyDuplicated(c(names(dea), names(dea_info(x))))) {
              stop("Names in dea must be unique!")
            }

            dea_contrasts <- dea_info(x)
            dde_ids <- rownames(x)

            # update rowData, naming them correctly
            for (i in names(dea)) {
              this_de <- dea[[i]]

              # do different things according to what these objects are
              if (is(this_de, "DESeqResults")) {


                matched_ids <- match(rownames(x), rownames(this_de)) # we align de res with se
                # only valid indices
                valid_matches <- !is.na(matched_ids)


                # Pre-fill rowData with NA
                rowData(x)[[paste0(i, "_log2FoldChange")]] <- NA
                rowData(x)[[paste0(i, "_pvalue")]]         <- NA
                rowData(x)[[paste0(i, "_padj")]]           <- NA


                # assign values only for matched indices, to have on both sides the
                # same length. we keep NA for unmatched genes
                rowData(x)[[paste0(i, "_log2FoldChange")]][valid_matches] <- this_de$log2FoldChange[matched_ids[valid_matches]]
                rowData(x)[[paste0(i, "_pvalue")]][valid_matches]         <- this_de$pvalue[matched_ids[valid_matches]]
                rowData(x)[[paste0(i, "_padj")]][valid_matches]           <- this_de$padj[matched_ids[valid_matches]]


                dea_contrasts[[i]] <- list(
                  alpha = metadata(this_de)$alpha,
                  lfcThreshold = metadata(this_de)$lfcThreshold,
                  metainfo_logFC = mcols(this_de)$description[colnames(this_de) == "log2FoldChange"],
                  metainfo_pvalue = mcols(this_de)$description[colnames(this_de) == "pvalue"],
                  original_object = this_de,
                  package = "DESeq2"
                )
              } else if (is(this_de, "DGEExact") || is(this_de, "DGELRT")) {
                res_tbl <- topTags(
                  this_de,
                  n = nrow(this_de),
                  sort.by = "none"
                )

                # p value different from NA respect the 0-1 interval
                stopifnot(all(na.omit(res_tbl$PValue <= 1)) &
                            all(na.omit(res_tbl$PValue > 0)))

                # identify the logFC cols
                logFC_cols <- grep("^logFC", colnames(res_tbl), value = TRUE)

                matched_ids <- match(rownames(x), rownames(res_tbl)) # we align de res with se
                # only valid indices
                valid_matches <- !is.na(matched_ids)

                # pre-fill rowData with NA the assign the corresponding values only for matched
                # indices for logFC, accounting for the fact that the logFC column name in edgeR
                # depends on whether we have 1 or multiple contrasts
                for (j in logFC_cols) {
                  rowData(x)[[paste0(i, "_log2FoldChange")]] <- NA
                  # assign corresponding values
                  rowData(x)[[paste0(i, "_log2FoldChange")]][valid_matches] <- res_tbl$table[[j]][matched_ids[valid_matches]]
                }

                # pre-fill rowData with NA the assign the corresponding values for matched indices for pval and padj
                rowData(x)[[paste0(i, "_pvalue")]]         <- NA
                rowData(x)[[paste0(i, "_padj")]]           <- NA


                # assign values only for matched indices, to have on both sides the
                # same length. we keep NA for unmatched genes
                rowData(x)[[paste0(i, "_pvalue")]][valid_matches]         <- res_tbl$table$PValue[matched_ids[valid_matches]]
                rowData(x)[[paste0(i, "_padj")]][valid_matches]           <- res_tbl$table$FDR[matched_ids[valid_matches]]


                #print(colnames(rowData(x)))

                # store metadata
                dea_contrasts[[i]] <- list(
                  alpha = NA,
                  lfcThreshold = NA,
                  metainfo_logFC = res_tbl$comparison,
                  metainfo_pvalue = NA,
                  original_object = this_de,
                  package = "edgeR"
                )
              } else if (is(this_de, "MArrayLM")) {
                res_tbl <- topTable(
                  this_de,
                  coef    = 2,
                  number  = nrow(this_de),
                  sort.by = "none"
                )

                # p value different from NA respect the 0-1 interval
                stopifnot(all(na.omit(res_tbl$P.Value <= 1)) &
                            all(na.omit(res_tbl$P.Value > 0)))

                matched_ids <- match(rownames(x), rownames(res_tbl)) # we align de res with se
                # only valid indices
                valid_matches <- !is.na(matched_ids)


                # Pre-fill rowData with NA
                rowData(x)[[paste0(i, "_log2FoldChange")]] <- NA
                rowData(x)[[paste0(i, "_pvalue")]]         <- NA
                rowData(x)[[paste0(i, "_padj")]]           <- NA


                # assign values only for matched indices, to have on both sides the
                # same length. we keep NA for unmatched genes
                rowData(x)[[paste0(i, "_log2FoldChange")]][valid_matches] <- res_tbl$logFC[matched_ids[valid_matches]]
                rowData(x)[[paste0(i, "_pvalue")]][valid_matches]         <- res_tbl$P.Value[matched_ids[valid_matches]]
                rowData(x)[[paste0(i, "_padj")]][valid_matches]           <- res_tbl$adj.P.Val[matched_ids[valid_matches]]


                # matched_ids <- match(rownames(x), rownames(res_tbl))
                #
                # # if not tested, add NA - everywhere? -> pre-fill?
                # rowData(x)[[paste0(i,"_log2FoldChange")]] <- NA
                # rowData(x)[[paste0(i,"_pvalue")]]         <- NA
                # rowData(x)[[paste0(i,"_padj")]]           <- NA
                #
                # # populate using limma columns
                # rowData(x)[[paste0(i,"_log2FoldChange")]][!is.na(matched_ids)] <- res_tbl$logFC
                # rowData(x)[[paste0(i,"_pvalue")]][!is.na(matched_ids)]         <- res_tbl$P.Value
                # rowData(x)[[paste0(i,"_padj")]][!is.na(matched_ids)]           <- res_tbl$adj.P.Val

                # store metadata
                dea_contrasts[[i]] <- list(
                  alpha = NA,
                  lfcThreshold = NA,
                  metainfo_logFC = NA,
                  metainfo_pvalue = NA,
                  original_object = this_de,
                  package = "limma"
                )
              }
              else {
                stop("The dea result class '",i,
                     "' is not recognized (supported classes: DESeqResults, MArrayLM, DGEExact and DGELRT)")
              }
            }

            # update the dea slot
            dea_info(x) <- dea_contrasts

            # check here the validity
            validObject(x)

            # return the object
            return(x)
          }
)


# TODO: might need one where I also simply add ONE single DE object, and that gets autoconverted to a named list (of length 1)
## this one was half addressed, dede accepts 1 single DE object now, the corresponding name is still to be generated
## also add_dea() need a way to handle adding again 1 entry


#' @rdname DeeDeeExperiment-methods
#' @export
setMethod("remove_dea",
          signature = c("DeeDeeExperiment", "character"),
          definition = function(x, dea_name) {
            # x must be a DeeDeeExp

            # dea must be char vector
            deas <- names(dea_info(x))

            deas_to_remove <- intersect(dea_name, deas)

            # warning() if nothing to remove
            if(length(deas_to_remove) == 0){
              warning("No matching dea entries found to remove.")
            }

            for (i in deas_to_remove) {
              cols_to_remove <- c(paste0(i, c("_log2FoldChange", "_pvalue", "_padj")))
              rowData(x) <- rowData(x)[, !(colnames(rowData(x)) %in% cols_to_remove)]
              # update the de slot
              dea_info(x)[[i]] <- NULL
            }

            # here check some validity?
            validObject(x)

            # return the object
            return(x)
          }
)




#' @rdname DeeDeeExperiment-methods
#' @export
setMethod("dea",
          signature = c("DeeDeeExperiment", "character"),
          definition = function(x,
                                dea_name,
                                verbose = TRUE) {
            deas <- dea_info(x)
            dea_names <- names(deas)

            if (!(dea_name %in% dea_names)) {
              stop("dea not found")
            }

            #

            rd_info <- paste0(dea_name,
                              c("_log2FoldChange", "_pvalue", "_padj"))
            #print(rd_info)

            # if (! all(rd_info %in% colnames(rowData(x)))) {
            #   stop("Columns not found")
            # }

            # check for missing columns, for a more precise feedback on the error
            missing_cols <- rd_info[!rd_info %in% colnames(rowData(x))]
            #print(missing_cols)


            # maybe check for rowname mismatches potential gene version issue?
            # maybe not interesting to print back all missmatches in casee all rownames
            # dont match
            rownames_x <- rownames(rowData(x))
            rownames_y <- rownames(dea_info(x)[[dea_name]][["original_object"]])
            mismatched_rows <- sum(!rownames_x %in% rownames_y)

            affected_deas <- character()
            if (mismatched_rows > 0) {
              affected_deas <- c(affected_deas, dea_name)
            }

            if (length(affected_deas) > 0) {
              if (verbose)
                warning(
                  "Mismatch detected between `rownames(rowData(x))` and rownames for the following dea element(s): ",
                  paste(unique(affected_deas), collapse = ", ")
                )
            }

            if (length(missing_cols) > 0) {
              stop("The following columns are missing: ",
                   paste(missing_cols, collapse = ", "))
            }



            out <- rowData(x)[, rd_info]

            return(out)
          }
)

#' @rdname DeeDeeExperiment-methods
#' @export
setMethod("get_dea_list",
          signature = c("DeeDeeExperiment"),
          definition = function(x, verbose = TRUE) {
            verbose <- verbose
            deas <- dea_info(x)
            dea_names <- names(deas)

            dea_list <- list()
            mismatch_info <- list()
            affected_deas <- character()

            for (i in dea_names) {
              dea_list[[i]] <- as.data.frame(dea(x, i, verbose))
              colnames(dea_list[[i]]) <- c("log2FoldChange", "pvalue", "padj")

              # maybe check for rowname mismatches potential gene version issue?
              # maybe not interesting to print back all missmatches in casee all rownames
              # dont match
              rownames_x <- rownames(rowData(x))
              rownames_y <- rownames(deas[[i]][["original_object"]])

              mismatched_rows <- sum(!rownames_x %in% rownames_y)

              if (mismatched_rows > 0) {
                affected_deas <- c(affected_deas, i)
              }
            }


            # if (length(affected_deas) > 0) {
            #   warning(
            #     "Mismatch detected between `rownames(rowData(x))` and rownames for dea element(s): ",
            #     paste(unique(affected_deas), collapse = ", ")
            #   )
            # } # not needed since the warnings will be triggered from dea

            return(dea_list)
          }
)



# fea slot - get & set ---------------------------------------------------------

#' @rdname DeeDeeExperiment-methods
#' @export
setMethod("fea_info",
          signature = "DeeDeeExperiment",
          definition = function(x) {
            x@fea
          })

#' @rdname DeeDeeExperiment-methods
#' @export
setReplaceMethod("fea_info",
                 signature = c("DeeDeeExperiment", "ANY"),
                 definition = function(x, value) {
                   x@fea <- value
                   validObject(x)
                   x
                 })


# fea info - add, remove, get --------------------------------------------------
#' @rdname DeeDeeExperiment-methods
#' @export
setMethod("fea_names",
          signature = "DeeDeeExperiment",
          definition = function(x){
            names(fea_info(x))
          }
)

#' @rdname DeeDeeExperiment-methods
#' @export
setMethod("add_fea",
          signature = c("DeeDeeExperiment",
                        "ANY"), # or should i force it to data.frame??
          definition = function(x,
                                fea_res,
                                de_name= NA_character_,
                                fe_name= NULL) {

            # x must be a DeeDeeExperiment
            if (!is(x, "DeeDeeExperiment")) {
              stop("x must be DeeDeeExperiment object!")
            }

            # capture name inside the env where the func is called
            entry_name <- deparse(substitute(fea_res, env = parent.frame()))

            # TODO
            # check and preocess fea
            fea_list <- .check_enrich_results(fea_res, entry_name)

            # fea must be named list
            if (is.null(names(fea_list))) {
              stop("All elements in 'fea_res' list must have names!")
            }
            # check that names are all unique, and do not overlap with the existing ones
            if (anyDuplicated(c(names(fea_list), names(fea_info(x))))) {
              stop("Names in fea must be unique!")
            }
            # get existing results in the fea slot
            fea_contrasts <- fea_info(x)

            for (fe in names(fea_list)) {
              res_enrich <- fea_list[[fe]]
              if (!is.null(dea_info(x)) &&
                  length(dea_info(x)) > 0) {
                matched_name <- .match_fe_to_de(fe, names(dea_info(x)))
                if (!is.na(matched_name) &&
                    matched_name %in% names(dea_info(x))) {
                  de_res_name <- matched_name
                  if (fe != matched_name) {
                    message("FEA '", fe, "' matched to DE contrast '", matched_name, "'")
                  }
                } else {
                  de_res_name <- NA_character_
                  warning(
                    "Could not match FEA '",
                    fe,
                    "' to a DE contrast.\n",
                    "Available DE results: ",
                    paste(names(dea_info(x)), collapse = ", "),
                    "\n",
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

              if ("GO.ID" %in% colnames(res_enrich)) {
                fe_type <- "TopGo"
              } else {
                fe_type <- "NULL" # placeholder for now
              }

              fea_contrast <- list(
                de_name = de_res_name, # links to de result
                fe_name = fe,
                original_object = res_enrich,
                GeneTonicList = NULL , # we'll put back a GT ready obj
                fe_type = fe_type
              )
              fea_contrasts[[fe]] <- fea_contrast
            }
            # update the fea slot
            fea_info(x) <- fea_contrasts
            # check here the validity
            validObject(x)
            # return the object
            return(x)
          }
)


#' @rdname DeeDeeExperiment-methods
#' @export
setMethod("remove_fea",
          signature = c("DeeDeeExperiment", "character"),
          definition = function(x, fea_name) {

            # x must be a DeeDeeExp
            if(!is(x,"DeeDeeExperiment")) {
              stop("x must be a DeeDeeExperiment object!")
            }
            # fea must be char vector
            if(!is(fea_name,"character")) {
              stop("'fea_name' must be a character vector!") # is it necessary??? since we define that in signature
            }
            feas <- fea_names(x)

            feas_to_remove <- intersect(fea_name, feas)

            # warning() if nothing to remove
            if(length(feas_to_remove) == 0){
              warning("No matching fea entries found to remove.")
            }

            for (i in feas_to_remove) {
              # update the de slot
              fea_info(x)[[i]] <- NULL
            }

            # here check some validity?
            validObject(x)

            # return the object
            return(x)
          }
)







# misc - show & more ------------------------------------------------------

#' @name DeeDeeExperiment-misc
#'
#' @title Miscellaneous DeeDeeExperiment methods
#'
#' @description
#' Miscellaneous methods for the \code{\link{DeeDeeExperiment}} class and its
#' descendants that do not fit into any other documentation category such as,
#' for example, show methods.
#'
#' @param object a \code{\link{DeeDeeExperiment}} object
#'
#' @return Returns NULL
NULL


#' @rdname DeeDeeExperiment-misc
#' @export
setMethod("show",
          signature = signature(object = "DeeDeeExperiment"),
          definition = function(object) {

            callNextMethod()
            cat(
              "dea(",length(object@dea), "): ",
              paste(names(object@dea), collapse = ", "), " \n",
              "fea(",length(object@fea), "): ",
              paste(names(object@fea), collapse = ", "),
              sep = ""
            )
          })


#' @rdname DeeDeeExperiment-misc
#' @param ... additional arguments passed to the summary method.
#' @export
setMethod("summary",
          signature = signature(object = "DeeDeeExperiment"),
          definition = function(object,
                                ...) {

            args <- list(...)
            FDR <- if (!is.null(args$FDR)) args$FDR else 0.05

            # dea summary
            dea <- dea_info(object)

            if (length(dea) > 0) {
              cat("DE Results Summary:\n")
              de_table <- data.frame(
                Contrast = names(dea),

                Up = sapply(names(dea), function(contrast) {
                  lfc_col <- paste0(contrast, "_log2FoldChange")
                  padj_col <- paste0(contrast, "_padj")
                  if (all(c(lfc_col, padj_col) %in% colnames(rowData(object)))) {
                    lfc <- rowData(object)[[lfc_col]]
                    padj <- rowData(object)[[padj_col]]
                    sum(lfc > 0 & padj < FDR, na.rm = TRUE)
                  } else {
                    NA_integer_
                  }
                }),

                Down = sapply(names(dea), function(contrast) {
                  lfc_col <- paste0(contrast, "_log2FoldChange")
                  padj_col <- paste0(contrast, "_padj")
                  if (all(c(lfc_col, padj_col) %in% colnames(rowData(object)))) {
                    lfc <- rowData(object)[[lfc_col]]
                    padj <- rowData(object)[[padj_col]]
                    sum(lfc < 0 & padj < FDR, na.rm = TRUE)
                  } else {
                    NA_integer_
                  }
                }),

                FDR_Cutoff = rep(FDR, length(dea))
              )
              print(de_table, row.names = FALSE)
              cat("\n")
            } else {
              cat("No DEA results stored.\n\n")
            }

            # fea summary

            fea <- fea_info(object)
            if (length(fea) > 0) {
              cat("FE Results Summary:\n")
              fea_table <- data.frame(
                FEA_Name = names(fea),
                Linked_DE = sapply(fea, function(object) ifelse(is.na(object$de_name),
                                                           ".", object$de_name)),
                Type = sapply(fea, function(object) object$fe_type),
                Terms = sapply(fea, function(object) nrow(object$original_object))
              )
              print(fea_table, row.names = FALSE)

            } else {
              cat("No FEA results stored.\n")
            }

          })


