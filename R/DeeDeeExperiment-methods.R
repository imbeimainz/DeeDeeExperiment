#' @name DeeDeeExperiment-methods
#'
#' @title Methods for [DeeDeeExperiment] objects
#'
#' @aliases
#' dea_info
#' dea_info<-
#' dea_names
#' dea_rename
#' add_dea
#' remove_dea
#' dea
#' get_dea_list
#' add_scenario_info
#' fea_info
#' fea_info<-
#' fea_names
#' fea_rename
#' add_fea
#' remove_fea
#' fea
#' get_fea_list
#' link_dea_and_fea
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
#' get or remove, or match against (e.g., to fetch associated FEA results), or to which
#' additional context and information can be attached
#' @param verbose Logical, whether or not to display warnings. If TRUE, warnings/messages
#' will be displayed. If FALSE, the function runs silently
#' @param extra_rd A character vector of additional columns from rowData(x) to include. It
#' defaults to c("gene_id", "SYMBOL").
#' @param old_name A character vector of existing DEA names to be renamed in a `DeeDeeExperiment` object
#' @param new_name A character vector with new names to assign to existing DEA names in a
#' `DeeDeeExperiment` object. It must be the same length of `old_name`, and contains unique values that
#' don't overlap with existing DEA names.
#' @param fea A named list of Functional Enrichment results. Each element can be
#' either a data.frame (currently supports results from `topGO`, `enrichR`, `gProfiler`,
#' `fgsea`, `gsea`, `DAVID`, and output of `GeneTonic` shakers), or an `enrichResult`/`gseaResult`
#' objects (currently supports `clusterProfiler`)
#' @param fea_name Character value, specifying the name of the functional enrichment
#' result to add or remove
#' @param de_name A character string to explicitly specify the name of the de result this fea should be linked to.
#' If not provided, the function will attempt to match fea names to de results automatically.
#' @param fe_name A character string giving a name to the FE results.
#' @param remove_linked_fea A logical, specifying whether to remove or not the linked FEA when
#' a DEA results is removed
#' @param fea_tool A character string indicating the FEA tool used. It can take
#' any of the following values : "topGO", "clusterProfiler", "GeneTonic", "DAVID", "gsea",
#' "fgsea", "enrichr", "gProfiler". When not specified, it defaults to "auto" and
#' the tool is inferred automatically based on the input.
#' @param force A logical, indicating whether to overwrite results when introducing the same
#' results name. It defaults to FALSE.
#' @param format A character string, specifying the DEA/FEAs output format.
#' It takes either "minimal" to return only essential columns
#' (e.g. log2FC, p-value, adjusted p-value for DEAs,
#' or gs_id, gs_description, gs_pvalue, gs_genes... for FEAs), or "original" to return the full
#' result object. It defaults to "minimal"
#' @param info A character vector, containing contextual information about the
#' specified DE analysis. It defaults to NULL
#'
#' @return Return value varies depending on the individual methods, as described
#' below.
#'
#' @details
#'
#' DEAs
#'
#' * `dea_info` and `dea_info<-` are the methods to get and set the `dea` information as a
#' whole. These methods return `DeeDeeExperiment` objects.
#' * `dea_names` returns the names of the available DE contrasts in `DeeDeeExperiment` objects.
#' * `dea_rename` is the method to rename one or multiple DEAs stored in a `DeeDeeExperiment` object.
#' * `add_dea` and `remove_dea` are used to respectively add or remove DE-results
#' items. These methods also return `DeeDeeExperiment` objects, with updated
#' content in the `dea` slot.
#' * `dea` and `get_dea_list` retrieve the DEA information, as well as some extra rowData information and provide
#' this as a `DataFrame` object (for a specific analysis) or as a list, with one
#' element for each reported analysis.
#' * `add_scenario_info` is the method to add user defined contextual information for a specific DE analysis.
#' It allows users to attach free-text notes to a specific DEA results that stored in a
#' `DeeDeeExperiment` object. This information can include any other relevant information to help document
#' that DEA scenario. This context is stored in the `dea` slot under the name `scenario_info`,
#' which is not a default element in `dea`.
#'
#' FEAs
#'
#' * `fea_info` and `fea_info<-` are the methods to get and set the `fea` information as a
#' whole. These methods return `DeeDeeExperiment` objects.
#' * `fea_names` returns the names of the available enrichment results in `DeeDeeExperiment` objects.
#' * `fea_rename` is the method to rename one or multiple FEAs stored in a `DeeDeeExperiment``` object.
#' * `add_fea` and `remove_fea` are used to respectively add or remove functional
#' enrichment results items. These methods also return `DeeDeeExperiment` objects, with updated
#' content in the `fea` slot.
#' * `fea` is the method to retrieve FE results stored in a `DeeDeeExperiment` object
#' for a specific contrast, as a standardized format similar to the output of `GeneTonic` shakers.
#' * `get_fea_list` is the method that retrieves FEA results as a list. if the `dea_name` is indicated, the method
#' will return only FEAs linked to that `dea_name`, otherwise it returns all FEAs in the `fea` slot.
#' * `link_dea_and_fea` is the method that allows the user to manually link a FEA result to a specific DEA result
#'
#' * `show` is the method to nicely print out the information of a `DeeDeeExperiment`
#' object.
#' * `summary` is the method to print a summary of the available DE and FE results in a `DeeDeeExperiment`
#' object.
#'
#' @examples
#' data("de_named_list", package = "DeeDeeExperiment")
#' data("topGO_results_list", package = "DeeDeeExperiment")
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
#'
#' # add a new (set of) FE result(s)
#' dde_new <- add_fea(dde, fea = topGO_results)
#'
#' # removing FEAs
#' dde_rem <- remove_fea(dde_new, "ifng_vs_naive")
#'
#' # display available DEAs
#' dea_names(dde)
#'
#' # display available FEAs
#' fea_names(dde)
#'
#' # print a summary of the available DEAs and FEAs
#' summary(dde, FDR= 0.01)
#'
#' # rename DEA
#' dde_new <- dea_rename(dde_new, old_name = "salmonella_vs_naive",
#'                       new_name = "Salmo_vs_Naive_renamed")
#'
#' # assign DEA to FEA
#'
#' dde_new <- link_dea_and_fea(dde_new,
#'                              dea_name = "ifngsalmo_vs_naive",
#'                              fea_name = "ifngsalmo_vs_naive")
#'
NULL


# dea slot - get & set ---------------------------------------------------------

## dea_info --------------------------------------------------------------------

#' @rdname DeeDeeExperiment-methods
#' @export
setMethod("dea_info",
          signature = "DeeDeeExperiment",
          definition = function(x) {
            x@dea
          })

## dea_info <- -----------------------------------------------------------------

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

## dea_names -------------------------------------------------------------------

#' @rdname DeeDeeExperiment-methods
#' @export
setMethod("dea_names",
          signature = "DeeDeeExperiment",
          definition = function(x){
            names(dea_info(x))
          }
          )

## dea_rename ------------------------------------------------------------------

### TODO: add a setter for dea_names, in case one wants to rename the de res in
### dde? or a new method
#' @rdname DeeDeeExperiment-methods
#' @export
setMethod("dea_rename",
          signature = c("DeeDeeExperiment"),
          definition = function(x,
                                old_name,
                                new_name){

            # check uniqueness of new names, and that they don't overlap with
            # existing ones
            # dont forget to handle the naming of the current columns in the rowdata!!

            if (!is.character(old_name) ||  length(old_name) == 0) {
              stop("'old_name' must be a non empty character vector!")
            }

            if (!is.character(new_name) ||  length(new_name) == 0) {
              stop("'new_name' must be a non empty character vector!")
            }

            deas <- dea_info(x)
            current_names <- dea_names(x)
            if (length(current_names) == 0) {
              stop("No DEA results found")
            }

            if (length(old_name) != length(new_name)) {
              stop("'old_name' and 'new_name' must be the same length!")
            }

            matching_index <- match(old_name, current_names)

            if (any(is.na(matching_index))) {
              missing_names <- old_name[is.na(matching_index)]
              stop("The following DEA names where not found in dea slot:",
                   paste(missing_names, collapse = ", "))
            }

            if (anyDuplicated(new_name)) {
              stop("New names must be unique!")
            }

            overlapping_names <- intersect(new_name, current_names)
            if (length(overlapping_names) > 0) {
              stop("New names overlap with existing DEA names: ",
                   paste(overlapping_names, collapse = ", "))
            }


            names(deas)[matching_index] <- new_name
            x@dea <- deas

            rd <- rowData(x)
            rd_colnames <- colnames(rd)
            suffix <- c("_log2FoldChange","_pvalue","_padj")

            for (i in seq_along(old_name)) {
              old_prefix <- old_name[i]
              new_prefix <- new_name[i]

              for (j in suffix) {
                old_col <- paste0(old_prefix, j)
                new_col <- paste0(new_prefix, j)

                if (old_col %in% rd_colnames) {
                  colnames(rd)[which(rd_colnames == old_col)] <- new_col
                }
              }
            }

            rowData(x) <- rd

            # also rename in fea slot in there is a linked fea

            fea_names <- fea_names(x)

            for (fea in fea_names) {

              current_link <- fea_info(x)[[fea]][["de_name"]]
              if (!is.null(current_link) && current_link %in% old_name) {
                new_index <- match(current_link, old_name)
                updated_name <- new_name[new_index]

                fea_info(x)[[fea]][["de_name"]]  <- updated_name
              }

            }

            cli::cli_alert_success("Renamed DEA entries: {.val {old_name}} to {.val {new_name}}")

            validObject(x)
            x

          }
)


## add_dea ---------------------------------------------------------------------

#' @rdname DeeDeeExperiment-methods
#' @export
setMethod("add_dea",
          signature = c("DeeDeeExperiment"),
          definition = function(x,
                                dea,
                                force = FALSE) {

            # dde must be a DeeDeeExp
            # if (!is(x, "DeeDeeExperiment")) {
            #   stop("x must be DeeDeeExperiment object!")
            # }

            # dea must be named list
            if (is.null(names(dea))) {
              stop("All elements in dea list must have names!")
            }

            #check that names are all unique
            if (anyDuplicated(names(dea))) {
              stop("Names in dea must be unique!")
            }

            # unless force is TRUE

            new_names <- names(dea)
            existing_names <- names(dea_info(x))

            overlapping_names <- intersect(new_names, existing_names)

            if (length(overlapping_names) > 0 && !force) {
              stop("Names in 'dea' overlap with existing DEA results: ",
                   paste(overlapping_names, collapse = ", "),
                   ". Set force = TRUE to overwrite.")
            }


            # capture name inside the env where the func is called
            entry_name <- deparse(substitute(dea))

            # check and preocess dea
            dea <- .check_de_results(dea, entry_name)
            # names(dea)
            # names(dea_info(x))

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



## remove_dea ------------------------------------------------------------------

#' @rdname DeeDeeExperiment-methods
#' @export
setMethod("remove_dea",
          signature = c("DeeDeeExperiment"),
          definition = function(x,
                                dea_name,
                                remove_linked_fea = FALSE) {
            # x must be a DeeDeeExp

            if (!is.character(dea_name) || length(dea_name) == 0) {
              stop("'dea_name' must be a non empty character vector!")
            }

            # dea must be char vector
            deas <- names(dea_info(x))

            deas_to_remove <- intersect(dea_name, deas)

            # warning() if nothing to remove
            if (length(deas_to_remove) == 0){
              warning("Some elements in 'dea_name' were not found among DEA results.\n",
                      "Available results: ", paste(deas,collapse = ","))
            }


            for (i in deas_to_remove) {
              cols_to_remove <- c(paste0(i, c("_log2FoldChange", "_pvalue", "_padj")))
              rowData(x) <- rowData(x)[, !(colnames(rowData(x)) %in% cols_to_remove)]
              # update the de slot
              dea_info(x)[[i]] <- NULL

              if (remove_linked_fea) { ## fea is not removed!!!!!
                feas <- fea_info(x)
                removed_fea <- character()
                for (fea_name in names(feas)) {
                  if (!is.null(feas[[fea_name]][["de_name"]]) &&
                      feas[[fea_name]][["de_name"]] %in% deas_to_remove) {
                    removed_fea <- c(removed_fea, fea_name)
                    feas[[fea_name]] <- NULL
                    fea_info(x) <- feas
                  }
                }
                if (length(removed_fea) > 0) {
                  # message("The following linked FEA entries were removed: ",
                  #         paste(removed_fea, collapse = ", "))

                  cli::cli_alert_success("The following linked FEA entries were removed: {.val {paste(removed_fea, collapse = ', ')}} ")
                }
              }


              # unlink
              fea_info(x)[[i]][["de_name"]] <- NULL

            }
            removed_fea <- character()





            # here check some validity?
            validObject(x)

            # return the object
            return(x)
          }
)



## dea -------------------------------------------------------------------------

#' @rdname DeeDeeExperiment-methods
#' @export
setMethod("dea",
          signature = c("DeeDeeExperiment"),
          definition = function(x,
                                dea_name = NULL,
                                format = "minimal",
                                extra_rd = NULL,
                                verbose = TRUE) {

            if (!is.null(extra_rd) && !is.character(extra_rd)) {
              stop("'extra_rd' must be a character vector!")
            }

            deas <- dea_info(x)
            dea_names <- names(deas)

            if (!(format %in% c("minimal", "original"))) {
              stop("'format' not supported. Please use 'minimal' to return the ",
                   "essential columns, or 'original' to return the original object")
            }

            if (is.null(dea_name)) {
              if (length(dea_names) == 0) {
                stop("No DEA results found")
              }

              warning("'dea_name' was not specified. Returning the 1st entry: ",
                      dea_names[1])

              dea_name <- dea_names[1]
            }


            if (!is.character(dea_name) || length(dea_name) != 1) {
              stop("'dea_name' must be a single character string!")
            }

            if (!(dea_name %in% dea_names)) {
              stop("Could not find '",dea_name,"' among DEA results.\n",
                   "Available results: ", paste(dea_names,collapse = ","))
            }

            #
            if (format == "minimal") {
              rd_info <- paste0(dea_name,
                                c("_log2FoldChange", "_pvalue", "_padj"))

              extra_info <- extra_rd
              extra_cols <- extra_info[extra_info %in% colnames(rowData(x))] # drop if missing
              all_cols <- c(extra_cols,rd_info)

              overlap <- intersect(extra_info, rd_info)

              if (length(overlap) > 0) {
                stop("The following `extra_rd` are already part of the core `dea` columns and should not be repeated: ",
                     paste(overlap, collapse = ", "))
              }


              if (verbose && length(setdiff(extra_info, extra_cols)) > 0) {
                warning("Some 'extra_rd' are not available in rowData: ",
                        paste(setdiff(extra_info, extra_cols), collapse = ", "))
              }

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



              out <- rowData(x)[, all_cols]


            } else if (format == "original") {
              out <- dea_info(x)[[dea_name]][["original_object"]]

            }
            return(out)

          }
)



## get_dea_list ----------------------------------------------------------------

#' @rdname DeeDeeExperiment-methods
#' @export
setMethod("get_dea_list",
          signature = c("DeeDeeExperiment"),
          definition = function(x,
                                format = "minimal",
                                verbose = TRUE) {

            if (!(format %in% c("minimal", "original"))) {
              stop("'format' not supported. Please use 'minimal' to return the ",
                   "essential columns, or 'original' to return the original object")
            }

            deas <- dea_info(x)
            dea_names <- names(deas)

            dea_list <- list()
            affected_deas <- character()

            for (i in dea_names) {
              # dea_list[[i]] <- as.data.frame(dea(x, i, verbose))
              dea_list[[i]] <- as.data.frame(
                dea(x, dea_name = i, format = format, verbose = verbose))

              if (format == "minimal") {
                # remove the first two columns
                # dea_list[[i]] <- dea_list[[i]][, -c(1,2)]
                colnames(dea_list[[i]]) <- c("log2FoldChange", "pvalue", "padj")

                # maybe check for rowname mismatches potential gene version issue?
                # maybe not interesting to print back all missmatches in casee all rownames
                # dont match
                # rownames_x <- rownames(rowData(x))
                # rownames_y <- rownames(deas[[i]][["original_object"]])
                #
                # mismatched_rows <- sum(!rownames_x %in% rownames_y)
                #
                # if (mismatched_rows > 0) {
                #   affected_deas <- c(affected_deas, i)
                # }
              }

            }

            # not needed to check mismatch since the warnings will be triggered from dea

            return(dea_list)
          }
)


## add_scenario_info -----------------------------------------------------------

#' @rdname DeeDeeExperiment-methods
#' @export
setMethod("add_scenario_info",
          signature = c("DeeDeeExperiment"),
          definition = function(x,
                                dea_name,
                                info = NULL,
                                force = FALSE){
            dea_names <- dea_names(x)
            existing_info <- dea_info(x)[[dea_name]][["scenario_info"]]

            # checks on dea_name
            if (!is.character(dea_name) || length(dea_name) != 1) {
              stop("'dea_name' must be a single character string!")
            }

            # checks on info
            if (!is.null(info) && !is.character(info)) {
              stop("'info' must be a character vector (e.g. one or more strings)")
            }

            if (!(dea_name %in% dea_names)) {
              stop("'dea_name'", dea_name,"not found among DEA results.\n",
                   "Available results: ", paste(dea_names,collapse = ","))
            }

            if (!is.null(existing_info) && !force) {
              stop("Existing scenario_info for '", dea_name, "' already exists.",
                   "Set force = TRUE to overwrite")
            }

            # when both info and existing_info are null -> do nothing

            dea_info(x)[[dea_name]][["scenario_info"]] <- info

            # update object
            validObject(x)
            x

          }
)




# fea slot - get & set ---------------------------------------------------------

## fea_info --------------------------------------------------------------------

#' @rdname DeeDeeExperiment-methods
#' @export
setMethod("fea_info",
          signature = "DeeDeeExperiment",
          definition = function(x) {
            x@fea
          })

## fea_info <- -----------------------------------------------------------------

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

## fea_names -------------------------------------------------------------------

#' @rdname DeeDeeExperiment-methods
#' @export
setMethod("fea_names",
          signature = "DeeDeeExperiment",
          definition = function(x){
            names(fea_info(x))
          }
)

## fea_rename ------------------------------------------------------------------

#' @rdname DeeDeeExperiment-methods
#' @export
setMethod("fea_rename",
          signature = c("DeeDeeExperiment"),
          definition = function(x,
                                old_name,
                                new_name){

            # check uniqueness of new names, and that they don't overlap with
            # existing ones

            if (!is.character(old_name) ||  length(old_name) == 0) {
              stop("'old_name' must be a non empty character vector!")
            }

            if (!is.character(new_name) ||  length(new_name) == 0) {
              stop("'new_name' must be a non empty character vector!")
            }

            feas <- fea_info(x)
            current_names <- fea_names(x)
            if (length(current_names) == 0) {
              stop("No FEA results found")
            }

            if (length(old_name) != length(new_name)) {
              stop("'old_name' and 'new_name' must be the same length!")
            }

            matching_index <- match(old_name, current_names)

            if(any(is.na(matching_index))) {
              missing_names <- old_name[is.na(matching_index)]
              stop("The following FEA names where not found in fea slot:",
                   paste(missing_names, collapse = ", "))
            }

            if(anyDuplicated(new_name)) {
              stop("New names must be unique!")
            }

            overlapping_names <- intersect(new_name, current_names)
            if(length(overlapping_names) > 0) {
              stop("New names overlap with existing FEA names: ",
                   paste(overlapping_names, collapse = ", "))
            }


            names(feas)[matching_index] <- new_name
            x@fea <- feas

            cli::cli_alert_success("Renamed FEA entries: {.val {old_name}} to {.val {new_name}}")

            validObject(x)
            x

          }
)

## add_fea ---------------------------------------------------------------------

#' @rdname DeeDeeExperiment-methods
#' @export
setMethod(
  "add_fea",
  signature = c("DeeDeeExperiment"),
  definition = function(x,
                        fea,
                        de_name = NA_character_,
                        fe_name = NULL,
                        fea_tool = "auto",
                        force = FALSE) {
    # x must be a DeeDeeExperiment
    # if (!is(x, "DeeDeeExperiment")) {
    #   stop("x must be DeeDeeExperiment object!")
    # }

    if (!is.character(de_name) || length(de_name) != 1) {
      stop("'de_name' must be a single character string or NA_character_")
    } # should it be a vector of different de_name???


    # allowed fea_tools

    allowed_fea_tools <- c("auto", "topGO", "clusterProfiler", "GeneTonic",
                           "DAVID", "gsea", "fgsea", "enrichr", "gProfiler")

    if (!is.character(fea_tool)) {
      stop("fea_tool should be a character vector!")
    }

    if (!all(fea_tool %in% allowed_fea_tools)) {
      stop("fea_tool should be one of the following: topGO, clusterProfiler, GeneTonic,
                          DAVID, gsea, fgsea, enrichr, gProfiler")
    }

    # capture name inside the env where the func is called
    entry_name <- deparse(substitute(fea))

    # check and preocess fea
    fea_list <- .check_enrich_results(fea, entry_name)

    # fea must be named list
    if (is.null(names(fea_list))) {
      stop("All elements in 'fea' list must have names!")
    }

    #check that names are all unique
    if (anyDuplicated(names(fea))) {
      stop("Names in dea must be unique!")
    }

    # check that names are all unique, and do not overlap with the existing ones
    # unless force is TRUE

    new_names <- names(fea_list)
    existing_names <- names(fea_info(x))

    overlapping_names <- intersect(new_names, existing_names)

    if (length(overlapping_names) > 0 && !force) {
      stop(
        "Names in 'fea' overlap with existing FEA results: ",
        paste(overlapping_names, collapse = ", "),
        ". Set force = TRUE to overwrite."
      )
    }


    # get existing results in the fea slot
    fea_contrasts <- fea_info(x)

    for (fe in names(fea_list)) {
      res_enrich <- fea_list[[fe]]
      if (!is.null(dea_info(x)) && length(dea_info(x)) > 0) {

        if (!is.na(de_name)) {
          if (de_name %in% names(dea_info(x))) {
            de_res_name <- de_name
          } else {
            warning("Provided 'de_name' ('", de_name,"') not found among DE results. Coercing into NA_character_")
            de_res_name <- NA_character_
          }
        } else {
          matched_name <- .match_fe_to_de(fe, names(dea_info(x)))
          if (!is.na(matched_name) && matched_name %in% names(dea_info(x))) {
            de_res_name <- matched_name
            if (fe != matched_name) {
              ### if the name is exactly the same do we need a msg or it s obvious???
              #message("FEA '", fe, "' matched to DE contrast '", matched_name,"'")
              if (verbose) {
                cli::cli_alert_info("FEA {.val {fe}} matched to DE contrast {.val {matched_name}}")
              }

            } else{
              # in case of the same name
              #message("FEA '", fe, "' matched **directly** to DE contrast '", matched_name,"'")
              if (verbose) {
                cli::cli_alert_info("FEA {.val {fe}} matched directly to DE contrast {.val {matched_name}}")
              }

            }
          } else {
            de_res_name <- NA_character_
            warning("Could not match FEA '", fe, "' to any DE contrast.\n",
              "Available DE results: ", paste(names(dea_info(x)), collapse = ", "), "\n",
              "Consider naming your enrich_results starting with one of the following prefixes:",
              " 'topGO_', 'clusterProfiler_','GeneTonic_', 'DAVID_','gsea_', 'fgsea_', 'enrichr_', 'gPro_',",
              "followed by the contrast name"
            )
          }

        }

      } else {
        de_res_name <- NA_character_
        warning("Could not match FEA '", fe, "' to a DE contrast because no DE results were provided.\n")
      }


      n_fea <- length(fea_list)

      if (length(fea_tool) == 1 && fea_tool == "auto") {
        fea_tool_vec <- rep("auto", n_fea)
      } else if (length(fea_tool) == 1 && fea_tool %in% allowed_fea_tools) {
        fea_tool_vec <- rep(fea_tool, n_fea)
      } else if (length(fea_tool) == n_fea && all(fea_tool %in% allowed_fea_tools)) {
        fea_tool_vec <- fea_tool
      } else {
        stop("'fea_tool' must be either: A single valid tool name (e.g.",
             paste(allowed_fea_tools, collapse = ", "), "or a character vector of length ",
             n_fea, " with tool names matching the FEA elements in order")
      }

      names(fea_tool_vec) <- names(fea_list)

      this_tool <- fea_tool_vec[[fe]]

      if (this_tool == "auto") {
        # auto detect
        fe_tool <- .detect_fea_tool(res_enrich)
      } else {
        fe_tool <- this_tool
      }
      res_enrich_shaken <- NULL # default

      if (fe_tool == "topGO") {
        # shake using shake_topGOtableResult
        res_enrich_shaken <- .DeeDeefy_topGOtableResult(res_enrich)

        } else if (fe_tool == "clusterProfiler") {
          #shake using shake_enrichResult
          res_enrich_shaken <- .DeeDeefy_enrichResult(res_enrich)

        } else if (fe_tool == "GeneTonic") {
          # shake based on specific columns or return original object
          res_enrich_shaken <- res_enrich # input already shaken

        } else if (fe_tool == "DAVID") {
          # we are not taking the output of the file!!  so we cannot
          # use genetonic shakers!!
          # create shakers for that
          res_enrich_shaken <- .DeeDeefy_david(res_enrich)

        } else if (fe_tool == "fgsea") {
          res_enrich_shaken <- .DeeDeefy_fgseaResult(res_enrich)

        } else if (fe_tool == "gsea") {
          res_enrich_shaken <- .DeeDeefy_gsenrichResult(res_enrich)

        } else if (fe_tool == "enrichr") {
          res_enrich_shaken <- .DeeDeefy_enrichr(res_enrich)

        } else if (fe_tool == "gProfiler") {
          res_enrich_shaken <- .DeeDeefy_gprofiler(res_enrich)
      }

      if (is.null(res_enrich_shaken)) {
        # message(
        #   "No shaking method available for this functional enrichment results.",
        #   " Returning only the original object."
        # )

        cli::cli_alert_info("No shaking method available for this functional enrichment results.
                            Returning only the original object.")
      }

      fea_contrast <- list(
        de_name = de_res_name,
        # links to de result
        fe_name = if (!is.null(fe_name)) fe_name else fe,
        shaken_results = res_enrich_shaken ,
        # return shaken results for later use in GeneTonic
        original_object = res_enrich,
        fe_tool = fe_tool
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



## remove_fea ------------------------------------------------------------------

#' @rdname DeeDeeExperiment-methods
#' @export
setMethod("remove_fea",
          signature = c("DeeDeeExperiment"),
          definition = function(x, fea_name) {

            # x must be a DeeDeeExp
            # if(!is(x,"DeeDeeExperiment")) {
            #   stop("x must be a DeeDeeExperiment object!")
            # }

            if (!is.character(fea_name) || length(fea_name) == 0) {
              stop("'fea_name' must be a non empty character vector!")
            }

            feas <- fea_names(x)

              if (!all(fea_name %in% feas)) {
                stop("Some elements in 'fea_name' were not found among FEA results.\n",
                     "Available results: ", paste(feas,collapse = ","))
              }

            feas_to_remove <- intersect(fea_name, feas)

            # warning() if nothing to remove
            if(length(feas_to_remove) == 0){
              warning("No matching fea entries found to remove.")
            }

            for (i in feas_to_remove) {
              # update the fea slot
              fea_info(x)[[i]] <- NULL
            }

            # here check some validity?
            validObject(x)

            # return the object
            return(x)
          }
)


## fea -------------------------------------------------------------------------

#' @rdname DeeDeeExperiment-methods
#' @export
setMethod("fea",
          signature = c("DeeDeeExperiment"),
          definition = function(x,
                                fea_name = NULL,
                                format = "minimal") {

            # get returns shaken table by default for a specific contrast
            # for now the only case where we won't have shaken results if the user
            # introduces fea that is not generate with {topGO,clusterProfiler,...}
            #we can handle the other types later

            # or should we give the user the freedom to choose which table to fetch??? using another arg

            # check
            # x must be a DeeDeeExperiment
            # if (!is(x, "DeeDeeExperiment")) {
            #   stop("x must be DeeDeeExperiment object!")
            # }

            if (!(format %in% c("minimal", "original"))) {
              stop("'format' not supported. Please use 'minimal' to return the ",
                   "essential columns, or 'original' to return the original object")
            }

            fea_names <- fea_names(x)

            if (is.null(fea_name)) {
              if (length(fea_names) == 0) {
                stop("No FEA results found")
              }

              warning("'fea_name' was not specified. Returning the 1st entry: ",
                      fea_names[1])

              fea_name <- fea_names[1]
            }

            if (!is.character(fea_name) || length(fea_name) != 1) {
              stop("'fea_name' must be a single character string!")
            }

            if (!(fea_name %in% fea_names)) {
              stop("Could not find '",fea_name,"' among FEA results.\n",
                   "Available results: ", paste(fea_names,collapse = ","))
            }


            if (format == "minimal") {

              fea <- fea_info(x)[[fea_name]][["shaken_results"]]

              if (is.null(fea)) {
                warning("No shaken results available for '", fea_name,
                        "'. Returning original enrichment results instead.")
                fea <- fea_info(x)[[fea_name]]$original_object
              }


            } else if (format == "original") {
              fea <- fea_info(x)[[fea_name]][["original_object"]]


            }

            return(fea)
          }

)



## get_fea_list ----------------------------------------------------------------

#' @rdname DeeDeeExperiment-methods
#' @export
setMethod("get_fea_list",
          signature = c("DeeDeeExperiment"),
          definition = function(x,
                                dea_name = NULL,
                                format = "minimal") {

    if (!(format %in% c("minimal", "original"))) {
      stop(
        "'format' not supported. Please use 'minimal' to return the ",
        "essential columns, or 'original' to return the original object"
      )
    }

    all_fea_names <- fea_names(x)

    if (length(all_fea_names) == 0) {
      stop("No FEA results found")
    }

    matched_feas <- list()

    if (!is.null(dea_name) && (!is.character(dea_name) || length(dea_name) != 1)) {
      stop("'dea_name' must be a single character string")
    }

    for (i in all_fea_names) {
      # catch the corresponding dea
      de_name <- fea_info(x)[[i]][["de_name"]]

      #if dea_name is not indicated, return all feas
      # otherwise return only the specific feas associated with that dea_name


      if (is.null(dea_name) || (!is.na(de_name) && de_name == dea_name)) {

        if (format == "minimal") {
          fe_res <- fea_info(x)[[i]][["shaken_results"]]

          if (!is.null(fe_res)) {
            matched_feas[[i]] <- fe_res

          } else {
            warning(
              "No shaken results available for '",
              i,
              "'. Returning original enrichment results instead."
            )

            matched_feas[[i]] <- fea_info(x)[[i]][["original_object"]]

          }
        } else if (format == "original") {
          matched_feas[[i]] <- fea_info(x)[[i]][["original_object"]]
        }
      }
    }

    if (length(matched_feas) == 0) {
      if (!is.null(dea_name)) {
      warning("No FEA results found for '", dea_name, "'")
      } else {
      warning("No FEA results returned")
      }
    }

    return(matched_feas)

  }
)



## link_dea_and_fea -----------------------------------------------------------

#' @rdname DeeDeeExperiment-methods
#' @export
setMethod("link_dea_and_fea",
          signature = c("DeeDeeExperiment"),
          definition = function(x,
                                dea_name,
                                fea_name,
                                force = FALSE) {

            # check fea_name & dea_name are character
            if (!is.character(dea_name) ||  length(dea_name) == 0) {
              stop("'dea_name' must be a single character string!")
            }

            if (!is.character(fea_name) ||  length(fea_name) == 0) {
              stop("'fea_name' must be a non empty character vector!")
            }


            dea_names <- dea_names(x)
            fea_names <- fea_names(x)


            if (!(dea_name %in% dea_names)) {
              stop("DEA result: '",dea_name,"' not found")
            }

            for (fea in fea_name) {
              if (!(fea %in% fea_names)) {
                stop("FEA result: '",fea,"' not found")
              }

              # do we have existing link?
              current_de_name <- fea_info(x)[[fea]][["de_name"]]

              if (!is.null(current_de_name) && !is.na(current_de_name) && current_de_name != dea_name) {
                if(!force) {
                  stop("FEA '", fea, "' is already linked to DEA '", current_de_name,
                       "'. Use `force = TRUE` to overwrite")
                } else {
                  warning("FEA '", fea, "' was linked to DEA '", current_de_name,
                          ". Now linked to '", dea_name, "'")
                }

              }

              # assign
              message("Assigning DEA '", dea_name, "' to FEA '", fea, "'")

              fea_info(x)[[fea]][["de_name"]] <- dea_name

            }

            validObject(x)
            x
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
#' for example, show and summary methods.
#'
#' @param object a \code{\link{DeeDeeExperiment}} object
#' @param ... additional argument passed to the summary method.
#' Currently supports:
#' \itemize{
#'   \item `FDR`: Numeric, sets the significance threshold for subsetting differentially
#'   expressed genes based on adjusted p-values. Defaults to 0.05
#'   \item `show_scenario_info`: Logical; if \code{TRUE}, displays the associated
#'   scenario info for each DE contrast, if available.
#'   }
#' @return Returns NULL
NULL


## show ------------------------------------------------------------------------

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


## summary ---------------------------------------------------------------------

#' @rdname DeeDeeExperiment-misc
#' @export
setMethod("summary",
          signature = signature(object = "DeeDeeExperiment"),
          definition = function(object, ...) {
            # using ellipsis because we can't change the summary method
            args <- list(...)
            FDR <- if (!is.null(args$FDR))
              args$FDR
            else 0.05
            show_scenario_info <- isTRUE(args$show_scenario_info)
            # dea summary
            dea <- dea_info(object)

            if (length(dea) > 0) {
              cat("DE Results Summary:\n")
              de_table <- data.frame(

                DEA_name = names(dea),

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

                FDR = rep(FDR, length(dea))
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
                Linked_DE = sapply(fea, function(object) {
                  if (!is.null(object$de_name) && !is.na(object$de_name)) {
                    object$de_name
                  } else {
                    "."
                  }
                  }),
                FE_Type = sapply(fea, function(object) {
                  if (!is.null(object$fe_tool)) {
                    object$fe_tool
                  } else {
                    "Not Specified"
                  }
                  }),
                Term_Number = sapply(fea, function(object) {
                  if (!is.null(object$original_object)) {
                    NROW(object$original_object)
                    } else {
                      NA_integer_
                      }
                  })
                )
              print(fea_table, row.names = FALSE)

              } else {
                cat("No FEA results stored.\n")
              }

            # scenario info (only if show_scenario_info is TRUE)
            if (show_scenario_info && length(dea) > 0) {
              cat("\nScenario Info:\n")
              missing <- character()
              for (de_name in names(dea)) {
                scenario_info <- dea[[de_name]][["scenario_info"]]
                if (!is.null(scenario_info)) {
                  cat(" -", de_name, ":\n")

                  wrapped_txt <- stringr::str_wrap(scenario_info,
                                          width = 80,
                                          indent = 1,
                                          exdent = 2)

                  cat(paste(wrapped_txt, "\n"), "\n")
                } else {
                  missing <- c(missing, de_name)
                }
              }

              if (length(missing) > 0) {
                cat("\nNo scenario info for:", paste(missing, collapse = ", "), "\n")
              }

              cat("\n")
            }

    }
)
