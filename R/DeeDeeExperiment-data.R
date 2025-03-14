#' A sample `DESeqResults` object
#'
#' A sample `DESeqResults` object, generated with `DESeq2`
#'
#' @details This `DESeqResults` object contains the results of a Differential
#' Expression Analysis performed on data from the `macrophage` package, more
#' precisely contrasting the counts from naive macrophage to those associated
#' with IFNg.
#'
#' The code to create said object can be found in the folder `/inst/scripts` in
#' the DeeDeeExperiment package, the file is called `generate_data.R`.
#'
#' @return A sample `DESeqResults` object, generated with `DESeq2`
#'
#' @references Alasoo, et al. "Shared genetic effects on chromatin and gene
#' expression indicate a role for enhancer priming in immune response",
#' Nature Genetics, January 2018 doi: 10.1038/s41588-018-0046-7.
#'
#' @name IFNg_naive
#' @docType data
NULL


#' A sample `DESeqResults` object
#'
#' A sample `DESeqResults` object, generated with `DESeq2`
#'
#' @details This `DESeqResults` object contains the results of a Differential
#' Expression Analysis performed on data from the `macrophage` package, more
#' precisely contrasting the counts from macrophage associated with IFNg to
#' those associated with IFNg and Salmonella.
#'
#' The code to create said object can be found in the folder `/inst/scripts` in
#' the DeeDeeExperiment package, the file is called `generate_data.R`.
#'
#' @return A sample `DESeqResults` object, generated with `DESeq2`
#'
#' @references Alasoo, et al. "Shared genetic effects on chromatin and gene
#' expression indicate a role for enhancer priming in immune response",
#' Nature Genetics, January 2018 doi: 10.1038/s41588-018-0046-7.
#'
#' @name IFNg_both
#' @docType data
NULL


#' A sample `DESeqResults` object
#'
#' A sample `DESeqResults` object, generated with `DESeq2`
#'
#' @details This `DESeqResults` object contains the results of a Differential
#' Expression Analysis performed on data from the `macrophage` package, more
#' precisely contrasting the counts from naive macrophage to those associated
#' with Salmonella.
#'
#' The code to create said object can be found in the folder `/inst/scripts` in
#' the DeeDeeExperiment package, the file is called `generate_data.R`.
#'
#' @return A sample `DESeqResults` object, generated with `DESeq2`
#'
#' @references Alasoo, et al. "Shared genetic effects on chromatin and gene
#' expression indicate a role for enhancer priming in immune response",
#' Nature Genetics, January 2018 doi: 10.1038/s41588-018-0046-7.
#'
#' @name Salm_naive
#' @docType data
NULL



#' A sample `DESeqResults` object
#'
#' A sample `DESeqResults` object, generated with `DESeq2`
#'
#' @details This `DESeqResults` object contains the results of a Differential
#' Expression Analysis performed on data from the `macrophage` package, more
#' precisely contrasting the counts from macrophage associated with Salmonella
#' to those associated with IFNg and Salmonella.
#'
#' The code to create said object can be found in the folder `/inst/scripts` in
#' the DeeDeeExperiment package, the file is called `generate_data.R`.
#'
#' @return A sample `DESeqResults` object, generated with `DESeq2`
#'
#' @references Alasoo, et al. "Shared genetic effects on chromatin and gene
#' expression indicate a role for enhancer priming in immune response",
#' Nature Genetics, January 2018 doi: 10.1038/s41588-018-0046-7.
#'
#' @name Salm_both
#' @docType data
NULL




#' dd_list_original
#'
#' A list of `deedee_prepare`d DE results.
#'
#' @return A list of DE results
#'
#' @details documented creation in ... script
#'
#' @name dd_list_original
#' @docType data
NULL



#' de_named_list
#'
#' A named list of DE results, in their original format (from DESeq2, edgeR or limma)
#'
#' @return A named list of DE results, in their original format
#'
#' @details documented creation in the `create_dataset_example.R` script in the
#' `scripts` package folder
#'
#' @name de_named_list
#' @docType data
NULL


#' A sample `MArrayLM` object
#'
#' A sample `MArrayLM` object, generated with `limma`
#'
#' @details This `MArrayLM` object contains the results of a Differential
#' Expression Analysis performed on data from the `macrophage` package, more
#' precisely contrasting the counts from macrophage associated with different
#' stimulation conditions, including IFNg treatment, Salmonella infection, and
#' their combined effects.
#'
#' The code to create said object can be found in the folder `/inst/scripts` in
#' the DeeDeeExperiment package, the file is called `generate_data.R`.
#'
#' @return A sample `MArrayLM` object, generated with `limma`
#'
#' @references Alasoo, et al. "Shared genetic effects on chromatin and gene
#' expression indicate a role for enhancer priming in immune response",
#' Nature Genetics, January 2018 doi: 10.1038/s41588-018-0046-7.
#'
#' @name de_limma
#' @docType data
NULL

#' A sample `DGELRT` object
#'
#' A sample `DGELRT` object, generated with `edgeR`
#'
#' @details This `DGELRT` object contains the results of a Differential
#' Expression Analysis performed on data from the `macrophage` package, more
#' precisely contrasting the counts from macrophage associated with different
#' stimulation conditions, including IFNg treatment, Salmonella infection, and
#' their combined effects.
#'
#' The code to create said object can be found in the folder `/inst/scripts` in
#' the DeeDeeExperiment package, the file is called `generate_data.R`.
#'
#' @return A sample `DGELRT` object, generated with `edgeR`
#'
#' @references Alasoo, et al. "Shared genetic effects on chromatin and gene
#' expression indicate a role for enhancer priming in immune response",
#' Nature Genetics, January 2018 doi: 10.1038/s41588-018-0046-7.
#'
#' @name dge_lrt
#' @docType data
NULL


#' A sample `DGEExact` object
#'
#' A sample `DGEExact` object, generated with `edgeR`
#'
#' @details This `DGEExact` object contains the results of a Differential
#' Expression Analysis performed on data from the `macrophage` package, more
#' precisely contrasting the counts from macrophage associated with IFNg to
#' those associated with IFNg and Salmonella.
#'
#' The code to create said object can be found in the folder `/inst/scripts` in
#' the DeeDeeExperiment package, the file is called `generate_data.R`.
#'
#' @return A sample `DGEExact` object, generated with `edgeR`
#'
#' @references Alasoo, et al. "Shared genetic effects on chromatin and gene
#' expression indicate a role for enhancer priming in immune response",
#' Nature Genetics, January 2018 doi: 10.1038/s41588-018-0046-7.
#'
#' @name dge_exact_IFNg_both
#' @docType data
NULL

#' A sample `DGEExact` object
#'
#' A sample `DGEExact` object, generated with `edgeR`
#'
#' @details This `DGEExact` object contains the results of a Differential
#' Expression Analysis performed on data from the `macrophage` package, more
#' precisely contrasting the counts from naive macrophage to those associated
#' with IFNg.
#'
#' The code to create said object can be found in the folder `/inst/scripts` in
#' the DeeDeeExperiment package, the file is called `generate_data.R`.
#'
#' @return A sample `DGEExact` object, generated with `edgeR`
#'
#' @references Alasoo, et al. "Shared genetic effects on chromatin and gene
#' expression indicate a role for enhancer priming in immune response",
#' Nature Genetics, January 2018 doi: 10.1038/s41588-018-0046-7.
#'
#' @name dge_exact_IFNg_naive
#' @docType data
NULL

#' A sample `DGEExact` object
#'
#' A sample `DGEExact` object, generated with `edgeR`
#'
#' @details This `DGEExact` object contains the results of a Differential
#' Expression Analysis performed on data from the `macrophage` package, more
#' precisely contrasting the counts from naive macrophage to those associated
#' with Salmonella.
#'
#' The code to create said object can be found in the folder `/inst/scripts` in
#' the DeeDeeExperiment package, the file is called `generate_data.R`.
#'
#' @return A sample `DGEExact` object, generated with `edgeR`
#'
#' @references Alasoo, et al. "Shared genetic effects on chromatin and gene
#' expression indicate a role for enhancer priming in immune response",
#' Nature Genetics, January 2018 doi: 10.1038/s41588-018-0046-7.
#'
#' @name dge_exact_Salm_naive
#' @docType data
NULL

#' A sample `DGEExact` object
#'
#' A sample `DGEExact` object, generated with `edgeR`
#'
#' @details This `DGEExact` object contains the results of a Differential
#' Expression Analysis performed on data from the `macrophage` package, more
#' precisely contrasting the counts from macrophage associated with Salmonella
#' to those associated with IFNg and Salmonella.
#'
#' The code to create said object can be found in the folder `/inst/scripts` in
#' the DeeDeeExperiment package, the file is called `generate_data.R`.
#'
#' @return A sample `DGEExact` object, generated with `edgeR`
#'
#' @references Alasoo, et al. "Shared genetic effects on chromatin and gene
#' expression indicate a role for enhancer priming in immune response",
#' Nature Genetics, January 2018 doi: 10.1038/s41588-018-0046-7.
#'
#' @name dge_exact_Salm_both
#' @docType data
NULL
