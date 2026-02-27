test_that("misc", {
  ## misc ----------------------------------------------------------------------
  ### testing DAVID's output ---------------------------------------------------

  dde <- DeeDeeExperiment(
    se_macrophage_noassays,
    de_results = de_named_list
  )

  df_david <- read.delim(
    file = system.file("extdata", "david_output_chart_BPonly_ifng_vs_naive.txt",
      package = "DeeDeeExperiment"
    ),
    sep = "\t"
  )

  dde_w_david <- addFEA(dde,
    fea = df_david, fea_tool = "DAVID",
    de_name = "ifng_vs_naive"
  )

  expect_true("df_david" %in% getFEANames(dde_w_david))

  expect_error(
    .DeeDeefy_david(matrix())
  )

  expect_error(
    {
      fail_david <- df_david
      fail_david$Term <- NULL
      .DeeDeefy_david(fail_david)
    },
    "I could not find some of the usual"
  )

  ### testing available fea formats --------------------------------------------

  expect_s3_class(supported_fea_formats(), "data.frame")
  expect_equal(nrow(supported_fea_formats()), 8)

  ### testing helper functions

 expect_true(is.list(limma_list_for_dde(de_limma)))

 expect_true(is.list(muscat_list_for_dde(list(`stim-ctrl` = muscat_res))))

 expect_length(limma_list_for_dde(de_limma), 4)

 de_limma_list <- limma_list_for_dde(de_limma)

 expect_true(all(vapply(de_limma_list, inherits, logical(1), "data.frame")))

 expected_cols <- c("log2FoldChange", "pvalue", "padj")

 expect_true(all(vapply(de_limma_list, function(x) {
   all(expected_cols %in% colnames(x))
 }, logical(1))))

 expect_identical(attr(de_limma_list$IFNgNaive, "package"), "limma")
 expect_true(!is.null(attr(de_limma_list$IFNgNaive, "package_version")))

 expect_length(muscat_list_for_dde(list(`stim-ctrl` = muscat_res)), 3)

 expect_error(
   {
     muscat_list_for_dde(muscat_res)
   },
   "No valid muscat results found to convert"
 )

 dde_muscat_list <- muscat_list_for_dde(list(`stim-ctrl` = muscat_res))

 expect_true(all(vapply(dde_muscat_list, function(x) {
   all(expected_cols %in% colnames(x))
 }, logical(1))))


 expect_true(all(vapply(dde_muscat_list, function(x){
   !("table" %in% names(x))
 }, logical(1)
 )))

 expect_identical(attr(dde_muscat_list$`stim-ctrl_B cells`,
                       "package"), "muscat")
 expect_true(!is.null(attr(dde_muscat_list$`stim-ctrl_B cells`,
                           "package_version")))


 expect_error(limma_list_for_dde(fit = topTable(de_limma,
                                                coef = 2,
                                                number = 10,
                                                sort.by = "none")))

 expect_identical(attr(de_limma_list$IFNgNaive,
                       "package"), "limma")
 expect_true(!is.null(attr(de_limma_list$IFNgNaive,
                           "package_version")))


 ### testing export function

 expect_error(export_result_for_dde(dde_w_david,res_type = "guiga"))

 expect_error(export_result_for_dde(dde_w_david,res_format = "guiga"))

 expect_error(export_result_for_dde(dde_w_david, output_dir = character()))

 expect_error(export_result_for_dde(dde_w_david, file_name = character()))

 expect_error(export_result_for_dde(dde_w_david,
                                    output_dir = c("chico", "guiga")))

 expect_error(export_result_for_dde(dde_w_david,
                                    force = "yes"))

 outdir <- tempdir()

 path <- export_result_for_dde(dde_w_david,
                               res_type = "dea",
                               output_dir = outdir,
                               force = TRUE)
 expect_true(file.exists(path))
 expect_length(path, 1)

 paths <- export_result_for_dde(dde_w_david,
                                res_type = "both",
                                res_format = "original",
                                output_dir = outdir,
                                force = TRUE)

 expect_length(paths, 2)


})
