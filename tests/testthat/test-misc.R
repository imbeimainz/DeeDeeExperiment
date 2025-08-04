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

  dde_w_david <- addFea(dde,
    fea = df_david, fea_tool = "DAVID",
    de_name = "ifng_vs_naive"
  )

  expect_true("df_david" %in% feaNames(dde_w_david))

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
})
