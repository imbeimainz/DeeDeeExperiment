test_that("validity and so", {
  dde2 <- DeeDeeExperiment(
    se_macrophage_noassays,
    de_results = de_named_list
  )

  expect_true(validObject(dde2))

  rowData(dde2)[["ifng_vs_naive_log2FoldChange"]] <- NULL
  expect_error(validObject(dde2))

  dde3 <- DeeDeeExperiment(
    se_macrophage_noassays,
    de_results = de_named_list
  )

  dde3@dea <- list("foo", "bar")
  expect_error(validObject(dde3))

  dea_not_list <- dde3

  expect_error({
    deaInfo(dea_not_list) <- data.frame()
  })

  expect_error({
    feaInfo(dea_not_list) <- data.frame()
  })

  dde4 <- DeeDeeExperiment(
    se_macrophage_noassays,
    enrich_results = topGO_results_list
  )
  dde4@fea <- list("foo", "bar")
  expect_error(validObject(dde4))

  expect_error({
    addFea(dde4, fea = "meow", fea_tool = "fujitsu")
  })


  expect_error({
    feaInfo(dde4) <- list(list(
      de_name = "c1",
      fe_name = "c1",
      shaken_results = NULL,
      original_object = "meow",
      fe_tool = "topGO"
    ))
  })

  expect_error({
    feaInfo(dde4) <- list(list(
      de_name = "c1",
      fe_name = "c1",
      shaken_results = NULL,
      original_object = NULL,
      fe_tool = "chico"
    ))
  })
})
