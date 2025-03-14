test_that("creating", {
  dde <- DeeDeeExperiment(
    se_macrophage_noassays,
    de_results = de_named_list
  )

  print(dde)

  expect_s4_class(dde, "DeeDeeExperiment")

  dde_only_de <- DeeDeeExperiment(
    de_results = de_named_list
  )
  expect_s4_class(dde_only_de, "DeeDeeExperiment")

  dde_nodd <- DeeDeeExperiment(
    se = se_macrophage_noassays,
  )
  expect_s4_class(dde_nodd, "DeeDeeExperiment")


  expect_s4_class(
    get_dea_df(dde, "ifng_vs_naive"), "DataFrame"
  )

  expect_error(
    get_dea_df(dde, "wrong_name")
  )

  dde_gone_wrong <- dde
  rowData(dde_gone_wrong)[["ifng_vs_naive_log2FoldChange"]] <- NULL
  expect_error(
    get_dea_df(dde_gone_wrong, "ifng_vs_naive")
  )

  expect_error(
    DeeDeeExperiment(
      rowData(se_macrophage_noassays),
      de_results = de_named_list
    )
  )

  expect_error(
    DeeDeeExperiment(
      assay(se_macrophage),
      de_results = de_named_list
    )
  )

  expect_error(
    DeeDeeExperiment()
  )

  salmo_both <- de_named_list$salmo_both

  dde_one <- DeeDeeExperiment(se = se_macrophage_noassays,
                              de_results = salmo_both)

  expect_false(is.list(salmo_both))

  expect_s4_class(salmo_both, "DESeqResults")

  expect_type(dea(dde_one), "list")

  expect_length(dea(dde_one), 1)

  expect_true("salmo_both" == names(dea(dde_one)))

  dea1 <- de_limma
  de_res_list <- list(de_deseq = salmo_both,
                      dge_lrt =dea1)

  dde_list <- DeeDeeExperiment(de_results = de_res_list)

  expect_warning(get_dea_list(dde_list))

  expect_warning(get_dea_df(dde_list, dea_name = "dge_lrt"))

}

)


test_that("adding and removing", {
  dde <- DeeDeeExperiment(
    se_macrophage_noassays,
    de_results = de_named_list
  )

  new_del <- list(
    ifng2 = de_named_list$ifng_vs_naive,
    ifngsalmo2 = de_named_list$ifngsalmo_vs_naive
  )
  # add a new (set of) DE result(s)
  dde_new <- add_dea(dde, new_del)
  expect_s4_class(dde_new, "DeeDeeExperiment")
  expect_equal(length(dea(dde)), 4)
  expect_equal(length(dea(dde_new)), 6)

  dde_removed <- remove_dea(dde, "ifngsalmo_vs_naive")
  expect_s4_class(dde_removed, "DeeDeeExperiment")
  expect_equal(length(dea(dde_removed)), 3)
})


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

  # invalid replacements
  ## actually, can enable it is an empty list (i.e. no DE (yet) inserted)
  # dde3@dea <- list()
  # expect_error(validObject(dde3))

  dde3@dea <- list("foo", "bar")
  expect_error(validObject(dde3))
})
