test_that("creating", {
  expect_error(DeeDeeExperiment())

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

  expect_type(dea_names(dde_only_de), "character")

  expect_equal(
    dea_names(dde_only_de),
    c("ifng_vs_naive", "ifngsalmo_vs_naive", "salmonella_vs_naive", "salmo_both")
  )

  dde_nodd <- DeeDeeExperiment(
    se = se_macrophage_noassays,
  )
  expect_s4_class(dde_nodd, "DeeDeeExperiment")


  expect_s4_class(
    dea(dde, "ifng_vs_naive"), "DataFrame"
  )

  expect_error(
    dea(dde, "wrong_name")
  )

  dde_gone_wrong <- dde
  rowData(dde_gone_wrong)[["ifng_vs_naive_log2FoldChange"]] <- NULL
  expect_error(
    dea(dde_gone_wrong, "ifng_vs_naive")
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

  salmo_both <- de_named_list$salmo_both

  dde_one <- DeeDeeExperiment(se = se_macrophage_noassays,
                              de_results = salmo_both)

  expect_false(is.list(salmo_both))

  expect_s4_class(salmo_both, "DESeqResults")

  expect_type(dea_info(dde_one), "list")

  expect_length(dea_info(dde_one), 1)

  expect_true("salmo_both" == names(dea_info(dde_one)))

  de_results_mismatch <- list(
    contrast = de_named_list$ifng_vs_naive
  )
  rownames(de_results_mismatch$contrast) <- paste0(
    "gene", 101:(100 + nrow(de_results_mismatch$contrast)))

  expect_warning(
    DeeDeeExperiment(
    se_macrophage_noassays,
    de_results = de_results_mismatch
  ))

  broken_de_res <- salmo_both
  rownames(broken_de_res) <- NULL

  expect_error(DeeDeeExperiment(de_results = list(salmo_both = salmo_both,
                                                  broken_salmo_both = broken_de_res)))


  dea1 <- de_limma

  de_res_list <- list(de_deseq = salmo_both,
                      dge_lrt = dea1)

  dde_list <- DeeDeeExperiment(de_results = de_res_list)

  expect_warning({get_dea_list(dde_list)})

  expect_warning({dea(dde_list, dea_name = "dge_lrt")})

  expect_warning({DeeDeeExperiment(se =se_macrophage_noassays,
                                   de_results = dea1)})

  dea2 <- dge_exact_IFNg_both
  expect_warning({DeeDeeExperiment(se =se_macrophage_noassays,
                                   de_results = dea2)})

  fea1 <- topGO_results

  expect_warning({DeeDeeExperiment(enrich_results = fea1)})

  contrast1 <- topGO_results$ifng_vs_naive

  expect_warning({DeeDeeExperiment(de_results = de_named_list,
                                  enrich_results = contrast1)})



  expect_message({DeeDeeExperiment(
    de_results = list(ifng_vs_naive = de_named_list$ifng_vs_naive,
                      salmonella_vs_naive = de_named_list$salmonella_vs_naive),
    enrich_results = list(
      topGO_ifng_vs_naive = topGO_results$ifng_vs_naive,
      salmonella_vs_naive = topGO_results$salmonella_vs_naive))}
    )

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
  expect_equal(length(dea_info(dde)), 4)
  expect_equal(length(dea_info(dde_new)), 6)

  expect_error({dde_add <- add_dea(x = "sthg else", dea = new_del)})

  expect_error({dde_add <- add_dea(x = de, dea = list(de_named_list$ifng_vs_naive))})

  expect_error({dde_add <- add_dea(x = dde, dea = list(
    ifng2 = de_named_list$ifng_vs_naive,
    ifng2 = de_named_list$ifngsalmo_vs_naive
  ))})

  dde_removed <- remove_dea(dde, "ifngsalmo_vs_naive")
  expect_s4_class(dde_removed, "DeeDeeExperiment")
  expect_equal(length(dea_info(dde_removed)), 3)

  topGO_Salm_naive <- topGO_results$salmonella_vs_naive
  topGO_IFNg_naive <- topGO_results$ifng_vs_naive


  dde2 <- add_fea(dde, fea_res = list(topGO_Salm_naive = topGO_Salm_naive,
                            topGO_IFNg_naive = topGO_IFNg_naive),
                  fea_type = "topGO")

  expect_s4_class(dde2, "DeeDeeExperiment")
  expect_equal(length(fea_info(dde2)), 2)

  dde2 <- remove_fea(dde2, "topGO_IFNg_naive")
  expect_s4_class(dde2, "DeeDeeExperiment")
  expect_equal(length(fea_info(dde2)), 1)


  expect_error(remove_fea(dde2,"IFNgVSnaive"))

  expect_s3_class(fea(dde2, "topGO_Salm_naive"), "data.frame")

  dde3 <- DeeDeeExperiment(se = se_macrophage_noassays)
  expect_warning(add_fea(dde3, fea_res = list(topGO_Salm_naive = topGO_Salm_naive,
                                      topGO_IFNg_naive = topGO_IFNg_naive)))



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

  dde3@dea <- list("foo", "bar")
  expect_error(validObject(dde3))

  dea_not_list <- dde3

  expect_error({
    dea_info(dea_not_list) <- data.frame()
  })

  expect_error({
    fea_info(dea_not_list) <- data.frame()
  })

  dde4 <- DeeDeeExperiment(
    se_macrophage_noassays,
    enrich_results = topGO_results
  )
  dde4@fea <- list("foo", "bar")
  expect_error(validObject(dde4))

  expect_error({
    add_fea(dde4, fea = "meow", fea_type = "fujitsu")
  })


  expect_error({
    fea_info(dde4) <- list(list(de_name = "c1",
                           fe_name = "c1",
                           shaken_results = NULL,
                           original_object = "meow",
                           fe_tool = "topGO"))
  }
  )

  expect_error({
    fea_info(dde4) <- list(list(de_name = "c1",
                                fe_name = "c1",
                                shaken_results = NULL,
                                original_object = NULL,
                                fe_tool = "chico"))
  }
  )


})
