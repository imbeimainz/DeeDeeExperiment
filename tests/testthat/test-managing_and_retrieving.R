test_that("managing and retrieving", {
  ## adding DEA ----------------------------------------------------------------
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

  expect_error({
    add_dea(x = dde, dea = list(de_named_list$ifng_vs_naive))
  })

  expect_error({
    dde_add <- add_dea(x = "sthg else", dea = new_del)
  })

  expect_error({
    dde_add <- add_dea(x = dde, dea = list(
      ifng2 = de_named_list$ifng_vs_naive,
      ifng2 = de_named_list$ifngsalmo_vs_naive
    ))
  })

  dde_edgeR <- add_dea(dde, dea = list(DGEExact_IFNg_both = dge_exact_IFNg_both))
  expect_s4_class(dde_edgeR, "DeeDeeExperiment")
  expect_equal(length(dea_info(dde_edgeR)), 5)

  dde_limma <- add_dea(dde, dea = list(de_limma = de_limma))
  expect_s4_class(dde_limma, "DeeDeeExperiment")
  expect_equal(length(dea_info(dde_limma)), 5)

  expect_no_error(add_dea(dde_limma,
    dea = list(de_limma = de_limma),
    force = TRUE
  ))

  genes <- rownames(de_limma)[1:20]

  de_custom <- data.frame(
    p_val = rep(0.5, 20),
    adj_pvalue = rep(0.5, 20),
    gene = genes
  )
  expect_error(add_dea(dde, dea = list(de_custom = de_custom)))

  expect_error(add_dea(dde, "ifng_vs_naive"))

  dde_overlap <- DeeDeeExperiment(se_macrophage_noassays,
    de_results = de_named_list,
    enrich_results = topGO_results_list
  )

  expect_error(add_dea(dde_overlap,
    dea = list(ifng_vs_naive = de_named_list$ifng_vs_naive)
  ))

  ## removing DEA --------------------------------------------------------------

  dde_removed <- remove_dea(dde, "ifngsalmo_vs_naive")
  expect_s4_class(dde_removed, "DeeDeeExperiment")
  expect_equal(length(dea_info(dde_removed)), 3)

  expect_warning(dde_removed <- remove_dea(dde, "lol"))

  expect_error(remove_dea(dde))

  dde_overlap_add <- add_fea(dde_overlap, fea = list(INFg_vs_Naive = topGO_results_list$ifng_vs_naive))
  dde_overlap_add <- add_fea(dde_overlap_add,
    fea = list(gPro_res = gost_res$result)
  )

  new_remove_dea <- remove_dea(dde_overlap_add,
    dea_name = "ifng_vs_naive",
    remove_linked_fea = TRUE
  )

  expect_equal(length(fea_info(new_remove_dea)), 5)

  expect_error(remove_dea(dde, dea_name = NULL))

  ## adding FEA ----------------------------------------------------------------

  topGO_Salm_naive <- topGO_results_list$salmonella_vs_naive
  topGO_IFNg_naive <- topGO_results_list$ifng_vs_naive

  dde2 <- add_fea(dde,
    fea = list(
      topGO_Salm_naive = topGO_Salm_naive,
      topGO_IFNg_naive = topGO_IFNg_naive
    ),
    fea_tool = "topGO"
  )

  expect_s4_class(dde2, "DeeDeeExperiment")
  expect_equal(length(fea_info(dde2)), 2)

  dde3 <- DeeDeeExperiment(se = se_macrophage_noassays)
  expect_warning(add_fea(dde3,
    fea = list(
      topGO_Salm_naive = topGO_Salm_naive,
      topGO_IFNg_naive = topGO_IFNg_naive
    ),
    verbose = TRUE
  ))

  expect_error(add_fea(dde3, fea = list(
    topGO_Salm_naive = topGO_Salm_naive,
    topGO_IFNg_naive
  )))

  dde3 <- add_fea(dde3,
    fea = list(
      topGO_Salm_naive = topGO_Salm_naive,
      topGO_IFNg_naive = topGO_IFNg_naive
    ),
    verbose = TRUE
  )
  expect_error({
    add_fea(dde3,
      fea = list(topGO_Salm_naive = topGO_Salm_naive),
      force = FALSE
    )
  })

  expect_error(add_fea(dde3, fea = list(
    FE1 = topGO_Salm_naive,
    FE1 = topGO_IFNg_naive
  )))

  dde3 <- add_fea(dde3, fea = list(gPro_salmonella_vs_naive = gost_res$result))

  expect_equal(fea_info(dde3)$gPro_salmonella_vs_naive$fe_tool, "gProfiler")

  expect_message(dde3 <- add_fea(dde3,
    fea =
      list(salmonella_vs_naive = enrichr_res$KEGG_2019_Human),
    verbose = TRUE
  ))
  dde_de_empty <- DeeDeeExperiment(se_macrophage_noassays)

  expect_error(add_fea(dde_de_empty, fea = list(topGO_results_list$salmo_both)))

  expect_error(add_fea(dde_overlap_add,
    fea = list(gPro_res = gost_res)
  ))

  expect_error(add_fea(dde_overlap_add,
    fea = topGO_results_list$ifng_vs_naive,
    de_name = NA
  ))

  expect_error(add_fea(dde_overlap_add,
    fea = topGO_results_list$ifng_vs_naive,
    fea_tool = NA
  ))

  dde_with_info <- add_scenario_info(dde_overlap,
    dea_name = "ifng_vs_naive",
    info = "here goes some txt about the contrast ifng_vs_naive"
  )

  expect_s4_class(dde_with_info, "DeeDeeExperiment")

  expect_no_error(add_fea(dde_with_info,
    fea = topGO_results_list$ifng_vs_naive,
    de_name = "ifng_vs_naive",
    verbose = TRUE
  ))

  expect_warning(add_fea(dde_with_info,
    fea = topGO_results_list$ifng_vs_naive,
    de_name = "IFNg_vs_naive"
  ))

  expect_message(add_fea(dde_with_info,
    fea = list(topGO_ifng_vs_naive = topGO_results_list$ifng_vs_naive),
    verbose = TRUE
  ))

  expect_message(add_fea(dde_with_info,
    fea = list(ifng_vs_naive = topGO_results_list$ifng_vs_naive),
    verbose = TRUE,
    force = TRUE
  ))

  expect_error(add_fea(dde_with_info,
    fea = topGO_results_list$ifng_vs_naive,
    fea_tool = c("topGO", "topGO")
  ))

  expect_message(add_fea(dde_with_info,
    fea = list(
      FEA1 = clusterPro_res$ifng_vs_naive,
      FEA2 = gost_res$result,
      FEA3 = gsea_res,
      FEA4 = fgseaRes
    ),
    fea_tool = c(
      "clusterProfiler",
      "gProfiler",
      "gsea",
      "fgsea"
    )
  ))


  ## removing FEA --------------------------------------------------------------

  dde2 <- remove_fea(dde2, "topGO_IFNg_naive")
  expect_s4_class(dde2, "DeeDeeExperiment")
  expect_equal(length(fea_info(dde2)), 1)

  expect_error(remove_fea(dde2, "IFNgVSnaive"))

  fea_name <- character(0)
  expect_error({
    remove_fea(dde3, fea_name)
  })

  ## retrieving DEA ------------------------------------------------------------

  expect_warning(dea(dde, verbose = TRUE))

  expect_error(dea(dde3, dea_name = "contrast1"))

  expect_error(dea(dde, dea_name = c("salmonella_vs_naive", "salmo_both")))

  expect_error(dea(dde, format = "simple"))

  expect_error(dea(dde3))

  expect_error(dea(dde_overlap,
    dea_name = "ifng_vs_naive",
    extra_rd = "ifng_vs_naive_pvalue"
  ))

  extract_dea <- dea(dde_overlap, dea_name = "ifng_vs_naive", format = "original")

  expect_s4_class(extract_dea, "DESeqResults")

  expect_error(get_dea_list(dde_overlap, format = "simple"))

  expect_error(dea(dde_overlap_add, extra_rd = NA))

  expect_warning(dea(dde_overlap_add,
    extra_rd = c("guiga", "other"),
    verbose = TRUE
  ))

  extract_dea <- dea(dde_overlap, dea_name = "ifng_vs_naive", type = "data.frame")

  expect_s3_class(extract_dea, "data.frame")

  extract_dea <- dea(dde_overlap, dea_name = "ifng_vs_naive")

  expect_s4_class(extract_dea, "DFrame")


  ## retrieving FEA ------------------------------------------------------------

  expect_s3_class(fea(dde2, "topGO_Salm_naive"), "data.frame")

  expect_error(get_fea_list(dde3, dea_name = c("salmonella_vs_naive", "salmo_both")))

  expect_error(get_fea_list(dde))

  expect_error(fea(dde, format = "simple"))

  expect_error(fea(dde))

  expect_warning(fea(dde_overlap, verbose = TRUE))

  expect_error(fea(dde_overlap, fea_name = c("salmonella_vs_naive", "ifng_vs_naive")))

  expect_error(fea(dde_overlap, fea_name = "sthg else"))

  expect_s3_class(fea(dde_overlap,
    fea_name = "salmonella_vs_naive",
    format = "original"
  ), "data.frame")

  expect_error(get_fea_list(dde_overlap,
    dea_name = "ifng_vs_naive",
    format = "simple"
  ))

  expect_warning(get_fea_list(dde_overlap_add, dea_name = "INFg_vs_Naive"))

  expect_warning(fea(dde_overlap_add, verbos = TRUE, format = "original"))

  expect_length(get_fea_list(dde_overlap_add, format = "minimal", verbose = TRUE), 6)

  expect_length(get_fea_list(dde_overlap_add, dea_name = "ifng_vs_naive", format = "original", verbose = TRUE), 1)

  expect_equal(length(fea_info(dde_overlap_add)), 6)

  ## adding scenario info ------------------------------------------------------

  expect_error(add_scenario_info(dde_overlap,
    dea_name = "i dont exist"
  ))

  expect_error(add_scenario_info(dde_overlap,
    dea_name = 2
  ))



  expect_error(add_scenario_info(dde_with_info,
    dea_name = "ifng_vs_naive",
    info = "sthg else about ifng_vs_naive"
  ))

  expect_error(add_scenario_info(dde_overlap,
    dea_name = c("i dont exist", "ifng_vs_naive")
  ))

  expect_error(add_scenario_info(dde_overlap_add,
    dea_name = "salmo_both",
    info = NA
  ))

  expect_error(add_scenario_info(dde_overlap_add,
    dea_name = "salmo_both",
    info = data.frame(Info = "here is some context")
  ))
})
