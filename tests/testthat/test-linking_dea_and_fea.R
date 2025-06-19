test_that("linking dea and fea", {

  ## linking DEA and FEA -------------------------------------------------------

  dde <- DeeDeeExperiment(
    se_macrophage_noassays,
    de_results = de_named_list
  )

  expect_error(link_dea_and_fea(dde,
                                dea_name = "contrast1",
                                fea_name = "topGO_IFNg_naive"))

  expect_error(link_dea_and_fea(dde,
                                dea_name = "salmonella_vs_naive",
                                fea_name = "topGO"))
  expect_error(link_dea_and_fea(dde,
                                "not there",
                                "not there"))

  dde_overlap <- DeeDeeExperiment(se_macrophage_noassays,
                                  de_results = de_named_list,
                                  enrich_results = topGO_results_list)

  expect_error(link_dea_and_fea(dde_overlap,
                                "new_name",
                                "salmonella_vs_naive"))

  expect_error(link_dea_and_fea(dde_overlap,
                                "ifng_vs_naive",
                                "salmonella_vs_naive",
                                force = FALSE))

  expect_error(link_dea_and_fea(dde_overlap,
                                dea_name = 2,
                                fea_name = "salmonella_vs_naive"))

  expect_error(link_dea_and_fea(dde_overlap,
                                dea_name = "salmonella_vs_naive",
                                fea_name = 2))

  dde_overlap_add <- add_fea(dde_overlap, fea = list(INFg_vs_Naive = topGO_results_list$ifng_vs_naive))
  dde_overlap_add <- add_fea(dde_overlap_add,
                             fea = list(gPro_res = gost_res$result))

  expect_warning(link_dea_and_fea(dde_overlap_add,
                                  "ifng_vs_naive",
                                  "salmonella_vs_naive",
                                  force = TRUE))
})
