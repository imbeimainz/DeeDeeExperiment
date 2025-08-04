test_that("summary", {
  ## summary -------------------------------------------------------------------

  dde_no_fea <- DeeDeeExperiment(
    sce = se_macrophage_noassays,
    de_results = de_named_list
  )
  expect_no_error(summary(dde_no_fea))

  dde_no_dea <- DeeDeeExperiment(
    sce = se_macrophage_noassays,
    enrich_results = topGO_results_list
  )
  expect_no_error(summary(dde_no_dea))

  dde_empty <- DeeDeeExperiment(sce = se_macrophage_noassays)
  expect_no_error(summary(dde_empty))

  dde_with_scenario <-
    addScenarioInfo(dde_no_fea,
      dea_name = "ifng_vs_naive",
      info = "this is the interferon gamma vs naive setting"
    )

  expect_no_error(summary(dde_with_scenario, show_scenario_info = TRUE))

  dde_with_info <-
    addScenarioInfo(dde_no_fea,
      dea_name = "salmo_both",
      info = "this is the interferon gamma vs salmonella setting"
    )

  expect_no_error(
    summary(dde_with_info, show_scenario_info = TRUE)
  )
})
