test_that("renaming", {
  ## renaming DEA --------------------------------------------------------------

  dde_de_empty <- DeeDeeExperiment(se_macrophage_noassays)

  expect_error(renameDEA(dde_de_empty,
    old_name = "ifng_vs_naive", new_name = "IFNgvsNaive"
  ))

  dde <- DeeDeeExperiment(
    se_macrophage_noassays,
    de_results = de_named_list
  )

  dde_rename <- renameDEA(dde,
    old_name = "salmonella_vs_naive",
    new_name = "SalmvsNaive"
  )

  expect_s4_class(dde_rename, "DeeDeeExperiment")

  expect_equal(names(DEAInfo(dde_rename)), c(
    "ifng_vs_naive", "ifngsalmo_vs_naive",
    "SalmvsNaive", "salmo_both"
  ))

  expect_error(renameDEA(dde,
    old_name = "salmonella_vs_naive",
    new_name = "ifng_vs_naive"
  ))

  expect_error(renameDEA(dde,
    old_name = NULL,
    new_name = "ifng_vs_naive"
  ))

  expect_error(renameDEA(dde,
    old_name = "salmonella_vs_naive",
    new_name = NULL
  ))

  expect_error(renameDEA(dde,
    old_name = "contrast1",
    new_name = "ifng_vs_naive"
  ))

  expect_error(renameDEA(dde,
    old_name = 1,
    new_name = "1"
  ))

  dde3 <- addFEA(dde_rename, fea = topGO_results_list)

  expect_error(renameDEA(dde3,
    old_name = 1,
    new_name = "1"
  ))

  expect_error(renameDEA(dde,
    old_name = "salmonella_vs_naive",
    new_name = c("ifng_vs_naive", "new_column")
  ))

  expect_error(renameDEA(dde,
    old_name = c("salmonella_vs_naive", "salmo_both"),
    new_name = c("salmonel_vs_naive", "salmonel_vs_naive")
  ))


  dde3_rename <- renameDEA(dde3,
    old_name = "ifngsalmo_vs_naive",
    new_name = "IFNg_SalmvsNaive"
  )

  expect_s4_class(dde3_rename, "DeeDeeExperiment")

  expect_equal(names(DEAInfo(dde3_rename)), c(
    "ifng_vs_naive", "IFNg_SalmvsNaive",
    "SalmvsNaive", "salmo_both"
  ))

  ## renaming FEA --------------------------------------------------------------

  expect_error(renameFEA(dde3,
    old_name = "contrast1",
    new_name = "topGO_IFNg_naive"
  ))

  expect_error(renameFEA(dde,
    old_name = "topGO_Salm_naive",
    new_name = c("salmonella_vs_naive", "new_column")
  ))

  expect_error(renameFEA(dde_overlap,
    old_name = c("ifng_vs_naive", "ifngsalmo_vs_naive"),
    new_name = "naive"
  ))

  expect_error(renameFEA(dde_overlap,
    old_name = c("ifng_vs_naive", "ifngsalmo_vs_naive"),
    new_name = c("naive", "naive")
  ))

  expect_error(renameFEA(dde_overlap,
    old_name = c("ifng_vs_naive"),
    new_name = c("ifngsalmo_vs_naive")
  ))

  expect_error(renameFEA(dde3,
    old_name = NULL,
    new_name = "ifng_vs_naive"
  ))

  expect_error(renameFEA(dde3,
    old_name = "salmonella_vs_naive",
    new_name = NULL
  ))

  expect_equal(
    FEAInfo(dde3_rename)[["ifngsalmo_vs_naive"]][["de_name"]],
    "IFNg_SalmvsNaive"
  )

  fea_rename <- renameFEA(dde3,
    old_name = "salmonella_vs_naive",
    new_name = "topGO_SalmonellavsNaive"
  )
  expect_s4_class(fea_rename, "DeeDeeExperiment")

  expect_equal(FEANames(fea_rename), c(
    "ifng_vs_naive",
    "ifngsalmo_vs_naive",
    "topGO_SalmonellavsNaive",
    "salmo_both"
  ))
})
