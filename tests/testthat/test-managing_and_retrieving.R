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
  dde_new <- addDEA(dde, new_del)
  expect_s4_class(dde_new, "DeeDeeExperiment")
  expect_equal(length(infoDEA(dde)), 4)
  expect_equal(length(infoDEA(dde_new)), 6)

  expect_error({
    addDEA(x = dde, dea = list(de_named_list$ifng_vs_naive))
  })

  expect_error({
    dde_add <- addDEA(x = "sthg else", dea = new_del)
  })

  expect_error({
    dde_add <- addDEA(x = dde, dea = list(
      ifng2 = de_named_list$ifng_vs_naive,
      ifng2 = de_named_list$ifngsalmo_vs_naive
    ))
  })

  dde_edgeR <- addDEA(dde, dea = list(DGEExact_IFNg_both = dge_exact_IFNg_both))
  expect_s4_class(dde_edgeR, "DeeDeeExperiment")
  expect_equal(length(infoDEA(dde_edgeR)), 5)

  dde_limma <- addDEA(dde, dea = list(de_limma = de_limma))
  expect_s4_class(dde_limma, "DeeDeeExperiment")
  expect_equal(length(infoDEA(dde_limma)), 5)

  expect_no_error(addDEA(dde_limma,
    dea = list(de_limma = de_limma),
    force = TRUE
  ))

  genes <- rownames(de_limma)[1:20]

  de_custom <- data.frame(
    p_val = rep(0.5, 20),
    adj_pvalue = rep(0.5, 20),
    gene = genes
  )
  expect_error(addDEA(dde, dea = list(de_custom = de_custom)))

  expect_error(addDEA(dde, "ifng_vs_naive"))

  dde_overlap <- DeeDeeExperiment(se_macrophage_noassays,
    de_results = de_named_list,
    enrich_results = topGO_results_list
  )

  expect_error(addDEA(dde_overlap,
    dea = list(ifng_vs_naive = de_named_list$ifng_vs_naive)
  ))

  gene_ids <- paste0("gene", 1:5)

  res_de_df <- data.frame(
    log2FoldChange = c(2.1, -1.3, 0.5, NA, -0.9),
    pvalue = c(0.01, 0.20, 0.05, NA, 0.80),
    padj = c(0.05, 0.25, 0.10, NA, 0.90),
    row.names = gene_ids
  )

  dde_overlap <- addDEA(dde_overlap, dea = res_de_df)
  expect_s4_class(dde_overlap, "DeeDeeExperiment")
  expect_equal(length(infoDEA(dde_overlap)), 5)
  expect_true(all(
    c("res_de_df_log2FoldChange", "res_de_df_pvalue", "res_de_df_padj")
    %in% colnames(rowData(dde_overlap))
  ))

  expect_false(all(rownames(res_de_df) %in% rownames(dde_overlap)))

  gene_ids <- c(gene_ids, rownames(de_named_list$ifng_vs_naive)[1:6])

  de_custom <- data.frame(
    log2FoldChange = rep(0.5, 11),
    pvalue = rep(0.5, 11),
    padj = rep(0.5, 11),
    row.names = gene_ids
  )

  expect_warning(DeeDeeExperiment(se_macrophage_noassays, de_results = de_custom))

  ## removing DEA --------------------------------------------------------------

  dde_removed <- removeDEA(dde, "ifngsalmo_vs_naive")
  expect_s4_class(dde_removed, "DeeDeeExperiment")
  expect_equal(length(infoDEA(dde_removed)), 3)

  expect_warning(dde_removed <- removeDEA(dde, "lol"))

  expect_error(removeDEA(dde))

  dde_overlap_add <- addFEA(dde_overlap,
                            fea = list(INFg_vs_Naive = topGO_results_list$ifng_vs_naive))
  dde_overlap_add <- addFEA(dde_overlap_add,
    fea = list(gPro_res = gost_res$result)
  )

  new_remove_dea <- removeDEA(dde_overlap_add,
    dea_name = "ifng_vs_naive",
    remove_linked_fea = TRUE
  )

  expect_equal(length(infoFEA(new_remove_dea)), 5)

  expect_error(removeDEA(dde, dea_name = NULL))

  expect_equal(length(infoDEA(removeDEA(dde_overlap,
                                        dea_name = "res_de_df"))), 4)

  ## adding FEA ----------------------------------------------------------------

  topGO_Salm_naive <- topGO_results_list$salmonella_vs_naive
  topGO_IFNg_naive <- topGO_results_list$ifng_vs_naive

  dde2 <- addFEA(dde,
    fea = list(
      topGO_Salm_naive = topGO_Salm_naive,
      topGO_IFNg_naive = topGO_IFNg_naive
    ),
    fea_tool = "topGO"
  )

  expect_s4_class(dde2, "DeeDeeExperiment")
  expect_equal(length(infoFEA(dde2)), 2)

  dde3 <- DeeDeeExperiment(sce = se_macrophage_noassays)
  expect_warning(addFEA(dde3,
    fea = list(
      topGO_Salm_naive = topGO_Salm_naive,
      topGO_IFNg_naive = topGO_IFNg_naive
    ),
    verbose = TRUE
  ))

  expect_error(addFEA(dde3, fea = list(
    topGO_Salm_naive = topGO_Salm_naive,
    topGO_IFNg_naive
  )))

  dde3 <- addFEA(dde3,
    fea = list(
      topGO_Salm_naive = topGO_Salm_naive,
      topGO_IFNg_naive = topGO_IFNg_naive
    ),
    verbose = TRUE
  )
  expect_error({
    addFEA(dde3,
      fea = list(topGO_Salm_naive = topGO_Salm_naive),
      force = FALSE
    )
  })

  expect_error(addFEA(dde3, fea = list(
    FE1 = topGO_Salm_naive,
    FE1 = topGO_IFNg_naive
  )))

  dde3 <- addFEA(dde3, fea = list(gPro_salmonella_vs_naive = gost_res$result))

  expect_equal(infoFEA(dde3)$gPro_salmonella_vs_naive$fe_tool, "gProfiler")

  expect_message(dde3 <- addFEA(dde3,
    fea =
      list(salmonella_vs_naive = enrichr_res$KEGG_2019_Human),
    verbose = TRUE
  ))
  dde_de_empty <- DeeDeeExperiment(se_macrophage_noassays)

  expect_error(addFEA(dde_de_empty, fea = list(topGO_results_list$salmo_both)))

  expect_error(addFEA(dde_overlap_add,
    fea = list(gPro_res = gost_res)
  ))

  expect_error(addFEA(dde_overlap_add,
    fea = topGO_results_list$ifng_vs_naive,
    de_name = NA
  ))

  expect_error(addFEA(dde_overlap_add,
    fea = topGO_results_list$ifng_vs_naive,
    fea_tool = NA
  ))

  dde_with_info <- addScenarioInfo(dde_overlap,
    dea_name = "ifng_vs_naive",
    info = "here goes some txt about the contrast ifng_vs_naive"
  )

  expect_s4_class(dde_with_info, "DeeDeeExperiment")

  expect_no_error(addFEA(dde_with_info,
    fea = topGO_results_list$ifng_vs_naive,
    de_name = "ifng_vs_naive",
    verbose = TRUE
  ))

  expect_warning(addFEA(dde_with_info,
    fea = topGO_results_list$ifng_vs_naive,
    de_name = "IFNg_vs_naive"
  ))

  expect_message(addFEA(dde_with_info,
    fea = list(topGO_ifng_vs_naive = topGO_results_list$ifng_vs_naive),
    verbose = TRUE
  ))

  expect_message(addFEA(dde_with_info,
    fea = list(ifng_vs_naive = topGO_results_list$ifng_vs_naive),
    verbose = TRUE,
    force = TRUE
  ))

  expect_error(addFEA(dde_with_info,
    fea = topGO_results_list$ifng_vs_naive,
    fea_tool = c("topGO", "topGO")
  ))

  expect_message(addFEA(dde_with_info,
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

  dde2 <- removeFEA(dde2, "topGO_IFNg_naive")
  expect_s4_class(dde2, "DeeDeeExperiment")
  expect_equal(length(infoFEA(dde2)), 1)

  expect_error(removeFEA(dde2, "IFNgVSnaive"))

  fea_name <- character(0)
  expect_error({
    removeFEA(dde3, fea_name)
  })

  ## retrieving DEA ------------------------------------------------------------

  expect_warning(getDEA(dde, verbose = TRUE))

  expect_error(getDEA(dde3, dea_name = "contrast1"))

  expect_error(getDEA(dde, dea_name = c("salmonella_vs_naive", "salmo_both")))

  expect_error(getDEA(dde, format = "simple"))

  expect_error(getDEA(dde, type = c("DFrame", "dataframe")))

  expect_error(getDEA(dde, type = c("table")))

  expect_error(getDEA(dde3))



  expect_error(getDEA(dde_overlap,
    dea_name = "ifng_vs_naive",
    extra_rd = "ifng_vs_naive_pvalue"
  ))

  extract_dea <- getDEA(dde_overlap, dea_name = "ifng_vs_naive", format = "original")

  expect_s4_class(extract_dea, "DESeqResults")

  expect_error(getDEAList(dde_overlap, format = "simple"))

  expect_error(getDEA(dde_overlap_add, extra_rd = NA))

  expect_warning(getDEA(dde_overlap_add,
    extra_rd = c("guiga", "other"),
    verbose = TRUE
  ))

  extract_dea <- getDEA(dde_overlap, dea_name = "ifng_vs_naive", type = "data.frame")

  expect_s3_class(extract_dea, "data.frame")

  extract_dea <- getDEA(dde_overlap, dea_name = "ifng_vs_naive")

  expect_s4_class(extract_dea, "DFrame")


  ## retrieving FEA ------------------------------------------------------------

  expect_s3_class(getFEA(dde2, "topGO_Salm_naive"), "data.frame")

  expect_error(getFEAList(dde3, dea_name = c("salmonella_vs_naive", "salmo_both")))

  expect_error(getFEAList(dde))

  expect_error(getFEA(dde, format = "simple"))

  expect_error(getFEA(dde))

  expect_warning(getFEA(dde_overlap, verbose = TRUE))

  expect_error(getFEA(dde_overlap, fea_name = c("salmonella_vs_naive", "ifng_vs_naive")))

  expect_error(getFEA(dde_overlap, fea_name = "sthg else"))

  expect_s3_class(getFEA(dde_overlap,
    fea_name = "salmonella_vs_naive",
    format = "original"
  ), "data.frame")

  expect_error(getFEAList(dde_overlap,
    dea_name = "ifng_vs_naive",
    format = "simple"
  ))

  expect_warning(getFEAList(dde_overlap_add, dea_name = "INFg_vs_Naive"))

  expect_warning(getFEA(dde_overlap_add, verbos = TRUE, format = "original"))

  expect_length(getFEAList(dde_overlap_add, format = "minimal", verbose = TRUE), 6)

  expect_length(getFEAList(dde_overlap_add, dea_name = "ifng_vs_naive", format = "original", verbose = TRUE), 1)

  expect_equal(length(infoFEA(dde_overlap_add)), 6)

  ## adding scenario info ------------------------------------------------------

  expect_error(addScenarioInfo(dde_overlap,
    dea_name = "i dont exist"
  ))

  expect_error(addScenarioInfo(dde_overlap,
    dea_name = 2
  ))

  expect_error(addScenarioInfo(dde_with_info,
    dea_name = "ifng_vs_naive",
    info = "sthg else about ifng_vs_naive"
  ))

  expect_error(addScenarioInfo(dde_overlap,
    dea_name = c("i dont exist", "ifng_vs_naive")
  ))

  expect_error(addScenarioInfo(dde_overlap_add,
    dea_name = "salmo_both",
    info = NA
  ))

  expect_error(addScenarioInfo(dde_overlap_add,
    dea_name = "salmo_both",
    info = data.frame(Info = "here is some context")
  ))
})
