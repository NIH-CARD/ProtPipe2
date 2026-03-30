dat <- data.table::fread(basic_example_path("iPSC.csv"))
ipsc_genes <- data.table::fread(basic_example_path("stem_cell_gene.csv"))$Gene
neuron_genes <- data.table::fread(basic_example_path("neuron_gene.csv"))$Gene
dat_pro <- create_se(dat)

test_that("basic heatmap", {
  p <- plot_proteomics_heatmap(dat_pro)
  expect_s3_class(p, "ggplot")
})

test_that("heatmap subset", {
  p <- plot_proteomics_heatmap(dat_pro, protmeta_col = "PG.Genes", genes = ipsc_genes)
  expect_s3_class(p, "ggplot")
})

test_that("heatmap subset with condition", {
  p <- plot_proteomics_heatmap(dat_pro, protmeta_col = "PG.Genes", genes = ipsc_genes, condition = 'base_condition')
  expect_s3_class(p, "ggplot")
})

test_that("heatmap subset with condition and row clustering", {
  p <- plot_proteomics_heatmap(dat_pro, protmeta_col = "PG.Genes", genes = ipsc_genes,
                               condition = 'base_condition', cluster_rows = T)
  expect_s3_class(p, "ggplot")
})

test_that("heatmap subset with condition and col clustering", {
  p <- plot_proteomics_heatmap(dat_pro, protmeta_col = "PG.Genes", genes = ipsc_genes,
                               condition = 'base_condition', cluster_cols = T)
  expect_s3_class(p, "ggplot")
})

test_that("heatmap subset with condition and col and row clustering", {
  p <- plot_proteomics_heatmap(dat_pro, protmeta_col = "PG.Genes", genes = ipsc_genes,
                               condition = 'base_condition', cluster_cols = T,
                               cluster_rows = T)
  expect_s3_class(p, "ggplot")
})
