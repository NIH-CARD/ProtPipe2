dat <- data.table::fread("EXAMPLES/basic_example_data/iPSC.csv")
ipsc_genes <- data.table::fread("EXAMPLES/basic_example_data/stem_cell_gene.csv")$Gene
neuron_genes <- data.table::fread("EXAMPLES/basic_example_data/neuron_gene.csv")$Gene
dat_pro <- create_protdata(dat)

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
