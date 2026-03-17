test_that("read_ontology parses gmt files", {
  gmt <- tempfile(fileext = ".gmt")
  writeLines(
    c(
      "TERM_A\tTerm A\t101\t102",
      "TERM_B\tTerm B\t103"
    ),
    gmt
  )

  ontology <- read_ontology(gmt)

  expect_equal(names(ontology), c("term2gene", "term2name"))
  expect_equal(ncol(ontology$term2gene), 2)
  expect_equal(ncol(ontology$term2name), 2)
  expect_true(all(c("TERM_A", "TERM_B") %in% ontology$term2name$term))
})

test_that("read_ontology parses delimited ontology files", {
  tsv <- tempfile(fileext = ".tsv")
  write.table(
    data.frame(
      term = c("TERM_A", "TERM_A", "TERM_B"),
      name = c("Term A", "Term A", "Term B"),
      gene = c("101", "102", "103")
    ),
    tsv,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  ontology <- read_ontology(tsv)

  expect_equal(nrow(ontology$term2gene), 3)
  expect_equal(nrow(ontology$term2name), 2)
})
