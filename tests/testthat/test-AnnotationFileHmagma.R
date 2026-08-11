test_that("AnnotationFileHmagma produces the expected annotation file", {
  data("test_dataset.microglia.hg19")
  out <- tempfile()

  # The bundled toy dataset mixes NCBI-style ("chr23") and UCSC-style
  # ("chrX"/"chrY") seqlevels, which triggers a harmless Seqinfo merge warning.
  suppressWarnings(
    AnnotationFileHmagma(
      hic = hic,
      promoterRegions = promoterRegions,
      enhancerRegions = enhancerRegions,
      snps = snps,
      annotated_genes = annotated_genes,
      snpgeneexon = snpgeneexon,
      AnnotationFile = out
    )
  )

  result_file <- paste0(out, ".transcript.annot")
  expect_true(file.exists(result_file))

  result <- readLines(result_file)
  expect_equal(length(result), 53)

  # Each line should be: ENSG id, "chr:start:end" index, then one or more rsids,
  # all tab-separated.
  fields <- strsplit(result, "\t")
  expect_true(all(vapply(fields, length, integer(1)) >= 3))
  expect_true(all(grepl("^ENSG", vapply(fields, `[[`, character(1), 1))))
  expect_true(all(grepl("^[0-9XYM_a-z]+:[0-9]+:[0-9]+$", vapply(fields, `[[`, character(1), 2))))
})

test_that("AnnotationFileHmagma requires regulatoryRegions or promoter+enhancer regions", {
  data("test_dataset.microglia.hg19")

  expect_error(
    AnnotationFileHmagma(
      hic = hic,
      snps = snps,
      annotated_genes = annotated_genes,
      snpgeneexon = snpgeneexon,
      AnnotationFile = tempfile()
    ),
    "regulatoryRegions"
  )
})
