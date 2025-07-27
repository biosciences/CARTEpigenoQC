test_that("FRiP and CAR input files are loadable", {
  frip_file <- system.file("extdata", "frip_results.tsv", package = "CARTEpigenoQC")
  insert_file <- system.file("extdata", "car_coverage.tsv", package = "CARTEpigenoQC")
  expect_true(file.exists(frip_file))
  expect_true(file.exists(insert_file))
})

test_that("calc_frip returns expected columns", {
  library(GenomicRanges)
  fragments <- GRanges(seqnames = "chr1", ranges = IRanges::IRanges(start = 100, end = 200), barcode = "cell_1")
  peaks <- GRanges(seqnames = "chr1", ranges = IRanges::IRanges(start = 150, end = 180))
  result <- CARTEpigenoQC::calc_frip(fragments, peaks)
  expect_true("frip" %in% colnames(result))
  expect_gte(result$frip[1], 0)
})
