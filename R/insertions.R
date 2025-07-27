# Functions for CAR insertion site coverage

#' Compute coverage of CAR insertion sites per cell
#'
#' @param fragments GRanges object with cell barcode metadata
#' @param car_sites GRanges object of CAR insertion regions
#' @return data.table with per-cell insertion site overlap counts
#'
compute_car_coverage <- function(fragments, car_sites) {
  require(GenomicRanges)
  require(data.table)

  if (is.null(mcols(fragments)$barcode)) {
    stop("Fragments must contain a 'barcode' metadata column.")
  }

  # Find overlapping fragments
  hits <- findOverlaps(fragments, car_sites)
  overlapping_barcodes <- mcols(fragments[queryHits(hits)])$barcode
  site_names <- if (!is.null(mcols(car_sites)$name)) {
    as.character(mcols(car_sites)$name[subjectHits(hits)])
  } else {
    as.character(seqnames(car_sites)[subjectHits(hits)])
  }

  # Create data.table for counting
  dt <- data.table(
    barcode = overlapping_barcodes,
    site = site_names
  )

  coverage_dt <- dt[, .N, by = .(barcode, site)]
  setnames(coverage_dt, "N", "insert_reads")

  return(coverage_dt)
}