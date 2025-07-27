# Functions for FRiP, TSS enrichment, etc.

#' Calculate FRiP (Fraction of Reads in Peaks)
#'
#' @param fragments GRanges object with fragment coordinates and cell barcode metadata
#' @param peaks GRanges object of peak regions (e.g., from MACS2 or ArchR)
#' @return data.table with cell barcodes and FRiP values
#'
calc_frip <- function(fragments, peaks) {
  require(GenomicRanges)
  require(data.table)

  if (is.null(mcols(fragments)$barcode)) {
    stop("Fragments must have a 'barcode' metadata column.")
  }

  # Total fragments per cell
  total_per_cell <- table(mcols(fragments)$barcode)

  # Count overlaps
  overlaps <- findOverlaps(fragments, peaks)
  peak_barcodes <- mcols(fragments[queryHits(overlaps)])$barcode
  peak_per_cell <- table(peak_barcodes)

  # Handle NA safely
  peak_counts <- as.integer(peak_per_cell[names(total_per_cell)])
  peak_counts[is.na(peak_counts)] <- 0

  frip_dt <- data.table(
    barcode = names(total_per_cell),
    total = as.integer(total_per_cell),
    in_peaks = peak_counts
  )
  frip_dt[, frip := in_peaks / total]

  return(frip_dt)
}