# Functions for loading fragments and metadata

#' Load fragment.tsv.gz file as GRanges
#'
#' @param path Path to the fragment.tsv.gz file
#' @param barcode_filter Optional vector of barcodes to keep
#' @return GRanges object with barcode metadata column
#'
load_fragments <- function(path, barcode_filter = NULL) {
  require(data.table)
  require(GenomicRanges)

  message("Reading fragments from: ", path)
  frag_dt <- fread(path, col.names = c("chr", "start", "end", "barcode", "count"))

  if (!is.null(barcode_filter)) {
    frag_dt <- frag_dt[barcode %in% barcode_filter]
    message("Filtered to ", length(unique(frag_dt$barcode)), " barcodes.")
  }

  gr <- GRanges(
    seqnames = frag_dt$chr,
    ranges = IRanges(start = frag_dt$start, end = frag_dt$end),
    barcode = frag_dt$barcode
  )

  return(gr)
}


#' Load CAR-T insertion sites from BED file
#'
#' @param bed_path Path to BED file
#' @return GRanges object with insertion site ranges
#'
load_car_sites <- function(bed_path) {
  require(rtracklayer)
  car_sites <- import(bed_path, format = "bed")
  message("Loaded ", length(car_sites), " CAR insertion sites.")
  return(car_sites)
}