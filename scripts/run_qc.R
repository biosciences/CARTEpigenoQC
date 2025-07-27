# Main QC pipeline entry point
# CARTEpigenoQC QC Pipeline using real fragment.tsv.gz file

library(data.table)
library(GenomicRanges)
library(rtracklayer)

# ---- INPUT ----
fragment_path <- "../inst/extdata/atac_v1_pbmc_10k_fragments.tsv.gz"
bed_path <- "../inst/extdata/CAR_standard_insertion_sites.bed"

# ---- Load CAR Insertion Sites ----
car_sites <- import(bed_path, format = "bed")
print("Loaded CAR insertion sites:")
print(car_sites)

# ---- Read 10x fragment.tsv.gz ----
# Format: chr, start, end, cell barcode, read count
fragments_dt <- fread(fragment_path, col.names = c("chr", "start", "end", "barcode", "count"))
cat("Read", nrow(fragments_dt), "fragments.\n")

# Convert to GRanges
frag_gr <- GRanges(
  seqnames = fragments_dt$chr,
  ranges = IRanges(fragments_dt$start, fragments_dt$end),
  barcode = fragments_dt$barcode
)

# ---- Overlap with CAR Insertion Sites ----
overlaps <- countOverlaps(frag_gr, car_sites)
cat("Number of fragments overlapping CAR sites:", sum(overlaps > 0), "\n")

# ---- Calculate FRiP ----
source("../R/qc_metrics.R")
frip_df <- calc_frip(fragments = frag_gr, peaks = car_sites)

# ---- Save FRiP Results ----
frip_output_path <- "../inst/extdata/frip_results.tsv"
fwrite(frip_df, file = frip_output_path, sep = "\t")
cat("FRiP results written to frip_results.tsv\n")

# ---- Compute CAR Insertion Site Coverage ----
source("../R/insertions.R")
insert_df <- compute_car_coverage(fragments = frag_gr, car_sites = car_sites)

# ---- Save CAR Insertion Coverage Results ----
insert_output_path <- "../inst/extdata/car_coverage.tsv"
fwrite(insert_df, file = insert_output_path, sep = "\t")
cat("CAR insertion site coverage results written to car_coverage.tsv\n")

# ---- Generate HTML Report (Optional) ----
source("../R/report.R")
frip_output_path_rmd <- "../extdata/frip_results.tsv"
insert_output_path_rmd <- "../extdata/car_coverage.tsv"
render_qc_report(
  input_rmd = "../inst/templates/report_template.Rmd",
  output_html = "../../docs/qc_report.html",
  params = list(
    sample_name = "PBMC_10k",
    frip_file = frip_output_path_rmd,
    insert_file = insert_output_path_rmd
  )
)
cat("HTML report generated at qc_report.html\n")