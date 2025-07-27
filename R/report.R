# Functions for rendering HTML QC reports

#' Render CARTEpigenoQC HTML report
#'
#' @param input_rmd Path to the R Markdown report template
#' @param output_html Path to output HTML file
#' @param params List of named parameters to pass to the report
#' @return Path to the rendered HTML file
#'
render_qc_report <- function(input_rmd, output_html, params = list()) {
  require(rmarkdown)

  message("Rendering HTML report to: ", output_html)

  rendered_path <- rmarkdown::render(
    input = input_rmd,
    output_file = output_html,
    params = params,
    envir = new.env(parent = globalenv()),
    clean = TRUE,
    quiet = TRUE
  )

  return(normalizePath(rendered_path))
}