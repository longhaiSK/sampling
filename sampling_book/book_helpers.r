## Helpers shared by the book chapters (not the estimating functions)

library(gt)

## A long data table. The HTML book shows every row in a paged, interactive gt
## table; a PDF has no pager, so it shows a snapshot of the first `n` rows with
## a note pointing to the full table in the HTML version.
book_table <- function(df, n = 10) {
  if (knitr::is_latex_output() && nrow(df) > n) {
    gt(head(df, n)) |>
      tab_source_note(sprintf(
        "Showing the first %d of %d rows; the full table is in the HTML version of the book.",
        n, nrow(df)))
  } else {
    gt(df) |>
      opt_interactive(use_pagination = TRUE, page_size_default = n)
  }
}
