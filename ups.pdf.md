# Unequal Probability Sampling

## Key Concepts and Formulas

Unequal Probability Sampling (UPS) improves estimation precision by sampling units with probabilities proportional to their size (PPS). If a cluster or county is larger, it has a higher chance of being included in the sample. This approach prevents extreme variance when cluster sizes vary drastically.

**1. UPS With Replacement (UPSWR)**
When sampling with replacement, we draw $n$ units. A single unit $i$ has a selection probability $\psi_i$ on each draw. 

*   **Total Estimator (Hansen-Hurwitz):** $\hat{t}_{pwr} = \frac{1}{n} \sum_{i=1}^n \frac{y_i}{\psi_i}$
*   **Ratio Estimator (Mean):** $\hat{\bar{y}}_{r} = \frac{\sum_{i=1}^n \frac{y_i}{\psi_i}}{\sum_{i=1}^n \frac{M_i}{\psi_i}}$

**2. UPS Without Replacement (UPSWO)**
When sampling without replacement, the inclusion probability for unit $i$ across all draws is $\pi_i$.

*   **Total Estimator (Horvitz-Thompson):** $\hat{t}_{HT} = \sum_{i=1}^n \frac{y_i}{\pi_i}$

---

## Functions and Packages

The estimating functions used throughout this book (including `ratio_est`, `reg_est`, `srs_est`, `upswr_total`, `upswr_ratio`, `cluster_upswr_ratio`, `cluster_upswo_ratio`, and the `format_est_gt` table-formatting helper used below) live in a single shared file, `samplingestimate.r`, which we source below.


::: {.cell}

:::


::: {.content-visible when-format="html"}

<details>
<summary>Show the shared estimating functions (`samplingestimate.r`)</summary>

```r
# =============================================================
# Shared estimation functions for the Elements of Sampling Survey
# book. Source this file from any chapter that needs point/SE/CI
# estimators or their gt-table formatting helpers, e.g.:
#
#   source("samplingestimate.r")
#
# =============================================================

suppressPackageStartupMessages({
    library(gt)
})

# ---------------------------------------------------------
# Shared gt-table helpers
# ---------------------------------------------------------

## Formats a gt table column as ordinary fixed-point, UNLESS its values are
## too small (< `small` in absolute value) or too large (>= `large`), in
## which case it switches that column to scientific notation (rendered as
## "m x 10^n" via gt's exp_style = "x10n") instead.
fmt_auto <- function (gt_tbl, data, columns, decimals = 3, small = 1e-3, large = 1e5)
{
    for (col in columns)
    {
        x <- data[[col]]
        x <- x[is.finite (x) & x != 0]
        if (length (x) == 0) next
        if (any (abs (x) < small) || any (abs (x) >= large))
            gt_tbl <- fmt_scientific (gt_tbl, columns = tidyselect::all_of (col), decimals = decimals, exp_style = "x10n")
        else
            gt_tbl <- fmt_number (gt_tbl, columns = tidyselect::all_of (col), decimals = decimals)
    }
    gt_tbl
}

## Formats a single scalar for inline formula display (fixed-point with
## thousands separators, switching to "m \times 10^{n}" LaTeX scientific
## notation for very small or very large magnitudes, mirroring fmt_auto()'s
## use of exp_style = "x10n") --- used to show worked-example steps such as
## \bar y = sum/n with the actual numbers plugged in.
fmt_num <- function (x, decimals = 3, small = 1e-3, large = 1e5)
{
    use_sci <- is.finite (x) & x != 0 & (abs (x) < small | abs (x) >= large)
    exps    <- ifelse (x == 0, 0, floor (log10 (abs (x))))
    mant    <- ifelse (x == 0, 0, x / 10^exps)
    ifelse (use_sci,
            sprintf ("%s \\times 10^{%d}", formatC (mant, format = "f", digits = decimals), exps),
            formatC (x, format = "f", digits = decimals, big.mark = ","))
}

## Truncates a row-level working data.frame whose LAST row is a summary row
## (e.g. "Sum"/"Total"): if there are more than head_n data rows (excluding
## that summary row), keeps only the first head_n data rows, inserts one
## "..." placeholder row, then the summary row --- so long working tables
## stay readable. id_col --- name of the row-label column.
truncate_working <- function (working, id_col, head_n = 6)
{
    n <- nrow (working) - 1
    if (n <= head_n + 1) return (working)

    dots <- working[1, ]
    dots[1, ] <- NA
    dots[1, id_col] <- "..."

    rbind (working[seq_len (head_n), ], dots, working[nrow (working), ])
}

# Helper function to format estimation outputs with appropriate mathematical headers
format_est_gt <- function(est_data, est_type = c("mean", "total", "ratio", "reg_mean", "reg_total", "ratio_mean", "domain_mean")) {
  est_type <- match.arg(est_type)

  sym_map <- list(
    "mean"        = c(est = "$\\bar{y}$",         se = "$\\mathrm{SE}(\\bar{y})$"),
    "total"       = c(est = "$\\hat{t}$",         se = "$\\mathrm{SE}(\\hat{t})$"),
    "ratio"       = c(est = "$\\hat{B}$",         se = "$\\mathrm{SE}(\\hat{B})$"),
    "reg_mean"    = c(est = "$\\bar{y}_{reg}$",   se = "$\\mathrm{SE}(\\bar{y}_{reg})$"),
    "reg_total"   = c(est = "$\\hat{t}_{reg}$",   se = "$\\mathrm{SE}(\\hat{t}_{reg})$"),
    "ratio_mean"  = c(est = "$\\bar{y}_{r}$",     se = "$\\mathrm{SE}(\\bar{y}_{r})$"),
    "domain_mean" = c(est = "$\\bar{y}_{d}$",     se = "$\\mathrm{SE}(\\bar{y}_{d})$")
  )

  est_sym <- sym_map[[est_type]]["est"]
  se_sym  <- sym_map[[est_type]]["se"]

  # Convert to data frame. Use unname() to strip hidden lm() coefficient names
  if (is.vector(est_data) && !is.list(est_data)) {
    df <- as.data.frame(t(unname(est_data)))
  } else {
    df <- as.data.frame(est_data)
  }

  # Force standard column names so gt() never gets confused
  colnames(df) <- c("Est.", "S.E.", "ci.low", "ci.upp")

  # Check if we should use rownames as a stub
  use_stub <- !is.null(rownames(df)) && !all(rownames(df) == as.character(1:nrow(df)))

  df |>
    gt(rownames_to_stub = use_stub) |>
    cols_label(
      Est.   = md(est_sym),
      S.E.   = md(se_sym),
      ci.low = md("95% CI Lower"),
      ci.upp = md("95% CI Upper")
    ) |>
    fmt_number(
      columns = everything(),
      decimals = 4
    ) |>
    tab_options(table.width = pct(70))
}

## Formats an lm object as two gt tables, replacing the plain-text
## summary(lmfit) output:
##   $coef --- coefficient table (Estimate, S.E., t value, p value)
##   $fit  --- one-row table of model fit statistics, including R^2
## returns list($coef, $fit)
format_lm_gt <- function (lmfit)
{
    s <- summary (lmfit)

    coef_df <- as.data.frame (s$coefficients)
    colnames (coef_df) <- c ("Estimate", "S.E.", "t value", "p value")

    coef_gt <- coef_df |>
        gt (rownames_to_stub = TRUE) |>
        fmt_number (columns = c ("Estimate", "S.E.", "t value"), decimals = 4) |>
        fmt_number (columns = "p value", decimals = 4) |>
        tab_header (title = "Coefficient Estimates") |>
        tab_options (table.width = pct(70))

    fstat <- s$fstatistic
    fit_df <- data.frame (
        R.squared     = s$r.squared,
        Adj.R.squared = s$adj.r.squared,
        Sigma         = s$sigma,
        F.statistic   = unname (fstat["value"]),
        p.value       = unname (pf (fstat["value"], fstat["numdf"], fstat["dendf"], lower.tail = FALSE))
    )

    fit_gt <- fit_df |>
        gt () |>
        cols_label (
            R.squared     = md ("$R^2$"),
            Adj.R.squared = md ("Adj. $R^2$"),
            Sigma         = md ("$\\hat\\sigma$"),
            F.statistic   = "F-statistic",
            p.value       = "p-value"
        ) |>
        fmt_number (columns = everything (), decimals = 4) |>
        tab_header (title = "Model Fit Summary") |>
        tab_options (table.width = pct(70))

    list (coef = coef_gt, fit = fit_gt)
}

# ---------------------------------------------------------
# Base Estimators
# ---------------------------------------------------------

## sdata --- a vector of original survey data
## N --- population size
## estimate --- "mean" (default) returns the population MEAN estimate;
##              "total" returns the population TOTAL estimate (N * mean,
##              with SE/CI scaled accordingly); requires finite N
## show.details --- if TRUE (default), also return a gt table of the
##                  row-level working values y_i, y_i-\bar y, (y_i-\bar y)^2,
##                  with a Sum row and a footnote showing s_y^2; NULL otherwise
## returns list ($estimate = c(Est., S.E., ci.low, ci.upp), $table)
srs_est <- function (sdata, N = Inf, estimate = c ("mean", "total"),
                      show.details = TRUE, col_labels = list (y = "y_i"))
{
    estimate <- match.arg (estimate)

    n <- length (sdata)
    ybar <- mean (sdata)
    dev <- sdata - ybar
    s2y <- sum (dev^2) / (n - 1)
    se.ybar <- sqrt((1 - n / N)) * sd (sdata) / sqrt(n)

    if (estimate == "total")
    {
        if (!is.finite (N))
            stop ("N must be finite to compute estimate = \"total\"")
        est <- N * ybar
        se.est <- N * se.ybar
    }
    else
    {
        est <- ybar
        se.est <- se.ybar
    }
    mem <- qt (0.975, df = n - 1) * se.est
    estimate_vec <- c (Est. = est, S.E. = se.est, ci.low = est - mem, ci.upp = est + mem)

    if (!show.details)
        return (list (estimate = estimate_vec, table = NULL))

    working <- data.frame (
        i    = c (seq_len (n), "Sum"),
        y    = c (sdata, sum (sdata)),
        dev  = c (dev, sum (dev)),
        dev2 = c (dev^2, sum (dev^2))
    )
    working <- truncate_working (working, id_col = "i")

    table <- gt (working) |>
        tab_header (title = "SRS Mean Estimation: Working Table", subtitle = md (sprintf ("$n = %d$", n))) |>
        cols_label (
            i    = md ("$i$"),
            y    = md (sprintf ("$%s$", col_labels$y)),
            dev  = md (sprintf ("$%s-\\bar y$", col_labels$y)),
            dev2 = md (sprintf ("$(%s-\\bar y)^2$", col_labels$y))
        ) |>
        cols_width (i ~ px (40), everything () ~ px (110)) |>
        fmt_auto (data = working, columns = c ("y", "dev", "dev2"), decimals = 3) |>
        sub_missing (missing_text = "...") |>
        tab_style (
            style     = cell_text (weight = "bold"),
            locations = cells_body (rows = i == "Sum")
        ) |>
        tab_footnote (
            footnote  = md (sprintf ("$s_y^2 = \\sum_i (%s-\\bar y)^2/(n-1) = %.4f$.", col_labels$y, s2y)),
            locations = cells_column_labels (columns = dev2)
        ) |>
        tab_source_note (
            source_note = md (if (estimate == "total")
                sprintf ("**Point estimate:** $\\hat t = N\\bar y = %s \\times %s = %s$",
                         fmt_num (N, 0), fmt_num (ybar), fmt_num (est))
            else
                sprintf ("**Point estimate:** $\\bar y = \\dfrac{\\sum_i %s}{n} = \\dfrac{%s}{%d} = %s$",
                         col_labels$y, fmt_num (sum (sdata)), n, fmt_num (ybar)))
        ) |>
        tab_source_note (
            source_note = md (if (estimate == "total")
                sprintf ("**Standard error:** $\\mathrm{SE}(\\hat t) = N \\times \\mathrm{SE}(\\bar y) = %s \\times %s = %s$",
                         fmt_num (N, 0), fmt_num (se.ybar), fmt_num (se.est))
            else if (is.finite (N))
                sprintf ("**Standard error:** $\\mathrm{SE}(\\bar y) = \\sqrt{1-\\dfrac{n}{N}}\\,\\dfrac{s}{\\sqrt n} = \\sqrt{1-\\dfrac{%d}{%d}}\\,\\dfrac{%s}{\\sqrt{%d}} = %s$",
                         n, N, fmt_num (sqrt (s2y)), n, fmt_num (se.ybar))
            else
                sprintf ("**Standard error:** $\\mathrm{SE}(\\bar y) = \\dfrac{s}{\\sqrt n} = \\dfrac{%s}{\\sqrt{%d}} = %s$",
                         fmt_num (sqrt (s2y)), n, fmt_num (se.ybar)))
        )

    list (estimate = estimate_vec, table = table)
}

## ydata --- observations of the variable of interest
## xdata --- observations of the auxiliary variable
## xbarU --- population mean of the auxiliary variable; required when
##           estimate is "mean" or "total"; ignored when estimate = "model"
## N --- population size
## estimate --- "mean" (default) returns the ratio estimate of the
##              population MEAN of y, Bhat*xbarU; "total" returns the ratio
##              estimate of the population TOTAL of y (requires finite N);
##              "model" returns the raw ratio coefficient Bhat = ybar/xbar
##              itself (with its own SE)
## show.details --- if TRUE (default), also return a gt table of the
##                  row-level working values y_i, x_i, fitted yhat_i =
##                  Bhat*x_i, residual e_i, and e_i^2, with a Sum row and a
##                  footnote showing Bhat and s_e^2; NULL otherwise
## col_labels --- named list (y, x, yhat) of RAW latex for those
##                working-table columns, so this function can be reused
##                (e.g. by cluster_ratio(), upswr_ratio()) with
##                context-appropriate labels
## extra_col --- optional list (label, values, formula) adding one extra
##               RAW-latex-labelled column of per-unit `values` right before
##               the y column of the working table (e.g. the within-cluster
##               means $\bar y_i$ used by cluster_ratio() to build $\hat
##               t_i = \bar y_i M_i$); `formula` (optional) is shown as a
##               footnote on the y column documenting that relationship
## returns list ($estimate = c(Est., S.E., ci.low, ci.upp), $table)
ratio_est <- function (ydata, xdata, xbarU = NULL, N = Inf,
                        estimate = c ("mean", "total", "model"), show.details = TRUE,
                        col_labels = list (y = "y_i", x = "x_i", yhat = "\\hat y_i"),
                        extra_col = NULL)
{
  estimate <- match.arg (estimate)

  n <- length (xdata)
  xbar <- mean (xdata)
  ybar <- mean (ydata)
  B_hat <- ybar / xbar
  yhat <- B_hat * xdata
  e <- ydata - yhat
  var_e <- sum (e^2) / (n - 1)
  sd_B_hat <- sqrt ((1 - n/N) * var_e / n) / xbar

  if (estimate == "model")
  {
      est <- B_hat
      sd_est <- sd_B_hat
  }
  else
  {
      if (is.null (xbarU))
          stop ("xbarU (population mean of x) is required when estimate is \"mean\" or \"total\"")
      if (estimate == "total")
      {
          if (!is.finite (N))
              stop ("N must be finite to compute estimate = \"total\"")
          est <- B_hat * xbarU * N
          sd_est <- sd_B_hat * xbarU * N
      }
      else
      {
          est <- B_hat * xbarU
          sd_est <- sd_B_hat * xbarU
      }
  }

  mem <- qt (0.975, df = n - 1) * sd_est
  estimate_vec <- c (Est. = est, S.E. = sd_est, ci.low = est - mem, ci.upp = est + mem)

  if (!show.details)
      return (list (estimate = estimate_vec, table = NULL))

  working_cols <- list (i = c (seq_len (n), "Sum"))
  if (!is.null (extra_col))
      working_cols$extra <- c (extra_col$values, NA)
  working_cols$y    <- c (ydata, sum (ydata))
  working_cols$x    <- c (xdata, sum (xdata))
  working_cols$yhat <- c (yhat, sum (yhat))
  working_cols$e    <- c (e, sum (e))
  working_cols$e2   <- c (e^2, sum (e^2))
  working <- as.data.frame (working_cols, check.names = FALSE)
  working <- truncate_working (working, id_col = "i")

  fmt_cols <- c ("y", "x", "yhat", "e", "e2")
  if (!is.null (extra_col)) fmt_cols <- c ("extra", fmt_cols)

  label_args <- list (i = md ("$i$"))
  if (!is.null (extra_col))
      label_args$extra <- md (sprintf ("$%s$", extra_col$label))
  label_args$y    <- md (sprintf ("$%s$", col_labels$y))
  label_args$x    <- md (sprintf ("$%s$", col_labels$x))
  label_args$yhat <- md (sprintf ("$%s$", col_labels$yhat))
  label_args$e    <- md ("$e_i$")
  label_args$e2   <- md ("$e_i^2$")

  table <- gt (working) |>
      tab_header (title = "Ratio Estimation: Working Table", subtitle = md (sprintf ("$n = %d$", n)))
  table <- do.call (cols_label, c (list (table), label_args))
  table <- table |>
      cols_width (i ~ px (40), everything () ~ px (100)) |>
      fmt_auto (data = working, columns = fmt_cols, decimals = 3) |>
      sub_missing (missing_text = "...") |>
      tab_style (
          style     = cell_text (weight = "bold"),
          locations = cells_body (rows = i == "Sum")
      ) |>
      tab_footnote (
          footnote  = md (sprintf ("$\\hat B = \\sum_i %s / \\sum_i %s = %.4f$, $%s = \\hat B\\, %s$, $s_e^2 = %.4f$.",
                                    col_labels$y, col_labels$x, B_hat, col_labels$yhat, col_labels$x, var_e)),
          locations = cells_column_labels (columns = yhat)
      )
  if (!is.null (extra_col) && !is.null (extra_col$formula))
      table <- table |>
          tab_footnote (
              footnote  = md (sprintf ("$%s$.", extra_col$formula)),
              locations = cells_column_labels (columns = y)
          )
  table <- table |>
      tab_source_note (
          source_note = md (if (estimate == "model")
              sprintf ("**Point estimate:** $\\hat B = \\dfrac{\\sum_i %s}{\\sum_i %s} = \\dfrac{%s}{%s} = %s$",
                       col_labels$y, col_labels$x, fmt_num (sum (ydata)), fmt_num (sum (xdata)), fmt_num (B_hat))
          else if (estimate == "total")
              sprintf ("**Point estimate:** $\\hat t = \\hat B\\,\\bar x_U\\,N = %s \\times %s \\times %s = %s$",
                       fmt_num (B_hat), fmt_num (xbarU), fmt_num (N, 0), fmt_num (est))
          else
              sprintf ("**Point estimate:** $\\bar y_r = \\hat B\\,\\bar x_U = %s \\times %s = %s$",
                       fmt_num (B_hat), fmt_num (xbarU), fmt_num (est)))
      ) |>
      tab_source_note (
          source_note = md (if (estimate == "model")
              (if (is.finite (N))
                  sprintf ("**Standard error:** $\\mathrm{SE}(\\hat B) = \\dfrac{1}{\\bar x}\\sqrt{\\left(1-\\dfrac{n}{N}\\right)\\dfrac{s_e^2}{n}} = \\dfrac{1}{%s}\\sqrt{\\left(1-\\dfrac{%d}{%d}\\right)\\dfrac{%s}{%d}} = %s$",
                          fmt_num (xbar), n, N, fmt_num (var_e), n, fmt_num (sd_B_hat))
              else
                  sprintf ("**Standard error:** $\\mathrm{SE}(\\hat B) = \\dfrac{1}{\\bar x}\\sqrt{\\dfrac{s_e^2}{n}} = \\dfrac{1}{%s}\\sqrt{\\dfrac{%s}{%d}} = %s$",
                          fmt_num (xbar), fmt_num (var_e), n, fmt_num (sd_B_hat)))
          else if (estimate == "total")
              sprintf ("**Standard error:** $\\mathrm{SE}(\\hat t) = \\mathrm{SE}(\\hat B)\\,\\bar x_U\\,N = %s \\times %s \\times %s = %s$",
                       fmt_num (sd_B_hat), fmt_num (xbarU), fmt_num (N, 0), fmt_num (sd_est))
          else
              sprintf ("**Standard error:** $\\mathrm{SE}(\\bar y_r) = \\mathrm{SE}(\\hat B)\\,\\bar x_U = %s \\times %s = %s$",
                       fmt_num (sd_B_hat), fmt_num (xbarU), fmt_num (sd_est)))
      )

  list (estimate = estimate_vec, table = table)
}

## ydata --- observations of the variable of interest
## xdata --- observations of the auxiliary variable
## xbarU --- population mean of the auxiliary variable; required when
##           estimate is "mean" or "total"; ignored when estimate = "model"
## N --- population size
## estimate --- "mean" (default) returns the regression estimate of the
##              population MEAN of y, Bhat0 + Bhat1*xbarU; "total" returns
##              the regression estimate of the population TOTAL of y
##              (requires finite N); "model" returns the raw fitted slope
##              Bhat1 itself (with its own SE from the lm fit)
## show.details --- if TRUE (default), also return a gt table of the
##                  row-level working values y_i, x_i, fitted yhat_i, e_i,
##                  e_i^2, with a Sum row and a footnote showing Bhat0,
##                  Bhat1, and s_e^2; NULL otherwise
## returns list ($estimate = c(Est., S.E., ci.low, ci.upp), $table)
reg_est <- function (ydata, xdata, xbarU = NULL, N = Inf,
                      estimate = c ("mean", "total", "model"), show.details = TRUE)
{
  estimate <- match.arg (estimate)

  n <- length (ydata)
  lmfit <- lm (ydata ~ xdata)
  Bhat <- lmfit$coefficients
  yhat <- lmfit$fitted.values
  e <- lmfit$residuals
  SSe <- sum (e^2) / (n - 2)

  if (estimate == "model")
  {
      est <- unname (Bhat[2])
      sd_est <- summary (lmfit)$coefficients[2, "Std. Error"]
  }
  else
  {
      if (is.null (xbarU))
          stop ("xbarU (population mean of x) is required when estimate is \"mean\" or \"total\"")
      yhat_reg <- unname (Bhat[1] + Bhat[2] * xbarU)
      se_yhat_reg <- sqrt ((1-n/N) * SSe / n)
      if (estimate == "total")
      {
          if (!is.finite (N))
              stop ("N must be finite to compute estimate = \"total\"")
          est <- yhat_reg * N
          sd_est <- se_yhat_reg * N
      }
      else
      {
          est <- yhat_reg
          sd_est <- se_yhat_reg
      }
  }

  mem <- qt (0.975, df = n - 2) * sd_est
  estimate_vec <- c (Est. = est, S.E. = sd_est, ci.low = est - mem, ci.upp = est + mem)

  if (!show.details)
      return (list (estimate = estimate_vec, table = NULL))

  Sxx <- sum ((xdata - mean (xdata))^2)
  Sxy <- sum ((xdata - mean (xdata)) * (ydata - mean (ydata)))

  working <- data.frame (
      i    = c (seq_len (n), "Sum"),
      y    = c (ydata, sum (ydata)),
      x    = c (xdata, sum (xdata)),
      yhat = c (yhat, sum (yhat)),
      e    = c (e, sum (e)),
      e2   = c (e^2, sum (e^2))
  )
  working <- truncate_working (working, id_col = "i")

  table <- gt (working) |>
      tab_header (title = "Regression Estimation: Working Table", subtitle = md (sprintf ("$n = %d$", n))) |>
      cols_label (
          i    = md ("$i$"),
          y    = md ("$y_i$"),
          x    = md ("$x_i$"),
          yhat = md ("$\\hat y_i$"),
          e    = md ("$e_i$"),
          e2   = md ("$e_i^2$")
      ) |>
      cols_width (i ~ px (40), everything () ~ px (100)) |>
      fmt_auto (data = working, columns = c ("y", "x", "yhat", "e", "e2"), decimals = 3) |>
      sub_missing (missing_text = "...") |>
      tab_style (
          style     = cell_text (weight = "bold"),
          locations = cells_body (rows = i == "Sum")
      ) |>
      tab_footnote (
          footnote  = md (sprintf ("$\\hat B_0 = %.4f$, $\\hat B_1 = %.4f$, $\\hat y_i = \\hat B_0 + \\hat B_1 x_i$, $s_e^2 = %.4f$.",
                                    Bhat[1], Bhat[2], SSe)),
          locations = cells_column_labels (columns = yhat)
      ) |>
      tab_source_note (
          source_note = md (if (estimate == "model")
              sprintf ("**Point estimate:** $\\hat B_1 = \\dfrac{\\sum_i(x_i-\\bar x)(y_i-\\bar y)}{\\sum_i(x_i-\\bar x)^2} = \\dfrac{%s}{%s} = %s$",
                       fmt_num (Sxy), fmt_num (Sxx), fmt_num (unname (Bhat[2])))
          else if (estimate == "total")
              sprintf ("**Point estimate:** $\\hat t_{reg} = N\\bar y_{reg} = %s \\times %s = %s$",
                       fmt_num (N, 0), fmt_num (yhat_reg), fmt_num (est))
          else
              sprintf ("**Point estimate:** $\\bar y_{reg} = \\hat B_0+\\hat B_1\\bar x_U = %s + %s \\times %s = %s$",
                       fmt_num (Bhat[1]), fmt_num (Bhat[2]), fmt_num (xbarU), fmt_num (yhat_reg)))
      ) |>
      tab_source_note (
          source_note = md (if (estimate == "model")
              sprintf ("**Standard error:** $\\mathrm{SE}(\\hat B_1) = \\sqrt{\\dfrac{s_e^2}{\\sum_i(x_i-\\bar x)^2}} = \\sqrt{\\dfrac{%s}{%s}} = %s$",
                       fmt_num (SSe), fmt_num (Sxx), fmt_num (sd_est))
          else if (estimate == "total")
              sprintf ("**Standard error:** $\\mathrm{SE}(\\hat t_{reg}) = N\\times\\mathrm{SE}(\\bar y_{reg}) = %s \\times %s = %s$",
                       fmt_num (N, 0), fmt_num (se_yhat_reg), fmt_num (sd_est))
          else if (is.finite (N))
              sprintf ("**Standard error:** $\\mathrm{SE}(\\bar y_{reg}) = \\sqrt{\\left(1-\\dfrac{n}{N}\\right)\\dfrac{s_e^2}{n}} = \\sqrt{\\left(1-\\dfrac{%d}{%d}\\right)\\dfrac{%s}{%d}} = %s$",
                       n, N, fmt_num (SSe), n, fmt_num (se_yhat_reg))
          else
              sprintf ("**Standard error:** $\\mathrm{SE}(\\bar y_{reg}) = \\sqrt{\\dfrac{s_e^2}{n}} = \\sqrt{\\dfrac{%s}{%d}} = %s$",
                       fmt_num (SSe), n, fmt_num (se_yhat_reg)))
      )

  list (estimate = estimate_vec, table = table)
}

## sdata --- a vector of original survey data in a domain
## N --- population size
## n --- total sample size (not the sample size in the domain)
## to find total, multiply domain size N_d to the estimate returned by this function
domain_est <- function (sdata, n, N = Inf)
{
    n_d <- length (sdata)
    ybar <- mean (sdata)
    se.ybar <- sqrt((1 - n / N)) * sd (sdata) / sqrt(n_d)
    mem <- qt (0.975, df = n_d - 1) * se.ybar
    c (Est. = ybar, S.E. = se.ybar, ci.low = ybar - mem, ci.upp = ybar + mem)
}

# ---------------------------------------------------------
# Stratified Estimator
# ---------------------------------------------------------

## Given per-stratum summaries, compute the stratified (or post-stratified)
## mean or total.
## estimate --- "mean" (default) returns the stratified MEAN estimate;
##              "total" returns the stratified TOTAL estimate
## show.details --- if TRUE (default), also return a gt table showing the
##                  stratum-by-stratum working, with a Total row
## for poststratification, use nh = n * Nh/N
## returns list ($estimate = c(Est., S.E., ci.low, ci.upp), $table)
str_est <- function (ybarh, sh, nh, Nh, estimate = c ("mean", "total"), show.details = TRUE)
{
    estimate <- match.arg (estimate)

    N    <- sum (Nh)
    Pi_h <- Nh / N
    v_h  <- (1 - nh / Nh) * Pi_h^2 * sh^2 / nh

    ybar   <- sum (ybarh * Pi_h)
    seybar <- sqrt (sum (v_h))

    if (estimate == "total")
    {
        est    <- N * ybar
        se.est <- N * seybar
    }
    else
    {
        est    <- ybar
        se.est <- seybar
    }
    mem          <- 1.96 * se.est
    estimate_vec <- c (Est. = est, S.E. = se.est, ci.low = est - mem, ci.upp = est + mem)

    if (!show.details)
        return (list (estimate = estimate_vec, table = NULL))

    stratum <- names (Nh)
    if (is.null (stratum)) stratum <- seq_along (Nh)

    working <- data.frame (
        Stratum       = c (stratum, "Total"),
        Nh            = c (Nh, N),
        nh            = c (nh, sum (nh)),
        Pi_h          = c (Pi_h, sum (Pi_h)),
        ybarh         = c (ybarh, NA),
        sh2           = c (sh^2, NA),
        weighted_mean = c (Pi_h * ybarh, ybar),
        v_h           = c (v_h, seybar^2)
    )

    table <- gt (working) |>
        tab_header (title = "Stratified Mean Estimation: Working Table") |>
        cols_label (
            Stratum       = md ("Stratum"),
            Nh            = md ("$N_h$"),
            nh            = md ("$n_h$"),
            Pi_h          = md ("$\\pi_h$"),
            ybarh         = md ("$\\bar y_h$"),
            sh2           = md ("$s_h^2$"),
            weighted_mean = md ("$\\pi_h\\bar y_h$"),
            v_h           = md ("$v_h$")
        ) |>
        cols_width (c (Nh, nh) ~ px (60), everything () ~ px (110)) |>
        fmt_number (columns = c (Nh, nh), decimals = 0) |>
        fmt_number (columns = Pi_h, decimals = 2) |>
        fmt_auto (data = working, columns = c ("ybarh", "sh2", "weighted_mean", "v_h"), decimals = 2) |>
        sub_missing (missing_text = "") |>
        tab_style (
            style     = cell_text (weight = "bold"),
            locations = cells_body (rows = Stratum == "Total")
        ) |>
        tab_footnote (
            footnote  = md ("$\\pi_h = N_h/N$ is the stratum weight."),
            locations = cells_column_labels (columns = Pi_h)
        ) |>
        tab_footnote (
            footnote  = md ("$v_h = (1-n_h/N_h)\\,\\pi_h^2\\, s_h^2/n_h$ is stratum $h$'s contribution to $V(\\bar y_{str}) = \\sum_h v_h$."),
            locations = cells_column_labels (columns = v_h)
        ) |>
        tab_source_note (
            source_note = md (if (estimate == "total")
                sprintf ("**Point estimate:** $\\hat t_{str} = N\\bar y_{str} = %s \\times %s = %s$",
                         fmt_num (N, 0), fmt_num (ybar), fmt_num (est))
            else
                sprintf ("**Point estimate:** $\\bar y_{str} = \\sum_h \\pi_h\\bar y_h = %s$", fmt_num (ybar)))
        ) |>
        tab_source_note (
            source_note = md (if (estimate == "total")
                sprintf ("**Standard error:** $\\mathrm{SE}(\\hat t_{str}) = N \\times \\mathrm{SE}(\\bar y_{str}) = %s \\times %s = %s$",
                         fmt_num (N, 0), fmt_num (seybar), fmt_num (se.est))
            else
                sprintf ("**Standard error:** $\\mathrm{SE}(\\bar y_{str}) = \\sqrt{\\sum_h v_h} = \\sqrt{%s} = %s$",
                         fmt_num (sum (v_h)), fmt_num (seybar)))
        )

    list (estimate = estimate_vec, table = table)
}

## this function finds statistical estimates given a dataset with sampling weight
# stratdata --- data.frame containing stratified sample
# y --- name of variable for which we want to estiamte population mean
# stratum --- name of variable that will be used as stratum variable
# weight --- name of variable indicating sampling weight
# note: from weights we can find Nh (see the code for formula)
# returns the same list ($estimate, $table) as str_est(), which
# this function calls after deriving ybarh/sh/nh/Nh from the raw data
str_est_data <- function (stratdata, y, stratum, weight, estimate = c ("mean", "total"), show.details = TRUE)
{
    ## compute stratum-wise data
    sh <- tapply (stratdata[, y], stratdata[,stratum], sd)
    ybarh <- tapply (stratdata[, y], stratdata[,stratum], mean)
    ## find population stratum size using sampling weight included in the data set
    Nh <- tapply (stratdata[, weight], stratdata[,stratum], sum)
    nh <- tapply (1:nrow(stratdata), stratdata[,stratum], length)

    str_est (ybarh, sh, nh, Nh, estimate = estimate, show.details = show.details)
}

# ---------------------------------------------------------
# Cluster (ratio-to-size) Estimator
# ---------------------------------------------------------

## Cluster ratio estimation (one-stage if Mi = mi, i.e. every element in the
## sampled cluster is measured; two-stage/ratio-to-size if only mi < Mi
## elements are subsampled, so t_hat_i = Mi * ybari is an ESTIMATED total).
## data --- original data frame for holding data
## cname --- variable recording cluster (psu) identity
## csize --- variable recording cluster (psu) population size (not sample size)
## yvar --- variable of interest
## N --- total number of clusters (psus) in the population
## estimate --- "mean" (default) or "model" both return the raw ratio
##              Bhat = sum(t_hat_i)/sum(Mi) --- for this ratio-to-size
##              cluster estimator, that raw B IS the population MEAN per
##              element, so the two coincide; "total" returns the
##              population TOTAL of y, Bhat*Mtotal_U (requires Mtotal_U,
##              the population total of the cluster-size variable)
## Mtotal_U --- population total of the cluster-size variable; required
##              when estimate = "total"
## show.details --- if TRUE, also return a gt table of the cluster-level
##                  working values (t_hat_i, M_i, fitted, e_i, e_i^2), via
##                  ratio_est()'s own working table
## returns list ($estimate = c(Est., S.E., ci.low, ci.upp), $table)
cluster_ratio <- function (data, cname, csize, yvar, N = Inf,
                            estimate = c ("mean", "total", "model"),
                            show.details = TRUE, Mtotal_U = NULL)
{
  estimate <- match.arg (estimate)

  clust <- data[, cname]
  ydata <- data[, yvar]

  ybari <- tapply (ydata, clust, mean)
  Mi    <- tapply (data[, csize], clust, function (x) x[1])
  ## same as the cluster total if all Mi elements were measured (mi = Mi)
  t_hat_cls <- ybari * Mi

  ybar_col <- list (label = "\\bar y_i", values = ybari, formula = "\\hat t_i = \\bar y_i \\, M_i")

  if (estimate == "total")
  {
      if (is.null (Mtotal_U))
          stop ("Mtotal_U (population total of the cluster-size variable) is required when estimate = \"total\"")
      ratio_est (t_hat_cls, Mi, xbarU = Mtotal_U, N = N, estimate = "mean", show.details = show.details,
                 col_labels = list (y = "\\hat t_i", x = "M_i", yhat = "\\hat B M_i"),
                 extra_col = ybar_col)
  }
  else
  {
      ratio_est (t_hat_cls, Mi, N = N, estimate = "model", show.details = show.details,
                 col_labels = list (y = "\\hat t_i", x = "M_i", yhat = "\\hat B M_i"),
                 extra_col = ybar_col)
  }
}

# ---------------------------------------------------------
# UPS Estimators
# ---------------------------------------------------------

## Total estimation for datasets collected with UPS with Replacement
## (psi-estimator / Hansen-Hurwitz).
## estimate --- "total" (default) is this function's whole purpose, computed
##              as the SAMPLE MEAN of the expanded pseudo-values t_i/psi_i
##              (a with-replacement HT-type property); "mean" instead
##              divides that total by N (requires finite N)
upswr_total <- function (total, psi, N = Inf, estimate = c ("total", "mean"), show.details = TRUE)
{
  estimate <- match.arg (estimate)
  res <- srs_est (total/psi, N = N, estimate = "mean", show.details = show.details,
                  col_labels = list (y = "t_i/\\psi_i"))
  if (estimate == "mean")
  {
      if (!is.finite (N)) stop ("N must be finite to compute estimate = \"mean\" (mean total per psu)")
      res$estimate <- res$estimate / N
  }
  res
}

## Ratio estimation for datasets collected with UPS with Replacement.
## estimate --- "mean" (default) or "model" both return the raw ratio
##              Bhat = sum(t_i/psi_i)/sum(M_i/psi_i) --- that raw B IS the
##              population MEAN per element; "total" returns the population
##              TOTAL of y, Bhat*Mtotal_U (requires Mtotal_U, the population
##              total of the size variable)
upswr_ratio <- function (total, M, psi, N = Inf, estimate = c ("mean", "total", "model"),
                          show.details = TRUE, Mtotal_U = NULL)
{
  estimate <- match.arg (estimate)

  if (estimate == "total")
  {
      if (is.null (Mtotal_U))
          stop ("Mtotal_U (population total of the size variable) is required when estimate = \"total\"")
      ratio_est (total/psi, M/psi, xbarU = Mtotal_U, N = N, estimate = "mean", show.details = show.details,
                 col_labels = list (y = "t_i/\\psi_i", x = "M_i/\\psi_i", yhat = "\\hat B\\,M_i/\\psi_i"))
  }
  else
  {
      ratio_est (total/psi, M/psi, N = N, estimate = "model", show.details = show.details,
                 col_labels = list (y = "t_i/\\psi_i", x = "M_i/\\psi_i", yhat = "\\hat B\\,M_i/\\psi_i"))
  }
}

## Cluster UPS with Replacement: aggregates to the cluster level (ybari, Mi,
## psi per sampled cluster), then delegates to upswr_ratio(), so the
## working table and estimate/total logic are shared with the base function.
## data --- data frame holding the sample
## sid --- variable recording cluster (psu) identity
## csize --- variable recording cluster (psu) population size
## cpsi --- variable recording the cluster's selection probability psi_i
## yvar --- variable of interest (element level)
## estimate, Mtotal_U --- see upswr_ratio()
cluster_upswr_ratio <- function (data, sid, csize, cpsi, yvar, N = Inf,
                                  estimate = c ("mean", "total", "model"),
                                  show.details = TRUE, Mtotal_U = NULL)
{
    estimate <- match.arg (estimate)

    clust <- data[, sid]
    ydata <- data[, yvar]

    ybari <- tapply (ydata, clust, mean)
    Mi <- tapply (data[,csize], clust, function (x) x[1])
    psi <- tapply (data[,cpsi], clust, function (x) x[1])
    t_hat_cls <- ybari * Mi

    upswr_ratio (t_hat_cls, Mi, psi, N = N, estimate = estimate,
                 show.details = show.details, Mtotal_U = Mtotal_U)
}

## Cluster UPS without Replacement: same idea as cluster_upswr_ratio(),
## but dividing by the cluster's inclusion probability pik instead of psi.
## data --- data frame holding the sample
## cname --- variable recording cluster identity
## csize --- variable recording cluster population size
## cpik --- variable recording the cluster's inclusion probability pi_k
## yvar --- variable of interest (element level)
cluster_upswo_ratio <- function (data, cname, csize, cpik, yvar, N = Inf,
                                  estimate = c ("mean", "total", "model"),
                                  show.details = TRUE, Mtotal_U = NULL)
{
  estimate <- match.arg (estimate)

  clust <- data[, cname]
  ydata <- data[, yvar]

  ybari <- tapply (ydata, clust, mean)
  Mi <- tapply (data[,csize], clust, function (x) x[1])
  pik <- tapply (data[,cpik], clust, function (x) x[1])
  t_hat_cls <- ybari * Mi

  if (estimate == "total")
  {
      if (is.null (Mtotal_U))
          stop ("Mtotal_U (population total of the size variable) is required when estimate = \"total\"")
      ratio_est (t_hat_cls/pik, Mi/pik, xbarU = Mtotal_U, N = N, estimate = "mean", show.details = show.details,
                 col_labels = list (y = "\\hat t_i/\\pi_i", x = "M_i/\\pi_i", yhat = "\\hat B\\,M_i/\\pi_i"))
  }
  else
  {
      ratio_est (t_hat_cls/pik, Mi/pik, N = N, estimate = "model", show.details = show.details,
                 col_labels = list (y = "\\hat t_i/\\pi_i", x = "M_i/\\pi_i", yhat = "\\hat B\\,M_i/\\pi_i"))
  }
}

```

</details>

:::

## Analysis of a dataset on student study hours collected by UPSWR

We first manually define a small dataset describing hours spent studying, where the probability of sampling a specific student group ($\psi_i$) was proportional to the number of students in that group ($M_i$).


::: {.cell}
::: {.cell-output-display}
\begin{table}
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}rrr}
\toprule
\(M_i\) & \(\psi_i\) & \textbf{Total Hours (\(y_i\))} \\ 
\midrule\addlinespace[2.5pt]
24 & 0.037 & 75 \\ 
100 & 0.155 & 203 \\ 
100 & 0.155 & 203 \\ 
76 & 0.117 & 191 \\ 
44 & 0.068 & 168 \\ 
\bottomrule
\end{tabular*}
\end{table}

:::

::: {.cell-output-display}
![](ups_files/figure-pdf/studyhrs-data-1.pdf)
:::
:::


**Estimate the total and mean study hours:**


::: {.cell}
::: {.cell-output-display}
\begin{table}
\caption*{
{\fontsize{20}{25}\selectfont  SRS Mean Estimation: Working Table\fontsize{12}{15}\selectfont } \\ 
{\fontsize{14}{17}\selectfont  \(n = 5\)\fontsize{12}{15}\selectfont }
} 
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}>{\raggedright\arraybackslash}p{\dimexpr 30.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 82.50pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 82.50pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 82.50pt -2\tabcolsep-1.5\arrayrulewidth}}
\toprule
\(i\) & \(t_i/\psi_i\) & \(t_i/\psi_i-\bar y\) & \((t_i/\psi_i-\bar y)^2\)\textsuperscript{\textit{1}} \\ 
\midrule\addlinespace[2.5pt]
1 & 2,021.875 & 2.729 $\times$ 10\textsuperscript{2} & 7.445 $\times$ 10\textsuperscript{4} \\ 
2 & 1,313.410 & -4.356 $\times$ 10\textsuperscript{2} & 1.898 $\times$ 10\textsuperscript{5} \\ 
3 & 1,313.410 & -4.356 $\times$ 10\textsuperscript{2} & 1.898 $\times$ 10\textsuperscript{5} \\ 
4 & 1,626.013 & -1.230 $\times$ 10\textsuperscript{2} & 1.513 $\times$ 10\textsuperscript{4} \\ 
5 & 2,470.364 & 7.213 $\times$ 10\textsuperscript{2} & 5.203 $\times$ 10\textsuperscript{5} \\ 
{\bfseries Sum} & {\bfseries 8,745.072} & {\bfseries -4.547 $\times$ 10\textsuperscript{-13}} & {\bfseries 9.894 $\times$ 10\textsuperscript{5}} \\ 
\bottomrule
\end{tabular*}
\begin{minipage}{\linewidth}
\vspace{.05em}
\textsuperscript{\textit{1}}\(s_y^2 = \sum_i (t_i/\psi_i-\bar y)^2/(n-1) = 247357.3300\).\\
\textbf{Point estimate:} \(\bar y = \dfrac{\sum_i t_i/\psi_i}{n} = \dfrac{8,745.072}{5} = 1,749.014\)\\
\textbf{Standard error:} \(\mathrm{SE}(\bar y) = \dfrac{s}{\sqrt n} = \dfrac{497.350}{\sqrt{5}} = 222.422\)\\
\end{minipage}
\end{table}

:::

::: {.cell-output-display}
\begin{table}
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{0.7\linewidth}{@{\extracolsep{\fill}}rrrr}
\toprule
\(\hat{t}\) & \(\mathrm{SE}(\hat{t})\) & 95\% CI Lower & 95\% CI Upper \\ 
\midrule\addlinespace[2.5pt]
1,749.0144 & 222.4218 & 1,131.4724 & 2,366.5563 \\ 
\bottomrule
\end{tabular*}
\end{table}

:::

::: {.cell-output-display}
\begin{table}
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{0.7\linewidth}{@{\extracolsep{\fill}}rrrr}
\toprule
\(\bar{y}\) & \(\mathrm{SE}(\bar{y})\) & 95\% CI Lower & 95\% CI Upper \\ 
\midrule\addlinespace[2.5pt]
2.7033 & 0.3438 & 1.7488 & 3.6577 \\ 
\bottomrule
\end{tabular*}
\end{table}

:::

::: {.cell-output-display}
\begin{table}
\caption*{
{\fontsize{20}{25}\selectfont  Ratio Estimation: Working Table\fontsize{12}{15}\selectfont } \\ 
{\fontsize{14}{17}\selectfont  \(n = 5\)\fontsize{12}{15}\selectfont }
} 
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}>{\raggedright\arraybackslash}p{\dimexpr 30.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}}
\toprule
\(i\) & \(t_i/\psi_i\) & \(M_i/\psi_i\) & \(\hat B\,M_i/\psi_i\)\textsuperscript{\textit{1}} & \(e_i\) & \(e_i^2\) \\ 
\midrule\addlinespace[2.5pt]
1 & 2,021.875 & 647.000 & 1,749.014 & 2.729 $\times$ 10\textsuperscript{2} & 7.445 $\times$ 10\textsuperscript{4} \\ 
2 & 1,313.410 & 647.000 & 1,749.014 & -4.356 $\times$ 10\textsuperscript{2} & 1.898 $\times$ 10\textsuperscript{5} \\ 
3 & 1,313.410 & 647.000 & 1,749.014 & -4.356 $\times$ 10\textsuperscript{2} & 1.898 $\times$ 10\textsuperscript{5} \\ 
4 & 1,626.013 & 647.000 & 1,749.014 & -1.230 $\times$ 10\textsuperscript{2} & 1.513 $\times$ 10\textsuperscript{4} \\ 
5 & 2,470.364 & 647.000 & 1,749.014 & 7.213 $\times$ 10\textsuperscript{2} & 5.203 $\times$ 10\textsuperscript{5} \\ 
{\bfseries Sum} & {\bfseries 8,745.072} & {\bfseries 3,235.000} & {\bfseries 8,745.072} & {\bfseries -4.547 $\times$ 10\textsuperscript{-13}} & {\bfseries 9.894 $\times$ 10\textsuperscript{5}} \\ 
\bottomrule
\end{tabular*}
\begin{minipage}{\linewidth}
\vspace{.05em}
\textsuperscript{\textit{1}}\(\hat B = \sum_i t_i/\psi_i / \sum_i M_i/\psi_i = 2.7033\), \(\hat B\,M_i/\psi_i = \hat B\, M_i/\psi_i\), \(s_e^2 = 247357.3300\).\\
\textbf{Point estimate:} \(\hat B = \dfrac{\sum_i t_i/\psi_i}{\sum_i M_i/\psi_i} = \dfrac{8,745.072}{3,235.000} = 2.703\)\\
\textbf{Standard error:} \(\mathrm{SE}(\hat B) = \dfrac{1}{\bar x}\sqrt{\dfrac{s_e^2}{n}} = \dfrac{1}{647.000}\sqrt{\dfrac{2.474 \times 10^{5}}{5}} = 0.344\)\\
\end{minipage}
\end{table}

:::

::: {.cell-output-display}
\begin{table}
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{0.7\linewidth}{@{\extracolsep{\fill}}rrrr}
\toprule
\(\bar{y}_{r}\) & \(\mathrm{SE}(\bar{y}_{r})\) & 95\% CI Lower & 95\% CI Upper \\ 
\midrule\addlinespace[2.5pt]
2.7033 & 0.3438 & 1.7488 & 3.6577 \\ 
\bottomrule
\end{tabular*}
\end{table}

:::
:::


## Analysis of `statepop.csv` collected with one-stage UPSWR

The file `statepop.csv` contains data from an unequal-probability sample of 100 U.S. counties. Counties were chosen with probabilities proportional to their populations. The total U.S. population is $M_0 = 255,077,536$. Because sampling was done **with replacement**, large counties like Los Angeles occur multiple times in the sample.


::: {.cell}
::: {.cell-output-display}
\begin{table}
\caption*{
{\fontsize{20}{25}\selectfont  U.S. County Data (Sample)\fontsize{12}{15}\selectfont }
} 
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}llrrr}
\toprule
state & county & popn & phys & psi \\ 
\midrule\addlinespace[2.5pt]
AL & Wilcox & 13672 & 4 & 5.359939e-05 \\ 
AZ & Maricopa & 2209567 & 4320 & 8.662335e-03 \\ 
AZ & Maricopa & 2209567 & 4320 & 8.662335e-03 \\ 
AZ & Pinal & 120786 & 61 & 4.735266e-04 \\ 
AR & Garland & 76100 & 131 & 2.983407e-04 \\ 
AR & Mississippi & 55060 & 48 & 2.158559e-04 \\ 
CA & Contra\_Costa & 840585 & 1761 & 3.295410e-03 \\ 
CA & Kern & 587680 & 682 & 2.303927e-03 \\ 
CA & Los\_Angeles & 9053645 & 23677 & 3.549370e-02 \\ 
CA & Los\_Angeles & 9053645 & 23677 & 3.549370e-02 \\ 
CA & Los\_Angeles & 9053645 & 23677 & 3.549370e-02 \\ 
CA & Los\_Angeles & 9053645 & 23677 & 3.549370e-02 \\ 
CA & Merced & 189107 & 185 & 7.413707e-04 \\ 
CA & Orange & 2484789 & 6062 & 9.741309e-03 \\ 
CA & Riverside & 1288435 & 1385 & 5.051150e-03 \\ 
CA & San\_Francisco & 728921 & 4761 & 2.857645e-03 \\ 
CA & Ventura & 686560 & 1168 & 2.691574e-03 \\ 
CO & Arapahoe & 420882 & 583 & 1.650016e-03 \\ 
FL & Alachua & 189409 & 1180 & 7.425546e-04 \\ 
FL & Lee & 352051 & 509 & 1.380172e-03 \\ 
FL & Orange & 714579 & 1367 & 2.801419e-03 \\ 
FL & Pinellas & 854976 & 1620 & 3.351828e-03 \\ 
GA & Forsyth & 49859 & 21 & 1.954661e-04 \\ 
IL & Cook & 5139341 & 15153 & 2.014815e-02 \\ 
IL & Cook & 5139341 & 15153 & 2.014815e-02 \\ 
IL & Cook & 5139341 & 15153 & 2.014815e-02 \\ 
IL & DuPage & 816116 & 2157 & 3.199482e-03 \\ 
IL & Kane & 333626 & 473 & 1.307940e-03 \\ 
IL & Kane & 333626 & 473 & 1.307940e-03 \\ 
IL & Lake & 541047 & 1093 & 2.121108e-03 \\ 
IL & St.\_Clair & 263124 & 329 & 1.031545e-03 \\ 
IL & Williamson & 58449 & 63 & 2.291421e-04 \\ 
IN & Bartholomew & 65498 & 102 & 2.567768e-04 \\ 
KS & Douglas & 84538 & 107 & 3.314208e-04 \\ 
KY & Boone & 63107 & 62 & 2.474032e-04 \\ 
KY & Jefferson & 670837 & 2171 & 2.629934e-03 \\ 
KY & Jefferson & 670837 & 2171 & 2.629934e-03 \\ 
LA & Jefferson & 457738 & 1237 & 1.794505e-03 \\ 
LA & St.\_Charles & 44372 & 18 & 1.739549e-04 \\ 
MD & Frederick & 159738 & 172 & 6.262331e-04 \\ 
MD & Worcester & 37306 & 24 & 1.462536e-04 \\ 
MI & Oakland & 1118611 & 4020 & 4.385376e-03 \\ 
MI & Oakland & 1118611 & 4020 & 4.385376e-03 \\ 
MN & Hennepin & 1041332 & 3706 & 4.082414e-03 \\ 
MN & Otter\_Tail & 51362 & 55 & 2.013584e-04 \\ 
MS & Leflore & 37477 & 42 & 1.469240e-04 \\ 
MS & Pearl\_River & 39879 & 22 & 1.563407e-04 \\ 
MO & Jasper & 91895 & 135 & 3.602630e-04 \\ 
MO & Nodaway & 21236 & 15 & 8.325312e-05 \\ 
NE & Lancaster & 219582 & 389 & 8.608441e-04 \\ 
NJ & Atlantic & 229430 & 379 & 8.994520e-04 \\ 
NJ & Burlington & 397631 & 719 & 1.558863e-03 \\ 
NJ & Camden & 507735 & 1255 & 1.990512e-03 \\ 
NM & Lea & 56659 & 40 & 2.221246e-04 \\ 
NM & Santa\_Fe & 105178 & 235 & 4.123374e-04 \\ 
NY & Kings & 2286167 & 4861 & 8.962636e-03 \\ 
NY & Kings & 2286167 & 4861 & 8.962636e-03 \\ 
NY & Kings & 2286167 & 4861 & 8.962636e-03 \\ 
NY & Monroe & 724418 & 2438 & 2.839991e-03 \\ 
NY & Nassau & 1302067 & 6147 & 5.104593e-03 \\ 
NY & New\_York & 1489066 & 14052 & 5.837699e-03 \\ 
NY & New\_York & 1489066 & 14052 & 5.837699e-03 \\ 
NY & Richmond & 391085 & 1242 & 1.533200e-03 \\ 
NY & Richmond & 391085 & 1242 & 1.533200e-03 \\ 
NY & Saratoga & 188387 & 215 & 7.385480e-04 \\ 
NY & Suffolk & 1338204 & 3021 & 5.246264e-03 \\ 
NY & Suffolk & 1338204 & 3021 & 5.246264e-03 \\ 
NC & Rutherford & 57951 & 46 & 2.271897e-04 \\ 
OH & Allen & 26405 & 168 & 1.035175e-04 \\ 
OH & Cuyahoga & 1411209 & 5620 & 5.532471e-03 \\ 
OH & Franklin & 992095 & 2675 & 3.889386e-03 \\ 
OH & Hamilton & 872026 & 3164 & 3.418670e-03 \\ 
OH & Harrison & 15946 & 12 & 6.251433e-05 \\ 
OH & Lucas & 461508 & 1391 & 1.809285e-03 \\ 
OH & Stark & 372125 & 595 & 1.458870e-03 \\ 
OK & Tulsa & 519847 & 1158 & 2.037996e-03 \\ 
OR & Multnomah & 600811 & 2571 & 2.355405e-03 \\ 
OR & Wasco & 22219 & 44 & 8.710685e-05 \\ 
PA & Allegheny & 1334396 & 5281 & 5.231335e-03 \\ 
PA & Allegheny & 1334396 & 5281 & 5.231335e-03 \\ 
PA & Bucks & 556279 & 879 & 2.180823e-03 \\ 
PA & Delaware & 549506 & 1374 & 2.154270e-03 \\ 
PA & Luzerne & 328927 & 594 & 1.289518e-03 \\ 
PA & Northampton & 252393 & 459 & 9.894756e-04 \\ 
SC & Greenville & 324558 & 650 & 1.272390e-03 \\ 
SD & Charles\_Mix & 9265 & 6 & 3.632229e-05 \\ 
TN & Coffee & 41641 & 50 & 1.632484e-04 \\ 
TN & Dickson & 36509 & 27 & 1.431290e-04 \\ 
TN & Hamblen & 51657 & 54 & 2.025149e-04 \\ 
TN & Shelby & 844847 & 2489 & 3.312118e-03 \\ 
TN & Sullivan & 146676 & 377 & 5.750252e-04 \\ 
TX & Cooke & 30948 & 21 & 1.213278e-04 \\ 
TX & Ector & 122338 & 153 & 4.796110e-04 \\ 
TX & Fort\_Bend & 255788 & 231 & 1.002785e-03 \\ 
TX & Wichita & 120386 & 243 & 4.719585e-04 \\ 
UT & Davis & 199953 & 166 & 7.838911e-04 \\ 
VA & Chesterfield & 225225 & 181 & 8.829668e-04 \\ 
WA & King & 1557537 & 5280 & 6.106132e-03 \\ 
WI & Lincoln & 27822 & 28 & 1.090727e-04 \\ 
WI & Waukesha & 320306 & 687 & 1.255720e-03 \\ 
\bottomrule
\end{tabular*}
\end{table}

:::

::: {.cell-output-display}
![](ups_files/figure-pdf/statepop-load-1.pdf)
:::
:::


**Estimate the Total and Per Capita Physicians:**


::: {.cell}
::: {.cell-output-display}
\begin{table}
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{0.7\linewidth}{@{\extracolsep{\fill}}rrrr}
\toprule
\(\hat{t}\) & \(\mathrm{SE}(\hat{t})\) & 95\% CI Lower & 95\% CI Upper \\ 
\midrule\addlinespace[2.5pt]
570,304.2951 & 41,401.2299 & 488,155.2729 & 652,453.3173 \\ 
\bottomrule
\end{tabular*}
\end{table}

:::

::: {.cell-output-display}
\begin{table}
\caption*{
{\fontsize{20}{25}\selectfont  Ratio Estimation: Working Table\fontsize{12}{15}\selectfont } \\ 
{\fontsize{14}{17}\selectfont  \(n = 100\)\fontsize{12}{15}\selectfont }
} 
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}>{\raggedright\arraybackslash}p{\dimexpr 30.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}}
\toprule
\(i\) & \(t_i/\psi_i\) & \(M_i/\psi_i\) & \(\hat B\,M_i/\psi_i\)\textsuperscript{\textit{1}} & \(e_i\) & \(e_i^2\) \\ 
\midrule\addlinespace[2.5pt]
1 & 7.463 $\times$ 10\textsuperscript{4} & 2.551 $\times$ 10\textsuperscript{8} & 5.703 $\times$ 10\textsuperscript{5} & -4.957 $\times$ 10\textsuperscript{5} & 2.457 $\times$ 10\textsuperscript{11} \\ 
2 & 4.987 $\times$ 10\textsuperscript{5} & 2.551 $\times$ 10\textsuperscript{8} & 5.703 $\times$ 10\textsuperscript{5} & -7.159 $\times$ 10\textsuperscript{4} & 5.126 $\times$ 10\textsuperscript{9} \\ 
3 & 4.987 $\times$ 10\textsuperscript{5} & 2.551 $\times$ 10\textsuperscript{8} & 5.703 $\times$ 10\textsuperscript{5} & -7.159 $\times$ 10\textsuperscript{4} & 5.126 $\times$ 10\textsuperscript{9} \\ 
4 & 1.288 $\times$ 10\textsuperscript{5} & 2.551 $\times$ 10\textsuperscript{8} & 5.703 $\times$ 10\textsuperscript{5} & -4.415 $\times$ 10\textsuperscript{5} & 1.949 $\times$ 10\textsuperscript{11} \\ 
5 & 4.391 $\times$ 10\textsuperscript{5} & 2.551 $\times$ 10\textsuperscript{8} & 5.703 $\times$ 10\textsuperscript{5} & -1.312 $\times$ 10\textsuperscript{5} & 1.722 $\times$ 10\textsuperscript{10} \\ 
6 & 2.224 $\times$ 10\textsuperscript{5} & 2.551 $\times$ 10\textsuperscript{8} & 5.703 $\times$ 10\textsuperscript{5} & -3.479 $\times$ 10\textsuperscript{5} & 1.211 $\times$ 10\textsuperscript{11} \\ 
... & ... & ... & ... & ... & ... \\ 
{\bfseries Sum} & {\bfseries 5.703 $\times$ 10\textsuperscript{7}} & {\bfseries 2.551 $\times$ 10\textsuperscript{10}} & {\bfseries 5.703 $\times$ 10\textsuperscript{7}} & {\bfseries -2.095 $\times$ 10\textsuperscript{-9}} & {\bfseries 1.697 $\times$ 10\textsuperscript{13}} \\ 
\bottomrule
\end{tabular*}
\begin{minipage}{\linewidth}
\vspace{.05em}
\textsuperscript{\textit{1}}\(\hat B = \sum_i t_i/\psi_i / \sum_i M_i/\psi_i = 0.0022\), \(\hat B\,M_i/\psi_i = \hat B\, M_i/\psi_i\), \(s_e^2 = 171406183906.6202\).\\
\textbf{Point estimate:} \(\hat B = \dfrac{\sum_i t_i/\psi_i}{\sum_i M_i/\psi_i} = \dfrac{5.703 \times 10^{7}}{2.551 \times 10^{10}} = 0.002\)\\
\textbf{Standard error:} \(\mathrm{SE}(\hat B) = \dfrac{1}{\bar x}\sqrt{\dfrac{s_e^2}{n}} = \dfrac{1}{2.551 \times 10^{8}}\sqrt{\dfrac{1.714 \times 10^{11}}{100}} = 1.623 \times 10^{-4}\)\\
\end{minipage}
\end{table}

:::

::: {.cell-output-display}
\begin{table}
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{0.7\linewidth}{@{\extracolsep{\fill}}rrrr}
\toprule
\(\bar{y}_{r}\) & \(\mathrm{SE}(\bar{y}_{r})\) & 95\% CI Lower & 95\% CI Upper \\ 
\midrule\addlinespace[2.5pt]
0.0022 & 0.0002 & 0.0019 & 0.0026 \\ 
\bottomrule
\end{tabular*}
\end{table}

:::
:::

**Output Note:** The total estimator heavily leverages the selection probabilities to estimate that there are roughly \ensuremath{5.70304\times 10^{5}} physicians in the US.

### Compare with an estimate based on an SRS sample

Let's compare our UPSWR results against a naive Simple Random Sample of 100 counties.


::: {.cell}
::: {.cell-output-display}
\begin{table}
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{0.7\linewidth}{@{\extracolsep{\fill}}rrrr}
\toprule
\(\hat{t}\) & \(\mathrm{SE}(\hat{t})\) & 95\% CI Lower & 95\% CI Upper \\ 
\midrule\addlinespace[2.5pt]
933,410.9700 & 491,982.7871 & -42,789.6161 & 1,909,611.5561 \\ 
\bottomrule
\end{tabular*}
\end{table}

:::

::: {.cell-output-display}
\begin{table}
\caption*{
{\fontsize{20}{25}\selectfont  Ratio Estimation: Working Table\fontsize{12}{15}\selectfont } \\ 
{\fontsize{14}{17}\selectfont  \(n = 100\)\fontsize{12}{15}\selectfont }
} 
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}>{\raggedright\arraybackslash}p{\dimexpr 30.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}}
\toprule
\(i\) & \(y_i\) & \(x_i\) & \(\hat y_i\)\textsuperscript{\textit{1}} & \(e_i\) & \(e_i^2\) \\ 
\midrule\addlinespace[2.5pt]
1 & 24.000 & 3.602 $\times$ 10\textsuperscript{4} & 90.313 & -6.631 $\times$ 10\textsuperscript{1} & 4.397 $\times$ 10\textsuperscript{3} \\ 
2 & 44.000 & 7.352 $\times$ 10\textsuperscript{4} & 184.332 & -1.403 $\times$ 10\textsuperscript{2} & 1.969 $\times$ 10\textsuperscript{4} \\ 
3 & 7.000 & 6.408 $\times$ 10\textsuperscript{3} & 16.066 & -9.066 & 8.218 $\times$ 10\textsuperscript{1} \\ 
4 & 11.000 & 1.926 $\times$ 10\textsuperscript{4} & 48.289 & -3.729 $\times$ 10\textsuperscript{1} & 1.390 $\times$ 10\textsuperscript{3} \\ 
5 & 3.000 & 7.649 $\times$ 10\textsuperscript{3} & 19.177 & -1.618 $\times$ 10\textsuperscript{1} & 2.617 $\times$ 10\textsuperscript{2} \\ 
6 & 327.000 & 1.884 $\times$ 10\textsuperscript{5} & 472.281 & -1.453 $\times$ 10\textsuperscript{2} & 2.111 $\times$ 10\textsuperscript{4} \\ 
... & ... & ... & ... & ... & ... \\ 
{\bfseries Sum} & {\bfseries 29,717.000} & {\bfseries 1.185 $\times$ 10\textsuperscript{7}} & {\bfseries 29,717.000} & {\bfseries -2.675 $\times$ 10\textsuperscript{-12}} & {\bfseries 1.705 $\times$ 10\textsuperscript{7}} \\ 
\bottomrule
\end{tabular*}
\begin{minipage}{\linewidth}
\vspace{.05em}
\textsuperscript{\textit{1}}\(\hat B = \sum_i y_i / \sum_i x_i = 0.0025\), \(\hat y_i = \hat B\, x_i\), \(s_e^2 = 172267.8740\).\\
\textbf{Point estimate:} \(\hat t = \hat B\,\bar x_U\,N = 0.003 \times 81,209.021 \times 3,141 = 6.395 \times 10^{5}\)\\
\textbf{Standard error:} \(\mathrm{SE}(\hat t) = \mathrm{SE}(\hat B)\,\bar x_U\,N = 3.445 \times 10^{-4} \times 81,209.021 \times 3,141 = 87,885.265\)\\
\end{minipage}
\end{table}

:::

::: {.cell-output-display}
\begin{table}
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{0.7\linewidth}{@{\extracolsep{\fill}}rrrr}
\toprule
\(\hat{t}\) & \(\mathrm{SE}(\hat{t})\) & 95\% CI Lower & 95\% CI Upper \\ 
\midrule\addlinespace[2.5pt]
639,506.0284 & 87,885.2651 & 465,122.5955 & 813,889.4613 \\ 
\bottomrule
\end{tabular*}
\end{table}

:::

::: {.cell-output-display}
\begin{table}
\caption*{
{\fontsize{20}{25}\selectfont  Regression Estimation: Working Table\fontsize{12}{15}\selectfont } \\ 
{\fontsize{14}{17}\selectfont  \(n = 100\)\fontsize{12}{15}\selectfont }
} 
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}>{\raggedright\arraybackslash}p{\dimexpr 30.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}}
\toprule
\(i\) & \(y_i\) & \(x_i\) & \(\hat y_i\)\textsuperscript{\textit{1}} & \(e_i\) & \(e_i^2\) \\ 
\midrule\addlinespace[2.5pt]
1 & 24.000 & 3.602 $\times$ 10\textsuperscript{4} & 52.564 & -2.856 $\times$ 10\textsuperscript{1} & 8.159 $\times$ 10\textsuperscript{2} \\ 
2 & 44.000 & 7.352 $\times$ 10\textsuperscript{4} & 163.740 & -1.197 $\times$ 10\textsuperscript{2} & 1.434 $\times$ 10\textsuperscript{4} \\ 
3 & 7.000 & 6.408 $\times$ 10\textsuperscript{3} & -35.234 & 4.223 $\times$ 10\textsuperscript{1} & 1.784 $\times$ 10\textsuperscript{3} \\ 
4 & 11.000 & 1.926 $\times$ 10\textsuperscript{4} & 2.871 & 8.129 & 6.609 $\times$ 10\textsuperscript{1} \\ 
5 & 3.000 & 7.649 $\times$ 10\textsuperscript{3} & -31.555 & 3.455 $\times$ 10\textsuperscript{1} & 1.194 $\times$ 10\textsuperscript{3} \\ 
6 & 327.000 & 1.884 $\times$ 10\textsuperscript{5} & 504.237 & -1.772 $\times$ 10\textsuperscript{2} & 3.141 $\times$ 10\textsuperscript{4} \\ 
... & ... & ... & ... & ... & ... \\ 
{\bfseries Sum} & {\bfseries 29,717.000} & {\bfseries 1.185 $\times$ 10\textsuperscript{7}} & {\bfseries 29,717.000} & {\bfseries 6.963 $\times$ 10\textsuperscript{-13}} & {\bfseries 1.135 $\times$ 10\textsuperscript{7}} \\ 
\bottomrule
\end{tabular*}
\begin{minipage}{\linewidth}
\vspace{.05em}
\textsuperscript{\textit{1}}\(\hat B_0 = -54.2313\), \(\hat B_1 = 0.0030\), \(\hat y_i = \hat B_0 + \hat B_1 x_i\), \(s_e^2 = 115813.8922\).\\
\textbf{Point estimate:} \(\hat t_{reg} = N\bar y_{reg} = 3,141 \times 186.524 = 5.859 \times 10^{5}\)\\
\textbf{Standard error:} \(\mathrm{SE}(\hat t_{reg}) = N\times\mathrm{SE}(\bar y_{reg}) = 3,141 \times 33.485 = 1.052 \times 10^{5}\)\\
\end{minipage}
\end{table}

:::

::: {.cell-output-display}
\begin{table}
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{0.7\linewidth}{@{\extracolsep{\fill}}rrrr}
\toprule
\(\hat{t}_{reg}\) & \(\mathrm{SE}(\hat{t}_{reg})\) & 95\% CI Lower & 95\% CI Upper \\ 
\midrule\addlinespace[2.5pt]
585,870.5992 & 105,177.4184 & 377,149.4353 & 794,591.7630 \\ 
\bottomrule
\end{tabular*}
\end{table}

:::
:::

**Output Note:** The true value is **532,638**. Notice how poorly the naive SRS estimate performs (predicting over 800,000) because it ignores county sizes. Both the SRS Ratio and SRS Regression methods perform much better by incorporating auxiliary population data.

## Analysis of `studenthrs.csv` collected with two-stage UPSWR

Here, we sample classes (Primary Sampling Units) with replacement based on class size, and then subsample students within them to ask about their study hours.


::: {.cell}
::: {.cell-output-display}
\begin{table}
\caption*{
{\fontsize{20}{25}\selectfont  Ratio Estimation: Working Table\fontsize{12}{15}\selectfont } \\ 
{\fontsize{14}{17}\selectfont  \(n = 5\)\fontsize{12}{15}\selectfont }
} 
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}>{\raggedright\arraybackslash}p{\dimexpr 30.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}}
\toprule
\(i\) & \(t_i/\psi_i\) & \(M_i/\psi_i\) & \(\hat B\,M_i/\psi_i\)\textsuperscript{\textit{1}} & \(e_i\) & \(e_i^2\) \\ 
\midrule\addlinespace[2.5pt]
1 & 2,393.900 & 647.000 & 1,617.500 & 776.400 & 6.028 $\times$ 10\textsuperscript{5} \\ 
2 & 1,811.600 & 647.000 & 1,617.500 & 194.100 & 3.767 $\times$ 10\textsuperscript{4} \\ 
3 & 1,552.800 & 647.000 & 1,617.500 & -64.700 & 4.186 $\times$ 10\textsuperscript{3} \\ 
4 & 1,035.200 & 647.000 & 1,617.500 & -582.300 & 3.391 $\times$ 10\textsuperscript{5} \\ 
5 & 1,294.000 & 647.000 & 1,617.500 & -323.500 & 1.047 $\times$ 10\textsuperscript{5} \\ 
{\bfseries Sum} & {\bfseries 8,087.500} & {\bfseries 3,235.000} & {\bfseries 8,087.500} & {\bfseries 0.000} & {\bfseries 1.088 $\times$ 10\textsuperscript{6}} \\ 
\bottomrule
\end{tabular*}
\begin{minipage}{\linewidth}
\vspace{.05em}
\textsuperscript{\textit{1}}\(\hat B = \sum_i t_i/\psi_i / \sum_i M_i/\psi_i = 2.5000\), \(\hat B\,M_i/\psi_i = \hat B\, M_i/\psi_i\), \(s_e^2 = 272095.8500\).\\
\textbf{Point estimate:} \(\hat B = \dfrac{\sum_i t_i/\psi_i}{\sum_i M_i/\psi_i} = \dfrac{8,087.500}{3,235.000} = 2.500\)\\
\textbf{Standard error:} \(\mathrm{SE}(\hat B) = \dfrac{1}{\bar x}\sqrt{\dfrac{s_e^2}{n}} = \dfrac{1}{647.000}\sqrt{\dfrac{2.721 \times 10^{5}}{5}} = 0.361\)\\
\end{minipage}
\end{table}

:::
:::


**Estimate the overall mean study hours:**


::: {.cell}
::: {.cell-output-display}
\begin{table}
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{0.7\linewidth}{@{\extracolsep{\fill}}rrrr}
\toprule
\(\bar{y}_{r}\) & \(\mathrm{SE}(\bar{y}_{r})\) & 95\% CI Lower & 95\% CI Upper \\ 
\midrule\addlinespace[2.5pt]
2.5000 & 0.3606 & 1.4989 & 3.5011 \\ 
\bottomrule
\end{tabular*}
\end{table}

:::
:::


## Simulation Studies for UPSWR using `agpop` data

Let's simulate 2000 random samples to see how UPSWR empirically compares to Simple Random Sampling. 


::: {.cell}

:::


### Visualizing Estimator Variance


::: {.cell}
::: {.cell-output-display}
![](ups_files/figure-pdf/upswr-sim-plot-1.pdf)
:::
:::

**Output Note:** The red line represents the true total acreage. The naive SRS estimator has massive variance. Both SRS Regression and UPSWR drastically reduce variance, centering tightly around the true mean.

## Analysis of a dataset on student study hours collected by UPSWO 

When sampling **Without Replacement** with unequal probabilities, we rely on the inclusion probabilities ($\pi_k$) rather than single-draw probabilities ($\psi_i$). 


::: {.cell}
::: {.cell-output-display}
\begin{table}
\caption*{
{\fontsize{20}{25}\selectfont  Ratio Estimation: Working Table\fontsize{12}{15}\selectfont } \\ 
{\fontsize{14}{17}\selectfont  \(n = 5\)\fontsize{12}{15}\selectfont }
} 
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}>{\raggedright\arraybackslash}p{\dimexpr 30.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}}
\toprule
\(i\) & \(\hat t_i/\pi_i\) & \(M_i/\pi_i\) & \(\hat B\,M_i/\pi_i\)\textsuperscript{\textit{1}} & \(e_i\) & \(e_i^2\) \\ 
\midrule\addlinespace[2.5pt]
1 & 3.500 & 1.000 & 3.450 & 5.000 $\times$ 10\textsuperscript{-2} & 0.002 \\ 
2 & 5.000 & 1.000 & 3.450 & 1.550 & 2.402 \\ 
3 & 3.625 & 1.000 & 3.450 & 1.750 $\times$ 10\textsuperscript{-1} & 0.031 \\ 
4 & 3.125 & 1.000 & 3.450 & -3.250 $\times$ 10\textsuperscript{-1} & 0.106 \\ 
5 & 2.000 & 1.000 & 3.450 & -1.450 & 2.103 \\ 
{\bfseries Sum} & {\bfseries 17.250} & {\bfseries 5.000} & {\bfseries 17.250} & {\bfseries -8.882 $\times$ 10\textsuperscript{-16}} & {\bfseries 4.644} \\ 
\bottomrule
\end{tabular*}
\begin{minipage}{\linewidth}
\vspace{.05em}
\textsuperscript{\textit{1}}\(\hat B = \sum_i \hat t_i/\pi_i / \sum_i M_i/\pi_i = 3.4500\), \(\hat B\,M_i/\pi_i = \hat B\, M_i/\pi_i\), \(s_e^2 = 1.1609\).\\
\textbf{Point estimate:} \(\hat B = \dfrac{\sum_i \hat t_i/\pi_i}{\sum_i M_i/\pi_i} = \dfrac{17.250}{5.000} = 3.450\)\\
\textbf{Standard error:} \(\mathrm{SE}(\hat B) = \dfrac{1}{\bar x}\sqrt{\dfrac{s_e^2}{n}} = \dfrac{1}{1.000}\sqrt{\dfrac{1.161}{5}} = 0.482\)\\
\end{minipage}
\end{table}

:::
:::


**Estimate the overall mean study hours (Horvitz-Thompson ratio style):**


::: {.cell}
::: {.cell-output-display}
\begin{table}
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{0.7\linewidth}{@{\extracolsep{\fill}}rrrr}
\toprule
\(\bar{y}_{r}\) & \(\mathrm{SE}(\bar{y}_{r})\) & 95\% CI Lower & 95\% CI Upper \\ 
\midrule\addlinespace[2.5pt]
3.4500 & 0.4819 & 2.1121 & 4.7879 \\ 
\bottomrule
\end{tabular*}
\end{table}

:::
:::
