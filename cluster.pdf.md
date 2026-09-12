# Cluster Sampling

## Key Formulas for Cluster Sampling

Cluster sampling divides the population into $N$ clusters, out of which a simple random sample of $n$ clusters is selected. This design is highly practical when a sampling frame of individual elements is unavailable, or when individuals are geographically grouped.

#### 1. One-Stage Cluster Sampling {-}
In a one-stage design, every element $M_i$ within the sampled cluster $i$ is surveyed. 

*   **Total in Cluster $i$:** $t_i = \sum_{j=1}^{M_i} y_{ij}$
*   **Ratio Estimator of the Population Mean:** The ratio estimator leverages cluster sizes to estimate the overall element mean:
    $$
\overline{y}_r = \frac{\sum_{i=1}^n t_i}{\sum_{i=1}^n M_i}.
$$

*   **Estimated Variance of the Mean:**
    $$
\hat{V}(\overline{y}_r) = \left(1 - \frac{n}{N}\right) \frac{1}{n\overline{M}^2} \frac{\sum_{i=1}^n (t_i - \overline{y}_r M_i)^2}{n-1}.
$$

#### 2. Two-Stage Cluster Sampling {-}
In a two-stage design, we first sample $n$ clusters, then take a subsample of $m_i$ elements from the $M_i$ elements within each selected cluster. Because we don't observe the true cluster total $t_i$, we estimate it:

*   **Estimated Total in Cluster $i$:** $\hat{t}_i = M_i \overline{y}_i$, where $\overline{y}_i = \frac{1}{m_i} \sum_{j=1}^{m_i} y_{ij}$
*   **Ratio Estimator of the Population Mean:**
    $$
\overline{y}_r = \frac{\sum_{i=1}^n \hat{t}_i}{\sum_{i=1}^n M_i}
$$

---

## Functions and packages for Analyzing Data

The estimating functions used throughout this book (including `cluster_ratio`, used below for one-stage and two-stage cluster samples) live in a single shared file, `samplingestimate.r`, which we source below.


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

## Estimating for One-stage cluster sampling: analysis of algebra.csv

Consider a population of 187 high school algebra classes in a city. An investigator takes an SRS of 12 of those classes (the clusters). Because this is a one-stage sample, every student in the 12 sampled classes takes the test to evaluate their function knowledge.

We read the raw data and display it using `gt`'s interactive pagination feature to keep the document clean.


::: {.cell}
::: {.cell-output-display}
\begin{table}
\caption*{
{\fontsize{20}{25}\selectfont  Raw Data: Algebra Classes\fontsize{12}{15}\selectfont }
} 
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}rrr}
\toprule
class & M<sub>i</sub> & score \\ 
\midrule\addlinespace[2.5pt]
23 & 20 & 57 \\ 
23 & 20 & 90 \\ 
23 & 20 & 56 \\ 
23 & 20 & 57 \\ 
23 & 20 & 46 \\ 
23 & 20 & 55 \\ 
23 & 20 & 62 \\ 
23 & 20 & 66 \\ 
23 & 20 & 78 \\ 
23 & 20 & 76 \\ 
23 & 20 & 57 \\ 
23 & 20 & 84 \\ 
23 & 20 & 27 \\ 
23 & 20 & 70 \\ 
23 & 20 & 49 \\ 
23 & 20 & 82 \\ 
23 & 20 & 52 \\ 
23 & 20 & 82 \\ 
23 & 20 & 59 \\ 
23 & 20 & 25 \\ 
37 & 26 & 57 \\ 
37 & 26 & 34 \\ 
37 & 26 & 68 \\ 
37 & 26 & 99 \\ 
37 & 26 & 24 \\ 
37 & 26 & 83 \\ 
37 & 26 & 40 \\ 
37 & 26 & 98 \\ 
37 & 26 & 90 \\ 
37 & 26 & 50 \\ 
37 & 26 & 60 \\ 
37 & 26 & 64 \\ 
37 & 26 & 58 \\ 
37 & 26 & 100 \\ 
37 & 26 & 68 \\ 
37 & 26 & 52 \\ 
37 & 26 & 77 \\ 
37 & 26 & 69 \\ 
37 & 26 & 64 \\ 
37 & 26 & 64 \\ 
37 & 26 & 73 \\ 
37 & 26 & 46 \\ 
37 & 26 & 44 \\ 
37 & 26 & 95 \\ 
37 & 26 & 28 \\ 
37 & 26 & 65 \\ 
38 & 24 & 41 \\ 
38 & 24 & 66 \\ 
38 & 24 & 79 \\ 
38 & 24 & 24 \\ 
38 & 24 & 77 \\ 
38 & 24 & 60 \\ 
38 & 24 & 37 \\ 
38 & 24 & 26 \\ 
38 & 24 & 74 \\ 
38 & 24 & 58 \\ 
38 & 24 & 65 \\ 
38 & 24 & 76 \\ 
38 & 24 & 92 \\ 
38 & 24 & 59 \\ 
38 & 24 & 51 \\ 
38 & 24 & 87 \\ 
38 & 24 & 45 \\ 
38 & 24 & 40 \\ 
38 & 24 & 51 \\ 
38 & 24 & 100 \\ 
38 & 24 & 46 \\ 
38 & 24 & 68 \\ 
38 & 24 & 40 \\ 
38 & 24 & 40 \\ 
39 & 34 & 54 \\ 
39 & 34 & 84 \\ 
39 & 34 & 51 \\ 
39 & 34 & 75 \\ 
39 & 34 & 34 \\ 
39 & 34 & 53 \\ 
39 & 34 & 46 \\ 
39 & 34 & 53 \\ 
39 & 34 & 82 \\ 
39 & 34 & 58 \\ 
39 & 34 & 63 \\ 
39 & 34 & 64 \\ 
39 & 34 & 49 \\ 
39 & 34 & 34 \\ 
39 & 34 & 54 \\ 
39 & 34 & 62 \\ 
39 & 34 & 78 \\ 
39 & 34 & 55 \\ 
39 & 34 & 54 \\ 
39 & 34 & 89 \\ 
39 & 34 & 46 \\ 
39 & 34 & 67 \\ 
39 & 34 & 49 \\ 
39 & 34 & 67 \\ 
39 & 34 & 32 \\ 
39 & 34 & 50 \\ 
39 & 34 & 55 \\ 
39 & 34 & 78 \\ 
39 & 34 & 51 \\ 
39 & 34 & 65 \\ 
39 & 34 & 52 \\ 
39 & 34 & 44 \\ 
39 & 34 & 72 \\ 
39 & 34 & 52 \\ 
41 & 26 & 82 \\ 
41 & 26 & 65 \\ 
41 & 26 & 55 \\ 
41 & 26 & 65 \\ 
41 & 26 & 57 \\ 
41 & 26 & 65 \\ 
41 & 26 & 50 \\ 
41 & 26 & 19 \\ 
41 & 26 & 60 \\ 
41 & 26 & 69 \\ 
41 & 26 & 40 \\ 
41 & 26 & 37 \\ 
41 & 26 & 49 \\ 
41 & 26 & 48 \\ 
41 & 26 & 65 \\ 
41 & 26 & 60 \\ 
41 & 26 & 76 \\ 
41 & 26 & 87 \\ 
41 & 26 & 57 \\ 
41 & 26 & 57 \\ 
41 & 26 & 62 \\ 
41 & 26 & 51 \\ 
41 & 26 & 55 \\ 
41 & 26 & 55 \\ 
41 & 26 & 62 \\ 
41 & 26 & 60 \\ 
44 & 28 & 59 \\ 
44 & 28 & 65 \\ 
44 & 28 & 60 \\ 
44 & 28 & 78 \\ 
44 & 28 & 64 \\ 
44 & 28 & 52 \\ 
44 & 28 & 59 \\ 
44 & 28 & 51 \\ 
44 & 28 & 67 \\ 
44 & 28 & 79 \\ 
44 & 28 & 75 \\ 
44 & 28 & 67 \\ 
44 & 28 & 62 \\ 
44 & 28 & 69 \\ 
44 & 28 & 63 \\ 
44 & 28 & 27 \\ 
44 & 28 & 70 \\ 
44 & 28 & 57 \\ 
44 & 28 & 86 \\ 
44 & 28 & 100 \\ 
44 & 28 & 41 \\ 
44 & 28 & 74 \\ 
44 & 28 & 66 \\ 
44 & 28 & 54 \\ 
44 & 28 & 75 \\ 
44 & 28 & 76 \\ 
44 & 28 & 62 \\ 
44 & 28 & 58 \\ 
46 & 19 & 34 \\ 
46 & 19 & 42 \\ 
46 & 19 & 49 \\ 
46 & 19 & 49 \\ 
46 & 19 & 71 \\ 
46 & 19 & 49 \\ 
46 & 19 & 88 \\ 
46 & 19 & 46 \\ 
46 & 19 & 82 \\ 
46 & 19 & 57 \\ 
46 & 19 & 34 \\ 
46 & 19 & 64 \\ 
46 & 19 & 54 \\ 
46 & 19 & 47 \\ 
46 & 19 & 69 \\ 
46 & 19 & 63 \\ 
46 & 19 & 72 \\ 
46 & 19 & 35 \\ 
46 & 19 & 43 \\ 
51 & 32 & 58 \\ 
51 & 32 & 22 \\ 
51 & 32 & 82 \\ 
51 & 32 & 99 \\ 
51 & 32 & 88 \\ 
51 & 32 & 52 \\ 
51 & 32 & 61 \\ 
51 & 32 & 44 \\ 
51 & 32 & 98 \\ 
51 & 32 & 70 \\ 
51 & 32 & 59 \\ 
51 & 32 & 92 \\ 
51 & 32 & 72 \\ 
51 & 32 & 49 \\ 
51 & 32 & 76 \\ 
51 & 32 & 77 \\ 
51 & 32 & 79 \\ 
51 & 32 & 84 \\ 
51 & 32 & 65 \\ 
51 & 32 & 88 \\ 
51 & 32 & 55 \\ 
51 & 32 & 100 \\ 
51 & 32 & 65 \\ 
51 & 32 & 84 \\ 
51 & 32 & 87 \\ 
51 & 32 & 71 \\ 
51 & 32 & 47 \\ 
51 & 32 & 70 \\ 
51 & 32 & 66 \\ 
51 & 32 & 69 \\ 
51 & 32 & 91 \\ 
51 & 32 & 88 \\ 
58 & 17 & 61 \\ 
58 & 17 & 44 \\ 
58 & 17 & 100 \\ 
58 & 17 & 28 \\ 
58 & 17 & 28 \\ 
58 & 17 & 63 \\ 
58 & 17 & 100 \\ 
58 & 17 & 51 \\ 
58 & 17 & 34 \\ 
58 & 17 & 49 \\ 
58 & 17 & 100 \\ 
58 & 17 & 43 \\ 
58 & 17 & 48 \\ 
58 & 17 & 52 \\ 
58 & 17 & 70 \\ 
58 & 17 & 69 \\ 
58 & 17 & 49 \\ 
62 & 21 & 31 \\ 
62 & 21 & 99 \\ 
62 & 21 & 43 \\ 
62 & 21 & 57 \\ 
62 & 21 & 69 \\ 
62 & 21 & 73 \\ 
62 & 21 & 43 \\ 
62 & 21 & 71 \\ 
62 & 21 & 81 \\ 
62 & 21 & 74 \\ 
62 & 21 & 81 \\ 
62 & 21 & 49 \\ 
62 & 21 & 63 \\ 
62 & 21 & 85 \\ 
62 & 21 & 75 \\ 
62 & 21 & 73 \\ 
62 & 21 & 85 \\ 
62 & 21 & 57 \\ 
62 & 21 & 77 \\ 
62 & 21 & 51 \\ 
62 & 21 & 61 \\ 
106 & 26 & 65 \\ 
106 & 26 & 72 \\ 
106 & 26 & 14 \\ 
106 & 26 & 82 \\ 
106 & 26 & 69 \\ 
106 & 26 & 86 \\ 
106 & 26 & 82 \\ 
106 & 26 & 73 \\ 
106 & 26 & 65 \\ 
106 & 26 & 74 \\ 
106 & 26 & 61 \\ 
106 & 26 & 73 \\ 
106 & 26 & 68 \\ 
106 & 26 & 77 \\ 
106 & 26 & 50 \\ 
106 & 26 & 55 \\ 
106 & 26 & 69 \\ 
106 & 26 & 42 \\ 
106 & 26 & 55 \\ 
106 & 26 & 67 \\ 
106 & 26 & 94 \\ 
106 & 26 & 76 \\ 
106 & 26 & 38 \\ 
106 & 26 & 28 \\ 
106 & 26 & 43 \\ 
106 & 26 & 43 \\ 
108 & 26 & 54 \\ 
108 & 26 & 72 \\ 
108 & 26 & 84 \\ 
108 & 26 & 68 \\ 
108 & 26 & 72 \\ 
108 & 26 & 31 \\ 
108 & 26 & 75 \\ 
108 & 26 & 81 \\ 
108 & 26 & 51 \\ 
108 & 26 & 57 \\ 
108 & 26 & 62 \\ 
108 & 26 & 57 \\ 
108 & 26 & 71 \\ 
108 & 26 & 55 \\ 
108 & 26 & 58 \\ 
108 & 26 & 52 \\ 
108 & 26 & 79 \\ 
108 & 26 & 43 \\ 
108 & 26 & 81 \\ 
108 & 26 & 70 \\ 
108 & 26 & 64 \\ 
108 & 26 & 69 \\ 
108 & 26 & 100 \\ 
108 & 26 & 74 \\ 
108 & 26 & 89 \\ 
108 & 26 & 77 \\ 
\bottomrule
\end{tabular*}
\end{table}

:::
:::

**Output Note:** The table above displays the raw, unaggregated survey data. The interactive pagination allows you to flip through the dataset compactly.

Next, we run our updated function to compute both the cluster-level working table and the overall ratio estimate for the average student score.


::: {.cell}
::: {.cell-output-display}
\begin{table}
\caption*{
{\fontsize{20}{25}\selectfont  Ratio Estimation: Working Table\fontsize{12}{15}\selectfont } \\ 
{\fontsize{14}{17}\selectfont  \(n = 12\)\fontsize{12}{15}\selectfont }
} 
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}>{\raggedright\arraybackslash}p{\dimexpr 30.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}}
\toprule
\(i\) & \(\bar y_i\) & \(\hat t_i\)\textsuperscript{\textit{1}} & \(M_i\) & \(\hat B M_i\)\textsuperscript{\textit{2}} & \(e_i\) & \(e_i^2\) \\ 
\midrule\addlinespace[2.5pt]
1 & 61.500 & 1,230.000 & 20.000 & 1,251.371 & -2.137 $\times$ 10\textsuperscript{1} & 4.567 $\times$ 10\textsuperscript{2} \\ 
2 & 64.231 & 1,670.000 & 26.000 & 1,626.783 & 4.322 $\times$ 10\textsuperscript{1} & 1.868 $\times$ 10\textsuperscript{3} \\ 
3 & 58.417 & 1,402.000 & 24.000 & 1,501.645 & -9.965 $\times$ 10\textsuperscript{1} & 9.929 $\times$ 10\textsuperscript{3} \\ 
4 & 58.000 & 1,972.000 & 34.000 & 2,127.331 & -1.553 $\times$ 10\textsuperscript{2} & 2.413 $\times$ 10\textsuperscript{4} \\ 
5 & 58.000 & 1,508.000 & 26.000 & 1,626.783 & -1.188 $\times$ 10\textsuperscript{2} & 1.411 $\times$ 10\textsuperscript{4} \\ 
6 & 64.857 & 1,816.000 & 28.000 & 1,751.920 & 6.408 $\times$ 10\textsuperscript{1} & 4.106 $\times$ 10\textsuperscript{3} \\ 
... & ... & ... & ... & ... & ... & ... \\ 
{\bfseries Sum} & {\bfseries ...} & {\bfseries 18,708.000} & {\bfseries 299.000} & {\bfseries 18,708.000} & {\bfseries 1.364 $\times$ 10\textsuperscript{-12}} & {\bfseries 1.948 $\times$ 10\textsuperscript{5}} \\ 
\bottomrule
\end{tabular*}
\begin{minipage}{\linewidth}
\vspace{.05em}
\textsuperscript{\textit{1}}\(\hat t_i = \bar y_i \, M_i\).\\
\textsuperscript{\textit{2}}\(\hat B = \sum_i \hat t_i / \sum_i M_i = 62.5686\), \(\hat B M_i = \hat B\, M_i\), \(s_e^2 = 17711.5489\).\\
\textbf{Point estimate:} \(\hat B = \dfrac{\sum_i \hat t_i}{\sum_i M_i} = \dfrac{18,708.000}{299.000} = 62.569\)\\
\textbf{Standard error:} \(\mathrm{SE}(\hat B) = \dfrac{1}{\bar x}\sqrt{\left(1-\dfrac{n}{N}\right)\dfrac{s_e^2}{n}} = \dfrac{1}{24.917}\sqrt{\left(1-\dfrac{12}{187}\right)\dfrac{17,711.549}{12}} = 1.492\)\\
\end{minipage}
\end{table}

:::
:::

**Output Note:** The working table above displays, for each of the 12 randomly selected classes, the class size ($M_i$), the total class score ($t_i$), the fitted value $\hat B M_i$, and the residual and squared residual used to estimate $s_e^2$.

Finally, we extract the ratio estimate from the function output list to assess the overall population mean.


::: {.cell}
::: {.cell-output-display}
\begin{table}
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{0.7\linewidth}{@{\extracolsep{\fill}}rrrr}
\toprule
\(\bar{y}_{r}\) & \(\mathrm{SE}(\bar{y}_{r})\) & 95\% CI Lower & 95\% CI Upper \\ 
\midrule\addlinespace[2.5pt]
62.5686 & 1.4916 & 59.2856 & 65.8515 \\ 
\bottomrule
\end{tabular*}
\end{table}

:::
:::

**Output Note:** Using the ratio estimator, the estimated mean score for a student across all 187 algebra classes is 62.57, with a 95% confidence interval of [59.29, 65.85]. The standard error accounts for the clustering effect, which typically increases variance compared to a pure SRS of students.

## Estimating for two-stage cluster sampling: analysis of coots.csv

In two-stage cluster sampling, we sample clusters, and then subsample elements within them. In this dataset, researcher Arnold investigated coot eggs. He first selected a random sample of clutches (the clusters), and then measured the volume of a random subsample of eggs within each selected clutch.



::: {.cell}
::: {.cell-output-display}
\begin{table}
\caption*{
{\fontsize{20}{25}\selectfont  Raw Data: Coots\fontsize{12}{15}\selectfont }
} 
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}rrrrrr}
\toprule
clutch & csize & length & breadth & volume & tmt \\ 
\midrule\addlinespace[2.5pt]
1 & 13 & 44.30 & 31.10 & 3.7957569 & 1 \\ 
1 & 13 & 45.90 & 32.70 & 3.9328497 & 1 \\ 
2 & 13 & 49.20 & 34.40 & 4.2156036 & 1 \\ 
2 & 13 & 48.70 & 32.70 & 4.1727621 & 1 \\ 
3 & 6 & 51.05 & 34.25 & 0.9317646 & 0 \\ 
3 & 6 & 49.35 & 34.40 & 0.9007362 & 0 \\ 
4 & 11 & 49.20 & 31.55 & 3.0182724 & 1 \\ 
4 & 11 & 48.55 & 33.10 & 2.9783969 & 1 \\ 
5 & 10 & 49.40 & 34.55 & 2.5045800 & 1 \\ 
5 & 10 & 49.05 & 34.95 & 2.4868350 & 1 \\ 
6 & 13 & 46.35 & 33.10 & 3.9714070 & 1 \\ 
6 & 13 & 46.65 & 33.30 & 3.9971119 & 1 \\ 
7 & 9 & 48.15 & 33.95 & 1.9773760 & 1 \\ 
7 & 9 & 45.70 & 32.65 & 1.8767619 & 1 \\ 
8 & 11 & 47.45 & 34.10 & 2.9109152 & 1 \\ 
8 & 11 & 49.10 & 33.80 & 3.0121377 & 1 \\ 
9 & 12 & 47.50 & 33.50 & 3.4678800 & 1 \\ 
9 & 12 & 47.30 & 33.30 & 3.4532784 & 1 \\ 
10 & 11 & 50.00 & 33.45 & 3.0673500 & 1 \\ 
10 & 11 & 46.55 & 31.60 & 2.8557029 & 1 \\ 
11 & 12 & 48.50 & 33.05 & 3.5408880 & 0 \\ 
11 & 12 & 47.35 & 34.00 & 3.4569288 & 0 \\ 
12 & 11 & 48.75 & 34.50 & 2.9906662 & 1 \\ 
12 & 11 & 49.05 & 34.35 & 3.0090704 & 1 \\ 
13 & 12 & 49.95 & 33.70 & 3.6467496 & 1 \\ 
13 & 12 & 47.75 & 33.05 & 3.4861320 & 1 \\ 
14 & 11 & 49.30 & 35.30 & 3.0244071 & 1 \\ 
14 & 11 & 48.05 & 35.80 & 2.9477233 & 1 \\ 
15 & 11 & 49.00 & 33.55 & 3.0060030 & 1 \\ 
15 & 11 & 48.25 & 33.30 & 2.9599928 & 1 \\ 
16 & 10 & 48.10 & 34.10 & 2.4386700 & 0 \\ 
16 & 10 & 46.85 & 33.15 & 2.3752950 & 0 \\ 
17 & 9 & 48.55 & 35.35 & 1.9938029 & 0 \\ 
17 & 9 & 49.35 & 34.85 & 2.0266564 & 0 \\ 
18 & 10 & 48.15 & 33.60 & 2.4412050 & 0 \\ 
18 & 10 & 48.00 & 32.85 & 2.4336000 & 0 \\ 
19 & 11 & 47.80 & 33.45 & 2.9323866 & 1 \\ 
19 & 11 & 47.60 & 33.50 & 2.9201172 & 1 \\ 
20 & 11 & 50.40 & 32.95 & 3.0918888 & 1 \\ 
20 & 11 & 45.70 & 33.25 & 2.8035579 & 1 \\ 
21 & 10 & 50.70 & 36.00 & 2.5704900 & 1 \\ 
21 & 10 & 49.40 & 34.70 & 2.5045800 & 1 \\ 
22 & 13 & 50.90 & 34.00 & 4.3612647 & 1 \\ 
22 & 13 & 48.75 & 34.55 & 4.1770463 & 1 \\ 
23 & 12 & 52.15 & 35.70 & 3.8073672 & 0 \\ 
23 & 12 & 51.00 & 35.75 & 3.7234080 & 0 \\ 
24 & 10 & 50.00 & 33.00 & 2.5350000 & 0 \\ 
24 & 10 & 51.15 & 33.20 & 2.5933050 & 0 \\ 
25 & 9 & 47.45 & 32.50 & 1.9486291 & 1 \\ 
25 & 9 & 48.10 & 32.90 & 1.9753227 & 1 \\ 
26 & 12 & 49.30 & 33.80 & 3.5992944 & 1 \\ 
26 & 12 & 46.10 & 31.80 & 3.3656688 & 1 \\ 
27 & 10 & 50.55 & 35.50 & 2.5628850 & 1 \\ 
27 & 10 & 51.00 & 35.70 & 2.5857000 & 1 \\ 
28 & 9 & 46.55 & 33.95 & 1.9116689 & 1 \\ 
28 & 9 & 48.70 & 33.05 & 1.9999629 & 1 \\ 
29 & 12 & 47.55 & 34.15 & 3.4715304 & 1 \\ 
29 & 12 & 45.80 & 32.70 & 3.3437664 & 1 \\ 
30 & 11 & 48.35 & 32.20 & 2.9661275 & 0 \\ 
30 & 11 & 48.00 & 32.25 & 2.9446560 & 0 \\ 
31 & 10 & 49.35 & 33.85 & 2.5020450 & 0 \\ 
31 & 10 & 49.45 & 33.00 & 2.5071150 & 0 \\ 
32 & 8 & 47.65 & 32.10 & 1.5461472 & 0 \\ 
32 & 8 & 48.50 & 33.70 & 1.5737280 & 0 \\ 
33 & 10 & 49.80 & 34.40 & 2.5248600 & 1 \\ 
33 & 10 & 48.50 & 34.00 & 2.4589500 & 1 \\ 
34 & 10 & 48.05 & 35.15 & 2.4361350 & 1 \\ 
34 & 10 & 50.05 & 35.25 & 2.5375350 & 1 \\ 
35 & 9 & 48.00 & 33.30 & 1.9712160 & 0 \\ 
35 & 9 & 48.75 & 33.30 & 2.0020163 & 0 \\ 
36 & 9 & 47.90 & 35.40 & 1.9671093 & 0 \\ 
36 & 9 & 45.85 & 33.80 & 1.8829220 & 0 \\ 
37 & 8 & 45.10 & 33.30 & 1.4634048 & 0 \\ 
37 & 8 & 47.20 & 33.45 & 1.5315456 & 0 \\ 
38 & 9 & 47.05 & 34.50 & 1.9322024 & 0 \\ 
38 & 9 & 47.15 & 33.30 & 1.9363090 & 0 \\ 
39 & 8 & 47.60 & 35.00 & 1.5445248 & 0 \\ 
39 & 8 & 48.85 & 35.10 & 1.5850848 & 0 \\ 
40 & 10 & 47.20 & 33.50 & 2.3930400 & 0 \\ 
40 & 10 & 48.00 & 34.65 & 2.4336000 & 0 \\ 
41 & 11 & 47.60 & 33.80 & 2.9201172 & 1 \\ 
41 & 11 & 44.30 & 31.60 & 2.7176721 & 1 \\ 
42 & 7 & 50.10 & 33.70 & 1.2446343 & 1 \\ 
42 & 7 & 49.25 & 33.15 & 1.2235178 & 1 \\ 
43 & 7 & 51.40 & 33.40 & 1.2769302 & 1 \\ 
43 & 7 & 50.10 & 33.35 & 1.2446343 & 1 \\ 
44 & 11 & 49.20 & 33.30 & 3.0182724 & 1 \\ 
44 & 11 & 49.60 & 33.70 & 3.0428112 & 1 \\ 
45 & 10 & 49.25 & 33.45 & 2.4969750 & 0 \\ 
45 & 10 & 48.25 & 33.90 & 2.4462750 & 0 \\ 
46 & 9 & 51.10 & 34.90 & 2.0985237 & 1 \\ 
46 & 9 & 50.00 & 34.80 & 2.0533500 & 1 \\ 
47 & 11 & 51.50 & 33.60 & 3.1593705 & 1 \\ 
47 & 11 & 51.50 & 34.40 & 3.1593705 & 1 \\ 
48 & 9 & 47.95 & 31.80 & 1.9691626 & 1 \\ 
48 & 9 & 46.25 & 31.55 & 1.8993487 & 1 \\ 
49 & 9 & 46.55 & 32.60 & 1.9116689 & 1 \\ 
49 & 9 & 46.70 & 32.50 & 1.9178289 & 1 \\ 
50 & 8 & 49.35 & 33.20 & 1.6013088 & 0 \\ 
50 & 8 & 48.80 & 33.75 & 1.5834624 & 0 \\ 
51 & 10 & 48.95 & 35.00 & 2.4817650 & 1 \\ 
51 & 10 & 48.50 & 34.20 & 2.4589500 & 1 \\ 
52 & 11 & 50.75 & 34.85 & 3.1133602 & 0 \\ 
52 & 11 & 48.45 & 33.50 & 2.9722622 & 0 \\ 
53 & 9 & 50.15 & 31.70 & 2.0595101 & 0 \\ 
53 & 9 & 50.25 & 31.85 & 2.0636167 & 0 \\ 
54 & 10 & 47.85 & 32.90 & 2.4259950 & 0 \\ 
54 & 10 & 47.85 & 33.50 & 2.4259950 & 0 \\ 
55 & 8 & 48.85 & 33.10 & 1.5850848 & 0 \\ 
55 & 8 & 48.35 & 32.60 & 1.5688608 & 0 \\ 
56 & 9 & 46.60 & 33.30 & 1.9137222 & 1 \\ 
56 & 9 & 46.00 & 31.40 & 1.8890820 & 1 \\ 
57 & 10 & 52.10 & 35.15 & 2.6414700 & 1 \\ 
57 & 10 & 51.90 & 35.15 & 2.6313300 & 1 \\ 
58 & 6 & 48.55 & 33.15 & 0.8861346 & 0 \\ 
58 & 6 & 47.65 & 33.50 & 0.8697078 & 0 \\ 
59 & 6 & 48.50 & 33.35 & 0.8852220 & 1 \\ 
59 & 6 & 47.75 & 32.80 & 0.8715330 & 1 \\ 
60 & 6 & 50.10 & 34.30 & 0.9144252 & 1 \\ 
60 & 6 & 48.80 & 34.10 & 0.8906976 & 1 \\ 
61 & 8 & 46.90 & 32.75 & 1.5218112 & 0 \\ 
61 & 8 & 45.20 & 33.00 & 1.4666496 & 0 \\ 
62 & 8 & 49.20 & 33.30 & 1.5964416 & 0 \\ 
62 & 8 & 49.25 & 33.20 & 1.5980640 & 0 \\ 
63 & 8 & 46.70 & 30.80 & 1.5153216 & 0 \\ 
63 & 8 & 45.60 & 30.40 & 1.4796288 & 0 \\ 
64 & 9 & 49.15 & 33.40 & 2.0184431 & 0 \\ 
64 & 9 & 50.30 & 33.00 & 2.0656701 & 0 \\ 
65 & 6 & 44.35 & 31.80 & 0.8094762 & 1 \\ 
65 & 6 & 44.80 & 32.50 & 0.8176896 & 1 \\ 
66 & 7 & 47.65 & 31.75 & 1.1837689 & 0 \\ 
66 & 7 & 49.30 & 32.85 & 1.2247599 & 0 \\ 
67 & 5 & 51.80 & 34.35 & 0.6565650 & 0 \\ 
67 & 5 & 52.35 & 34.05 & 0.6635363 & 0 \\ 
68 & 9 & 51.05 & 34.15 & 2.0964704 & 0 \\ 
68 & 9 & 49.50 & 33.65 & 2.0328165 & 0 \\ 
69 & 8 & 46.90 & 33.85 & 1.5218112 & 0 \\ 
69 & 8 & 48.50 & 33.45 & 1.5737280 & 0 \\ 
70 & 12 & 47.95 & 34.15 & 3.5007336 & 1 \\ 
70 & 12 & 48.30 & 34.65 & 3.5262864 & 1 \\ 
71 & 10 & 48.80 & 34.70 & 2.4741600 & 0 \\ 
71 & 10 & 50.75 & 34.45 & 2.5730250 & 0 \\ 
72 & 9 & 49.45 & 35.25 & 2.0307631 & 1 \\ 
72 & 9 & 46.90 & 32.65 & 1.9260423 & 1 \\ 
73 & 9 & 50.50 & 33.95 & 2.0738835 & 0 \\ 
73 & 9 & 51.65 & 34.25 & 2.1211106 & 0 \\ 
74 & 8 & 47.10 & 34.20 & 1.5283008 & 1 \\ 
74 & 8 & 45.95 & 33.90 & 1.4909856 & 1 \\ 
75 & 13 & 45.75 & 34.05 & 3.9199973 & 0 \\ 
75 & 13 & 45.68 & 33.12 & 3.9139994 & 0 \\ 
76 & 9 & 51.40 & 35.35 & 2.1108438 & 1 \\ 
76 & 9 & 50.35 & 34.70 & 2.0677234 & 1 \\ 
77 & 10 & 49.30 & 33.80 & 2.4995100 & 0 \\ 
77 & 10 & 48.50 & 33.50 & 2.4589500 & 0 \\ 
78 & 8 & 48.95 & 33.30 & 1.5883296 & 0 \\ 
78 & 8 & 48.80 & 33.70 & 1.5834624 & 0 \\ 
79 & 11 & 49.45 & 33.60 & 3.0336092 & 1 \\ 
79 & 11 & 49.52 & 33.10 & 3.0379034 & 1 \\ 
80 & 8 & 50.00 & 34.85 & 1.6224000 & 0 \\ 
80 & 8 & 47.30 & 32.20 & 1.5347904 & 0 \\ 
81 & 7 & 47.90 & 35.00 & 1.1899797 & 0 \\ 
81 & 7 & 50.10 & 34.50 & 1.2446343 & 0 \\ 
82 & 7 & 47.75 & 32.35 & 1.1862533 & 0 \\ 
82 & 7 & 47.95 & 32.80 & 1.1912219 & 0 \\ 
83 & 7 & 46.35 & 31.80 & 1.1514731 & 1 \\ 
83 & 7 & 46.55 & 31.95 & 1.1564417 & 1 \\ 
84 & 5 & 45.91 & 33.19 & 0.5819093 & 1 \\ 
84 & 5 & 45.21 & 32.42 & 0.5730367 & 1 \\ 
85 & 5 & 50.60 & 33.90 & 0.6413550 & 0 \\ 
85 & 5 & 51.85 & 34.10 & 0.6571988 & 0 \\ 
86 & 9 & 52.04 & 34.32 & 2.1371267 & 1 \\ 
86 & 9 & 51.05 & 34.47 & 2.0964704 & 1 \\ 
87 & 10 & 49.70 & 33.10 & 2.5197900 & 0 \\ 
87 & 10 & 48.05 & 33.25 & 2.4361350 & 0 \\ 
88 & 9 & 45.17 & 32.69 & 1.8549964 & 0 \\ 
88 & 9 & 46.32 & 32.14 & 2.8415930 & 0 \\ 
89 & 11 & 49.35 & 33.50 & 3.0274745 & 0 \\ 
89 & 11 & 48.95 & 32.19 & 3.0029356 & 0 \\ 
90 & 9 & 46.72 & 32.19 & 1.9186502 & 0 \\ 
90 & 9 & 49.92 & 33.81 & 2.0500646 & 0 \\ 
91 & 9 & 50.81 & 33.52 & 2.0866143 & 0 \\ 
91 & 9 & 48.25 & 33.59 & 1.9814828 & 0 \\ 
92 & 10 & 47.05 & 32.48 & 2.3854350 & 0 \\ 
92 & 10 & 48.18 & 33.01 & 2.4427260 & 0 \\ 
93 & 7 & 47.95 & 34.07 & 1.1912219 & 0 \\ 
93 & 7 & 47.33 & 34.79 & 1.1758192 & 0 \\ 
94 & 10 & 50.73 & 33.98 & 2.5720110 & 0 \\ 
94 & 10 & 49.69 & 32.78 & 2.5192830 & 0 \\ 
95 & 7 & 49.23 & 33.82 & 1.2230209 & 1 \\ 
95 & 7 & 49.52 & 32.22 & 1.2302254 & 1 \\ 
96 & 6 & 49.22 & 33.71 & 0.8983634 & 1 \\ 
96 & 6 & 47.25 & 33.40 & 0.8624070 & 1 \\ 
97 & 5 & 46.75 & 34.05 & 0.5925563 & 1 \\ 
97 & 5 & 47.00 & 34.70 & 0.5957250 & 1 \\ 
98 & 9 & 48.09 & 33.82 & 1.9749120 & 1 \\ 
98 & 9 & 48.27 & 32.98 & 1.9823041 & 1 \\ 
99 & 11 & 48.70 & 33.70 & 2.9875989 & 1 \\ 
99 & 11 & 48.71 & 33.89 & 2.9882124 & 1 \\ 
100 & 10 & 47.70 & 32.65 & 2.4183900 & 1 \\ 
100 & 10 & 45.67 & 31.61 & 2.3154690 & 1 \\ 
101 & 12 & 46.78 & 33.49 & 3.4153142 & 0 \\ 
101 & 12 & 47.19 & 33.18 & 3.4452475 & 0 \\ 
102 & 5 & 49.67 & 33.65 & 0.6295673 & 0 \\ 
102 & 5 & 47.65 & 32.66 & 0.6039638 & 0 \\ 
103 & 7 & 48.03 & 32.85 & 1.1932093 & 0 \\ 
103 & 7 & 43.33 & 30.94 & 1.0764472 & 0 \\ 
104 & 7 & 46.22 & 32.45 & 1.1482435 & 1 \\ 
104 & 7 & 43.95 & 31.25 & 1.0918499 & 1 \\ 
105 & 9 & 48.12 & 33.18 & 1.9761440 & 0 \\ 
105 & 9 & 48.67 & 33.42 & 1.9987309 & 0 \\ 
106 & 11 & 48.36 & 34.35 & 2.9667409 & 0 \\ 
106 & 11 & 49.83 & 34.35 & 3.0569210 & 0 \\ 
107 & 9 & 51.96 & 35.27 & 2.1338413 & 0 \\ 
107 & 9 & 49.96 & 34.55 & 2.0517073 & 0 \\ 
108 & 9 & 47.13 & 32.58 & 1.9354877 & 0 \\ 
108 & 9 & 48.31 & 33.31 & 1.9839468 & 0 \\ 
109 & 10 & 46.65 & 31.56 & 2.3651550 & 0 \\ 
109 & 10 & 48.46 & 33.90 & 2.4569220 & 0 \\ 
110 & 9 & 46.02 & 33.58 & 1.8899033 & 0 \\ 
110 & 9 & 48.07 & 33.38 & 1.9740907 & 0 \\ 
111 & 10 & 51.42 & 35.49 & 2.6069940 & 0 \\ 
111 & 10 & 51.16 & 35.52 & 2.5938120 & 0 \\ 
112 & 12 & 49.32 & 34.43 & 3.6007546 & 0 \\ 
112 & 12 & 50.08 & 34.39 & 3.6562406 & 0 \\ 
113 & 8 & 50.92 & 34.55 & 1.6522522 & 0 \\ 
113 & 8 & 50.70 & 35.84 & 1.6451136 & 0 \\ 
114 & 8 & 50.14 & 35.23 & 1.6269427 & 0 \\ 
114 & 8 & 51.42 & 35.02 & 1.6684762 & 0 \\ 
115 & 11 & 46.82 & 33.50 & 2.8722665 & 0 \\ 
115 & 11 & 46.75 & 33.55 & 2.8679722 & 0 \\ 
116 & 9 & 51.88 & 34.77 & 2.1305560 & 0 \\ 
116 & 9 & 51.28 & 35.68 & 2.1059158 & 0 \\ 
117 & 10 & 51.36 & 36.19 & 2.6039520 & 0 \\ 
117 & 10 & 50.49 & 35.75 & 2.5598430 & 0 \\ 
118 & 8 & 49.32 & 34.10 & 1.6003354 & 0 \\ 
118 & 8 & 49.65 & 34.17 & 1.6110432 & 0 \\ 
119 & 13 & 46.52 & 33.59 & 3.9859732 & 0 \\ 
119 & 13 & 45.92 & 32.69 & 3.9345634 & 0 \\ 
120 & 11 & 51.58 & 35.22 & 3.1642783 & 0 \\ 
120 & 11 & 51.28 & 35.12 & 3.1458742 & 0 \\ 
121 & 9 & 50.69 & 34.18 & 2.0816862 & 0 \\ 
121 & 9 & 49.84 & 34.52 & 2.0467793 & 0 \\ 
122 & 12 & 49.53 & 34.22 & 3.6160862 & 0 \\ 
122 & 12 & 50.53 & 30.49 & 3.6890942 & 0 \\ 
123 & 10 & 50.74 & 35.28 & 2.5725180 & 0 \\ 
123 & 10 & 51.82 & 36.17 & 2.6272740 & 0 \\ 
124 & 11 & 47.16 & 34.45 & 2.8931245 & 0 \\ 
124 & 11 & 48.36 & 34.02 & 2.9667409 & 0 \\ 
125 & 11 & 49.65 & 35.52 & 3.0458785 & 0 \\ 
125 & 11 & 50.02 & 35.23 & 3.0685769 & 0 \\ 
126 & 10 & 49.24 & 32.18 & 2.4964680 & 0 \\ 
126 & 10 & 50.22 & 32.59 & 2.5461540 & 0 \\ 
127 & 12 & 47.98 & 33.98 & 3.5029238 & 0 \\ 
127 & 12 & 47.82 & 33.22 & 3.4912426 & 0 \\ 
128 & 11 & 49.82 & 32.76 & 3.0563075 & 0 \\ 
128 & 11 & 47.50 & 32.05 & 2.9139825 & 0 \\ 
129 & 9 & 46.49 & 32.20 & 1.9092048 & 0 \\ 
129 & 9 & 47.62 & 33.67 & 1.9556105 & 0 \\ 
130 & 9 & 52.78 & 35.49 & 2.1675163 & 0 \\ 
130 & 9 & 52.10 & 35.69 & 2.1395907 & 0 \\ 
131 & 12 & 48.34 & 32.15 & 3.5292067 & 0 \\ 
131 & 12 & 48.68 & 33.98 & 3.5540294 & 0 \\ 
132 & 10 & 52.73 & 34.52 & 2.6734110 & 0 \\ 
132 & 10 & 51.92 & 33.65 & 2.6323440 & 0 \\ 
133 & 11 & 48.35 & 33.76 & 2.9661275 & 0 \\ 
133 & 11 & 46.54 & 33.22 & 2.8550894 & 0 \\ 
134 & 8 & 49.16 & 33.13 & 1.5951437 & 0 \\ 
134 & 8 & 48.63 & 32.88 & 1.5779462 & 0 \\ 
135 & 11 & 46.69 & 33.55 & 2.8642914 & 0 \\ 
135 & 11 & 47.05 & 32.81 & 2.8863763 & 0 \\ 
136 & 9 & 51.09 & 34.18 & 2.0981130 & 0 \\ 
136 & 9 & 51.07 & 33.88 & 2.0972917 & 0 \\ 
137 & 10 & 46.54 & 33.27 & 2.3595780 & 0 \\ 
137 & 10 & 46.62 & 32.70 & 2.3636340 & 0 \\ 
138 & 10 & 52.52 & 34.42 & 2.6627640 & 0 \\ 
138 & 10 & 50.98 & 35.08 & 2.5846860 & 0 \\ 
139 & 8 & 48.92 & 33.91 & 1.5873562 & 0 \\ 
139 & 8 & 47.92 & 34.07 & 1.5549082 & 0 \\ 
140 & 9 & 50.40 & 36.08 & 2.0697768 & 0 \\ 
140 & 9 & 49.28 & 35.52 & 2.0237818 & 0 \\ 
141 & 11 & 46.24 & 32.42 & 2.8366853 & 0 \\ 
141 & 11 & 47.41 & 32.75 & 2.9084613 & 0 \\ 
142 & 11 & 46.72 & 35.32 & 2.8661318 & 0 \\ 
142 & 11 & 47.32 & 34.81 & 2.9029400 & 0 \\ 
143 & 11 & 47.63 & 34.02 & 2.9219576 & 0 \\ 
143 & 11 & 48.17 & 33.95 & 2.9550850 & 0 \\ 
144 & 12 & 50.82 & 35.14 & 3.7102666 & 0 \\ 
144 & 12 & 51.41 & 35.46 & 3.7533413 & 0 \\ 
145 & 10 & 50.83 & 35.97 & 2.5770810 & 0 \\ 
145 & 10 & 50.92 & 35.72 & 2.5816440 & 0 \\ 
146 & 10 & 49.90 & 34.92 & 2.5299300 & 0 \\ 
146 & 10 & 53.07 & 33.58 & 2.6906490 & 0 \\ 
147 & 9 & 49.42 & 36.25 & 2.0295311 & 0 \\ 
147 & 9 & 49.84 & 34.88 & 2.0467793 & 0 \\ 
148 & 10 & 51.13 & 33.93 & 2.5922910 & 0 \\ 
148 & 10 & 51.75 & 30.49 & 2.6237250 & 0 \\ 
149 & 11 & 51.14 & 33.65 & 3.1372856 & 0 \\ 
149 & 11 & 50.08 & 33.35 & 3.0722578 & 0 \\ 
150 & 11 & 48.94 & 33.37 & 3.0023222 & 0 \\ 
150 & 11 & 46.45 & 32.49 & 2.8495682 & 0 \\ 
151 & 7 & 47.09 & 33.18 & 1.1698569 & 0 \\ 
151 & 7 & 47.41 & 33.09 & 1.1778066 & 0 \\ 
152 & 8 & 47.96 & 33.82 & 1.5562061 & 0 \\ 
152 & 8 & 48.92 & 34.13 & 1.5873562 & 0 \\ 
153 & 9 & 47.80 & 33.80 & 1.9630026 & 0 \\ 
153 & 9 & 48.21 & 33.38 & 1.9798401 & 0 \\ 
154 & 10 & 46.72 & 33.61 & 2.3687040 & 0 \\ 
154 & 10 & 48.99 & 34.80 & 2.4837930 & 0 \\ 
155 & 9 & 41.00 & 28.32 & 1.6837470 & 0 \\ 
155 & 9 & 47.87 & 33.05 & 1.9658773 & 0 \\ 
156 & 10 & 49.01 & 35.12 & 2.4848070 & 0 \\ 
156 & 10 & 52.11 & 34.54 & 2.6419770 & 0 \\ 
157 & 10 & 44.02 & 33.06 & 2.2318140 & 0 \\ 
157 & 10 & 48.29 & 33.89 & 2.4483030 & 0 \\ 
158 & 9 & 50.23 & 33.25 & 2.0627954 & 0 \\ 
158 & 9 & 47.13 & 32.55 & 1.9354877 & 0 \\ 
159 & 9 & 41.44 & 29.92 & 1.7018165 & 0 \\ 
159 & 9 & 43.92 & 31.90 & 1.8036626 & 0 \\ 
160 & 11 & 47.17 & 33.82 & 2.8937380 & 0 \\ 
160 & 11 & 50.32 & 32.45 & 3.0869810 & 0 \\ 
161 & 9 & 46.68 & 32.64 & 1.9170076 & 0 \\ 
161 & 9 & 47.92 & 34.55 & 1.9679306 & 0 \\ 
162 & 8 & 46.36 & 33.34 & 1.5042893 & 0 \\ 
162 & 8 & 49.66 & 32.96 & 1.6113677 & 0 \\ 
163 & 8 & 50.48 & 35.02 & 1.6379750 & 0 \\ 
163 & 8 & 47.96 & 33.71 & 1.5562061 & 0 \\ 
164 & 7 & 47.50 & 33.08 & 1.1800425 & 0 \\ 
164 & 7 & 47.61 & 33.48 & 1.1827752 & 0 \\ 
165 & 10 & 45.02 & 33.06 & 2.2825140 & 0 \\ 
165 & 10 & 51.10 & 33.43 & 2.5907700 & 0 \\ 
166 & 8 & 50.62 & 34.35 & 1.6425178 & 0 \\ 
166 & 8 & 51.20 & 34.45 & 1.6613376 & 0 \\ 
167 & 8 & 49.56 & 31.26 & 1.6081229 & 0 \\ 
167 & 8 & 49.37 & 31.71 & 1.6019578 & 0 \\ 
168 & 6 & 48.48 & 33.67 & 0.8848570 & 0 \\ 
168 & 6 & 48.22 & 33.48 & 0.8801114 & 0 \\ 
169 & 10 & 47.25 & 33.10 & 2.3955750 & 0 \\ 
169 & 10 & 46.85 & 32.80 & 2.3752950 & 0 \\ 
170 & 10 & 47.85 & 34.05 & 2.4259950 & 0 \\ 
170 & 10 & 50.10 & 33.50 & 2.5400700 & 0 \\ 
171 & 10 & 50.90 & 33.55 & 2.5806300 & 0 \\ 
171 & 10 & 48.30 & 31.90 & 2.4488100 & 0 \\ 
172 & 11 & 48.75 & 33.95 & 2.9906662 & 0 \\ 
172 & 11 & 48.45 & 33.45 & 2.9722622 & 0 \\ 
173 & 10 & 46.25 & 33.10 & 2.3448750 & 0 \\ 
173 & 10 & 46.30 & 33.30 & 2.3474100 & 0 \\ 
174 & 12 & 48.75 & 32.95 & 3.5591400 & 0 \\ 
174 & 12 & 50.40 & 32.80 & 3.6796032 & 0 \\ 
175 & 10 & 48.40 & 32.60 & 2.4538800 & 0 \\ 
175 & 10 & 47.60 & 31.90 & 2.4133200 & 0 \\ 
176 & 12 & 51.90 & 34.30 & 3.7891152 & 0 \\ 
176 & 12 & 50.90 & 34.50 & 3.7161072 & 0 \\ 
177 & 10 & 49.60 & 34.30 & 2.5147200 & 0 \\ 
177 & 10 & 50.10 & 35.70 & 2.5400700 & 0 \\ 
178 & 11 & 51.85 & 35.70 & 3.1808420 & 0 \\ 
178 & 11 & 49.25 & 34.25 & 3.0213398 & 0 \\ 
179 & 9 & 47.90 & 33.00 & 1.9671093 & 0 \\ 
179 & 9 & 46.60 & 32.10 & 1.9137222 & 0 \\ 
180 & 9 & 47.25 & 33.50 & 1.9404158 & 0 \\ 
180 & 9 & 47.55 & 34.10 & 1.9527359 & 0 \\ 
181 & 12 & 47.70 & 33.55 & 3.4824816 & 0 \\ 
181 & 12 & 46.90 & 33.65 & 3.4240752 & 0 \\ 
182 & 13 & 49.30 & 33.95 & 4.2241719 & 0 \\ 
182 & 13 & 49.20 & 33.50 & 4.2156036 & 0 \\ 
183 & 13 & 52.30 & 33.95 & 4.4812209 & 0 \\ 
183 & 13 & 50.75 & 34.00 & 4.3484122 & 0 \\ 
184 & 12 & 47.75 & 33.30 & 3.4861320 & 0 \\ 
184 & 12 & 47.70 & 32.65 & 3.4824816 & 0 \\ 
\bottomrule
\end{tabular*}
\end{table}

:::
:::



We apply the exact same ratio estimator formulation, though mathematically, the cluster total $t_i$ is now an *estimated* total $\hat{t}_i$.


::: {.cell}
::: {.cell-output-display}
\begin{table}
\caption*{
{\fontsize{20}{25}\selectfont  Ratio Estimation: Working Table\fontsize{12}{15}\selectfont } \\ 
{\fontsize{14}{17}\selectfont  \(n = 184\)\fontsize{12}{15}\selectfont }
} 
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}>{\raggedright\arraybackslash}p{\dimexpr 30.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 75.00pt -2\tabcolsep-1.5\arrayrulewidth}}
\toprule
\(i\) & \(\bar y_i\) & \(\hat t_i\)\textsuperscript{\textit{1}} & \(M_i\) & \(\hat B M_i\)\textsuperscript{\textit{2}} & \(e_i\) & \(e_i^2\) \\ 
\midrule\addlinespace[2.5pt]
1 & 3.864 & 50.236 & 13.000 & 32.378 & 1.786 $\times$ 10\textsuperscript{1} & 318.923 \\ 
2 & 4.194 & 54.524 & 13.000 & 32.378 & 2.215 $\times$ 10\textsuperscript{1} & 490.483 \\ 
3 & 0.916 & 5.498 & 6.000 & 14.943 & -9.446 & 89.226 \\ 
4 & 2.998 & 32.982 & 11.000 & 27.396 & 5.585 & 31.196 \\ 
5 & 2.496 & 24.957 & 10.000 & 24.906 & 5.129 $\times$ 10\textsuperscript{-2} & 0.003 \\ 
6 & 3.984 & 51.795 & 13.000 & 32.378 & 1.942 $\times$ 10\textsuperscript{1} & 377.053 \\ 
... & ... & ... & ... & ... & ... & ... \\ 
{\bfseries Sum} & {\bfseries ...} & {\bfseries 4,375.947} & {\bfseries 1,757.000} & {\bfseries 4,375.947} & {\bfseries 4.476 $\times$ 10\textsuperscript{-13}} & {\bfseries 11,439.579} \\ 
\bottomrule
\end{tabular*}
\begin{minipage}{\linewidth}
\vspace{.05em}
\textsuperscript{\textit{1}}\(\hat t_i = \bar y_i \, M_i\).\\
\textsuperscript{\textit{2}}\(\hat B = \sum_i \hat t_i / \sum_i M_i = 2.4906\), \(\hat B M_i = \hat B\, M_i\), \(s_e^2 = 62.5114\).\\
\textbf{Point estimate:} \(\hat B = \dfrac{\sum_i \hat t_i}{\sum_i M_i} = \dfrac{4,375.947}{1,757.000} = 2.491\)\\
\textbf{Standard error:} \(\mathrm{SE}(\hat B) = \dfrac{1}{\bar x}\sqrt{\dfrac{s_e^2}{n}} = \dfrac{1}{9.549}\sqrt{\dfrac{62.511}{184}} = 0.061\)\\
\end{minipage}
\end{table}

:::
:::

**Output Note:** Unlike the previous example, we only measured a subset of eggs per clutch. The function automatically estimates the total volume of each clutch ($\hat{t}_i$) by multiplying the true clutch size ($M_i$) by the observed subsample mean ($\overline{y}_i$); the working table (truncated to the first 6 rows + Sum) shows these estimated totals alongside the ratio fit and residuals.


::: {.cell}
::: {.cell-output-display}
\begin{table}
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{0.7\linewidth}{@{\extracolsep{\fill}}rrrr}
\toprule
\(\bar{y}_{r}\) & \(\mathrm{SE}(\bar{y}_{r})\) & 95\% CI Lower & 95\% CI Upper \\ 
\midrule\addlinespace[2.5pt]
2.4906 & 0.0610 & 2.3701 & 2.6110 \\ 
\bottomrule
\end{tabular*}
\end{table}

:::
:::

**Output Note:** The estimated overall average egg volume is 2.491 units. Because the total number of clutches in the population ($N$) is not explicitly provided (set to `Inf` by default), the finite population correction factor is safely omitted from the standard error calculation.