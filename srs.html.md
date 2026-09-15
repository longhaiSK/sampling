---
title: "Simple Random Sampling"
format: html
filters:
  - shinylive
execute:
  warning: false
  message: false
  cache: false
  keep-md: true
---


## Analysis of agsrs.csv Data

### Step by step calculation without using a function

To evaluate a Simple Random Sample (SRS), we calculate the point estimate, the standard error (incorporating the Finite Population Correction factor), and the margin of error to construct a 95% confidence interval.

**Key Formulas:**

* **Sample Mean:** $\bar{y} = \frac{1}{n} \sum_{i=1}^n y_i$
* **Standard Error:** $\mathrm{SE}(\bar{y}) = \sqrt{1 - \frac{n}{N}} \left( \frac{s}{\sqrt{n}} \right)$
* **Margin of Error:** $\mathrm{ME} = t_{\alpha/2, n-1} \times \mathrm{SE}(\bar{y})$
* **Population Total:** $\hat{t} = N\bar{y}$


::: {.cell}

```{.r .cell-code}
## read survey data
agsrs <- read.csv ("data/agsrs.csv")
head(agsrs)
```

::: {.cell-output-display}
`````{=html}
<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["county"],"name":[1],"type":["chr"],"align":["left"]},{"label":["state"],"name":[2],"type":["chr"],"align":["left"]},{"label":["acres92"],"name":[3],"type":["int"],"align":["right"]},{"label":["acres87"],"name":[4],"type":["int"],"align":["right"]},{"label":["acres82"],"name":[5],"type":["int"],"align":["right"]},{"label":["farms92"],"name":[6],"type":["int"],"align":["right"]},{"label":["farms87"],"name":[7],"type":["int"],"align":["right"]},{"label":["farms82"],"name":[8],"type":["int"],"align":["right"]},{"label":["largef92"],"name":[9],"type":["int"],"align":["right"]},{"label":["largef87"],"name":[10],"type":["int"],"align":["right"]},{"label":["largef82"],"name":[11],"type":["int"],"align":["right"]},{"label":["smallf92"],"name":[12],"type":["int"],"align":["right"]},{"label":["smallf87"],"name":[13],"type":["int"],"align":["right"]},{"label":["smallf82"],"name":[14],"type":["int"],"align":["right"]},{"label":["region"],"name":[15],"type":["chr"],"align":["left"]}],"data":[{"1":"COFFEE COUNTY","2":"AL","3":"175209","4":"179311","5":"194509","6":"760","7":"842","8":"944","9":"29","10":"28","11":"21","12":"57","13":"47","14":"66","15":"S","_rn_":"1"},{"1":"COLBERT COUNTY","2":"AL","3":"138135","4":"145104","5":"161360","6":"488","7":"563","8":"686","9":"37","10":"41","11":"42","12":"12","13":"44","14":"47","15":"S","_rn_":"2"},{"1":"LAMAR COUNTY","2":"AL","3":"56102","4":"59861","5":"72334","6":"299","7":"362","8":"447","9":"4","10":"4","11":"3","12":"16","13":"20","14":"30","15":"S","_rn_":"3"},{"1":"MARENGO COUNTY","2":"AL","3":"199117","4":"220526","5":"231207","6":"434","7":"471","8":"622","9":"48","10":"66","11":"62","12":"14","13":"11","14":"28","15":"S","_rn_":"4"},{"1":"MARION COUNTY","2":"AL","3":"89228","4":"105586","5":"113618","6":"566","7":"658","8":"748","9":"7","10":"9","11":"9","12":"11","13":"23","14":"27","15":"S","_rn_":"5"},{"1":"TUSCALOOSA COUNTY","2":"AL","3":"96194","4":"120542","5":"134616","6":"436","7":"521","8":"650","9":"20","10":"17","11":"23","12":"18","13":"32","14":"29","15":"S","_rn_":"6"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>
`````
:::

```{.r .cell-code}
## extract the variable of interest
sdata <- agsrs$acres92
N <- 3078

## do calculation
n <- length (sdata)
ybar <- mean (sdata)
se.ybar <- sqrt((1 - n / N)) * sd (sdata) / sqrt(n)  
mem <- qt (0.975, df = n - 1) * se.ybar

## return estimate vector for pop mean
c (Est. = ybar, S.E. = se.ybar, ci.low = ybar - mem, ci.upp = ybar + mem)
```

::: {.cell-output .cell-output-stdout}

```
     Est.      S.E.    ci.low    ci.upp 
297897.05  18898.43 260706.26 335087.84 
```


:::

```{.r .cell-code}
## return estimate vector for pop total
c (Est. = ybar, S.E. = se.ybar, ci.low = ybar - mem, ci.upp = ybar + mem) * N
```

::: {.cell-output .cell-output-stdout}

```
      Est.       S.E.     ci.low     ci.upp 
 916927110   58169381  802453859 1031400361 
```


:::
:::


### Write a function for repeated use

#### A function for doing data analysis for srs sample

Wrap the previous manual calculations into a custom R function to quickly estimate means, standard errors, and confidence intervals for any variable of interest. To avoid repeating this code in every chapter, all estimating functions (and their `gt`-table formatting helpers) live in a single shared file, `samplingestimate.r`, which we source below.


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
                      show.details = TRUE, col_labels = list (y = "y_i", ybar="\\bar y", dev="e_i"))
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
        ybar = rep (ybar, n + 1),
        dev  = c (dev, sum (dev)),
        dev2 = c (dev^2, sum (dev^2))
    )
    working <- truncate_working (working, id_col = "i")

    table <- gt (working) |>
        tab_header (title = "SRS Mean Estimation: Working Table", subtitle = md (sprintf ("$n = %d$", n))) |>
        cols_label (
            i    = md ("$i$"),
            y    = md (sprintf ("$%s$", col_labels$y)),
            ybar    = md (sprintf ("$%s$", col_labels$ybar)),
            dev  = md (sprintf ("$%s$", col_labels$dev)),
            dev2 = md (sprintf ("$(%s-\\bar y)^2$", col_labels$y))
        ) |>
        cols_width (i ~ px (40), everything () ~ px (110)) |>
        fmt_auto (data = working, columns = c ("y", "ybar", "dev", "dev2"), decimals = 3) |>
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
                sprintf ("$\\hat t = N\\bar y = %s \\times %s = %s$",
                         fmt_num (N, 0), fmt_num (ybar), fmt_num (est))
            else
                sprintf ("$\\bar y = \\dfrac{\\sum_i %s}{n} = \\dfrac{%s}{%d} = %s$",
                         col_labels$y, fmt_num (sum (sdata)), n, fmt_num (ybar)))
        ) |>
        tab_source_note (
            source_note = md (if (estimate == "total")
                sprintf ("$\\mathrm{SE}(\\hat t) = N \\times \\mathrm{SE}(\\bar y) = %s \\times %s = %s$",
                         fmt_num (N, 0), fmt_num (se.ybar), fmt_num (se.est))
            else if (is.finite (N))
                sprintf ("$\\mathrm{SE}(\\bar y) = \\sqrt{1-\\dfrac{n}{N}}\\,\\dfrac{s}{\\sqrt n} = \\sqrt{1-\\dfrac{%d}{%d}}\\,\\dfrac{%s}{\\sqrt{%d}} = %s$",
                         n, N, fmt_num (sqrt (s2y)), n, fmt_num (se.ybar))
            else
                sprintf ("$\\mathrm{SE}(\\bar y) = \\dfrac{s}{\\sqrt n} = \\dfrac{%s}{\\sqrt{%d}} = %s$",
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
                        extra_col = NULL, B_label = "\\hat B")
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
              sprintf ("$%s = \\dfrac{\\sum_i %s}{\\sum_i %s} = \\dfrac{%s}{%s} = %s$",
                       B_label, col_labels$y, col_labels$x, fmt_num (sum (ydata)), fmt_num (sum (xdata)), fmt_num (B_hat))
          else if (estimate == "total")
              sprintf ("$\\hat t = \\hat B\\,\\bar x_U\\,N = %s \\times %s \\times %s = %s$",
                       fmt_num (B_hat), fmt_num (xbarU), fmt_num (N, 0), fmt_num (est))
          else
              sprintf ("$\\bar y_r = \\hat B\\,\\bar x_U = %s \\times %s = %s$",
                       fmt_num (B_hat), fmt_num (xbarU), fmt_num (est)))
      ) |>
      tab_source_note (
          source_note = md (if (estimate == "model")
              (if (is.finite (N))
                  sprintf ("$\\mathrm{SE}(%s) = \\dfrac{1}{\\bar x}\\sqrt{\\left(1-\\dfrac{n}{N}\\right)\\dfrac{s_e^2}{n}} = \\dfrac{1}{%s}\\sqrt{\\left(1-\\dfrac{%d}{%d}\\right)\\dfrac{%s}{%d}} = %s$",
                          B_label, fmt_num (xbar), n, N, fmt_num (var_e), n, fmt_num (sd_B_hat))
              else
                  sprintf ("$\\mathrm{SE}(%s) = \\dfrac{1}{\\bar x}\\sqrt{\\dfrac{s_e^2}{n}} = \\dfrac{1}{%s}\\sqrt{\\dfrac{%s}{%d}} = %s$",
                          B_label, fmt_num (xbar), fmt_num (var_e), n, fmt_num (sd_B_hat)))
          else if (estimate == "total")
              sprintf ("$\\mathrm{SE}(\\hat t) = \\mathrm{SE}(\\hat B)\\,\\bar x_U\\,N = %s \\times %s \\times %s = %s$",
                       fmt_num (sd_B_hat), fmt_num (xbarU), fmt_num (N, 0), fmt_num (sd_est))
          else
              sprintf ("$\\mathrm{SE}(\\bar y_r) = \\mathrm{SE}(\\hat B)\\,\\bar x_U = %s \\times %s = %s$",
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
              sprintf ("$\\hat B_1 = \\dfrac{\\sum_i(x_i-\\bar x)(y_i-\\bar y)}{\\sum_i(x_i-\\bar x)^2} = \\dfrac{%s}{%s} = %s$",
                       fmt_num (Sxy), fmt_num (Sxx), fmt_num (unname (Bhat[2])))
          else if (estimate == "total")
              sprintf ("$\\hat t_{reg} = N\\bar y_{reg} = %s \\times %s = %s$",
                       fmt_num (N, 0), fmt_num (yhat_reg), fmt_num (est))
          else
              sprintf ("$\\bar y_{reg} = \\hat B_0+\\hat B_1\\bar x_U = %s + %s \\times %s = %s$",
                       fmt_num (Bhat[1]), fmt_num (Bhat[2]), fmt_num (xbarU), fmt_num (yhat_reg)))
      ) |>
      tab_source_note (
          source_note = md (if (estimate == "model")
              sprintf ("$\\mathrm{SE}(\\hat B_1) = \\sqrt{\\dfrac{s_e^2}{\\sum_i(x_i-\\bar x)^2}} = \\sqrt{\\dfrac{%s}{%s}} = %s$",
                       fmt_num (SSe), fmt_num (Sxx), fmt_num (sd_est))
          else if (estimate == "total")
              sprintf ("$\\mathrm{SE}(\\hat t_{reg}) = N\\times\\mathrm{SE}(\\bar y_{reg}) = %s \\times %s = %s$",
                       fmt_num (N, 0), fmt_num (se_yhat_reg), fmt_num (sd_est))
          else if (is.finite (N))
              sprintf ("$\\mathrm{SE}(\\bar y_{reg}) = \\sqrt{\\left(1-\\dfrac{n}{N}\\right)\\dfrac{s_e^2}{n}} = \\sqrt{\\left(1-\\dfrac{%d}{%d}\\right)\\dfrac{%s}{%d}} = %s$",
                       n, N, fmt_num (SSe), n, fmt_num (se_yhat_reg))
          else
              sprintf ("$\\mathrm{SE}(\\bar y_{reg}) = \\sqrt{\\dfrac{s_e^2}{n}} = \\sqrt{\\dfrac{%s}{%d}} = %s$",
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
                sprintf ("$\\hat t_{str} = N\\bar y_{str} = %s \\times %s = %s$",
                         fmt_num (N, 0), fmt_num (ybar), fmt_num (est))
            else
                sprintf ("$\\bar y_{str} = \\sum_h \\pi_h\\bar y_h = %s$", fmt_num (ybar)))
        ) |>
        tab_source_note (
            source_note = md (if (estimate == "total")
                sprintf ("$\\mathrm{SE}(\\hat t_{str}) = N \\times \\mathrm{SE}(\\bar y_{str}) = %s \\times %s = %s$",
                         fmt_num (N, 0), fmt_num (seybar), fmt_num (se.est))
            else
                sprintf ("$\\mathrm{SE}(\\bar y_{str}) = \\sqrt{\\sum_h v_h} = \\sqrt{%s} = %s$",
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
                 extra_col = ybar_col, B_label = "\\bar y_r")
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

#### Apply srs_est to agsrs.csv data

*Import Data*

Load the survey sample data representing US counties.


::: {.cell}

```{.r .cell-code}
agsrs <- read.csv ("data/agsrs.csv")
```
:::


*Estimating the mean of acre92*

Calculate the estimated average farm acreage per county using the custom function. The function's `$table` shows the row-level working values (truncated to the first 6 rows + Sum), and `$estimate` holds the final point estimate, SE, and CI, formatted with `format_est_gt()`.


::: {.cell}

```{.r .cell-code}
res_mean <- srs_est(agsrs[,"acres92"], N = 3078, estimate = "mean")
res_mean$table
```

::: {.cell-output-display}

```{=html}
<div id="rlnyzttjpf" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#rlnyzttjpf table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#rlnyzttjpf thead, #rlnyzttjpf tbody, #rlnyzttjpf tfoot, #rlnyzttjpf tr, #rlnyzttjpf td, #rlnyzttjpf th {
  border-style: none;
}

#rlnyzttjpf p {
  margin: 0;
  padding: 0;
}

#rlnyzttjpf .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#rlnyzttjpf .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#rlnyzttjpf .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#rlnyzttjpf .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#rlnyzttjpf .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#rlnyzttjpf .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#rlnyzttjpf .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#rlnyzttjpf .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#rlnyzttjpf .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#rlnyzttjpf .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#rlnyzttjpf .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#rlnyzttjpf .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#rlnyzttjpf .gt_spanner_row {
  border-bottom-style: hidden;
}

#rlnyzttjpf .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#rlnyzttjpf .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#rlnyzttjpf .gt_from_md > :first-child {
  margin-top: 0;
}

#rlnyzttjpf .gt_from_md > :last-child {
  margin-bottom: 0;
}

#rlnyzttjpf .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#rlnyzttjpf .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#rlnyzttjpf .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#rlnyzttjpf .gt_row_group_first td {
  border-top-width: 2px;
}

#rlnyzttjpf .gt_row_group_first th {
  border-top-width: 2px;
}

#rlnyzttjpf .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#rlnyzttjpf .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#rlnyzttjpf .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#rlnyzttjpf .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#rlnyzttjpf .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#rlnyzttjpf .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#rlnyzttjpf .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#rlnyzttjpf .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#rlnyzttjpf .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#rlnyzttjpf .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#rlnyzttjpf .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#rlnyzttjpf .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#rlnyzttjpf .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#rlnyzttjpf .gt_left {
  text-align: left;
}

#rlnyzttjpf .gt_center {
  text-align: center;
}

#rlnyzttjpf .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#rlnyzttjpf .gt_font_normal {
  font-weight: normal;
}

#rlnyzttjpf .gt_font_bold {
  font-weight: bold;
}

#rlnyzttjpf .gt_font_italic {
  font-style: italic;
}

#rlnyzttjpf .gt_super {
  font-size: 65%;
}

#rlnyzttjpf .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#rlnyzttjpf .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#rlnyzttjpf .gt_indent_1 {
  text-indent: 5px;
}

#rlnyzttjpf .gt_indent_2 {
  text-indent: 10px;
}

#rlnyzttjpf .gt_indent_3 {
  text-indent: 15px;
}

#rlnyzttjpf .gt_indent_4 {
  text-indent: 20px;
}

#rlnyzttjpf .gt_indent_5 {
  text-indent: 25px;
}

#rlnyzttjpf .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#rlnyzttjpf div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" style="table-layout:fixed;width:0px;" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <colgroup>
    <col style="width:40px;"/>
    <col style="width:110px;"/>
    <col style="width:110px;"/>
    <col style="width:110px;"/>
    <col style="width:110px;"/>
  </colgroup>
  <thead>
    <tr class="gt_heading">
      <td colspan="5" class="gt_heading gt_title gt_font_normal" style>SRS Mean Estimation: Working Table</td>
    </tr>
    <tr class="gt_heading">
      <td colspan="5" class="gt_heading gt_subtitle gt_font_normal gt_bottom_border" style><span data-qmd-base64="JG4gPSAzMDAk"><span class='gt_from_md'>\(n = 300\)</span></span></td>
    </tr>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1" scope="col" id="i"><span data-qmd-base64="JGkk"><span class='gt_from_md'>\(i\)</span></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="y"><span data-qmd-base64="JHlfaSQ="><span class='gt_from_md'>\(y_i\)</span></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="ybar"><span data-qmd-base64="JFxiYXIgeSQ="><span class='gt_from_md'>\(\bar y\)</span></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="dev"><span data-qmd-base64="JGVfaSQ="><span class='gt_from_md'>\(e_i\)</span></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="dev2"><span data-qmd-base64="JCh5X2ktXGJhciB5KV4yJA=="><span class='gt_from_md'>\((y_i-\bar y)^2\)</span></span><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="i" class="gt_row gt_left">1</td>
<td headers="y" class="gt_row gt_right">1.752&nbsp;×&nbsp;10<sup style='font-size: 65%;'>5</sup></td>
<td headers="ybar" class="gt_row gt_right">2.979&nbsp;×&nbsp;10<sup style='font-size: 65%;'>5</sup></td>
<td headers="dev" class="gt_row gt_right">−1.227&nbsp;×&nbsp;10<sup style='font-size: 65%;'>5</sup></td>
<td headers="dev2" class="gt_row gt_right">1.505&nbsp;×&nbsp;10<sup style='font-size: 65%;'>10</sup></td></tr>
    <tr><td headers="i" class="gt_row gt_left">2</td>
<td headers="y" class="gt_row gt_right">1.381&nbsp;×&nbsp;10<sup style='font-size: 65%;'>5</sup></td>
<td headers="ybar" class="gt_row gt_right">2.979&nbsp;×&nbsp;10<sup style='font-size: 65%;'>5</sup></td>
<td headers="dev" class="gt_row gt_right">−1.598&nbsp;×&nbsp;10<sup style='font-size: 65%;'>5</sup></td>
<td headers="dev2" class="gt_row gt_right">2.552&nbsp;×&nbsp;10<sup style='font-size: 65%;'>10</sup></td></tr>
    <tr><td headers="i" class="gt_row gt_left">3</td>
<td headers="y" class="gt_row gt_right">5.610&nbsp;×&nbsp;10<sup style='font-size: 65%;'>4</sup></td>
<td headers="ybar" class="gt_row gt_right">2.979&nbsp;×&nbsp;10<sup style='font-size: 65%;'>5</sup></td>
<td headers="dev" class="gt_row gt_right">−2.418&nbsp;×&nbsp;10<sup style='font-size: 65%;'>5</sup></td>
<td headers="dev2" class="gt_row gt_right">5.846&nbsp;×&nbsp;10<sup style='font-size: 65%;'>10</sup></td></tr>
    <tr><td headers="i" class="gt_row gt_left">4</td>
<td headers="y" class="gt_row gt_right">1.991&nbsp;×&nbsp;10<sup style='font-size: 65%;'>5</sup></td>
<td headers="ybar" class="gt_row gt_right">2.979&nbsp;×&nbsp;10<sup style='font-size: 65%;'>5</sup></td>
<td headers="dev" class="gt_row gt_right">−9.878&nbsp;×&nbsp;10<sup style='font-size: 65%;'>4</sup></td>
<td headers="dev2" class="gt_row gt_right">9.757&nbsp;×&nbsp;10<sup style='font-size: 65%;'>9</sup></td></tr>
    <tr><td headers="i" class="gt_row gt_left">5</td>
<td headers="y" class="gt_row gt_right">8.923&nbsp;×&nbsp;10<sup style='font-size: 65%;'>4</sup></td>
<td headers="ybar" class="gt_row gt_right">2.979&nbsp;×&nbsp;10<sup style='font-size: 65%;'>5</sup></td>
<td headers="dev" class="gt_row gt_right">−2.087&nbsp;×&nbsp;10<sup style='font-size: 65%;'>5</sup></td>
<td headers="dev2" class="gt_row gt_right">4.354&nbsp;×&nbsp;10<sup style='font-size: 65%;'>10</sup></td></tr>
    <tr><td headers="i" class="gt_row gt_left">6</td>
<td headers="y" class="gt_row gt_right">9.619&nbsp;×&nbsp;10<sup style='font-size: 65%;'>4</sup></td>
<td headers="ybar" class="gt_row gt_right">2.979&nbsp;×&nbsp;10<sup style='font-size: 65%;'>5</sup></td>
<td headers="dev" class="gt_row gt_right">−2.017&nbsp;×&nbsp;10<sup style='font-size: 65%;'>5</sup></td>
<td headers="dev2" class="gt_row gt_right">4.068&nbsp;×&nbsp;10<sup style='font-size: 65%;'>10</sup></td></tr>
    <tr><td headers="i" class="gt_row gt_left">...</td>
<td headers="y" class="gt_row gt_right">...</td>
<td headers="ybar" class="gt_row gt_right">...</td>
<td headers="dev" class="gt_row gt_right">...</td>
<td headers="dev2" class="gt_row gt_right">...</td></tr>
    <tr><td headers="i" class="gt_row gt_left" style="font-weight: bold;">Sum</td>
<td headers="y" class="gt_row gt_right" style="font-weight: bold;">8.937&nbsp;×&nbsp;10<sup style='font-size: 65%;'>7</sup></td>
<td headers="ybar" class="gt_row gt_right" style="font-weight: bold;">2.979&nbsp;×&nbsp;10<sup style='font-size: 65%;'>5</sup></td>
<td headers="dev" class="gt_row gt_right" style="font-weight: bold;">−3.900&nbsp;×&nbsp;10<sup style='font-size: 65%;'>−8</sup></td>
<td headers="dev2" class="gt_row gt_right" style="font-weight: bold;">3.550&nbsp;×&nbsp;10<sup style='font-size: 65%;'>13</sup></td></tr>
  </tbody>
  <tfoot>
    <tr class="gt_footnotes">
      <td class="gt_footnote" colspan="5"><span class="gt_footnote_marks" style="white-space:nowrap;font-style:italic;font-weight:normal;line-height:0;"><sup>1</sup></span> <span data-qmd-base64="JHNfeV4yID0gXHN1bV9pICh5X2ktXGJhciB5KV4yLyhuLTEpID0gMTE4NzE2MDA4MTgwLjQwNTkkLg=="><span class='gt_from_md'>\(s_y^2 = \sum_i (y_i-\bar y)^2/(n-1) = 118716008180.4059\).</span></span></td>
    </tr>
    <tr class="gt_sourcenotes">
      <td class="gt_sourcenote" colspan="5"><span data-qmd-base64="JFxiYXIgeSA9IFxkZnJhY3tcc3VtX2kgeV9pfXtufSA9IFxkZnJhY3s4LjkzNyBcdGltZXMgMTBeezd9fXszMDB9ID0gMi45NzkgXHRpbWVzIDEwXns1fSQ="><span class='gt_from_md'>\(\bar y = \dfrac{\sum_i y_i}{n} = \dfrac{8.937 \times 10^{7}}{300} = 2.979 \times 10^{5}\)</span></span></td>
    </tr>
    <tr class="gt_sourcenotes">
      <td class="gt_sourcenote" colspan="5"><span data-qmd-base64="JFxtYXRocm17U0V9KFxiYXIgeSkgPSBcc3FydHsxLVxkZnJhY3tufXtOfX1cLFxkZnJhY3tzfXtcc3FydCBufSA9IFxzcXJ0ezEtXGRmcmFjezMwMH17MzA3OH19XCxcZGZyYWN7My40NDYgXHRpbWVzIDEwXns1fX17XHNxcnR7MzAwfX0gPSAxOCw4OTguNDM0JA=="><span class='gt_from_md'>\(\mathrm{SE}(\bar y) = \sqrt{1-\dfrac{n}{N}}\,\dfrac{s}{\sqrt n} = \sqrt{1-\dfrac{300}{3078}}\,\dfrac{3.446 \times 10^{5}}{\sqrt{300}} = 18,898.434\)</span></span></td>
    </tr>
  </tfoot>
</table>
</div>
```

:::

```{.r .cell-code}
res_mean$estimate |> format_est_gt("mean")
```

::: {.cell-output-display}

```{=html}
<div id="mdmlucupjt" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#mdmlucupjt table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#mdmlucupjt thead, #mdmlucupjt tbody, #mdmlucupjt tfoot, #mdmlucupjt tr, #mdmlucupjt td, #mdmlucupjt th {
  border-style: none;
}

#mdmlucupjt p {
  margin: 0;
  padding: 0;
}

#mdmlucupjt .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: 70%;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#mdmlucupjt .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#mdmlucupjt .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#mdmlucupjt .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#mdmlucupjt .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#mdmlucupjt .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#mdmlucupjt .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#mdmlucupjt .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#mdmlucupjt .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#mdmlucupjt .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#mdmlucupjt .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#mdmlucupjt .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#mdmlucupjt .gt_spanner_row {
  border-bottom-style: hidden;
}

#mdmlucupjt .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#mdmlucupjt .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#mdmlucupjt .gt_from_md > :first-child {
  margin-top: 0;
}

#mdmlucupjt .gt_from_md > :last-child {
  margin-bottom: 0;
}

#mdmlucupjt .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#mdmlucupjt .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#mdmlucupjt .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#mdmlucupjt .gt_row_group_first td {
  border-top-width: 2px;
}

#mdmlucupjt .gt_row_group_first th {
  border-top-width: 2px;
}

#mdmlucupjt .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#mdmlucupjt .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#mdmlucupjt .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#mdmlucupjt .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#mdmlucupjt .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#mdmlucupjt .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#mdmlucupjt .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#mdmlucupjt .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#mdmlucupjt .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#mdmlucupjt .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#mdmlucupjt .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#mdmlucupjt .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#mdmlucupjt .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#mdmlucupjt .gt_left {
  text-align: left;
}

#mdmlucupjt .gt_center {
  text-align: center;
}

#mdmlucupjt .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#mdmlucupjt .gt_font_normal {
  font-weight: normal;
}

#mdmlucupjt .gt_font_bold {
  font-weight: bold;
}

#mdmlucupjt .gt_font_italic {
  font-style: italic;
}

#mdmlucupjt .gt_super {
  font-size: 65%;
}

#mdmlucupjt .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#mdmlucupjt .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#mdmlucupjt .gt_indent_1 {
  text-indent: 5px;
}

#mdmlucupjt .gt_indent_2 {
  text-indent: 10px;
}

#mdmlucupjt .gt_indent_3 {
  text-indent: 15px;
}

#mdmlucupjt .gt_indent_4 {
  text-indent: 20px;
}

#mdmlucupjt .gt_indent_5 {
  text-indent: 25px;
}

#mdmlucupjt .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#mdmlucupjt div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Est."><span data-qmd-base64="JFxiYXJ7eX0k"><span class='gt_from_md'>\(\bar{y}\)</span></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="S.E."><span data-qmd-base64="JFxtYXRocm17U0V9KFxiYXJ7eX0pJA=="><span class='gt_from_md'>\(\mathrm{SE}(\bar{y})\)</span></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="ci.low"><span data-qmd-base64="OTUlIENJIExvd2Vy"><span class='gt_from_md'>95% CI Lower</span></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="ci.upp"><span data-qmd-base64="OTUlIENJIFVwcGVy"><span class='gt_from_md'>95% CI Upper</span></span></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="Est." class="gt_row gt_right">297,897.0467</td>
<td headers="S.E." class="gt_row gt_right">18,898.4344</td>
<td headers="ci.low" class="gt_row gt_right">260,706.2569</td>
<td headers="ci.upp" class="gt_row gt_right">335,087.8365</td></tr>
  </tbody>
  
</table>
</div>
```

:::
:::


*Estimating the total of acre92*

Set `estimate = "total"` to have the function itself scale the mean and confidence interval by the population size ($N = 3078$) to estimate the total farm acreage across all counties.


::: {.cell}

```{.r .cell-code}
res_total <- srs_est(agsrs[,"acres92"], N = 3078, estimate = "total", show.details = FALSE)
res_total$estimate |> format_est_gt("total")
```

::: {.cell-output-display}

```{=html}
<div id="rovffllter" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#rovffllter table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#rovffllter thead, #rovffllter tbody, #rovffllter tfoot, #rovffllter tr, #rovffllter td, #rovffllter th {
  border-style: none;
}

#rovffllter p {
  margin: 0;
  padding: 0;
}

#rovffllter .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: 70%;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#rovffllter .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#rovffllter .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#rovffllter .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#rovffllter .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#rovffllter .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#rovffllter .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#rovffllter .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#rovffllter .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#rovffllter .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#rovffllter .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#rovffllter .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#rovffllter .gt_spanner_row {
  border-bottom-style: hidden;
}

#rovffllter .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#rovffllter .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#rovffllter .gt_from_md > :first-child {
  margin-top: 0;
}

#rovffllter .gt_from_md > :last-child {
  margin-bottom: 0;
}

#rovffllter .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#rovffllter .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#rovffllter .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#rovffllter .gt_row_group_first td {
  border-top-width: 2px;
}

#rovffllter .gt_row_group_first th {
  border-top-width: 2px;
}

#rovffllter .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#rovffllter .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#rovffllter .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#rovffllter .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#rovffllter .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#rovffllter .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#rovffllter .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#rovffllter .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#rovffllter .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#rovffllter .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#rovffllter .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#rovffllter .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#rovffllter .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#rovffllter .gt_left {
  text-align: left;
}

#rovffllter .gt_center {
  text-align: center;
}

#rovffllter .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#rovffllter .gt_font_normal {
  font-weight: normal;
}

#rovffllter .gt_font_bold {
  font-weight: bold;
}

#rovffllter .gt_font_italic {
  font-style: italic;
}

#rovffllter .gt_super {
  font-size: 65%;
}

#rovffllter .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#rovffllter .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#rovffllter .gt_indent_1 {
  text-indent: 5px;
}

#rovffllter .gt_indent_2 {
  text-indent: 10px;
}

#rovffllter .gt_indent_3 {
  text-indent: 15px;
}

#rovffllter .gt_indent_4 {
  text-indent: 20px;
}

#rovffllter .gt_indent_5 {
  text-indent: 25px;
}

#rovffllter .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#rovffllter div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Est."><span data-qmd-base64="JFxoYXR7dH0k"><span class='gt_from_md'>\(\hat{t}\)</span></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="S.E."><span data-qmd-base64="JFxtYXRocm17U0V9KFxoYXR7dH0pJA=="><span class='gt_from_md'>\(\mathrm{SE}(\hat{t})\)</span></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="ci.low"><span data-qmd-base64="OTUlIENJIExvd2Vy"><span class='gt_from_md'>95% CI Lower</span></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="ci.upp"><span data-qmd-base64="OTUlIENJIFVwcGVy"><span class='gt_from_md'>95% CI Upper</span></span></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="Est." class="gt_row gt_right">916,927,109.6400</td>
<td headers="S.E." class="gt_row gt_right">58,169,381.1695</td>
<td headers="ci.low" class="gt_row gt_right">802,453,858.6054</td>
<td headers="ci.upp" class="gt_row gt_right">1,031,400,360.6746</td></tr>
  </tbody>
  
</table>
</div>
```

:::
:::


*Estimating the proportion of counties with fewer than 200K acres for farming in 1992*

To estimate a proportion under SRS, create an indicator variable ($y_i = 1$ if true, $0$ if false). The mean of this indicator variable serves as the sample proportion ($\hat{p} = \bar{y}_{indicator}$).


::: {.cell}

```{.r .cell-code}
acres92.is.fewer.200k <- as.numeric (agsrs[,"acres92"] < 200000)
head(acres92.is.fewer.200k)
```

::: {.cell-output .cell-output-stdout}

```
[1] 1 1 1 1 1 1
```


:::

```{.r .cell-code}
res_prop <- srs_est(acres92.is.fewer.200k, N = 3078, estimate = "mean", show.details = FALSE)
res_prop$estimate |> format_est_gt("mean")
```

::: {.cell-output-display}

```{=html}
<div id="ktwkutgoys" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#ktwkutgoys table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#ktwkutgoys thead, #ktwkutgoys tbody, #ktwkutgoys tfoot, #ktwkutgoys tr, #ktwkutgoys td, #ktwkutgoys th {
  border-style: none;
}

#ktwkutgoys p {
  margin: 0;
  padding: 0;
}

#ktwkutgoys .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: 70%;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#ktwkutgoys .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#ktwkutgoys .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#ktwkutgoys .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#ktwkutgoys .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#ktwkutgoys .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#ktwkutgoys .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#ktwkutgoys .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#ktwkutgoys .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#ktwkutgoys .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#ktwkutgoys .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#ktwkutgoys .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#ktwkutgoys .gt_spanner_row {
  border-bottom-style: hidden;
}

#ktwkutgoys .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#ktwkutgoys .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#ktwkutgoys .gt_from_md > :first-child {
  margin-top: 0;
}

#ktwkutgoys .gt_from_md > :last-child {
  margin-bottom: 0;
}

#ktwkutgoys .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#ktwkutgoys .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#ktwkutgoys .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#ktwkutgoys .gt_row_group_first td {
  border-top-width: 2px;
}

#ktwkutgoys .gt_row_group_first th {
  border-top-width: 2px;
}

#ktwkutgoys .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#ktwkutgoys .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#ktwkutgoys .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#ktwkutgoys .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#ktwkutgoys .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#ktwkutgoys .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#ktwkutgoys .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#ktwkutgoys .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#ktwkutgoys .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#ktwkutgoys .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#ktwkutgoys .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#ktwkutgoys .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#ktwkutgoys .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#ktwkutgoys .gt_left {
  text-align: left;
}

#ktwkutgoys .gt_center {
  text-align: center;
}

#ktwkutgoys .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#ktwkutgoys .gt_font_normal {
  font-weight: normal;
}

#ktwkutgoys .gt_font_bold {
  font-weight: bold;
}

#ktwkutgoys .gt_font_italic {
  font-style: italic;
}

#ktwkutgoys .gt_super {
  font-size: 65%;
}

#ktwkutgoys .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#ktwkutgoys .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#ktwkutgoys .gt_indent_1 {
  text-indent: 5px;
}

#ktwkutgoys .gt_indent_2 {
  text-indent: 10px;
}

#ktwkutgoys .gt_indent_3 {
  text-indent: 15px;
}

#ktwkutgoys .gt_indent_4 {
  text-indent: 20px;
}

#ktwkutgoys .gt_indent_5 {
  text-indent: 25px;
}

#ktwkutgoys .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#ktwkutgoys div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Est."><span data-qmd-base64="JFxiYXJ7eX0k"><span class='gt_from_md'>\(\bar{y}\)</span></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="S.E."><span data-qmd-base64="JFxtYXRocm17U0V9KFxiYXJ7eX0pJA=="><span class='gt_from_md'>\(\mathrm{SE}(\bar{y})\)</span></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="ci.low"><span data-qmd-base64="OTUlIENJIExvd2Vy"><span class='gt_from_md'>95% CI Lower</span></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="ci.upp"><span data-qmd-base64="OTUlIENJIFVwcGVy"><span class='gt_from_md'>95% CI Upper</span></span></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="Est." class="gt_row gt_right">0.5100</td>
<td headers="S.E." class="gt_row gt_right">0.0275</td>
<td headers="ci.low" class="gt_row gt_right">0.4560</td>
<td headers="ci.upp" class="gt_row gt_right">0.5640</td></tr>
  </tbody>
  
</table>
</div>
```

:::
:::


*Estimating the total number of counties with fewer than 200K acres for farming in 1992*

Set `estimate = "total"` to estimate the absolute count of counties meeting the criteria.


::: {.cell}

```{.r .cell-code}
res_prop_total <- srs_est(acres92.is.fewer.200k, N = 3078, estimate = "total", show.details = FALSE)
res_prop_total$estimate |> format_est_gt("total")
```

::: {.cell-output-display}

```{=html}
<div id="nskealpzbv" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#nskealpzbv table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#nskealpzbv thead, #nskealpzbv tbody, #nskealpzbv tfoot, #nskealpzbv tr, #nskealpzbv td, #nskealpzbv th {
  border-style: none;
}

#nskealpzbv p {
  margin: 0;
  padding: 0;
}

#nskealpzbv .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: 70%;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#nskealpzbv .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#nskealpzbv .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#nskealpzbv .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#nskealpzbv .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#nskealpzbv .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#nskealpzbv .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#nskealpzbv .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#nskealpzbv .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#nskealpzbv .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#nskealpzbv .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#nskealpzbv .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#nskealpzbv .gt_spanner_row {
  border-bottom-style: hidden;
}

#nskealpzbv .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#nskealpzbv .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#nskealpzbv .gt_from_md > :first-child {
  margin-top: 0;
}

#nskealpzbv .gt_from_md > :last-child {
  margin-bottom: 0;
}

#nskealpzbv .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#nskealpzbv .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#nskealpzbv .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#nskealpzbv .gt_row_group_first td {
  border-top-width: 2px;
}

#nskealpzbv .gt_row_group_first th {
  border-top-width: 2px;
}

#nskealpzbv .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#nskealpzbv .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#nskealpzbv .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#nskealpzbv .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#nskealpzbv .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#nskealpzbv .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#nskealpzbv .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#nskealpzbv .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#nskealpzbv .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#nskealpzbv .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#nskealpzbv .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#nskealpzbv .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#nskealpzbv .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#nskealpzbv .gt_left {
  text-align: left;
}

#nskealpzbv .gt_center {
  text-align: center;
}

#nskealpzbv .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#nskealpzbv .gt_font_normal {
  font-weight: normal;
}

#nskealpzbv .gt_font_bold {
  font-weight: bold;
}

#nskealpzbv .gt_font_italic {
  font-style: italic;
}

#nskealpzbv .gt_super {
  font-size: 65%;
}

#nskealpzbv .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#nskealpzbv .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#nskealpzbv .gt_indent_1 {
  text-indent: 5px;
}

#nskealpzbv .gt_indent_2 {
  text-indent: 10px;
}

#nskealpzbv .gt_indent_3 {
  text-indent: 15px;
}

#nskealpzbv .gt_indent_4 {
  text-indent: 20px;
}

#nskealpzbv .gt_indent_5 {
  text-indent: 25px;
}

#nskealpzbv .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#nskealpzbv div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Est."><span data-qmd-base64="JFxoYXR7dH0k"><span class='gt_from_md'>\(\hat{t}\)</span></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="S.E."><span data-qmd-base64="JFxtYXRocm17U0V9KFxoYXR7dH0pJA=="><span class='gt_from_md'>\(\mathrm{SE}(\hat{t})\)</span></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="ci.low"><span data-qmd-base64="OTUlIENJIExvd2Vy"><span class='gt_from_md'>95% CI Lower</span></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="ci.upp"><span data-qmd-base64="OTUlIENJIFVwcGVy"><span class='gt_from_md'>95% CI Upper</span></span></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="Est." class="gt_row gt_right">1,569.7800</td>
<td headers="S.E." class="gt_row gt_right">84.5372</td>
<td headers="ci.low" class="gt_row gt_right">1,403.4167</td>
<td headers="ci.upp" class="gt_row gt_right">1,736.1433</td></tr>
  </tbody>
  
</table>
</div>
```

:::
:::


### Comparing with true value

Load the complete census (population) dataset to calculate the true parameters. This allows us to assess the accuracy of our previous SRS estimates against the exact population values.


::: {.cell}

```{.r .cell-code}
agpop <- read.csv ("data/agpop.csv", na = "-99")
#true mean
mean (agpop[, "acres92"], na.rm = T)
```

::: {.cell-output .cell-output-stdout}

```
[1] 308582.4
```


:::

```{.r .cell-code}
# true total
sum (agpop[, "acres92"], na.rm = T)
```

::: {.cell-output .cell-output-stdout}

```
[1] 943953599
```


:::

```{.r .cell-code}
# true proportion of counties with less than 200K acres for farming
mean (agpop[, "acres92"] < 200000, na.rm = T)
```

::: {.cell-output .cell-output-stdout}

```
[1] 0.5145472
```


:::

```{.r .cell-code}
# true number of counties with less than 200K acres for farming
sum (agpop[, "acres92"] < 200000, na.rm = T)
```

::: {.cell-output .cell-output-stdout}

```
[1] 1574
```


:::
:::


## A Simulation Demonstration of SRS Inference

Load the full population data and filter out missing values to serve as a complete, known finite population for resampling experiments.


::: {.cell}

```{.r .cell-code}
# read population data 
agpop <- read.csv ("data/agpop.csv")
# remove those counties with na
agpop <- subset( agpop, acres92 != -99)
```
:::


### True Values

Define the target sample size and compute the true population mean ($\mu$) and theoretical standard error ($\mathrm{SE}_{true}$) using the known population standard deviation ($S$).

**Theoretical Standard Error:**


$$\mathrm{SE}_{true}(\bar{y}) = \sqrt{1 - \frac{n}{N}} \left(\frac{S}{\sqrt{n}}\right)$$


::: {.cell}

```{.r .cell-code}
# sample size
n <- 300
# population size
N <- nrow (agpop); N
```

::: {.cell-output .cell-output-stdout}

```
[1] 3059
```


:::

```{.r .cell-code}
# true value of population mean
ybarU <- mean (agpop[,"acres92"]); ybarU
```

::: {.cell-output .cell-output-stdout}

```
[1] 308582.4
```


:::

```{.r .cell-code}
# true value of deviation of sample mean
true.se.ybar <- sqrt (1- n/N) * sd (agpop[,"acres92"]) / sqrt (n); true.se.ybar
```

::: {.cell-output .cell-output-stdout}

```
[1] 23320.29
```


:::
:::


### One SRS sampling

Draw a single random sample of size $n=300$ from the population and run our estimator function to see a single instance of survey inference.


::: {.cell}

```{.r .cell-code}
##
# srs sampling
srs <- sample (1:N,n)
head(agpop [srs, ])
```

::: {.cell-output-display}
`````{=html}
<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["county"],"name":[1],"type":["chr"],"align":["left"]},{"label":["state"],"name":[2],"type":["chr"],"align":["left"]},{"label":["acres92"],"name":[3],"type":["int"],"align":["right"]},{"label":["acres87"],"name":[4],"type":["int"],"align":["right"]},{"label":["acres82"],"name":[5],"type":["int"],"align":["right"]},{"label":["farms92"],"name":[6],"type":["int"],"align":["right"]},{"label":["farms87"],"name":[7],"type":["int"],"align":["right"]},{"label":["farms82"],"name":[8],"type":["int"],"align":["right"]},{"label":["largef92"],"name":[9],"type":["int"],"align":["right"]},{"label":["largef87"],"name":[10],"type":["int"],"align":["right"]},{"label":["largef82"],"name":[11],"type":["int"],"align":["right"]},{"label":["smallf92"],"name":[12],"type":["int"],"align":["right"]},{"label":["smallf87"],"name":[13],"type":["int"],"align":["right"]},{"label":["smallf82"],"name":[14],"type":["int"],"align":["right"]},{"label":["region"],"name":[15],"type":["chr"],"align":["left"]}],"data":[{"1":"DAVIS COUNTY","2":"UT","3":"50357","4":"63244","5":"111721","6":"582","7":"647","8":"660","9":"4","10":"7","11":"6","12":"192","13":"205","14":"219","15":"W","_rn_":"2754"},{"1":"BUTTS COUNTY","2":"GA","3":"29213","4":"29157","5":"30014","6":"139","7":"146","8":"171","9":"4","10":"5","11":"2","12":"4","13":"7","14":"6","15":"S","_rn_":"379"},{"1":"MONTGOMERY COUNTY","2":"IA","3":"240100","4":"247891","5":"253540","6":"617","7":"701","8":"769","9":"40","10":"29","11":"21","12":"26","13":"39","14":"51","15":"NC","_rn_":"593"},{"1":"ELLIS COUNTY","2":"TX","3":"426189","4":"392585","5":"445380","6":"1521","7":"1612","8":"1608","9":"96","10":"94","11":"99","12":"115","13":"121","14":"110","15":"S","_rn_":"2564"},{"1":"LEE COUNTY","2":"AL","3":"67962","4":"79836","5":"100949","6":"336","7":"402","8":"407","9":"10","10":"10","11":"20","12":"15","13":"22","14":"20","15":"S","_rn_":"46"},{"1":"MOTLEY COUNTY","2":"TX","3":"479889","4":"467309","5":"512458","6":"190","7":"202","8":"233","9":"86","10":"83","11":"93","12":"4","13":"6","14":"7","15":"S","_rn_":"2667"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>
`````
:::

```{.r .cell-code}
# get data of variable "acres92"
sdata <- agpop [srs, "acres92"]
# analysis
res_one <- srs_est(sdata, N, estimate = "mean", show.details = FALSE)
res_one$estimate |> format_est_gt("mean")
```

::: {.cell-output-display}

```{=html}
<div id="bvcqwrbvhf" style="padding-left:0px;padding-right:0px;padding-top:10px;padding-bottom:10px;overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>#bvcqwrbvhf table {
  font-family: system-ui, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji', 'Segoe UI Symbol', 'Noto Color Emoji';
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

#bvcqwrbvhf thead, #bvcqwrbvhf tbody, #bvcqwrbvhf tfoot, #bvcqwrbvhf tr, #bvcqwrbvhf td, #bvcqwrbvhf th {
  border-style: none;
}

#bvcqwrbvhf p {
  margin: 0;
  padding: 0;
}

#bvcqwrbvhf .gt_table {
  display: table;
  border-collapse: collapse;
  line-height: normal;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: 70%;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#bvcqwrbvhf .gt_caption {
  padding-top: 4px;
  padding-bottom: 4px;
}

#bvcqwrbvhf .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#bvcqwrbvhf .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 3px;
  padding-bottom: 5px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#bvcqwrbvhf .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#bvcqwrbvhf .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#bvcqwrbvhf .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#bvcqwrbvhf .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#bvcqwrbvhf .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#bvcqwrbvhf .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#bvcqwrbvhf .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#bvcqwrbvhf .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 5px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#bvcqwrbvhf .gt_spanner_row {
  border-bottom-style: hidden;
}

#bvcqwrbvhf .gt_group_heading {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  text-align: left;
}

#bvcqwrbvhf .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#bvcqwrbvhf .gt_from_md > :first-child {
  margin-top: 0;
}

#bvcqwrbvhf .gt_from_md > :last-child {
  margin-bottom: 0;
}

#bvcqwrbvhf .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#bvcqwrbvhf .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
}

#bvcqwrbvhf .gt_stub_row_group {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 5px;
  padding-right: 5px;
  vertical-align: top;
}

#bvcqwrbvhf .gt_row_group_first td {
  border-top-width: 2px;
}

#bvcqwrbvhf .gt_row_group_first th {
  border-top-width: 2px;
}

#bvcqwrbvhf .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#bvcqwrbvhf .gt_first_summary_row {
  border-top-style: solid;
  border-top-color: #D3D3D3;
}

#bvcqwrbvhf .gt_first_summary_row.thick {
  border-top-width: 2px;
}

#bvcqwrbvhf .gt_last_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#bvcqwrbvhf .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#bvcqwrbvhf .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#bvcqwrbvhf .gt_last_grand_summary_row_top {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-bottom-style: double;
  border-bottom-width: 6px;
  border-bottom-color: #D3D3D3;
}

#bvcqwrbvhf .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#bvcqwrbvhf .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#bvcqwrbvhf .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#bvcqwrbvhf .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#bvcqwrbvhf .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#bvcqwrbvhf .gt_sourcenote {
  font-size: 90%;
  padding-top: 4px;
  padding-bottom: 4px;
  padding-left: 5px;
  padding-right: 5px;
}

#bvcqwrbvhf .gt_left {
  text-align: left;
}

#bvcqwrbvhf .gt_center {
  text-align: center;
}

#bvcqwrbvhf .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#bvcqwrbvhf .gt_font_normal {
  font-weight: normal;
}

#bvcqwrbvhf .gt_font_bold {
  font-weight: bold;
}

#bvcqwrbvhf .gt_font_italic {
  font-style: italic;
}

#bvcqwrbvhf .gt_super {
  font-size: 65%;
}

#bvcqwrbvhf .gt_footnote_marks {
  font-size: 75%;
  vertical-align: 0.4em;
  position: initial;
}

#bvcqwrbvhf .gt_asterisk {
  font-size: 100%;
  vertical-align: 0;
}

#bvcqwrbvhf .gt_indent_1 {
  text-indent: 5px;
}

#bvcqwrbvhf .gt_indent_2 {
  text-indent: 10px;
}

#bvcqwrbvhf .gt_indent_3 {
  text-indent: 15px;
}

#bvcqwrbvhf .gt_indent_4 {
  text-indent: 20px;
}

#bvcqwrbvhf .gt_indent_5 {
  text-indent: 25px;
}

#bvcqwrbvhf .katex-display {
  display: inline-flex !important;
  margin-bottom: 0.75em !important;
}

#bvcqwrbvhf div.Reactable > div.rt-table > div.rt-thead > div.rt-tr.rt-tr-group-header > div.rt-th-group:after {
  height: 0px !important;
}
</style>
<table class="gt_table" data-quarto-disable-processing="false" data-quarto-bootstrap="false">
  <thead>
    <tr class="gt_col_headings">
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="Est."><span data-qmd-base64="JFxiYXJ7eX0k"><span class='gt_from_md'>\(\bar{y}\)</span></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="S.E."><span data-qmd-base64="JFxtYXRocm17U0V9KFxiYXJ7eX0pJA=="><span class='gt_from_md'>\(\mathrm{SE}(\bar{y})\)</span></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="ci.low"><span data-qmd-base64="OTUlIENJIExvd2Vy"><span class='gt_from_md'>95% CI Lower</span></span></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_right" rowspan="1" colspan="1" scope="col" id="ci.upp"><span data-qmd-base64="OTUlIENJIFVwcGVy"><span class='gt_from_md'>95% CI Upper</span></span></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td headers="Est." class="gt_row gt_right">300,559.5133</td>
<td headers="S.E." class="gt_row gt_right">19,458.1741</td>
<td headers="ci.low" class="gt_row gt_right">262,267.1952</td>
<td headers="ci.upp" class="gt_row gt_right">338,851.8315</td></tr>
  </tbody>
  
</table>
</div>
```

:::
:::


### Repeating SRS sampling 5000 times

To demonstrate the Central Limit Theorem and unbiasedness of the estimator, draw 5,000 independent samples, storing the estimates from each iteration. We plot the resulting sampling distribution against the true population mean.

*(Note: We use `show.details = FALSE` and pull out just `$estimate` here, so the loop skips building 5,000 working tables and stays fast.)*


::: {.cell}

```{.r .cell-code}
nres <- 5000 # number of repeated sampling
simulation.results <- matrix (0, nres, 4) # matrix recording repeated results
colnames(simulation.results) <- c( "Est.",   "S.E.",   "ci.low", "ci.upp")

for (i in 1:nres)
{
    srs <- sample (N, n)
   sdata <- agpop [srs, "acres92"]
   simulation.results [i,] <- srs_est (sdata, N, estimate = "mean", show.details = FALSE)$estimate
}

head(simulation.results)
```

::: {.cell-output .cell-output-stdout}

```
         Est.     S.E.   ci.low   ci.upp
[1,] 308193.3 25761.98 257495.5 358891.1
[2,] 335715.1 24079.57 288328.2 383102.0
[3,] 270887.9 19534.37 232445.7 309330.2
[4,] 353224.0 26049.61 301960.2 404487.8
[5,] 283973.8 17035.49 250449.2 317498.5
[6,] 288223.7 16911.20 254943.6 321503.8
```


:::

```{.r .cell-code}
# look at the distribution of sample mean
par (mfrow= c(2,2))
hist (agpop$acres92,main = "Population Distribution of acre92")
hist (simulation.results[,1], main = "Sampling Distribution of Sample Mean for acre92")
abline (v = ybarU, col = "red")
qqnorm (simulation.results[,1], main="QQ plot of Sample Mean"); qqline(simulation.results[,1])
boxplot (simulation.results[,1], main = "Boxplot of Sample Mean")
abline (h = ybarU, col = "red")
```

::: {.cell-output-display}
![](srs_files/figure-html/srs-simulation-1.png){width=672}
:::

```{.r .cell-code}
mean (simulation.results[,1])
```

::: {.cell-output .cell-output-stdout}

```
[1] 308320.3
```


:::

```{.r .cell-code}
ybarU
```

::: {.cell-output .cell-output-stdout}

```
[1] 308582.4
```


:::

```{.r .cell-code}
sd (simulation.results [,1])
```

::: {.cell-output .cell-output-stdout}

```
[1] 23234.34
```


:::

```{.r .cell-code}
true.se.ybar
```

::: {.cell-output .cell-output-stdout}

```
[1] 23320.29
```


:::
:::


### Empirical Coverage Rate of CIs

Verify the reliability of the 95% confidence intervals by checking whether the true population mean falls within the calculated bounds for each of the 5,000 samples. The empirical coverage rate should closely approximate 0.95. The plot highlights the first 100 intervals, flagging those that fail to capture the true mean in a distinct color.


::: {.cell}

```{.r .cell-code}
library(dplyr)
library(plotrix) 

# Convert to data frame (if it was a matrix) and use mutate to evaluate coverage
simulation.results <- as.data.frame(simulation.results) |>
  mutate(
    `Covered?` = as.numeric(ci.low < ybarU & ybarU < ci.upp)
  )

head(simulation.results)
```

::: {.cell-output-display}
`````{=html}
<div data-pagedtable="false">
  <script data-pagedtable-source type="application/json">
{"columns":[{"label":[""],"name":["_rn_"],"type":[""],"align":["left"]},{"label":["Est."],"name":[1],"type":["dbl"],"align":["right"]},{"label":["S.E."],"name":[2],"type":["dbl"],"align":["right"]},{"label":["ci.low"],"name":[3],"type":["dbl"],"align":["right"]},{"label":["ci.upp"],"name":[4],"type":["dbl"],"align":["right"]},{"label":["Covered?"],"name":[5],"type":["dbl"],"align":["right"]}],"data":[{"1":"308193.3","2":"25761.98","3":"257495.5","4":"358891.1","5":"1","_rn_":"1"},{"1":"335715.1","2":"24079.57","3":"288328.2","4":"383102.0","5":"1","_rn_":"2"},{"1":"270887.9","2":"19534.37","3":"232445.7","4":"309330.2","5":"1","_rn_":"3"},{"1":"353224.0","2":"26049.61","3":"301960.2","4":"404487.8","5":"1","_rn_":"4"},{"1":"283973.8","2":"17035.49","3":"250449.2","4":"317498.5","5":"1","_rn_":"5"},{"1":"288223.7","2":"16911.20","3":"254943.6","4":"321503.8","5":"1","_rn_":"6"}],"options":{"columns":{"min":{},"max":[10]},"rows":{"min":[10],"max":[10]},"pages":{}}}
  </script>
</div>
`````
:::

```{.r .cell-code}
par(mfrow=c(1,1))
plotCI(
  x = 1:100,
  y = simulation.results$Est.[1:100], 
  li = simulation.results$ci.low[1:100],
  ui = simulation.results$ci.upp[1:100],
  col = 2 - simulation.results$`Covered?`[1:100]
)
abline(h = ybarU, col = "blue")
```

::: {.cell-output-display}
![](srs_files/figure-html/srs-coverage-1.png){width=672}
:::

```{.r .cell-code}
# Empirical coverage rate
mean(simulation.results$`Covered?`)
```

::: {.cell-output .cell-output-stdout}

```
[1] 0.9366
```


:::
:::


### Interactive Demonstration

The static plots above show a snapshot of 5,000 repeated samples. The app below
draws SRS samples one at a time (or automatically).

The top panel shows the population of `acres92`: the grey curve is a kernel
density estimate and the grey points below the axis are the raw values, jittered
vertically. The axis is broken at twice the population mean; counties above the
cut are compressed, in rank order, into the narrow lane at the right, so that a
tail unit drawn into the sample is still visible. Points selected into the
current sample turn red.

The bottom panel accumulates the sample means. Its window is fixed at five
theoretical standard errors for $n = 100$ on either side of $\mu$, so the spread
is comparable across sample sizes; for $n < 100$ some means fall outside the
window and are counted in the legend. The interval above the histogram belongs to
the current sample, drawn green when it covers $\mu$ and red when it does not,
and replaced at every draw. Both panels are centred on $\mu$, so the two dashed
lines are vertically aligned.



### Interactive Demonstration

The static plots above show a snapshot of 5,000 repeated samples. The app below
draws SRS samples one at a time (or automatically).

The top panel shows the population of `acres92`: the grey curve is a kernel
density estimate and the grey points below the axis are the raw values, jittered
vertically. The axis is broken at twice the population mean; counties above the
cut are compressed, in rank order, into the narrow lane at the right, so that a
tail unit drawn into the sample is still visible. Points selected into the
current sample turn red.

The bottom panel accumulates the sample means. Its window is fixed at five
theoretical standard errors for $n = 100$ on either side of $\mu$, so the spread
is comparable across sample sizes; for $n < 100$ some means fall outside the
window and are counted in the legend. The interval above the histogram belongs to
the current sample, drawn green when it covers $\mu$ and red when it does not,
and replaced at every draw. Both panels are centred on $\mu$, so the two dashed
lines are vertically aligned.

```{shinylive-r}
#| standalone: true
#| viewerHeight: 830

## file: app.R
library(shiny)

pop   <- read.csv("acres92.csv")$acres92
N     <- length(pop)
ybarU <- mean(pop)
S     <- sd(pop)

## ---- broken x-axis for the population panel --------------------------------
## The drawn region is [0, 2 * ybarU]: the visible scale [0, cut] plus a
## compressed lane holding the right tail. That total width puts the population
## mean exactly at the midpoint, matching the histogram panel below.
lane_frac <- 0.10
cut       <- 2 * ybarU / (1 + lane_frac)
lane_w    <- lane_frac * cut
xmax_top  <- cut + lane_w                      # == 2 * ybarU

tail_id <- which(pop > cut)
n_tail  <- length(tail_id)
px      <- pop                                 # plotting positions
if (n_tail > 0) {
  r <- rank(pop[tail_id], ties.method = "first") / (n_tail + 1)
  px[tail_id] <- cut + lane_w * (0.10 + 0.80 * r)
}

# fixed jitter band for the rug of raw population values
set.seed(42)
rug_y <- runif(N, -0.42, -0.08)

dens   <- density(pop, from = 0, to = cut)
dens_y <- dens$y / max(dens$y)

## ---- fixed window for the sampling distribution ----------------------------
se100  <- sqrt(1 - 100 / N) * S / sqrt(100)
half   <- 5 * se100
xlim_b <- ybarU + c(-1, 1) * half

MAR <- c(4.5, 4.5, 3, 1.5)   # identical in both panels, so the scales align

ui <- fluidPage(
  titlePanel("Simple random sampling: sampling distribution and CI coverage"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      sliderInput("n", "Sample size (n)", min = 10, max = 1000, value = 300, step = 10),
      sliderInput("conf", "Confidence level", min = 0.80, max = 0.99, value = 0.95, step = 0.01),
      actionButton("draw", "Draw new SRS sample", class = "btn-primary"),
      br(), br(),
      checkboxInput("auto", "Auto-run", value = FALSE),
      sliderInput("speed", "Samples per second", min = 1, max = 10, value = 3, step = 1),
      actionButton("reset", "Reset"),
      hr(),
      verbatimTextOutput("stats")
    ),
    mainPanel(
      width = 9,
      plotOutput("popPlot",  height = "340px"),
      plotOutput("histPlot", height = "360px")
    )
  )
)

server <- function(input, output, session) {

  rv <- reactiveValues(
    idx = NULL, mean = NULL, ci = NULL, covered = NULL,
    means = numeric(0), hits = logical(0)
  )

  draw_one <- function() {
    n    <- isolate(input$n)
    conf <- isolate(input$conf)
    idx  <- sample(N, n)
    s    <- pop[idx]
    m    <- mean(s)
    se   <- sqrt(1 - n / N) * sd(s) / sqrt(n)
    me   <- qt(1 - (1 - conf) / 2, df = n - 1) * se

    rv$idx     <- idx
    rv$mean    <- m
    rv$ci      <- c(m - me, m + me)
    rv$covered <- (m - me) <= ybarU && ybarU <= (m + me)
    rv$means   <- c(rv$means, m)
    rv$hits    <- c(rv$hits, rv$covered)
  }

  observeEvent(input$draw, draw_one())

  observe({
    if (isTRUE(input$auto)) {
      invalidateLater(1000 / input$speed, session)
      isolate(draw_one())
    }
  })

  reset_history <- function() {
    rv$idx <- NULL; rv$mean <- NULL; rv$ci <- NULL; rv$covered <- NULL
    rv$means <- numeric(0); rv$hits <- logical(0)
  }
  observeEvent(input$reset, reset_history())
  observeEvent(input$n,     reset_history())
  observeEvent(input$conf,  reset_history())

  fmt <- function(x) format(round(x), big.mark = ",")

  output$popPlot <- renderPlot({
    par(mar = MAR, xaxs = "i")
    plot(dens$x, dens_y, type = "l", col = "grey40", lwd = 2,
         xlim = c(0, xmax_top), ylim = c(-0.62, 1.15),
         xlab = "acres92", ylab = "", xaxt = "n", yaxt = "n", bty = "n",
         main = "Population distribution with the current SRS sample")

    # raw values, jittered; tail units sit in the compressed lane
    points(px, rug_y, pch = 16, cex = 0.5, col = adjustcolor("grey60", 0.5))
    if (!is.null(rv$idx))
      points(px[rv$idx], rug_y[rv$idx], pch = 16, cex = 0.6,
             col = adjustcolor("red", 0.85))

    # axis up to the cut, then the break marker and the lane
    at <- pretty(c(0, cut), 6); at <- at[at <= cut]
    axis(1, at = at, labels = format(at, big.mark = ",", trim = TRUE))
    segments(cut, -0.62, cut, 1.05, col = "grey75", lty = 3)
    text(cut, -0.70, "//", col = "grey40", cex = 1.3, xpd = NA)
    text(cut + lane_w / 2, -0.70,
         sprintf("> %s\n(%d counties)", fmt(cut), n_tail),
         col = "grey40", cex = 0.75, xpd = NA)

    # current sample mean and interval
    if (!is.null(rv$mean)) {
      segments(rv$ci[1], -0.55, rv$ci[2], -0.55, col = "blue", lwd = 3)
      points(rv$mean, -0.55, pch = 19, cex = 1.3, col = "blue")
    }
    abline(v = ybarU, col = "black", lwd = 2, lty = 2)

    legend("topright", bty = "n", inset = c(0.12, 0),
           legend = c("Population", "Current sample", "Sample mean and CI", "True mean"),
           col = c("grey60", "red", "blue", "black"),
           pch = c(16, 16, 19, NA), lty = c(NA, NA, NA, 2))
  })

  output$histPlot <- renderPlot({
    par(mar = MAR, xaxs = "i")

    if (length(rv$means) == 0) {
      plot(1, type = "n", xlim = xlim_b, ylim = c(0, 1), bty = "n", yaxt = "n",
           xlab = "Sample mean of acres92", ylab = "",
           main = "Sampling distribution of the sample mean")
      abline(v = ybarU, col = "black", lwd = 2, lty = 2)
      text(ybarU, 0.5, "Click 'Draw new SRS sample' to begin", col = "grey50")
      return(invisible())
    }

    # fixed bin width; breaks extended to cover means outside the fixed window
    bw <- 2 * half / 40
    lo <- min(xlim_b[1], min(rv$means)) - bw
    hi <- max(xlim_b[2], max(rv$means)) + bw
    h  <- hist(rv$means, breaks = seq(lo, hi, by = bw), plot = FALSE)

    vis <- h$counts[h$mids >= xlim_b[1] & h$mids <= xlim_b[2]]
    top <- max(c(vis, 1))

    plot(h, col = "skyblue", border = "white",
         xlim = xlim_b, ylim = c(0, top * 1.30),
         xlab = "Sample mean of acres92", ylab = "Frequency",
         main = "Sampling distribution of the sample mean")
    abline(v = ybarU, col = "black", lwd = 2, lty = 2)

    # transient interval, for the current sample only
    ci_col <- if (rv$covered) "forestgreen" else "red"
    y      <- top * 1.17
    segments(rv$ci[1], y, rv$ci[2], y, col = ci_col, lwd = 3)
    segments(rv$ci, y - top * 0.045, rv$ci, y + top * 0.045, col = ci_col, lwd = 3)
    points(rv$mean, y, pch = 19, cex = 1.3, col = ci_col)

    outside <- sum(rv$means < xlim_b[1] | rv$means > xlim_b[2])
    leg <- c(
      sprintf("True mean = %s", fmt(ybarU)),
      sprintf("Samples drawn = %d", length(rv$means)),
      sprintf("Empirical coverage = %.1f%% (nominal %.0f%%)",
              100 * mean(rv$hits), 100 * input$conf)
    )
    if (outside > 0)
      leg <- c(leg, sprintf("%d mean(s) beyond the axis", outside))
    legend("topleft", bty = "n", cex = 0.95, legend = leg)
  })

  output$stats <- renderText({
    if (is.null(rv$idx)) return("No sample drawn yet.")
    sprintf(
      "Sample mean: %s\n%.0f%% CI: [%s, %s]\nContains true mean? %s\nMisses so far: %d of %d",
      fmt(rv$mean), 100 * input$conf, fmt(rv$ci[1]), fmt(rv$ci[2]),
      if (rv$covered) "Yes" else "No", sum(!rv$hits), length(rv$hits)
    )
  })
}

shinyApp(ui, server)

## file: acres92.csv
acres92
683533
47146
141338
210
50810
107259
167832
177189
48022
137426
144799
96427
73841
109555
121504
99466
67950
61426
68478
47200
175209
138135
82466
40832
166490
111315
196859
134555
233422
210733
104364
85872
85821
64755
130063
195536
128357
167583
166949
191810
204487
35748
56102
201892
173468
67962
207226
199714
138437
224370
199117
89228
142873
104342
110066
231243
155914
144193
106206
179319
96435
112620
78176
71697
167923
104199
78889
96194
50257
85086
141260
56680
411473
151325
92708
293745
250819
30196
18818
246184
269122
98919
313573
108046
34115
57253
167572
350402
145744
326808
324539
20589
262021
110260
210692
168755
223889
42794
37606
251710
168848
78498
105721
263182
183895
367969
281864
108913
107841
281895
298547
186685
143104
186829
382714
268075
142856
173861
484751
219444
79803
69422
102560
32003
67044
357416
70872
404585
122871
156363
313232
111895
253948
305401
45609
114762
195510
115019
131353
159013
136309
31190
119930
352322
358904
274843
190363
5785707
1891644
5989961
1151284
1846497
137834
246038
729947
1981938
7229585
3472248
1902452
334284
2108834
229365
286288
4768
236222
452347
246077
450236
163036
12594
102028
1774664
473920
597766
532866
247550
2839531
775829
164130
487499
183569
749465
168879
206138
725118
978831
686876
103294
1372778
235290
72471
60740
137723
119514
423602
379044
600073
1287057
517860
7
783715
1324403
57418
836989
342653
52905
388084
55446
647446
340328
517114
759649
318156
1016851
116083
1354262
137530
320597
518907
234781
685813
207448
322823
155465
1257229
796892
157493
84172
914094
7129
304592
330826
423785
156801
260728
167106
231364
213004
1105614
857404
331639
440581
13296
299142
177333
9021
641755
472018
103470
878447
1341738
14411
587339
540412
2286947
1660146
1066453
420233
15539
1159813
834018
447412
751517
633279
119287
388902
459659
32072
1004360
896994
546538
219612
576397
462086
200674
310394
38467
104010
1333577
2086292
1433111
9975
56510
86581
19830
25882
65987
38715
55263
197375
87134
304680
191140
24489
9135
36230
199724
23735
43314
227202
70672
86026
301977
96968
83681
334623
31693
40039
57179
52259
57853
70987
369965
14203
69405
327611
529835
61019
483835
265443
86706
174673
244185
118352
95833
199098
106721
100764
190553
11738
132208
299699
296242
190788
32
44962
56704
351885
138418
716542
637934
221232
4123
611336
105621
48839
300622
79270
151242
59642
253330
161936
48280
138208
8679
96730
45214
105538
77659
78739
108840
32976
49397
35851
85075
47071
129216
17105
62983
27561
168861
15948
213943
166511
29213
113861
17944
57074
82549
29451
21697
8518
5901
52651
33641
11559
42678
4519
13563
10192
178861
198184
26984
72636
41972
37973
109923
25802
19060
168593
3046
97215
156805
71135
8151
184137
16362
43775
54233
123702
40608
15577
22212
73659
36260
74641
21975
25376
28535
9681
73869
137637
46748
24239
36074
53944
35387
31529
31037
58529
24242
45624
73417
135247
83074
60811
72626
136082
77532
71379
31394
39712
40955
168051
104768
15583
32657
11969
73023
23284
33785
8003
120839
61757
45448
68729
121588
205573
44599
64901
93061
32950
4870
45845
51836
55310
18644
44470
18254
80905
45450
46014
80396
34746
11559
12733
95876
15974
12836
37923
138803
108967
24086
15521
49043
169989
38313
19314
119873
54356
71097
142824
174020
114487
88811
9910
32800
40783
98824
31161
21973
32865
88829
55779
53895
47000
111801
54445
53291
48755
24127
38691
115516
93078
31838
200061
926607
91998
214452
355786
328970
239800
321728
238609
268506
427215
299502
330080
236668
333115
341923
315448
345567
359755
347353
338801
308497
336254
274905
236409
314812
456954
368114
415104
312173
275319
261494
336131
192467
202249
343870
224811
401625
287586
343367
302352
366927
317205
328885
332377
329151
332358
399155
225835
260781
280797
272831
321285
346569
431185
227073
284537
321950
322401
615034
266083
349252
191291
219370
347599
305685
314887
268520
312858
237862
263047
392835
223638
240100
219832
362109
260780
318778
338730
518247
359442
229818
542855
340982
293266
364172
233217
353570
495769
331211
401858
278922
236265
241422
195021
302487
309508
282723
408462
231977
357684
442247
224632
353683
232879
221209
325338
269435
111510
1371605
266293
80333
150021
453647
72664
159358
129490
391050
587693
666342
286711
103246
140701
353528
230086
380928
197176
227114
744295
311296
207552
131281
347293
193908
211039
132429
224369
208161
477839
271143
752032
148776
435069
4428
134788
489993
78813
556131
464834
69354
182572
135163
144435
482169
99675
238906
209437
571807
390149
259923
223764
229120
263425
40917
223561
176012
377512
206271
259498
18206
354480
116312
257761
341274
300127
160533
431415
171938
303715
225506
201567
433246
37976
203974
453944
662629
186425
258014
217191
180675
290454
94681
203590
358920
178222
385560
73142
612112
169292
414442
637551
369952
344649
249240
709106
310518
402310
299709
253916
203749
282222
98838
164158
312128
187039
371936
311266
184599
392639
261482
167602
251277
443475
67998
82426
78081
270598
188999
175847
264140
141703
446750
207388
128867
402212
169622
314886
336450
119370
488215
115517
317467
297003
333238
234973
399312
325227
89591
203428
295844
197724
285730
165091
270618
87329
223429
22555
220057
227711
105658
162433
236073
59734
222435
86236
202429
153213
169265
193381
192311
111500
29837
229097
148662
194312
241049
196537
207766
162670
163248
161745
187079
190798
148609
187955
202896
301962
182836
130826
124694
139638
305634
251603
187549
144305
267695
158788
223328
38783
219402
71596
188843
59282
282764
139523
206885
184118
32318
116068
113129
181653
80078
85366
142482
220959
242777
204165
236436
164025
233183
172348
63332
217288
175124
134960
121710
181020
79235
257351
160930
80069
80958
119318
144722
197947
201739
96219
189136
189467
198680
285169
162244
282862
378517
245099
639327
580199
337300
339138
765688
351941
386881
271015
592207
565274
380969
407464
353371
486997
627612
302849
526064
514436
201798
222028
403375
324063
547483
442362
745371
671223
316317
164081
671506
512728
341608
517623
424104
603755
532890
499112
319686
366764
479903
340035
271713
484823
141386
517376
544071
399835
346519
419423
206530
482434
273841
603177
485656
537914
588061
572989
596103
286989
479310
323769
409839
427403
441417
326716
668420
465527
349293
547369
380403
449151
582053
451362
432326
641109
700869
443290
432701
228178
578283
427459
463690
403276
484415
510319
328094
227349
535359
620144
537457
436242
411785
450829
687593
702549
484093
423064
471658
521110
443802
312717
265978
22553
177858
156590
90033
111913
248634
132979
5419
80864
206881
27836
108291
99009
42602
266730
60911
140810
128807
137337
43447
78966
60812
112831
192189
299321
144904
68373
75409
125133
108777
250128
92487
60294
69310
147154
193859
10919
86074
96829
41352
138061
127161
210275
206090
134811
100468
69711
226206
4854
178217
200455
197826
159966
99066
144828
80692
44709
98545
23062
44188
3224
46321
120959
99527
48509
20803
1501
3383
159710
173892
119218
278675
49509
62766
13887
135179
247266
45133
175541
76141
5256
144254
119533
42642
133173
136869
165917
113383
105068
117868
191002
112409
159794
84434
176828
35712
127403
4469
6158
33155
218145
53022
92782
49712
91365
154082
229838
117768
93887
128719
165015
111866
71324
196701
252817
165391
135850
140432
44548
61145
123655
267511
117599
63446
67928
256118
136534
51401
110637
170353
327511
66380
258035
250607
63674
227698
147826
79150
194044
136812
176952
264321
44490
110173
81401
17086
4127
291526
87574
132678
27469
58730
36059
246536
241731
182605
100
74678
46110
193137
210570
97643
247106
57789
6166
23185
50319
42922
17347
283135
70936
81747
40181
126886
245986
44146
61883
316691
51799
116221
58790
38566
126839
86856
23089
5340
60980
34235
5757
25470
74484
37477
53459
31583
9882
72247
114805
37802
43320
83232
37320
126981
157505
80241
59389
123762
222768
110699
97312
44623
131283
82470
54459
165349
77491
55657
109108
123932
91254
107519
62242
334040
53893
38853
50076
95402
27622
24350
63473
118152
35988
18793
106971
71890
94755
61797
42572
16099
246403
77493
51517
82418
14104
165371
181052
19844
166886
227665
244927
186431
41037
40871
92809
64084
256236
1402
72777
28415
233921
40365
137082
61535
5965
66789
277400
231557
29161
438914
193688
254793
47308
29763
199733
210638
154482
16076
190706
264
18047
193956
64973
336273
118764
9391
22178
70306
48029
23290
73437
121153
110014
89173
88322
217095
224030
22056
73661
115338
48236
129083
75345
32980
108726
14081
36272
176305
79921
3786
318125
181569
234823
444407
13908
236799
324111
206781
188958
22488
31427
168073
61832
377693
224923
183760
262207
383373
347420
113422
165961
200199
326804
138594
566981
210897
1249
374920
130683
221193
241148
260125
414710
443496
366534
379603
269147
79183
272049
112412
131563
107810
401039
145545
360500
482991
68778
405029
5262
103665
205031
255453
395023
250507
186573
744710
412660
301111
142432
422916
392615
375628
241930
416570
457670
305831
821073
280089
263274
252658
1042850
310135
2142
183208
491726
600114
227519
270332
536299
153188
131753
117701
311849
643762
231610
286337
389897
395071
310184
245686
171412
237239
100774
249731
420778
290627
272540
407953
267066
227156
304032
377059
292426
311161
430451
239298
197530
271914
181292
255498
232189
339372
163076
252890
377000
54767
325796
187856
403597
211148
244810
130358
207611
187239
298709
201670
252783
227783
278841
210829
219440
300970
288810
296281
196959
245827
285496
226336
399193
321181
174314
232592
238909
372292
68596
134196
281327
119595
371022
268447
304560
356164
332910
249046
253281
337893
270576
199292
381934
111549
233304
219894
210253
242018
265245
217116
306175
224716
201714
368849
256023
507875
252074
316809
250475
291846
209452
359434
200766
323465
188595
345673
138986
254493
228936
221122
277322
89683
152529
204171
257217
168586
116910
54082
414394
165225
216694
219042
120036
273393
438142
137747
329999
160576
459671
402202
126474
111850
93053
289729
134028
316617
79981
78653
112896
107526
91367
426584
127351
151743
149027
42712
88522
68663
126352
294547
126613
79962
139591
37147
46532
43498
49282
99726
30050
16665
230838
223406
180102
113734
76673
24845
89168
66257
80902
96919
93352
98816
53401
80683
62833
95736
140209
262371
98914
125713
198955
89816
182009
175231
80272
137267
96474
201759
80761
218154
93180
31587
80342
124202
86096
186297
118651
110124
181946
96540
95121
32666
361003
273117
141245
108314
40676
230524
100433
108236
114083
342237
72515
75551
89807
83445
78230
361634
1342484
3002378
2338866
449970
598694
1619482
1424228
2277936
2085181
761459
1334041
135126
944497
2232575
277050
699409
2000266
1730537
636514
349938
1644001
367482
868064
631377
883479
951780
50220
1290134
1271160
912154
19158
248215
1031872
777803
669516
1968857
893872
1629363
675569
683088
241655
1197028
1414415
2585834
381104
962450
99746
889294
837904
1178885
1063086
599014
1688070
850599
490988
1454669
101073
52974
72621
70697
104426
19712
144529
170006
127760
39667
93584
31671
63067
31184
43056
64031
125428
62854
108363
23929
53902
16405
93970
162634
88386
98531
41750
7046
92192
68718
248518
19676
180400
47837
118291
34717
64532
8882
156027
112291
113654
204443
127663
69961
52281
75496
56693
93728
148135
13310
230402
68577
37434
141761
58384
21218
22089
93320
131767
27901
23007
36975
86982
179051
3486
155213
63992
67491
44000
83218
65136
68736
115854
194015
23140
144858
51916
292152
130879
104874
55309
266067
53690
89063
104733
121404
5897
12388
67698
167379
67716
119855
87478
102944
46726
179554
115106
142312
103773
38394
594587
858267
777675
818893
950031
677945
559385
877382
1070528
855458
627774
725974
1352738
369140
834293
366292
505461
769225
1019300
396154
688468
723816
669049
598832
1048701
544767
1165695
1128346
531643
1233663
1000679
552707
384213
600845
586244
639709
485012
503575
799606
522536
495509
521343
745815
785713
841736
439846
1269572
591185
501057
737273
1160916
750913
1182658
335465
489384
460156
407678
460513
437826
649612
296164
649634
587595
270005
335849
296016
428769
521389
3887635
772453
357067
228988
345739
1425338
138022
842110
658572
265048
242898
298854
96093
528731
340471
323315
526476
430972
508754
1069778
338136
229703
545799
304180
316551
321080
305724
401978
403584
1387740
375188
325330
326831
186806
310042
668957
446007
505150
612694
414763
1450481
335820
330369
464534
322120
290608
724458
236950
226042
333488
325801
223949
532901
375771
297326
409715
250086
439475
301513
657906
312079
105085
437274
417698
314949
1481503
298115
1005877
217228
347598
360203
193556
339358
228167
248502
307527
263976
345509
20910
25439
33935
46056
75733
39844
46610
34292
24716
38297
29606
2636
97186
7799
11644
68627
613
61748
0
106324
35786
25011
58758
23915
10365
1838
98256
43989
75531
325
87638
415263
1553328
3112271
2080760
2085387
924678
1343237
526407
1138681
1209335
1532887
1289733
843401
2149450
1881764
10
797117
3224090
905235
1166009
1769177
1552865
1646707
770155
1896131
2579730
517952
1233794
1868074
324476
1797466
2364443
349231
268043
82100
79635
3149788
1949420
235826
738041
494304
48968
189223
140380
624606
710618
231758
5361
57889
161643
0
97869
203704
254002
259540
58963
188008
158392
111974
138620
192116
109692
145679
54986
138299
35343
171722
45820
163072
300559
4
169313
205105
195626
110150
138822
1890
0
135494
242637
145329
181624
102733
133854
112334
218306
3803
92683
831
396721
70404
19198
117799
65323
115071
363293
35353
56002
114859
91822
69643
5811
205954
174627
5709
209889
102024
196365
194022
171129
155895
80507
205105
126195
202188
138297
122480
215796
179280
99214
227382
143773
161724
223216
4060
335575
196759
169017
88899
197736
235971
97020
205633
109820
65266
184280
127867
28983
275644
248400
113329
245049
230988
48050
177194
219023
74461
68344
210601
17138
61999
227327
202927
142624
74037
261320
73484
217397
103857
87036
269163
200405
109986
107157
113892
164607
187718
104197
106573
219037
96060
259390
87954
96293
204079
286698
160734
253383
202324
95704
296368
203498
136612
19088
120519
148479
222957
241787
41666
129416
139655
246938
187175
302456
216318
207333
487638
375683
986621
493631
513789
412673
726481
499876
372901
218803
300829
1034980
157105
264890
405831
358446
447212
336285
633874
242097
583098
699515
662121
419760
566152
599536
338529
286379
607180
268038
346572
469883
405271
323534
477655
521820
558313
194253
380529
397909
344280
255538
253652
314987
236766
494277
162690
277981
230832
347480
390957
282659
257577
156748
280533
1116090
207118
282211
328976
481244
353045
299263
239971
660214
309614
250958
215222
421233
1051384
480948
134096
216638
216268
577693
743563
687299
818736
118818
148848
24740
71839
174872
894853
74375
139483
402023
766373
1154399
1457339
27201
262251
530960
31249
720153
833025
242121
34292
380464
1318447
302462
1119004
31294
167880
487534
39559
1466580
473316
694304
1152965
139820
728131
179787
172366
32526
119566
57960
199056
221981
76466
310672
76790
129323
76997
19026
139918
176643
94817
55023
39412
101816
211037
141919
90298
5095
16371
167863
106390
4702
234391
88982
125707
129503
143036
79310
85113
36963
388368
86402
104519
82982
49850
132999
39561
160802
81426
20777
44425
41347
81479
109438
104292
6197
90065
89045
87253
219933
30613
177215
212477
63159
52760
67388
203026
121907
153897
62740
252052
1468
5636
9631
12408
20458
89935
136565
81257
160659
87355
74733
44800
51020
90995
32392
66165
94193
109652
135766
126370
156853
108848
62108
68858
55712
194822
37287
66809
70277
97241
195697
72500
52978
58323
128124
136151
82634
19486
78382
104862
93970
69897
262093
44200
66355
110679
107058
138573
55992
173188
120755
379973
724776
787857
322432
444440
1026353
496799
279202
1243168
417697
688081
534829
236608
392935
1701991
462238
270665
561312
341131
1859161
252419
641911
974811
560057
373787
601034
1204465
276744
861129
245174
1657305
390720
502469
545064
1361106
334057
584231
460063
297819
195077
322802
846435
325998
661474
485748
2076199
701352
313435
425288
284888
1066060
1726299
507101
604219
322784
1417516
890711
903980
615479
1079266
1006831
367239
259517
448834
271200
1406379
41899
213603
62989
93098
96181
91858
30299
96550
165547
36633
58289
71513
142729
70483
84029
132388
144918
96874
47319
86858
95818
144267
230906
258265
70457
135469
272121
256272
104457
236912
42735
57216
62542
80348
159927
110215
155600
224247
146868
191486
130167
44415
119419
87298
98669
54518
94254
91343
182754
196733
36978
275219
73654
123557
121848
138618
141357
50767
161902
245681
56253
100176
174807
48098
43202
257000
105519
53026
37550
31368
116696
52462
52433
233312
200097
32714
24832
74107
144953
150309
53795
92773
177522
183178
55097
11292
49452
32892
165309
117608
125092
204146
123792
204391
214497
352488
962576
103063
19131
612718
500665
765139
337351
432939
396508
394923
357933
442173
416631
408710
371257
629681
547829
263077
563993
295601
2405018
408824
566400
513533
317187
548351
263925
208073
490739
329288
68448
622130
166939
518316
251249
269146
449972
670459
370572
523049
680567
275638
462533
549167
207350
543750
571684
430377
605252
471498
394805
2001152
451584
1584367
780925
123756
553047
856707
107565
365618
569212
561521
677308
599637
801159
493227
518788
1117134
426189
581776
378003
412632
496742
545666
629148
321752
422464
128533
370140
749504
686578
102229
573827
688330
477515
465365
665421
576013
409501
47196
352915
347313
560355
443027
461249
576468
313952
32436
308344
201952
686578
494177
463450
545664
356170
660412
470096
490578
225852
386546
417187
489398
2234262
345138
402897
656961
519043
461829
70165
1524636
322324
736407
517671
330173
513418
383573
386991
354917
553226
595420
531206
774804
436040
698832
576893
424701
518028
432379
766037
527837
318658
482165
342213
428068
588500
558553
507135
346653
481539
491015
675927
472332
515960
243989
49579
512473
547428
562612
678590
658204
487573
724706
551148
428243
587316
495364
193885
587264
73948
479889
220355
525885
30268
536300
443224
593819
847608
56975
517272
195147
415694
526276
2891640
141215
402011
1695484
98449
497106
617851
362642
425281
1555905
667177
510079
391842
47446
572607
268058
33544
54580
82721
358211
743638
764723
518371
564382
187728
507449
247626
62850
632622
536507
835337
512247
926093
505140
167569
509017
1396275
459120
581511
180332
1020756
332826
109635
58690
194008
694909
917186
1806639
378217
430736
213923
242901
455873
328367
1712044
644730
501692
307783
863384
260892
545670
476493
432887
461127
203667
344667
563183
484907
723018
192288
1449976
267924
291860
21958
50357
399011
240535
137530
63116
434183
332686
209819
484156
234576
58522
493073
107663
324921
447463
158189
373582
437238
1294703
450315
139347
167374
105576
256522
91568
188567
25810
70411
90296
78691
287442
47010
200507
81768
96833
84677
8627
66034
134474
51604
113165
28827
112944
17392
68326
45451
115295
61669
85954
56280
15714
235533
116509
58326
166477
98142
73097
24478
51442
136320
37044
51374
232852
96282
24201
48968
96910
86247
9559
52508
37777
59326
19657
129296
195476
81427
85600
100602
6152
167858
21507
98914
72681
18367
52469
41290
63991
107700
64856
79449
297064
43180
68584
48889
32973
71803
78977
38940
24924
141766
236074
160973
134191
125394
119973
178469
52904
20107
52770
82736
138531
38967
190062
56289
13247
131366
2358
54290
83047
43332
209677
33682
96704
82849
17710
203503
24848
41348
93364
149503
132674
58891
43987
89785
996742
274546
640370
112085
24253
82967
304928
35678
918033
748088
670149
325472
1086045
44742
19526
9603
42290
10302
355360
689639
112263
1465788
10965
1291118
32637
55360
58750
20529
92074
4043
74153
625769
546303
59890
12611
710546
118136
1404289
1639965
119354
51208
350866
97521
205363
323482
84091
162205
386857
426884
327185
248956
538582
414240
130051
70547
366593
189905
21008
351633
26456
620951
293134
163145
361918
10263
218145
232591
195287
92761
170228
182339
356651
120383
85500
248862
529966
145980
135538
8763
346398
208888
31777
263514
78772
113548
272876
282405
265731
94596
133197
270930
343115
168304
308460
335517
47376
297984
207128
231427
348602
365511
226096
86091
147207
114184
241778
167191
169876
221357
76080
73430
2531
73407
12175
35587
34919
14868
59184
20455
52748
106325
179736
135577
7710
141742
88571
101214
74268
19956
81096
30015
3113
1088
40837
64332
117168
56555
74760
258
54622
148842
21871
33085
21164
178160
15650
115487
138688
55827
32633
104194
70960
82154
57717
41830
32093
47366
58678
28622
9335
37130
35836
59846
5693
1868333
441321
2704163
2720903
2363248
1542262
2415873
1234542
908320
2055522
1700356
558966
2508749
1344561
797270
1364948
1208776
592754
1720737
62307
879694
397883
1484583

```

