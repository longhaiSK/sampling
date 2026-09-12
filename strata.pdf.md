# Stratified Random Sampling

## Key Formulas for Stratified Sampling

Stratified random sampling divides a population of size $N$ into $H$ mutually exclusive strata. A simple random sample of size $n_h$ is then drawn from each stratum of size $N_h$. This design often yields higher precision than Simple Random Sampling (SRS).

#### 1. Stratified Population Mean Estimator {-}
The estimated population mean is a weighted average of the sample means from each stratum:
$$
\overline{y}_{str} = \sum_{h=1}^H \pi_h \overline{y}_h.
$$
Where $\pi_h = \frac{N_h}{N}$ is the stratum weight, and $\overline{y}_h$ is the sample mean of stratum $h$.

#### 2. Variance of the Stratified Mean {-}
The estimated variance accounts for the finite population correction in each stratum:
$$
\hat{V}(\overline{y}_{str}) = \sum_{h=1}^H \pi_h^2 \left(1 - \frac{n_h}{N_h}\right) \frac{s_h^2}{n_h}.
$$
Where $s_h^2$ is the sample variance in stratum $h$.

#### 3. Proportional Allocation {-}
Sample sizes are allocated proportional to the size of the stratum:
$$
n_h = n \left( \frac{N_h}{N} \right).
$$

#### 4. Neyman (Optimal) Allocation {-}
If surveying costs are equal across strata, Neyman allocation minimizes the variance of the estimator by sampling more heavily from larger and more highly variable strata:
$$
n_h = n \frac{N_h S_h}{\sum_{h=1}^H N_h S_h}
$$

---

## Functions and packages for Analyzing Data

The estimating functions used throughout this book (including `str_est`, `str_est_data`, and `srs_est`, used below to compute the stratified mean, standard error, and confidence intervals) live in a single shared file, `samplingestimate.r`, which we source below. `str_est_data` computes these directly from a dataset, while `str_est` computes them from pre-calculated summary statistics.


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

## Analysis of "agstrat.csv" data

We will estimate the average 1992 farm acreage (`acres92`) across US counties, using geographic region (`region`) as the stratification variable.

### Importing agstrat.csv data


::: {.cell}
::: {.cell-output-display}
\begin{table}
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}llrrrrrrrrrrrrlrr}
\toprule
county & state & acres92 & acres87 & acres82 & farms92 & farms87 & farms82 & largef92 & largef87 & largef82 & smallf92 & smallf87 & smallf82 & region & rn & weight \\ 
\midrule\addlinespace[2.5pt]
PIERCE COUNTY & NE & 297326 & 332862 & 319619 & 725 & 857 & 865 & 54 & 54 & 42 & 58 & 67 & 48 & NC & 805 & 10.23301 \\ 
JENNINGS COUNTY & IN & 124694 & 131481 & 139111 & 658 & 671 & 751 & 14 & 13 & 14 & 42 & 36 & 38 & NC & 241 & 10.23301 \\ 
WAYNE COUNTY & OH & 246938 & 263457 & 268434 & 1582 & 1734 & 1866 & 20 & 19 & 16 & 175 & 186 & 184 & NC & 913 & 10.23301 \\ 
VAN BUREN COUNTY & MI & 206781 & 190251 & 197055 & 1164 & 1278 & 1464 & 23 & 17 & 9 & 56 & 66 & 55 & NC & 478 & 10.23301 \\ 
OZAUKEE COUNTY & WI & 78772 & 85201 & 89331 & 448 & 483 & 527 & 6 & 5 & 5 & 56 & 49 & 48 & NC & 1028 & 10.23301 \\ 
CLEARWATER COUNTY & MN & 210897 & 229537 & 213105 & 583 & 699 & 693 & 34 & 32 & 23 & 8 & 19 & 13 & NC & 496 & 10.23301 \\ 
\bottomrule
\end{tabular*}
\end{table}

:::
:::

**Output Note:** The table above shows the sample survey data, including the stratum (region) and weights for each county.

### Spreadsheet calculation

#### Summarizing acre92 in each stratum

We first manually compute the sample size ($n_h$), standard deviation ($s_h$), and mean ($\overline{y}_h$) for each region, matching them with the known total counties in each region ($N_h$).


::: {.cell}
::: {.cell-output-display}
\begin{table}
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}l|rrrr}
\toprule
 & \(N_h\) & \(n_h\) & \(\bar{y}_h\) & \(s_h\) \\ 
\midrule\addlinespace[2.5pt]
NC & 1054 & 103 & 300,504.16 & 172,099.34 \\ 
NE & 220 & 21 & 97,629.81 & 87,449.83 \\ 
S & 1382 & 135 & 211,315.04 & 231,489.71 \\ 
W & 422 & 41 & 662,295.51 & 629,433.04 \\ 
\bottomrule
\end{tabular*}
\end{table}

:::
:::

**Output Note:** This displays the summary statistics per stratum. Note the significantly high variance ($s_h$) in the W region compared to the others.

#### Estimates

Rather than building the "spreadsheet-style" working table by hand, we call `str_est()`, whose `$table` shows the same intermediate stratum-level values ($N_h$, $n_h$, $\pi_h$, $\bar y_h$, $s_h^2$, $\pi_h\bar y_h$, $v_h$) needed for the final stratified mean and variance, with a Total row and footnotes for $\pi_h$ and $v_h$.


::: {.cell}
::: {.cell-output-display}
\begin{table}
\caption*{
{\fontsize{20}{25}\selectfont  Stratified Mean Estimation: Working Table\fontsize{12}{15}\selectfont }
} 
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}>{\raggedright\arraybackslash}p{\dimexpr 82.50pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 45.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 45.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 82.50pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 82.50pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 82.50pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 82.50pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 82.50pt -2\tabcolsep-1.5\arrayrulewidth}}
\toprule
Stratum & \(N_h\) & \(n_h\) & \(\pi_h\)\textsuperscript{\textit{1}} & \(\bar y_h\) & \(s_h^2\) & \(\pi_h\bar y_h\) & \(v_h\)\textsuperscript{\textit{2}} \\ 
\midrule\addlinespace[2.5pt]
NC & 1,054 & 103 & 0.34 & 3.01 $\times$ 10\textsuperscript{5} & 2.96 $\times$ 10\textsuperscript{10} & 1.03 $\times$ 10\textsuperscript{5} & 3.04 $\times$ 10\textsuperscript{7} \\ 
NE & 220 & 21 & 0.07 & 9.76 $\times$ 10\textsuperscript{4} & 7.65 $\times$ 10\textsuperscript{9} & 6.98 $\times$ 10\textsuperscript{3} & 1.68 $\times$ 10\textsuperscript{6} \\ 
S & 1,382 & 135 & 0.45 & 2.11 $\times$ 10\textsuperscript{5} & 5.36 $\times$ 10\textsuperscript{10} & 9.49 $\times$ 10\textsuperscript{4} & 7.22 $\times$ 10\textsuperscript{7} \\ 
W & 422 & 41 & 0.14 & 6.62 $\times$ 10\textsuperscript{5} & 3.96 $\times$ 10\textsuperscript{11} & 9.08 $\times$ 10\textsuperscript{4} & 1.64 $\times$ 10\textsuperscript{8} \\ 
{\bfseries Total} & {\bfseries 3,078} & {\bfseries 300} & {\bfseries 1.00} & {\bfseries } & {\bfseries } & {\bfseries 2.96 $\times$ 10\textsuperscript{5}} & {\bfseries 2.68 $\times$ 10\textsuperscript{8}} \\ 
\bottomrule
\end{tabular*}
\begin{minipage}{\linewidth}
\vspace{.05em}
\textsuperscript{\textit{1}}\(\pi_h = N_h/N\) is the stratum weight.\\
\textsuperscript{\textit{2}}\(v_h = (1-n_h/N_h)\,\pi_h^2\, s_h^2/n_h\) is stratum \(h\)'s contribution to \(V(\bar y_{str}) = \sum_h v_h\).\\
\textbf{Point estimate:} \(\bar y_{str} = \sum_h \pi_h\bar y_h = 2.956 \times 10^{5}\)\\
\textbf{Standard error:} \(\mathrm{SE}(\bar y_{str}) = \sqrt{\sum_h v_h} = \sqrt{2.683 \times 10^{8}} = 16,379.873\)\\
\end{minipage}
\end{table}

:::
:::

**Output Note:** The table provides a detailed breakdown of how each individual stratum contributes to the final mean and variance.

Extracting `$estimate` gives our overall result.


::: {.cell}
::: {.cell-output-display}
\begin{table}
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{0.7\linewidth}{@{\extracolsep{\fill}}rrrr}
\toprule
\(\bar{y}\) & \(\mathrm{SE}(\bar{y})\) & 95\% CI Lower & 95\% CI Upper \\ 
\midrule\addlinespace[2.5pt]
295,560.7652 & 16,379.8727 & 263,456.2147 & 327,665.3158 \\ 
\bottomrule
\end{tabular*}
\end{table}

:::
:::

**Output Note:** The overall result indicates an estimated mean of 295,561 acres, with a 95% confidence interval of [263,456, 327,665].

To find the population total (total acreage across all US counties), set `estimate = "total"`.


::: {.cell}
::: {.cell-output-display}
\begin{table}
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{0.7\linewidth}{@{\extracolsep{\fill}}rrrr}
\toprule
\(\hat{t}\) & \(\mathrm{SE}(\hat{t})\) & 95\% CI Lower & 95\% CI Upper \\ 
\midrule\addlinespace[2.5pt]
909,736,035.3920 & 50,417,248.2519 & 810,918,228.8183 & 1,008,553,841.9656 \\ 
\bottomrule
\end{tabular*}
\end{table}

:::
:::

**Output Note:** This scales the estimates up to the total area of 909,736,035 acres.

### Using the function "str_est_data"

Instead of pre-calculating the summary stats, this function extracts everything directly from the dataset using the pre-assigned sampling weights.


::: {.cell}
::: {.cell-output .cell-output-stdout}

```
  NC   NE    S    W 
1054  220 1382  422 
```


:::
:::

**Output Note:** This reconstructs the total population sizes ($N_h$) per region.


::: {.cell}
::: {.cell-output-display}
\begin{table}
\caption*{
{\fontsize{20}{25}\selectfont  Stratified Mean Estimation: Working Table\fontsize{12}{15}\selectfont }
} 
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}>{\raggedright\arraybackslash}p{\dimexpr 82.50pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 45.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 45.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 82.50pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 82.50pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 82.50pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 82.50pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 82.50pt -2\tabcolsep-1.5\arrayrulewidth}}
\toprule
Stratum & \(N_h\) & \(n_h\) & \(\pi_h\)\textsuperscript{\textit{1}} & \(\bar y_h\) & \(s_h^2\) & \(\pi_h\bar y_h\) & \(v_h\)\textsuperscript{\textit{2}} \\ 
\midrule\addlinespace[2.5pt]
NC & 1,054 & 103 & 0.34 & 3.01 $\times$ 10\textsuperscript{5} & 2.96 $\times$ 10\textsuperscript{10} & 1.03 $\times$ 10\textsuperscript{5} & 3.04 $\times$ 10\textsuperscript{7} \\ 
NE & 220 & 21 & 0.07 & 9.76 $\times$ 10\textsuperscript{4} & 7.65 $\times$ 10\textsuperscript{9} & 6.98 $\times$ 10\textsuperscript{3} & 1.68 $\times$ 10\textsuperscript{6} \\ 
S & 1,382 & 135 & 0.45 & 2.11 $\times$ 10\textsuperscript{5} & 5.36 $\times$ 10\textsuperscript{10} & 9.49 $\times$ 10\textsuperscript{4} & 7.22 $\times$ 10\textsuperscript{7} \\ 
W & 422 & 41 & 0.14 & 6.62 $\times$ 10\textsuperscript{5} & 3.96 $\times$ 10\textsuperscript{11} & 9.08 $\times$ 10\textsuperscript{4} & 1.64 $\times$ 10\textsuperscript{8} \\ 
{\bfseries Total} & {\bfseries 3,078} & {\bfseries 300} & {\bfseries 1.00} & {\bfseries } & {\bfseries } & {\bfseries 2.96 $\times$ 10\textsuperscript{5}} & {\bfseries 2.68 $\times$ 10\textsuperscript{8}} \\ 
\bottomrule
\end{tabular*}
\begin{minipage}{\linewidth}
\vspace{.05em}
\textsuperscript{\textit{1}}\(\pi_h = N_h/N\) is the stratum weight.\\
\textsuperscript{\textit{2}}\(v_h = (1-n_h/N_h)\,\pi_h^2\, s_h^2/n_h\) is stratum \(h\)'s contribution to \(V(\bar y_{str}) = \sum_h v_h\).\\
\textbf{Point estimate:} \(\bar y_{str} = \sum_h \pi_h\bar y_h = 2.956 \times 10^{5}\)\\
\textbf{Standard error:} \(\mathrm{SE}(\bar y_{str}) = \sqrt{\sum_h v_h} = \sqrt{2.683 \times 10^{8}} = 16,379.873\)\\
\end{minipage}
\end{table}

:::

::: {.cell-output-display}
\begin{table}
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{0.7\linewidth}{@{\extracolsep{\fill}}rrrr}
\toprule
\(\bar{y}\) & \(\mathrm{SE}(\bar{y})\) & 95\% CI Lower & 95\% CI Upper \\ 
\midrule\addlinespace[2.5pt]
295,560.7652 & 16,379.8727 & 263,456.2147 & 327,665.3158 \\ 
\bottomrule
\end{tabular*}
\end{table}

:::
:::

**Output Note:** Once again, this perfectly matches the manual calculation.


::: {.cell}
::: {.cell-output-display}
\begin{table}
\caption*{
{\fontsize{20}{25}\selectfont  Stratified Mean Estimation: Working Table\fontsize{12}{15}\selectfont }
} 
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}>{\raggedright\arraybackslash}p{\dimexpr 82.50pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 45.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 45.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 82.50pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 82.50pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 82.50pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 82.50pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 82.50pt -2\tabcolsep-1.5\arrayrulewidth}}
\toprule
Stratum & \(N_h\) & \(n_h\) & \(\pi_h\)\textsuperscript{\textit{1}} & \(\bar y_h\) & \(s_h^2\) & \(\pi_h\bar y_h\) & \(v_h\)\textsuperscript{\textit{2}} \\ 
\midrule\addlinespace[2.5pt]
NC & 1,054 & 103 & 0.34 & 44.26 & 1.29 $\times$ 10\textsuperscript{3} & 15.16 & 1.32 \\ 
NE & 220 & 21 & 0.07 & 47.24 & 2.36 $\times$ 10\textsuperscript{3} & 3.38 & 0.52 \\ 
S & 1,382 & 135 & 0.45 & 47.39 & 6.21 $\times$ 10\textsuperscript{3} & 21.28 & 8.36 \\ 
W & 422 & 41 & 0.14 & 124.39 & 1.01 $\times$ 10\textsuperscript{5} & 17.05 & 41.66 \\ 
{\bfseries Total} & {\bfseries 3,078} & {\bfseries 300} & {\bfseries 1.00} & {\bfseries } & {\bfseries } & {\bfseries 56.86} & {\bfseries 51.86} \\ 
\bottomrule
\end{tabular*}
\begin{minipage}{\linewidth}
\vspace{.05em}
\textsuperscript{\textit{1}}\(\pi_h = N_h/N\) is the stratum weight.\\
\textsuperscript{\textit{2}}\(v_h = (1-n_h/N_h)\,\pi_h^2\, s_h^2/n_h\) is stratum \(h\)'s contribution to \(V(\bar y_{str}) = \sum_h v_h\).\\
\textbf{Point estimate:} \(\bar y_{str} = \sum_h \pi_h\bar y_h = 56.863\)\\
\textbf{Standard error:} \(\mathrm{SE}(\bar y_{str}) = \sqrt{\sum_h v_h} = \sqrt{51.860} = 7.201\)\\
\end{minipage}
\end{table}

:::

::: {.cell-output-display}
\begin{table}
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{0.7\linewidth}{@{\extracolsep{\fill}}rrrr}
\toprule
\(\bar{y}\) & \(\mathrm{SE}(\bar{y})\) & 95\% CI Lower & 95\% CI Upper \\ 
\midrule\addlinespace[2.5pt]
56.8628 & 7.2014 & 42.7480 & 70.9776 \\ 
\bottomrule
\end{tabular*}
\end{table}

:::
:::

**Output Note:** This demonstrates the function's flexibility: we can easily estimate a different variable (an average of 56.9 small farms per county).

### Comparing with SRS estimate

To see the benefit of stratifying, we compare our stratified estimates against a Simple Random Sample approach on the same population.


::: {.cell}
::: {.cell-output-display}
\begin{table}
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{0.7\linewidth}{@{\extracolsep{\fill}}rrrr}
\toprule
\(\bar{y}\) & \(\mathrm{SE}(\bar{y})\) & 95\% CI Lower & 95\% CI Upper \\ 
\midrule\addlinespace[2.5pt]
297,897.0467 & 18,898.4344 & 260,706.2569 & 335,087.8365 \\ 
\bottomrule
\end{tabular*}
\end{table}

:::
:::

**Output Note:** The Simple Random Sample standard error is noticeably higher at 18,898.


::: {.cell}
::: {.cell-output .cell-output-stdout}

```
[1] 0.7512239
```


:::
:::

**Output Note:** The stratified variance is approximately 75.1% of the SRS variance.


::: {.cell}
::: {.cell-output .cell-output-stdout}

```
[1] 0.2487761
```


:::
:::

**Output Note:** Stratification by region reduced the variance by 24.9%. This is a substantial gain in precision for no extra sampling cost.

## Allocation of stratum sample size

Deciding how to distribute the total sample $n$ among the strata significantly impacts the precision of the estimator.

### Analyzing seals.csv collected with stratified sampling

This dataset investigates ringed seal holes in different zones. The researchers used Proportional Allocation, taking 20% of the areas in each zone.


::: {.cell}
::: {.cell-output-display}
\begin{table}
\caption*{
{\fontsize{20}{25}\selectfont  Stratified Mean Estimation: Working Table\fontsize{12}{15}\selectfont }
} 
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}>{\raggedright\arraybackslash}p{\dimexpr 82.50pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 45.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 45.00pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 82.50pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 82.50pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 82.50pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 82.50pt -2\tabcolsep-1.5\arrayrulewidth}>{\raggedleft\arraybackslash}p{\dimexpr 82.50pt -2\tabcolsep-1.5\arrayrulewidth}}
\toprule
Stratum & \(N_h\) & \(n_h\) & \(\pi_h\)\textsuperscript{\textit{1}} & \(\bar y_h\) & \(s_h^2\) & \(\pi_h\bar y_h\) & \(v_h\)\textsuperscript{\textit{2}} \\ 
\midrule\addlinespace[2.5pt]
1 & 68 & 17 & 0.34 & 1.76 & 3.32 & 0.60 & 0.02 \\ 
2 & 84 & 12 & 0.42 & 4.42 & 11.54 & 1.85 & 0.15 \\ 
3 & 48 & 11 & 0.24 & 10.55 & 46.07 & 2.53 & 0.19 \\ 
{\bfseries Total} & {\bfseries 200} & {\bfseries 40} & {\bfseries 1.00} & {\bfseries } & {\bfseries } & {\bfseries 4.99} & {\bfseries 0.35} \\ 
\bottomrule
\end{tabular*}
\begin{minipage}{\linewidth}
\vspace{.05em}
\textsuperscript{\textit{1}}\(\pi_h = N_h/N\) is the stratum weight.\\
\textsuperscript{\textit{2}}\(v_h = (1-n_h/N_h)\,\pi_h^2\, s_h^2/n_h\) is stratum \(h\)'s contribution to \(V(\bar y_{str}) = \sum_h v_h\).\\
\textbf{Point estimate:} \(\bar y_{str} = \sum_h \pi_h\bar y_h = 4.986\)\\
\textbf{Standard error:} \(\mathrm{SE}(\bar y_{str}) = \sqrt{\sum_h v_h} = \sqrt{0.348} = 0.590\)\\
\end{minipage}
\end{table}

:::
:::

**Output Note:** Zone 3 has the highest standard deviation ($s_h \approx 6.8$), indicating a high degree of variability in seal density in that region.


::: {.cell}
::: {.cell-output-display}
\begin{table}
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{0.7\linewidth}{@{\extracolsep{\fill}}rrrr}
\toprule
\(\bar{y}\) & \(\mathrm{SE}(\bar{y})\) & 95\% CI Lower & 95\% CI Upper \\ 
\midrule\addlinespace[2.5pt]
4.9859 & 0.5901 & 3.8292 & 6.1426 \\ 
\bottomrule
\end{tabular*}
\end{table}

:::
:::

**Output Note:** This provides the estimated mean number of seal holes per square kilometer (5).

### Neyman allocation of stratum sample size

Because Zone 3 was highly variable, we could have achieved a lower overall variance if we had sampled it more heavily. Neyman allocation formally calculates this optimal distribution.


::: {.cell}
::: {.cell-output-display}
\begin{table}
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}l|rrrrr}
\toprule
 & \(C_h\) & \(N_h\) & \(S_h\) & \(N_h S_h / \sqrt{C_h}\) & \(L_h\) (Prop) \\ 
\midrule\addlinespace[2.5pt]
1 & 1 & 68 & 1.821 & 123.831 & 0.168 \\ 
2 & 1 & 84 & 3.397 & 285.327 & 0.388 \\ 
3 & 1 & 48 & 6.788 & 325.809 & 0.443 \\ 
\bottomrule
\end{tabular*}
\end{table}

:::
:::

**Output Note:** This displays the optimal allocation proportions ($L_h$). Even though Zone 2 is a larger area ($N_h=84$), Zone 3's high variance means we should allocate the largest share (44%) of our sampling effort there.

## Simulation to Study the Efficiency of Stratified Sampling with Different Allocation

We will use the full `agpop.csv` dataset (acting as our true population) to run 2000 repeated samples. We will compare SRS, Proportional (P2S), and Optimal (Neyman) allocations.

### Using "Region" to Stratify

#### Read Population Data


::: {.cell}

:::


#### Define Stratum Variable


::: {.cell}

:::


A stratification variable is most effective when it strongly explains the variance in the target variable. We check this via ANOVA and R-squared.


::: {.cell}
::: {.cell-output-display}
![](strata_files/figure-pdf/agpop-region-anova-1.pdf)
:::

::: {.cell-output .cell-output-stdout}

```
Analysis of Variance Table

Response: agpop$acres92
                Df     Sum Sq    Mean Sq F value    Pr(>F)    
agpop$stratum    3 1.0073e+14 3.3578e+13  226.73 < 2.2e-16 ***
Residuals     3055 4.5243e+14 1.4810e+11                      
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```


:::

::: {.cell-output-display}
\begin{table}
\caption*{
{\fontsize{20}{25}\selectfont  Model Fit Summary\fontsize{12}{15}\selectfont }
} 
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{0.7\linewidth}{@{\extracolsep{\fill}}rrrrr}
\toprule
\(R^2\) & Adj. \(R^2\) & \(\hat\sigma\) & F-statistic & p-value \\ 
\midrule\addlinespace[2.5pt]
0.1821 & 0.1813 & 384,831.7694 & 226.7299 & 0.0000 \\ 
\bottomrule
\end{tabular*}
\end{table}

:::
:::

**Output Note:** The R-squared is 0.182. The region variable only explains about 18.2% of the variance in acreage, which is helpful but not extremely strong.

#### Stratified Sampling with P2S allocation (Single Run)


::: {.cell}
::: {.cell-output-display}
\begin{table}
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{0.7\linewidth}{@{\extracolsep{\fill}}rrrr}
\toprule
\(\bar{y}\) & \(\mathrm{SE}(\bar{y})\) & 95\% CI Lower & 95\% CI Upper \\ 
\midrule\addlinespace[2.5pt]
295,589.9412 & 18,918.1011 & 258,510.4630 & 332,669.4194 \\ 
\bottomrule
\end{tabular*}
\end{table}

:::
:::

**Output Note:** This provides a single point estimate for demonstration.

#### Repeat stratified sampling with P2S allocation 2000 times


::: {.cell}

:::


#### Repeat stratified sampling with optimal allocation 2000 times


::: {.cell}
::: {.cell-output-display}
\begin{table}
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}l|rrr}
\toprule
 & Nh & Sh & nh\_opt \\ 
\midrule\addlinespace[2.5pt]
NC & 1052 & 271,187.98 & 87 \\ 
NE & 213 & 78,906.20 & 5 \\ 
S & 1376 & 244,131.98 & 102 \\ 
W & 418 & 836,613.56 & 106 \\ 
\bottomrule
\end{tabular*}
\end{table}

:::
:::

**Output Note:** This table compares the population sizes to the optimal sample sizes. Notice that the West region (high variance) gets a disproportionately larger sample under Neyman allocation.


::: {.cell}

:::


#### Repeat simple random sampling 2000 times


::: {.cell}

:::


#### Compare the efficiency of different methods


::: {.cell}
::: {.cell-output-display}
![](strata_files/figure-pdf/agpop-region-compare-plot-1.pdf)
:::
:::

**Output Note:** The plot visually demonstrates that Neyman allocation has the tightest spread around the true mean (the red line).


::: {.cell}
::: {.cell-output-display}
\begin{table}
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}l|rrrr}
\toprule
 & Mean & Variance & Relative Variance & Percentage of Variance Reduction \\ 
\midrule\addlinespace[2.5pt]
SRS & 308,818.2 & 547,160,924.7 & 1.000 & 0.00\% \\ 
Prop2size & 308,934.9 & 452,690,789.8 & 0.827 & 17.27\% \\ 
Neyman & 308,340.2 & 294,427,723.7 & 0.538 & 46.19\% \\ 
\bottomrule
\end{tabular*}
\end{table}

:::
:::

**Output Note:** The table summarizes the simulation. Proportional allocation reduces variance by 17.3%, while Neyman allocation reduces it by 46.2% compared to Simple Random Sampling.

### Using "acres82" to Define Strata

What if we stratify by a variable highly correlated with our target, like the acreage from a decade prior (`acres82`)?

#### Read Population Data


::: {.cell}

:::


#### Define Stratum Variable with Quantiles of "acres82"
We convert the continuous `acres82` variable into four distinct categorical strata based on quantiles.


::: {.cell}

:::


Let's check the explanatory power of this new stratification variable.


::: {.cell}
::: {.cell-output-display}
![](strata_files/figure-pdf/agpop-acres82-anova-1.pdf)
:::

::: {.cell-output-display}
\begin{table}
\caption*{
{\fontsize{20}{25}\selectfont  Model Fit Summary\fontsize{12}{15}\selectfont }
} 
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{0.7\linewidth}{@{\extracolsep{\fill}}rrrrr}
\toprule
\(R^2\) & Adj. \(R^2\) & \(\hat\sigma\) & F-statistic & p-value \\ 
\midrule\addlinespace[2.5pt]
0.7393 & 0.7390 & 217,264.2290 & 2,887.8873 & 0.0000 \\ 
\bottomrule
\end{tabular*}
\end{table}

:::
:::

**Output Note:** The R-squared is 0.739. The past acreage explains 73.9% of the variance in 1992 acreage, making this a vastly superior stratification variable.

#### Stratified sampling with P2S allocation


::: {.cell}

:::


#### Repeat stratified sampling with optimal allocation 2000 times


::: {.cell}

:::


#### Repeat simple random sampling 2000 times


::: {.cell}

:::


#### Compare the efficiency of different methods


::: {.cell}
::: {.cell-output-display}
![](strata_files/figure-pdf/agpop-acres82-compare-plot-1.pdf)
:::
:::

**Output Note:** The boxplot visually highlights a massive reduction in variance compared to SRS, especially for the Neyman allocation.


::: {.cell}
::: {.cell-output-display}
\begin{table}
\fontsize{12.0pt}{14.0pt}\selectfont
\begin{tabular*}{\linewidth}{@{\extracolsep{\fill}}l|rrrr}
\toprule
 & Mean & Variance & Relative Variance & Percentage of Variance Reduction \\ 
\midrule\addlinespace[2.5pt]
SRS & 309,511.5 & 549,602,264.4 & 1.000 & 0.00\% \\ 
Prop2size & 308,267.4 & 133,504,327.8 & 0.243 & 75.71\% \\ 
Neyman & 308,690.0 & 38,783,485.8 & 0.071 & 92.94\% \\ 
\bottomrule
\end{tabular*}
\end{table}

:::
:::

**Output Note:** The table shows dramatic results: Proportional allocation reduces variance by 75.7%, and Neyman allocation reduces it by 92.9%. Stratifying by a highly correlated variable is immensely powerful.