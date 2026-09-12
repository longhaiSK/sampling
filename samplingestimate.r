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
