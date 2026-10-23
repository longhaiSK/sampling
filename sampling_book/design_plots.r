## Illustrations of the three basic sampling designs (SRS, stratified, cluster)
## on synthetic populations of N = 1000 units. Each design selects n = 60 units.
## Used at the start of the introduction, SRS, stratified and cluster chapters.

## stratum and cluster sizes (both sum to N = 1000)
design_strat_sizes <- c (400, 300, 200, 100)
design_clus_sizes  <- c (60, 150, 80, 120, 40, 100, 130, 70, 160, 90)

## run `expr` with a fixed seed, leaving the caller's random number stream untouched
with_seed <- function (seed, expr)
{
    old <- if (exists (".Random.seed", envir = globalenv ())) get (".Random.seed", envir = globalenv ())
    on.exit (if (is.null (old)) rm (".Random.seed", envir = globalenv ())
             else assign (".Random.seed", old, envir = globalenv ()))
    set.seed (seed)
    expr
}

## sort the groups by their means, and lay the units out with gaps between groups
sort_groups <- function (pop, gap)
{
    m <- tapply (pop$y, pop$group, mean)
    pop$group <- match (pop$group, order (m))
    pop <- pop[order (pop$group), ]
    pop$x <- seq_len (nrow (pop)) + gap * (pop$group - 1)
    pop
}

## four strata of sizes 400, 300, 200 and 100: exponential with mean 1.5,
## shifted by 0, 4, 8 and 12
design_strat_pop <- with_seed (348, {
    pop <- data.frame (group = rep (1:4, times = design_strat_sizes))
    pop$y <- c (0, 4, 8, 12)[pop$group] + rexp (1000, rate = 1 / 1.5)
    sort_groups (pop, gap = 30)
})

## ten clusters of sizes 40 to 160 (total 1000): exponential with means from 2 to 8
design_clus_pop <- with_seed (349, {
    pop <- data.frame (group = rep (1:10, times = design_clus_sizes))
    pop$y <- rexp (1000, rate = 1 / seq (2, 8, length.out = 10)[pop$group])
    sort_groups (pop, gap = 12)
})

## points: unselected units in grey, selected units in red
draw_units <- function (x, y, sel)
{
    points (x[!sel], y[!sel], pch = 16, cex = 0.45, col = adjustcolor ("grey55", 0.6))
    points (x[sel],  y[sel],  pch = 19, cex = 0.8, col = "red")
}

## plot a population by group; `sel` marks the sampled units, `shade` the groups to highlight
plot_groups <- function (pop, sel, label, shade = integer (0))
{
    H <- max (pop$group)
    par (mar = c (3.6, 4, 1, 0.5), mgp = c (2.2, 0.7, 0))
    plot (pop$x, pop$y, type = "n", ylim = c (0, max (pop$y) * 1.03), xaxt = "n",
          xlab = label, ylab = "y", las = 1)
    ## each group gets a background panel, leaving a white strip between neighbours
    pad <- (min (pop$x[pop$group == 2]) - max (pop$x[pop$group == 1])) / 2 - 2
    for (h in 1:H)
    {
        xs <- range (pop$x[pop$group == h])
        rect (xs[1] - pad, par ("usr")[3], xs[2] + pad, par ("usr")[4], border = NA,
              col = if (h %in% shade) adjustcolor ("cyan", 0.3) else "grey92")
    }
    draw_units (pop$x, pop$y, sel)
    for (h in 1:H)
    {
        xs <- range (pop$x[pop$group == h])
        m  <- mean (pop$y[pop$group == h])
        segments (xs[1], m, xs[2], m, col = "black", lwd = 2, lty = 2)
        axis (1, at = mean (xs), labels = h, tick = FALSE)
    }
}

## draw an SRS within each group listed in `groups`; k is one size, or one size per group
srs_within <- function (pop, groups, k)
{
    k <- rep_len (k, length (groups))
    sel <- logical (nrow (pop))
    for (j in seq_along (groups)) sel[sample (which (pop$group == groups[j]), k[j])] <- TRUE
    sel
}

## SRS of 60 from the 1000 units of the stratified population, listed in an arbitrary order
plot_design_srs <- function (seed = 1) with_seed (seed, {
    y <- design_strat_pop$y[sample (nrow (design_strat_pop))]
    sel <- seq_along (y) %in% sample (length (y), 60)
    par (mar = c (3.6, 4, 1, 0.5), mgp = c (2.2, 0.7, 0))
    plot (seq_along (y), y, type = "n", ylim = c (0, max (y) * 1.03),
          xlab = "unit", ylab = "y", las = 1)
    draw_units (seq_along (y), y, sel)
    abline (h = mean (y), col = "black", lwd = 2, lty = 2)
})

## stratified sample with proportional allocation: n_h = 60 N_h / N = 24, 18, 12, 6
plot_design_strata <- function (seed = 2) with_seed (seed, {
    Nh <- tabulate (design_strat_pop$group)
    plot_groups (design_strat_pop, srs_within (design_strat_pop, 1:4, 60 * Nh / sum (Nh)), "Stratum")
})

## two-stage cluster sample: 3 of the 10 clusters, then an SRS of 20 in each
plot_design_cluster <- function (seed = 4) with_seed (seed, {
    chosen <- sort (sample (10, 3))
    plot_groups (design_clus_pop, srs_within (design_clus_pop, chosen, 20), "Cluster",
                 shade = chosen)
})

## the cell diagram of the notation: clusters are columns of cells y_ij, with
## t_i on top and M_i at the bottom; `sampled = TRUE` shades the sampled
## clusters and highlights the subsampled cells
cluster_cells <- function (sampled = FALSE)
{
    vd <- "vdots"
    cols <- list (
        list (rows = list (quote (y[11]), quote (y[12]), quote (y[13]), vd, quote (y[1*M[1]])),
              top = quote (t[1]), bot = quote (M[1]), m = quote (m[1]), sel = FALSE, pick = integer (0)),
        list (rows = list (quote (y[21]), quote (y[22]), vd, quote (y[2*M[2]])),
              top = quote (t[2]), bot = quote (M[2]), m = quote (m[2]), sel = TRUE, pick = c (1, 4)),
        "hdots",
        list (rows = list (quote (y[i*1]), quote (y[i*2]), vd, quote (y[i*j]), vd, quote (y[i*M[i]])),
              top = quote (t[i]), bot = quote (M[i]), m = quote (m[i]), sel = TRUE, pick = c (2, 4)),
        "hdots",
        list (rows = list (quote (y[N*1]), quote (y[N*2]), quote (y[N*3]), vd, quote (y[N*M[N]])),
              top = quote (t[N]), bot = quote (M[N]), m = quote (m[N]), sel = FALSE, pick = integer (0))
    )
    w <- 1.5; gap <- 0.5; hw <- 0.8
    widths <- sapply (cols, function (cl) if (is.character (cl)) hw else w)
    x0 <- cumsum (c (0, head (widths, -1) + gap))
    par (mar = c (0, 0, 0, 0))
    plot (NA, xlim = c (-0.2, max (x0 + widths) + 0.2), ylim = c (-7.6, 1.1),
          axes = FALSE, xlab = "", ylab = "", xaxs = "i", yaxs = "i")
    vdots <- function (x, y) points (rep (x, 3), y + c (-0.22, 0, 0.22), pch = 16, cex = 0.45)
    for (k in seq_along (cols))
    {
        cl <- cols[[k]]; x <- x0[k]
        if (is.character (cl))
        {
            points (x + hw / 2 + c (-0.22, 0, 0.22), rep (-3, 3), pch = 16, cex = 0.6)
            next
        }
        R <- length (cl$rows)
        shade <- sampled && cl$sel
        if (shade)
            rect (x - 0.12, -R - 1.3, x + w + 0.12, 1.05, col = adjustcolor ("cyan", 0.3), border = NA)
        for (r in seq_len (R))
        {
            picked <- shade && r %in% cl$pick
            rect (x, -r, x + w, -(r - 1), border = "grey30",
                  col = if (picked) "#FFD0D0" else "white")
            if (identical (cl$rows[[r]], vd)) vdots (x + w / 2, -r + 0.5)
            else text (x + w / 2, -r + 0.5, as.expression (cl$rows[r]), cex = 1.05,
                       col = if (picked) "red3" else "black")
        }
        top <- if (shade) substitute (hat (T), list (T = cl$top)) else cl$top
        text (x + w / 2, 0.55, as.expression (list (top)), cex = 1.25, col = "blue4")
        text (x + w / 2, -R - 0.5, as.expression (list (cl$bot)), cex = 1.25, col = "darkgreen")
        if (shade) text (x + w / 2, -R - 1.0, as.expression (list (cl$m)), cex = 1.1, col = "red3")
    }
}

## a population of rectangular grids of cells; sampled cells are filled red
grid_panel <- function (sizes, fill, title)
{
    plot (NA, xlim = c (0, 20), ylim = c (-13.2, 1.2), asp = 1, axes = FALSE, xlab = "", ylab = "")
    text (0, 0.6, title, adj = 0, font = 2)
    x <- 0; y <- 0; rowh <- 0
    for (k in seq_len (nrow (sizes)))
    {
        nr <- sizes[k, 1]; nc <- sizes[k, 2]
        if (x + nc > 20) { x <- 0; y <- y - rowh - 1; rowh <- 0 }
        for (a in 1:nr) for (b in 1:nc)
            rect (x + b - 1, y - a, x + b, y - a + 1, border = "grey30",
                  col = if (fill[[k]][(a - 1) * nc + b]) "red" else "white")
        x <- x + nc + 1; rowh <- max (rowh, nr)
    }
}

## stratified versus one-stage cluster sampling on the same population of 12 groups
plot_cluster_vs_strat <- function () with_seed (7, {
    sizes <- matrix (c (4,5, 3,4, 4,3, 3,3,
                        3,5, 3,3, 3,5, 2,2,
                        3,3, 4,5, 3,3, 4,4), ncol = 2, byrow = TRUE)
    strat <- lapply (seq_len (nrow (sizes)), function (k) {
        K <- prod (sizes[k, ]); f <- logical (K); f[sample (K, max (1, round (K / 4)))] <- TRUE; f
    })
    chosen <- c (3, 6, 12)
    clus <- lapply (seq_len (nrow (sizes)), function (k) rep (k %in% chosen, prod (sizes[k, ])))
    par (mfrow = c (1, 2), mar = c (0.2, 0.4, 0.2, 0.4))
    grid_panel (sizes, strat, "Stratified: an SRS from every stratum")
    grid_panel (sizes, clus,  "Cluster: an SRS of clusters")
})

## clusters whose means differ a lot (R_a^2 large) versus random-looking clusters (R_a^2 small)
plot_cluster_r2 <- function () with_seed (5, {
    n_each <- 10
    par (mfrow = c (1, 2), mar = c (3, 3.5, 2.5, 0.5))
    panel <- function (means, sd, main)
    {
        x <- rep (1:4, each = n_each) + runif (4 * n_each, -0.15, 0.15)
        y <- rep (means, each = n_each) + rnorm (4 * n_each, 0, sd)
        plot (x, y, pch = 19, col = "blue4", xaxt = "n", xlim = c (0.5, 4.5),
              ylim = c (0, 11), xlab = "", ylab = "y", main = main)
        axis (1, at = 1:4, labels = paste ("cluster", 1:4))
        for (h in 1:4) segments (h - 0.25, means[h], h + 0.25, means[h], col = "firebrick", lwd = 3)
    }
    panel (c (1.5, 4, 7, 9.5), 0.3, expression (R[a]^2~"large: cluster sampling loses more"))
    panel (c (4.7, 5, 5.3, 5.1), 1.6, expression (R[a]^2~"small: cluster sampling loses less"))
})

## cluster totals t_i against sizes M_i: the horizontal line at tbar (what the
## unbiased estimator uses) and the ratio line t = ybar_r M through the origin
plot_cluster_ratio <- function () with_seed (2, {
    n <- 30
    Mi <- runif (n, 1, 10)
    ti <- 0.6 * Mi + rnorm (n, 0, 0.4)
    Bhat <- sum (ti) / sum (Mi)
    par (mar = c (4, 4, 0.5, 0.5))
    plot (Mi, ti, pch = 19, xlab = expression (M[i]), ylab = expression (t[i]),
          xlim = c (0, 10.5), ylim = c (0, 7), las = 1)
    abline (h = mean (ti), col = "blue4", lwd = 2, lty = 2)
    abline (0, Bhat, col = "firebrick", lwd = 2)
    segments (Mi, Bhat * Mi, Mi, ti, col = "gray50", lty = 3)
    legend ("bottomright", bty = "n", lwd = 2, lty = c (2, 1), col = c ("blue4", "firebrick"),
            legend = c (expression (bar (t) ~~ "(unbiased estimator)"),
                        expression (t == hat (bar (y))[r] * M ~~ "(ratio estimator)")))
})

## coot egg volumes: by clutch, and with clutches ranked by their means
plot_coots <- function (coots)
{
    par (mfrow = c (1, 2), mar = c (4, 4, 0.5, 0.5))
    plot (coots$clutch, coots$volume, pch = 19, cex = 0.6,
          xlab = "Clutch number", ylab = "Egg volume", las = 1)
    rng <- tapply (coots$volume, coots$clutch, range)
    lo  <- sapply (rng, `[`, 1); hi <- sapply (rng, `[`, 2)
    ord <- order ((lo + hi) / 2)
    plot (NA, xlim = c (1, length (ord)), ylim = range (coots$volume),
          xlab = "Clutch ranked by mean", ylab = "Egg volume", las = 1)
    segments (seq_along (ord), lo[ord], seq_along (ord), hi[ord], lwd = 1.5)
}

## ratio estimation: an SRS (green), the fitted line yhat = Bhat x (red) and the fitted
## values for every unit; the ratio estimate is their mean (copied from ratio.qmd)
plot_ratio_concept <- function () with_seed (2046, {
  N <- 40
  B_true <- 0.55
  x <- runif(N, 0.5, 10)
  y <- B_true * x + rnorm(N, 0, 0.20 * sqrt(x))

  xbar_U <- mean(x)
  ybar_U <- mean(y)

  n <- 8
  sel <- sample(N, n)               # an ordinary SRS, spread across x
  Bhat <- sum(y[sel]) / sum(x[sel]) # ratio estimate, computed from the sample only

  yhat <- Bhat * x                  # fitted value for every unit in the population
  ybar_ratio <- mean(yhat)          # = Bhat * xbar_U (the ratio estimate)
  ybar_srs <- mean(y[sel])          # naive sample mean, ignoring x

  par(mar = c(4, 4, 0.5, 0.5))
  plot(x, y, pch = 19, col = "grey50", xlab = "x", ylab = "y",
       xlim = c(0, 10.5), ylim = c(0, 7))

  abline(0, Bhat, col = "firebrick", lwd = 2)
  text(7.0, Bhat * 7.0 + 0.55, expression(y == hat(B) * x), col = "firebrick")

  points(x[sel], y[sel], pch = 21, col = "darkgreen", lwd = 2, cex = 2.2)
  points(x, yhat, pch = 19, col = "firebrick", cex = 0.9)

  abline(h = ybar_U, col = "black", lty = 3)
  abline(h = ybar_ratio, col = "firebrick", lty = 3)
  abline(h = ybar_srs, col = "darkgreen", lty = 3)

  points(0, ybar_U, pch = 19, col = "black", cex = 2)
  points(0, ybar_ratio, pch = 19, col = "firebrick", cex = 1.6)
  points(0, ybar_srs, pch = 19, col = "darkgreen", cex = 1.6)

  text(0.2, ybar_U, expression(bar(y)[U]), col = "black", adj = c(0, -0.7), font = 2)
  text(0.2, ybar_ratio, expression(hat(bar(y))[r]), col = "firebrick", adj = c(0, 1.4), font = 2)
  text(0.2, ybar_srs, expression(bar(y)), col = "darkgreen", adj = c(0, -0.7), font = 2)
})

## unequal probability sampling on the SRS population: each unit's inclusion
## probability is proportional to a size measure x_i = y_i times a lognormal error,
## shown by the point size; n = 60 units are drawn by systematic PPS sampling
plot_design_ups <- function (seed = 3) with_seed (seed, {
    y <- with_seed (1, design_strat_pop$y[sample (nrow (design_strat_pop))])  # same listing as the SRS plot
    x <- y * exp (rnorm (length (y), 0, 0.4))
    pik <- 60 * x / sum (x)
    ## systematic PPS: a random start, then every unit whose cumulative pik crosses an integer
    cum <- cumsum (pik)
    sel <- logical (length (y))
    sel[findInterval (runif (1) + 0:59, c (0, cum), left.open = TRUE)] <- TRUE
    cex <- 0.25 + 1.6 * sqrt (pik / max (pik))
    par (mar = c (3.6, 4, 1, 0.5), mgp = c (2.2, 0.7, 0))
    plot (seq_along (y), y, type = "n", ylim = c (0, max (y) * 1.03),
          xlab = "unit", ylab = "y", las = 1)
    points (which (!sel), y[!sel], pch = 16, cex = cex[!sel], col = adjustcolor ("grey55", 0.5))
    points (which (sel),  y[sel],  pch = 19, cex = cex[sel],  col = "red")
    abline (h = mean (y), col = "black", lwd = 2, lty = 2)
})

## regression estimation: an SRS (green), the least-squares line (red) and the fitted
## values for every unit; the regression estimate is their mean
plot_reg_concept <- function () with_seed (498, {
    N <- 40
    B0 <- 1.5
    B1 <- 0.35
    x <- runif(N, 1, 10)
    y <- B0 + B1 * x + rnorm(N, 0, 0.30)

    xbar_U <- mean(x)
    ybar_U <- mean(y)

    n <- 8
    sel <- sample(N, n)               # an ordinary SRS, spread across x
    xs <- x[sel]
    ys <- y[sel]
    B1hat <- sum((xs - mean(xs)) * (ys - mean(ys))) / sum((xs - mean(xs))^2)
    B0hat <- mean(ys) - B1hat * mean(xs)  # least squares fit, from the sample only

    yhat <- B0hat + B1hat * x         # fitted value for every unit in the population
    ybar_reg <- mean(yhat)            # = B0hat + B1hat * xbar_U (the regression estimate)
    ybar_srs <- mean(ys)              # naive sample mean, ignoring x

    par(mar = c(4, 4, 0.5, 0.5))
    plot(x, y, pch = 19, col = "grey50", xlab = "x", ylab = "y",
         xlim = c(0, 10.5), ylim = c(0, 6.5))

    abline(B0hat, B1hat, col = "firebrick", lwd = 2)
    text(10.4, B0hat + B1hat * 10.4 - 1.25, expression(y == hat(B)[0] + hat(B)[1] * x),
         col = "firebrick", adj = c(1, 0))

    points(x[sel], y[sel], pch = 21, col = "darkgreen", lwd = 2, cex = 2.2)
    points(x, yhat, pch = 19, col = "firebrick", cex = 0.9)

    abline(h = ybar_U, col = "black", lty = 3)
    abline(h = ybar_reg, col = "firebrick", lty = 3)
    abline(h = ybar_srs, col = "darkgreen", lty = 3)

    points(0, ybar_U, pch = 19, col = "black", cex = 2)
    points(0, ybar_reg, pch = 19, col = "firebrick", cex = 1.6)
    points(0, ybar_srs, pch = 19, col = "darkgreen", cex = 1.6)

    text(0.2, ybar_U, expression(bar(y)[U]), col = "black", adj = c(0, -0.7), font = 2)
    text(0.2, ybar_reg, expression(hat(bar(y))[reg]), col = "firebrick", adj = c(0, 1.4), font = 2)
    text(0.2, ybar_srs, expression(bar(y)), col = "darkgreen", adj = c(0, -0.7), font = 2)
})

## post-stratification: an SRS (green), post-strata A-D of widths N_h, and the
## stratum sample means (red) whose pi_h-weighted average is the estimate
plot_poststrat_concept <- function () with_seed (3268, {
    Nh <- c(50, 20, 30, 40)      # population sizes per post-stratum
    mu_h <- c(1, 4, 7, 10)       # true stratum means, exaggerated and ascending
    labs <- c("A", "B", "C", "D")
    N <- sum(Nh)
    H <- length(Nh)
    n <- 25
    sdy <- 0.9

    # index plot: sequential index across the population, grouped contiguously by
    # stratum -- so each stratum's width on the x-axis is exactly its Nh
    y_pop <- unlist(lapply(1:H, function(h) mu_h[h] + rnorm(Nh[h], 0, sdy)))
    strat_id <- rep(1:H, Nh)
    boundaries <- cumsum(Nh)
    starts <- c(0, boundaries[-H])
    x_pop <- unlist(lapply(1:H, function(h) starts[h] + seq_len(Nh[h])))
    ybar_U <- mean(y_pop)                    # true population mean

    sel <- sample(N, n)                      # an ordinary SRS from the whole population
    nh <- tabulate(strat_id[sel], H)         # post-stratum sample sizes: random, not proportional
    ybar_h <- sapply(1:H, function(h) mean(y_pop[sel][strat_id[sel] == h]))

    ybar_naive <- mean(y_pop[sel])           # pooled sample mean, ignoring strata
    ybar_post <- sum(ybar_h * Nh / N)        # post-stratified estimate

    par(mar = c(3, 4, 0.5, 0.5))
    plot(x_pop, y_pop, pch = 19, col = "grey50", xlab = "", ylab = "y",
         xlim = c(-12, N), ylim = c(-2, 13), xaxt = "n")

    abline(v = boundaries[-H], col = "grey70")           # separate the groups
    mid <- (starts + boundaries) / 2
    axis(1, at = mid, labels = labs, tick = FALSE)

    points(x_pop[sel], y_pop[sel], pch = 21, col = "darkgreen", lwd = 2, cex = 2.2)

    for (h in 1:H)                                       # sample mean of each post-stratum
      segments(starts[h], ybar_h[h], boundaries[h], ybar_h[h], col = "firebrick", lwd = 3)

    abline(h = ybar_U, col = "black", lty = 3)
    abline(h = ybar_post, col = "firebrick", lty = 3)
    abline(h = ybar_naive, col = "darkgreen", lty = 3)

    points(-11, ybar_U, pch = 19, col = "black", cex = 2, xpd = NA)
    points(-11, ybar_post, pch = 19, col = "firebrick", cex = 1.6, xpd = NA)
    points(-11, ybar_naive, pch = 19, col = "darkgreen", cex = 1.6, xpd = NA)

    text(-10, ybar_U, expression(bar(y)[U]), col = "black", adj = c(0, -0.7), font = 2)
    text(-10, ybar_post, expression(hat(bar(y))[post]), col = "firebrick", adj = c(0, 1.4), font = 2)
    text(-10, ybar_naive, expression(bar(y)), col = "darkgreen", adj = c(0, -0.7), font = 2)
})

## data suited to ratio estimation (line through the origin, spread growing with x)
## versus regression estimation (line with an intercept, constant spread)
plot_ratio_vs_reg <- function () with_seed (9, {
    n <- 40
    x <- runif(n, 1, 10)

    par(mfrow = c(1, 2), mar = c(4, 4, 2, 0.5))

    B <- 0.5
    y1 <- B * x + rnorm(n, 0, 0.4 * sqrt(x))
    plot(x, y1, pch = 19, xlab = "x", ylab = "y", xlim = c(0, 10.5), ylim = c(-1, 6),
         main = "Ratio Estimation")
    abline(0, B, col = "firebrick", lwd = 2)
    text(8.5, B * 8.5 + 0.6, expression(y == B*x), col = "firebrick")

    B0 <- 1.5
    B1 <- 0.35
    y2 <- B0 + B1 * x + rnorm(n, 0, 0.6)
    plot(x, y2, pch = 19, xlab = "x", ylab = "y", xlim = c(0, 10.5), ylim = c(-1, 6),
         main = "Regression Estimation")
    abline(B0, B1, col = "firebrick", lwd = 2)
    text(8.5, B0 + B1 * 8.5 + 0.6, expression(y == B[0]+B[1]*x), col = "firebrick")
})

## a prediction estimator with a flexible (loess) curve, against a straight-line fit
plot_ml_concept <- function () with_seed (37, {
    N <- 40
    f <- function(x) 2 + 2.5 * sin(0.8 * x) + 0.35 * x   # true, nonlinear E[y|x]
    x <- runif(N, 0.5, 10)
    y <- f(x) + rnorm(N, 0, 0.3)

    ybar_U <- mean(y)

    n <- 15
    sel <- sample(N, n)                     # an ordinary SRS, spread across x
    df <- data.frame(xx = x[sel], yy = y[sel])

    fit_curve <- loess(yy ~ xx, data = df, span = 0.9,
                        control = loess.control(surface = "direct")) # "ML" model
    fit_lm <- lm(yy ~ xx, data = df)                                  # linear model

    xg <- seq(0.3, 10.2, length.out = 200)
    curve_g <- predict(fit_curve, newdata = data.frame(xx = xg))
    lm_g <- predict(fit_lm, newdata = data.frame(xx = xg))

    yhat <- predict(fit_curve, newdata = data.frame(xx = x)) # fitted value, every unit
    ybar_ml <- mean(yhat)                   # model-based (ML) estimate
    ybar_srs <- mean(y[sel])                # naive sample mean

    par(mar = c(4, 4, 0.5, 0.5))
    plot(x, y, pch = 19, col = "grey50", xlab = "x", ylab = "y",
         xlim = c(0, 10.5), ylim = c(0, 8))

    lines(xg, lm_g, col = "steelblue", lwd = 2, lty = 2)
    lines(xg, curve_g, col = "firebrick", lwd = 2)
    text(9.0, tail(curve_g, 1) + 0.5, "ML curve", col = "firebrick")
    text(9.0, tail(lm_g, 1) - 0.5, "linear fit", col = "steelblue")

    points(x[sel], y[sel], pch = 21, col = "darkgreen", lwd = 2, cex = 2.2)
    points(x, yhat, pch = 19, col = "firebrick", cex = 0.9)

    abline(h = ybar_U, col = "black", lty = 3)
    abline(h = ybar_ml, col = "firebrick", lty = 3)
    abline(h = ybar_srs, col = "darkgreen", lty = 3)

    points(0, ybar_U, pch = 19, col = "black", cex = 2)
    points(0, ybar_ml, pch = 19, col = "firebrick", cex = 1.6)
    points(0, ybar_srs, pch = 19, col = "darkgreen", cex = 1.6)

    text(0.2, ybar_U, expression(bar(y)[U]), col = "black", adj = c(0, -0.7), font = 2)
    text(0.2, ybar_ml, expression(hat(bar(y))[ML]), col = "firebrick", adj = c(0, 1.4), font = 2)
    text(0.2, ybar_srs, expression(bar(y)), col = "darkgreen", adj = c(0, -0.7), font = 2)
})

## the toy France example: n = 5 communes, the ratio line and the residuals
plot_france_toy <- function ()
{
    x <- c(1, 2, 3, 4, 5)
    y <- c(18, 78, 62, 132, 160)
    Bhat <- sum(y) / sum(x)
    yhat <- Bhat * x
    e <- y - yhat

    par(mar = c(4, 4, 0.5, 0.5))
    plot(x, y, pch = 19, cex = 1.4, xlab = "x: registered births (000s)", ylab = "y: residents (000s)",
         xlim = c(0, 6), ylim = c(0, 180))
    abline(0, Bhat, col = "firebrick", lwd = 2)
    segments(x, y, x, yhat, col = "grey40", lty = 2)
    points(x, yhat, pch = 21, col = "firebrick", bg = "white", cex = 1.2)
    text(x + 0.15, (y + yhat) / 2, sprintf("e=%d", e), cex = 0.8, adj = 0)
    text(5.3, Bhat * 5.3, expression(y == hat(B) * x), col = "firebrick", adj = 0)
}
