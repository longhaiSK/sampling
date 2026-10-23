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
