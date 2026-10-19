## Draws a genuine simple random sample of n = 300 counties, without replacement,
## from agpop.csv and saves it as agsrs_syn.csv (same columns as agsrs.csv).
##
## The population is the counties with no missing (-99) value in acres92,
## acres87, acres82 or farms92 (the same population as the post-stratification
## app), so that every sampled county has all of the variables.
pop <- read.csv("agpop.csv")
vars <- c("acres92", "acres87", "acres82", "farms92")
for (v in vars) pop[[v]][pop[[v]] == -99] <- NA
pop <- pop[complete.cases(pop[, vars]), ]

## The seed was chosen for a teaching example: it is the first seed in 1, 2, 3, ...
## for which both the SRS sample mean of acres92 and the mean post-stratified by
## acres87 quantile groups (cuts at the 25%, 75%, 95% population percentiles)
## are within 5 thousand acres of the population mean. An arbitrary SRS can be
## much less typical, because acres92 is heavily skewed.
set.seed(38)
idx <- sort(sample(nrow(pop), 300))
write.csv(pop[idx, ], "agsrs_syn.csv", row.names = FALSE)
