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

set.seed(348)
idx <- sort(sample(nrow(pop), 300))
write.csv(pop[idx, ], "agsrs_syn.csv", row.names = FALSE)
