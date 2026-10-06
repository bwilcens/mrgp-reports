# crosstab.R -- crosstab of each category against damage, following Quick-R "Frequencies and Crosstabs"
#   (https://www.datacamp.com/doc/r/frequencies): table, margin.table, prop.table, chisq.test, then
#   gmodels::CrossTable in its SPSS format.
# Set STORM to 2023 or 2024 and run from mrgp-roads-core:
#   "C:/Program Files/R/R-4.5.3/bin/Rscript.exe" twopart_simple_2026-09-17/R/crosstab.R
# Graded segments, Precip >= 3.51 in, Averill excluded; damage = cost_dedup > 0.

STORM <- 2023

library(gmodels)

d <- read.csv(paste0("twopart_simple_2026-09-17/data/analysis_", STORM, "_2026-09-27.csv"))
d <- d[d$Precip >= 3.51 & d$Town != "Averill", ]
d$damaged <- as.integer(d$cost_dedup > 0)

cat("CROSSTABS", STORM, "\n")

mytable <- table(d$PARENT, d$damaged)
print(mytable)
print(margin.table(mytable, 1))
print(margin.table(mytable, 2))
print(round(prop.table(mytable, 1), 4))
print(chisq.test(mytable))
CrossTable(d$PARENT, d$damaged, expected = TRUE, prop.r = TRUE, prop.c = FALSE, prop.t = FALSE,
           prop.chisq = FALSE, chisq = TRUE, asresid = TRUE, format = "SPSS")

mytable <- table(d$HYDROGROUP, d$damaged)
print(mytable)
print(margin.table(mytable, 1))
print(margin.table(mytable, 2))
print(round(prop.table(mytable, 1), 4))
print(chisq.test(mytable))
CrossTable(d$HYDROGROUP, d$damaged, expected = TRUE, prop.r = TRUE, prop.c = FALSE, prop.t = FALSE,
           prop.chisq = FALSE, chisq = TRUE, asresid = TRUE, format = "SPSS")

mytable <- table(d$vtrans_district, d$damaged)
print(mytable)
print(margin.table(mytable, 1))
print(margin.table(mytable, 2))
print(round(prop.table(mytable, 1), 4))
print(chisq.test(mytable))
CrossTable(d$vtrans_district, d$damaged, expected = TRUE, prop.r = TRUE, prop.c = FALSE, prop.t = FALSE,
           prop.chisq = FALSE, chisq = TRUE, asresid = TRUE, format = "SPSS")
