# crosstab.R -- crosstab of each category against damage: the standard R workflow (table, margin.table,
#   prop.table, chisq.test; Quick-R "Frequencies and Crosstabs", https://www.datacamp.com/doc/r/frequencies),
#   printed as the SPSS tables: Crosstabulation (Count), Crosstabulation (Expected Count), Row Percent with
#   Adjusted Residual, and Chi-Square Tests.
# Set STORM to 2023 or 2024 and run from mrgp-roads-core:
#   "C:/Program Files/R/R-4.5.3/bin/Rscript.exe" twopart_simple_2026-09-17/R/crosstab.R
# Graded segments, Precip >= 3.51 in, Averill excluded; damage = cost_dedup > 0. Output is markdown.

STORM <- 2023

d <- read.csv(paste0("twopart_simple_2026-09-17/data/analysis_", STORM, "_2026-09-27.csv"))
d <- d[d$Precip >= 3.51 & d$Town != "Averill", ]
d$damaged <- as.integer(d$cost_dedup > 0)

for (v in c("PARENT", "HYDROGROUP", "vtrans_district")) {
  mytable  <- table(d[[v]], d$damaged)
  rowtotal <- margin.table(mytable, 1)
  coltotal <- margin.table(mytable, 2)
  rowpct   <- prop.table(mytable, 1)
  chi      <- chisq.test(mytable)
  expected <- chi$expected
  LR       <- 2 * sum(mytable[mytable > 0] * log(mytable[mytable > 0] / expected[mytable > 0]))

  cat("\n**", v, " \\* damaged Crosstabulation** -- Count, ", STORM, " storm\n\n", sep = "")
  cat("| ", v, " | 0 = not damaged | 1 = damaged | Total |\n|---|---:|---:|---:|\n", sep = "")
  for (r in rownames(mytable)) cat("| ", r, " | ", mytable[r, "0"], " | ", mytable[r, "1"], " | ", rowtotal[r], " |\n", sep = "")
  cat("| Total | ", coltotal["0"], " | ", coltotal["1"], " | ", sum(mytable), " |\n", sep = "")

  cat("\n**", v, " \\* damaged Crosstabulation** -- Expected Count\n\n", sep = "")
  cat("| ", v, " | 0 = not damaged | 1 = damaged | Total |\n|---|---:|---:|---:|\n", sep = "")
  for (r in rownames(mytable)) cat("| ", r, " | ", sprintf("%.1f", expected[r, "0"]), " | ", sprintf("%.1f", expected[r, "1"]), " | ", sprintf("%.1f", rowtotal[r]), " |\n", sep = "")
  cat("| Total | ", sprintf("%.1f", coltotal["0"]), " | ", sprintf("%.1f", coltotal["1"]), " | ", sprintf("%.1f", sum(mytable)), " |\n", sep = "")

  cat("\n**", v, " \\* damaged Crosstabulation** -- Row Percent and Adjusted Residual\n\n", sep = "")
  cat("| ", v, " | % damaged | Adjusted Residual (damaged) |\n|---|---:|---:|\n", sep = "")
  for (r in rownames(mytable)) cat("| ", r, " | ", sprintf("%.1f%%", 100 * rowpct[r, "1"]), " | ", sprintf("%.1f", chi$stdres[r, "1"]), " |\n", sep = "")

  cat("\n**Chi-Square Tests**\n\n")
  cat("| | Value | df | Asymptotic Significance (2-sided) |\n|---|---:|---:|---:|\n")
  cat("| Pearson Chi-Square | ", sprintf("%.3f", chi$statistic), "<sup>a</sup> | ", chi$parameter, " | ",
      sub("^0", "", sprintf("%.3f", chi$p.value)), " |\n", sep = "")
  cat("| Likelihood Ratio | ", sprintf("%.3f", LR), " | ", chi$parameter, " | ",
      sub("^0", "", sprintf("%.3f", pchisq(LR, chi$parameter, lower.tail = FALSE))), " |\n", sep = "")
  cat("| N of Valid Cases | ", sum(mytable), " | | |\n", sep = "")
  cat("\na. ", sum(expected < 5), " cells (", sprintf("%.1f%%", 100 * mean(expected < 5)), ") have expected count less than 5. The minimum expected count is ",
      sprintf("%.2f", min(expected)), ".\n", sep = "")
}
