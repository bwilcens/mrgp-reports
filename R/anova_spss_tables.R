# anova_spss_tables.R -- one-way ANOVA of ln(cost) by category, printed as the two SPSS tables: ANOVA
#   (Between Groups / Within Groups / Total) and Multiple Comparisons (Bonferroni, every ordered pair).
# Set STORM to 2023 or 2024 and run from mrgp-roads-core:
#   "C:/Program Files/R/R-4.5.3/bin/Rscript.exe" twopart_simple_2026-09-17/R/anova_spss_tables.R
# Damaged graded segments, Precip >= 3.51 in, Averill excluded; cost = cost_dedup, uncapped. Output is markdown.

STORM <- 2023

d <- read.csv(paste0("twopart_simple_2026-09-17/data/analysis_", STORM, "_2026-09-27.csv"))
k <- d[d$Precip >= 3.51 & d$Town != "Averill" & d$cost_dedup > 0, ]
k$ln_cost <- log(k$cost_dedup)

for (v in c("PARENT", "HYDROGROUP", "vtrans_district")) {
  g <- factor(k[[v]])
  a <- summary(aov(k$ln_cost ~ g))[[1]]          # the ANOVA table R computes

  cat("\n**ANOVA** -- ln(cost) by ", v, ", ", STORM, " storm\n\n", sep = "")
  cat("| | Sum of Squares | df | Mean Square | F | Sig. |\n|---|---:|---:|---:|---:|---:|\n")
  cat("| Between Groups | ", sprintf("%.3f", a$`Sum Sq`[1]), " | ", a$Df[1], " | ", sprintf("%.3f", a$`Mean Sq`[1]),
      " | ", sprintf("%.3f", a$`F value`[1]), " | ", sub("^0", "", sprintf("%.3f", a$`Pr(>F)`[1])), " |\n", sep = "")
  cat("| Within Groups | ", sprintf("%.3f", a$`Sum Sq`[2]), " | ", a$Df[2], " | ", sprintf("%.3f", a$`Mean Sq`[2]), " | | |\n", sep = "")
  cat("| Total | ", sprintf("%.3f", sum(a$`Sum Sq`)), " | ", sum(a$Df), " | | | |\n", sep = "")

  # Bonferroni: p multiplied by the number of pairs; the interval uses the matching critical t
  n_g    <- as.numeric(table(g))
  mean_g <- as.numeric(tapply(k$ln_cost, g, mean))
  mse    <- a$`Mean Sq`[2]
  df_w   <- a$Df[2]
  pairs  <- nlevels(g) * (nlevels(g) - 1) / 2
  t_crit <- qt(1 - 0.05 / (2 * pairs), df_w)

  cat("\n**Multiple Comparisons** -- Dependent Variable: ln(cost); Bonferroni\n\n")
  cat("| (I) ", v, " | (J) ", v, " | Mean Difference (I-J) | Std. Error | Sig. | 95% CI Lower Bound | 95% CI Upper Bound |\n", sep = "")
  cat("|---|---|---:|---:|---:|---:|---:|\n")
  for (i in 1:nlevels(g)) for (j in 1:nlevels(g)) {
    if (i == j) next
    diff <- mean_g[i] - mean_g[j]
    se   <- sqrt(mse * (1 / n_g[i] + 1 / n_g[j]))
    p    <- min(1, 2 * pt(-abs(diff / se), df_w) * pairs)
    cat("| ", if (j == 1 || (i == 1 && j == 2)) levels(g)[i] else "", " | ", levels(g)[j], " | ",
        sprintf("%.3f", diff), if (p < 0.05) "\\*" else "", " | ", sprintf("%.3f", se), " | ",
        sub("^0", "", sprintf("%.3f", p)), " | ", sprintf("%.2f", diff - t_crit * se), " | ", sprintf("%.2f", diff + t_crit * se), " |\n", sep = "")
  }
  cat("\n\\*. The mean difference is significant at the 0.05 level.\n")
}
