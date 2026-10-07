# Crosstab and ANOVA on Categoricals

The approach mirrors Dr. Wemple's homework assignments on crosstab and ANOVA, with the tables laid out as SPSS prints them.

## Script: crosstab_anova.R

```r
# Data:  ../mrgp-roads-core/twopart_simple_2026-09-17/data/analysis_<STORM>_2026-09-27.csv
#        one row per road segment; compliance-graded segments with at least 3.51 in of rain,
#        Averill excluded; damaged = 1 if the de-duplicated repair cost (cost_dedup) is above zero
# Hydrogroup, Parent, VTRANS District
# run from projects/reports: "C:/Program Files/R/R-4.5.3/bin/Rscript.exe" R/crosstab_anova.R

STORM <- 2023                                                   # which storm to run: 2023 or 2024

library(gmodels)                                                # CrossTable: the crosstab laid out the way SPSS prints it
library(car)                                                    # leveneTest: the equal-variances test SPSS prints with a one-way anova


# step 0  read the data
roads <- read.csv(paste0("../mrgp-roads-core/twopart_simple_2026-09-17/data/analysis_", STORM, "_2026-09-27.csv"))   # the storm's file
roads <- roads[roads$Precip >= 3.51 & roads$Town != "Averill", ]   # keep segments at or above the rain cutoff; drop Averill
roads$damaged <- ifelse(roads$cost_dedup > 0, 1, 0)             # damaged = 1 if any de-duplicated cost was recorded, else 0

roads$PARENT          <- factor(roads$PARENT)                   # the three categoricals as factors (district arrives as a number in the file)
roads$HYDROGROUP      <- factor(roads$HYDROGROUP)
roads$vtrans_district <- factor(roads$vtrans_district)

cat("storm", STORM, "\n")
cat("segments", nrow(roads), "\n")
cat("damaged", sum(roads$damaged), "\n")
cat("overall percent damaged", round(100 * mean(roads$damaged), 2), "\n")


# step 1  crosstabs: each cell shows count, expected count, row percent, residual (observed - expected)
cat("\nstep 1 parent material by damage\n")
CrossTable(roads$PARENT, roads$damaged, expected = TRUE, prop.r = TRUE, prop.c = FALSE, prop.t = FALSE,
           prop.chisq = FALSE, resid = TRUE, sresid = FALSE, asresid = FALSE, format = "SPSS")   # crosstab + pearson chi-square test

cat("\nstep 2 hydrologic group by damage\n")
CrossTable(roads$HYDROGROUP, roads$damaged, expected = TRUE, prop.r = TRUE, prop.c = FALSE, prop.t = FALSE,
           prop.chisq = FALSE, resid = TRUE, sresid = FALSE, asresid = FALSE, format = "SPSS")

cat("\nstep 3 vtrans district by damage\n")
CrossTable(roads$vtrans_district, roads$damaged, expected = TRUE, prop.r = TRUE, prop.c = FALSE, prop.t = FALSE,
           prop.chisq = FALSE, resid = TRUE, sresid = FALSE, asresid = FALSE, format = "SPSS")


# step 4  the cost side: damaged segments only, cost on the natural-log scale
dmg <- roads[roads$damaged == 1, ]                              # the damaged segments
dmg$log_cost <- log(dmg$cost_dedup)                             # ln(cost): the dollar amounts are too skewed for an anova as they are
dmg$cost_k   <- dmg$cost_dedup / 1000                           # cost in thousands of dollars, for the box plots
dmg$PARENT          <- factor(dmg$PARENT)                       # re-make the factors so categories with no damaged segments drop out
dmg$HYDROGROUP      <- factor(dmg$HYDROGROUP)
dmg$vtrans_district <- factor(dmg$vtrans_district)
cat("\nstep 4 damaged segments for the anova\n")
cat("damaged segments", nrow(dmg), "\n")
cat("mean ln(cost), all damaged", round(mean(dmg$log_cost), 3), "\n")


# step 5  one-way anova of ln(cost) by parent material, in SPSS's order: descriptives, levene, anova, bonferroni
cat("\nstep 5 ln(cost) by parent material\n")
groups <- levels(dmg$PARENT)                                    # the categories with damaged segments
desc <- data.frame(group = groups,
                   N             = as.integer(table(dmg$PARENT)),                              # damaged segments in each group
                   Mean          = aggregate(log_cost ~ PARENT, data = dmg, mean)$log_cost,    # mean ln(cost) in each group
                   Std_Deviation = aggregate(log_cost ~ PARENT, data = dmg, sd)$log_cost,
                   Minimum       = aggregate(log_cost ~ PARENT, data = dmg, min)$log_cost,
                   Maximum       = aggregate(log_cost ~ PARENT, data = dmg, max)$log_cost)
desc$Std_Error <- desc$Std_Deviation / sqrt(desc$N)              # standard error of each group mean
desc$Lower_95  <- desc$Mean - qt(0.975, pmax(desc$N - 1, 1)) * desc$Std_Error   # (a group of one segment has no sd, so its interval comes out NA)   # 95% confidence interval for the mean, lower
desc$Upper_95  <- desc$Mean + qt(0.975, pmax(desc$N - 1, 1)) * desc$Std_Error   # and upper
total <- data.frame(group = "Total", N = nrow(dmg), Mean = mean(dmg$log_cost), Std_Deviation = sd(dmg$log_cost),   # the all-groups row
                    Minimum = min(dmg$log_cost), Maximum = max(dmg$log_cost))
total$Std_Error <- total$Std_Deviation / sqrt(total$N)
total$Lower_95  <- total$Mean - qt(0.975, total$N - 1) * total$Std_Error
total$Upper_95  <- total$Mean + qt(0.975, total$N - 1) * total$Std_Error
shown <- rbind(desc, total)[, c("group", "N", "Mean", "Std_Deviation", "Std_Error", "Lower_95", "Upper_95", "Minimum", "Maximum")]   # SPSS column order
shown[, -1] <- round(shown[, -1], 3)
cat("\ndescriptives\n"); print(shown, row.names = FALSE)

cat("\ntest of homogeneity of variances\n"); print(leveneTest(log_cost ~ PARENT, data = dmg, center = mean))   # the anova assumes equal spread in every group

a <- anova(lm(log_cost ~ PARENT, data = dmg))                   # R's anova: row 1 the factor (between groups), row 2 residuals (within groups)
anova_table <- data.frame(Source         = c("Between Groups", "Within Groups", "Total"),
                          Sum_of_Squares = format(round(c(a[1, "Sum Sq"], a[2, "Sum Sq"], a[1, "Sum Sq"] + a[2, "Sum Sq"]), 3), nsmall = 3),
                          df             = c(a[1, "Df"], a[2, "Df"], a[1, "Df"] + a[2, "Df"]),
                          Mean_Square    = c(format(round(c(a[1, "Mean Sq"], a[2, "Mean Sq"]), 3), nsmall = 3), ""),   # blank where SPSS leaves a blank
                          F              = c(format(round(a[1, "F value"], 3), nsmall = 3), "", ""),
                          Sig            = c(format(round(a[1, "Pr(>F)"], 3), nsmall = 3), "", ""))
cat("\nanova\n"); print(anova_table, row.names = FALSE)

mse       <- a[2, "Mean Sq"]                                    # within-groups mean square: the pooled variance the comparisons use
df_within <- a[2, "Df"]                                         # its degrees of freedom
n_pairs   <- length(groups) * (length(groups) - 1) / 2          # number of pairs: the bonferroni multiplier
t_crit    <- qt(1 - 0.025 / n_pairs, df_within)                 # t value for the bonferroni-adjusted 95% interval
comparisons <- data.frame()                                     # one row per pair, filled in the loop (SPSS lists each pair twice; once here)
for (i in 1:(length(groups) - 1)) {
  for (j in (i + 1):length(groups)) {
    diff <- desc$Mean[i] - desc$Mean[j]                         # mean difference (I - J)
    se   <- sqrt(mse * (1 / desc$N[i] + 1 / desc$N[j]))         # its standard error from the pooled variance
    p    <- min(1, n_pairs * 2 * pt(-abs(diff / se), df_within))   # two-sided p times the number of pairs, capped at 1
    comparisons <- rbind(comparisons, data.frame(I = groups[i], J = groups[j], Mean_Difference = round(diff, 3), Std_Error = round(se, 3),
                                                 Sig = round(p, 3), Lower_95 = round(diff - t_crit * se, 3), Upper_95 = round(diff + t_crit * se, 3),
                                                 flag = ifelse(p < 0.05, "*", "")))   # * = the mean difference is significant at the .05 level
  }
}
cat("\nmultiple comparisons, bonferroni\n"); print(comparisons, row.names = FALSE)

png(paste0("figures/crosstab_anova_", STORM, "_fig1.png"), width = 900, height = 600)   # figure 1 to a file
boxplot(cost_k ~ PARENT, data = dmg, xlab = "Parent material", ylab = "Repair cost ($ thousands)",
        main = paste("Repair cost by parent material,", STORM, "storm"),
        names = paste0(groups, "\nn=", desc$N), cex.axis = 0.8)   # box = middle half of the values, line = median, whiskers = the rest, circles = outliers
invisible(dev.off())                                            # close the file without printing "null device"


# step 6  one-way anova of ln(cost) by hydrologic group
cat("\nstep 6 ln(cost) by hydrologic group\n")
groups <- levels(dmg$HYDROGROUP)
desc <- data.frame(group = groups,
                   N             = as.integer(table(dmg$HYDROGROUP)),
                   Mean          = aggregate(log_cost ~ HYDROGROUP, data = dmg, mean)$log_cost,
                   Std_Deviation = aggregate(log_cost ~ HYDROGROUP, data = dmg, sd)$log_cost,
                   Minimum       = aggregate(log_cost ~ HYDROGROUP, data = dmg, min)$log_cost,
                   Maximum       = aggregate(log_cost ~ HYDROGROUP, data = dmg, max)$log_cost)
desc$Std_Error <- desc$Std_Deviation / sqrt(desc$N)
desc$Lower_95  <- desc$Mean - qt(0.975, pmax(desc$N - 1, 1)) * desc$Std_Error   # (a group of one segment has no sd, so its interval comes out NA)
desc$Upper_95  <- desc$Mean + qt(0.975, pmax(desc$N - 1, 1)) * desc$Std_Error
shown <- rbind(desc, total)[, c("group", "N", "Mean", "Std_Deviation", "Std_Error", "Lower_95", "Upper_95", "Minimum", "Maximum")]   # the Total row is the same for every grouping
shown[, -1] <- round(shown[, -1], 3)
cat("\ndescriptives\n"); print(shown, row.names = FALSE)

cat("\ntest of homogeneity of variances\n"); print(leveneTest(log_cost ~ HYDROGROUP, data = dmg, center = mean))

a <- anova(lm(log_cost ~ HYDROGROUP, data = dmg))
anova_table <- data.frame(Source         = c("Between Groups", "Within Groups", "Total"),
                          Sum_of_Squares = format(round(c(a[1, "Sum Sq"], a[2, "Sum Sq"], a[1, "Sum Sq"] + a[2, "Sum Sq"]), 3), nsmall = 3),
                          df             = c(a[1, "Df"], a[2, "Df"], a[1, "Df"] + a[2, "Df"]),
                          Mean_Square    = c(format(round(c(a[1, "Mean Sq"], a[2, "Mean Sq"]), 3), nsmall = 3), ""),   # blank where SPSS leaves a blank
                          F              = c(format(round(a[1, "F value"], 3), nsmall = 3), "", ""),
                          Sig            = c(format(round(a[1, "Pr(>F)"], 3), nsmall = 3), "", ""))
cat("\nanova\n"); print(anova_table, row.names = FALSE)

mse       <- a[2, "Mean Sq"]
df_within <- a[2, "Df"]
n_pairs   <- length(groups) * (length(groups) - 1) / 2
t_crit    <- qt(1 - 0.025 / n_pairs, df_within)
comparisons <- data.frame()
for (i in 1:(length(groups) - 1)) {
  for (j in (i + 1):length(groups)) {
    diff <- desc$Mean[i] - desc$Mean[j]
    se   <- sqrt(mse * (1 / desc$N[i] + 1 / desc$N[j]))
    p    <- min(1, n_pairs * 2 * pt(-abs(diff / se), df_within))
    comparisons <- rbind(comparisons, data.frame(I = groups[i], J = groups[j], Mean_Difference = round(diff, 3), Std_Error = round(se, 3),
                                                 Sig = round(p, 3), Lower_95 = round(diff - t_crit * se, 3), Upper_95 = round(diff + t_crit * se, 3),
                                                 flag = ifelse(p < 0.05, "*", "")))
  }
}
cat("\nmultiple comparisons, bonferroni\n"); print(comparisons, row.names = FALSE)

png(paste0("figures/crosstab_anova_", STORM, "_fig2.png"), width = 900, height = 600)   # figure 2 to a file
boxplot(cost_k ~ HYDROGROUP, data = dmg, xlab = "Hydrologic group", ylab = "Repair cost ($ thousands)",
        main = paste("Repair cost by hydrologic group,", STORM, "storm"),
        names = paste0(groups, "\nn=", desc$N), cex.axis = 0.8)
invisible(dev.off())                                            # close the file without printing "null device"


# step 7  one-way anova of ln(cost) by vtrans district
cat("\nstep 7 ln(cost) by vtrans district\n")
groups <- levels(dmg$vtrans_district)
desc <- data.frame(group = groups,
                   N             = as.integer(table(dmg$vtrans_district)),
                   Mean          = aggregate(log_cost ~ vtrans_district, data = dmg, mean)$log_cost,
                   Std_Deviation = aggregate(log_cost ~ vtrans_district, data = dmg, sd)$log_cost,
                   Minimum       = aggregate(log_cost ~ vtrans_district, data = dmg, min)$log_cost,
                   Maximum       = aggregate(log_cost ~ vtrans_district, data = dmg, max)$log_cost)
desc$Std_Error <- desc$Std_Deviation / sqrt(desc$N)
desc$Lower_95  <- desc$Mean - qt(0.975, pmax(desc$N - 1, 1)) * desc$Std_Error   # (a group of one segment has no sd, so its interval comes out NA)
desc$Upper_95  <- desc$Mean + qt(0.975, pmax(desc$N - 1, 1)) * desc$Std_Error
shown <- rbind(desc, total)[, c("group", "N", "Mean", "Std_Deviation", "Std_Error", "Lower_95", "Upper_95", "Minimum", "Maximum")]
shown[, -1] <- round(shown[, -1], 3)
cat("\ndescriptives\n"); print(shown, row.names = FALSE)

cat("\ntest of homogeneity of variances\n"); print(leveneTest(log_cost ~ vtrans_district, data = dmg, center = mean))

a <- anova(lm(log_cost ~ vtrans_district, data = dmg))
anova_table <- data.frame(Source         = c("Between Groups", "Within Groups", "Total"),
                          Sum_of_Squares = format(round(c(a[1, "Sum Sq"], a[2, "Sum Sq"], a[1, "Sum Sq"] + a[2, "Sum Sq"]), 3), nsmall = 3),
                          df             = c(a[1, "Df"], a[2, "Df"], a[1, "Df"] + a[2, "Df"]),
                          Mean_Square    = c(format(round(c(a[1, "Mean Sq"], a[2, "Mean Sq"]), 3), nsmall = 3), ""),   # blank where SPSS leaves a blank
                          F              = c(format(round(a[1, "F value"], 3), nsmall = 3), "", ""),
                          Sig            = c(format(round(a[1, "Pr(>F)"], 3), nsmall = 3), "", ""))
cat("\nanova\n"); print(anova_table, row.names = FALSE)

mse       <- a[2, "Mean Sq"]
df_within <- a[2, "Df"]
n_pairs   <- length(groups) * (length(groups) - 1) / 2
t_crit    <- qt(1 - 0.025 / n_pairs, df_within)
comparisons <- data.frame()
for (i in 1:(length(groups) - 1)) {
  for (j in (i + 1):length(groups)) {
    diff <- desc$Mean[i] - desc$Mean[j]
    se   <- sqrt(mse * (1 / desc$N[i] + 1 / desc$N[j]))
    p    <- min(1, n_pairs * 2 * pt(-abs(diff / se), df_within))
    comparisons <- rbind(comparisons, data.frame(I = groups[i], J = groups[j], Mean_Difference = round(diff, 3), Std_Error = round(se, 3),
                                                 Sig = round(p, 3), Lower_95 = round(diff - t_crit * se, 3), Upper_95 = round(diff + t_crit * se, 3),
                                                 flag = ifelse(p < 0.05, "*", "")))
  }
}
cat("\nmultiple comparisons, bonferroni\n"); print(comparisons, row.names = FALSE)

png(paste0("figures/crosstab_anova_", STORM, "_fig3.png"), width = 900, height = 600)   # figure 3 to a file
boxplot(cost_k ~ vtrans_district, data = dmg, xlab = "VTrans district", ylab = "Repair cost ($ thousands)",
        main = paste("Repair cost by VTrans district,", STORM, "storm"),
        names = paste0(groups, "\nn=", desc$N), cex.axis = 0.8)
invisible(dev.off())                                            # close the file without printing "null device"
```

## Output 2023

```
storm 2023 
segments 59695 
damaged 1457 
overall percent damaged 2.44 

step 1 parent material by damage

   Cell Contents
|-------------------------|
|                   Count |
|         Expected Values |
|             Row Percent |
|                Residual |
|-------------------------|

Total Observations in Table:  59695 

             | roads$damaged 
roads$PARENT |        0  |        1  | Row Total | 
-------------|-----------|-----------|-----------|
           A |     2389  |      115  |     2504  | 
             | 2442.884  |   61.116  |           | 
             |   95.407% |    4.593% |    4.195% | 
             |  -53.884  |   53.884  |           | 
-------------|-----------|-----------|-----------|
          DT |    23838  |      682  |    24520  | 
             | 23921.530  |  598.470  |           | 
             |   97.219% |    2.781% |   41.075% | 
             |  -83.530  |   83.530  |           | 
-------------|-----------|-----------|-----------|
       DT/GT |      169  |        8  |      177  | 
             |  172.680  |    4.320  |           | 
             |   95.480% |    4.520% |    0.297% | 
             |   -3.680  |    3.680  |           | 
-------------|-----------|-----------|-----------|
          GF |     9947  |      213  |    10160  | 
             | 9912.021  |  247.979  |           | 
             |   97.904% |    2.096% |   17.020% | 
             |   34.979  |  -34.979  |           | 
-------------|-----------|-----------|-----------|
          GL |     2456  |       43  |     2499  | 
             | 2438.006  |   60.994  |           | 
             |   98.279% |    1.721% |    4.186% | 
             |   17.994  |  -17.994  |           | 
-------------|-----------|-----------|-----------|
          GT |    17767  |      386  |    18153  | 
             | 17709.932  |  443.068  |           | 
             |   97.874% |    2.126% |   30.410% | 
             |   57.068  |  -57.068  |           | 
-------------|-----------|-----------|-----------|
       GT/DT |     1509  |        9  |     1518  | 
             | 1480.950  |   37.050  |           | 
             |   99.407% |    0.593% |    2.543% | 
             |   28.050  |  -28.050  |           | 
-------------|-----------|-----------|-----------|
           M |      163  |        1  |      164  | 
             |  159.997  |    4.003  |           | 
             |   99.390% |    0.610% |    0.275% | 
             |    3.003  |   -3.003  |           | 
-------------|-----------|-----------|-----------|
Column Total |    58238  |     1457  |    59695  | 
-------------|-----------|-----------|-----------|

 
Statistics for All Table Factors


Pearson's Chi-squared test 
------------------------------------------------------------
Chi^2 =  105.9693     d.f. =  7     p =  6.287026e-20 


 
       Minimum expected frequency: 4.002814 
Cells with Expected Frequency < 5: 2 of 16 (12.5%)

Warning message:
In chisq.test(t, correct = FALSE, ...) :
  Chi-squared approximation may be incorrect

step 2 hydrologic group by damage

   Cell Contents
|-------------------------|
|                   Count |
|         Expected Values |
|             Row Percent |
|                Residual |
|-------------------------|

Total Observations in Table:  59695 

                 | roads$damaged 
roads$HYDROGROUP |        0  |        1  | Row Total | 
-----------------|-----------|-----------|-----------|
               A |     9350  |      206  |     9556  | 
                 | 9322.763  |  233.237  |           | 
                 |   97.844% |    2.156% |   16.008% | 
                 |   27.237  |  -27.237  |           | 
-----------------|-----------|-----------|-----------|
             A/D |      745  |        8  |      753  | 
                 |  734.621  |   18.379  |           | 
                 |   98.938% |    1.062% |    1.261% | 
                 |   10.379  |  -10.379  |           | 
-----------------|-----------|-----------|-----------|
               B |     8310  |      254  |     8564  | 
                 | 8354.975  |  209.025  |           | 
                 |   97.034% |    2.966% |   14.346% | 
                 |  -44.975  |   44.975  |           | 
-----------------|-----------|-----------|-----------|
             B/D |     1976  |       75  |     2051  | 
                 | 2000.940  |   50.060  |           | 
                 |   96.343% |    3.657% |    3.436% | 
                 |  -24.940  |   24.940  |           | 
-----------------|-----------|-----------|-----------|
               C |    14074  |      341  |    14415  | 
                 | 14063.167  |  351.833  |           | 
                 |   97.634% |    2.366% |   24.148% | 
                 |   10.833  |  -10.833  |           | 
-----------------|-----------|-----------|-----------|
             C/D |     9331  |      193  |     9524  | 
                 | 9291.544  |  232.456  |           | 
                 |   97.974% |    2.026% |   15.954% | 
                 |   39.456  |  -39.456  |           | 
-----------------|-----------|-----------|-----------|
               D |    14173  |      370  |    14543  | 
                 | 14188.043  |  354.957  |           | 
                 |   97.456% |    2.544% |   24.362% | 
                 |  -15.043  |   15.043  |           | 
-----------------|-----------|-----------|-----------|
       not rated |      279  |       10  |      289  | 
                 |  281.946  |    7.054  |           | 
                 |   96.540% |    3.460% |    0.484% | 
                 |   -2.946  |    2.946  |           | 
-----------------|-----------|-----------|-----------|
    Column Total |    58238  |     1457  |    59695  | 
-----------------|-----------|-----------|-----------|

 
Statistics for All Table Factors


Pearson's Chi-squared test 
------------------------------------------------------------
Chi^2 =  41.04512     d.f. =  7     p =  7.936052e-07 


 
       Minimum expected frequency: 7.05374 


step 3 vtrans district by damage

   Cell Contents
|-------------------------|
|                   Count |
|         Expected Values |
|             Row Percent |
|                Residual |
|-------------------------|

Total Observations in Table:  59695 

                      | roads$damaged 
roads$vtrans_district |        0  |        1  | Row Total | 
----------------------|-----------|-----------|-----------|
                    1 |     6046  |       33  |     6079  | 
                      | 5930.627  |  148.373  |           | 
                      |   99.457% |    0.543% |   10.183% | 
                      |  115.373  | -115.373  |           | 
----------------------|-----------|-----------|-----------|
                    2 |     6069  |      520  |     6589  | 
                      | 6428.180  |  160.820  |           | 
                      |   92.108% |    7.892% |   11.038% | 
                      | -359.180  |  359.180  |           | 
----------------------|-----------|-----------|-----------|
                    3 |    10360  |      147  |    10507  | 
                      | 10250.551  |  256.449  |           | 
                      |   98.601% |    1.399% |   17.601% | 
                      |  109.449  | -109.449  |           | 
----------------------|-----------|-----------|-----------|
                    4 |    10425  |      137  |    10562  | 
                      | 10304.209  |  257.791  |           | 
                      |   98.703% |    1.297% |   17.693% | 
                      |  120.791  | -120.791  |           | 
----------------------|-----------|-----------|-----------|
                    5 |     2645  |       35  |     2680  | 
                      | 2614.588  |   65.412  |           | 
                      |   98.694% |    1.306% |    4.489% | 
                      |   30.412  |  -30.412  |           | 
----------------------|-----------|-----------|-----------|
                    6 |     9650  |      505  |    10155  | 
                      | 9907.143  |  247.857  |           | 
                      |   95.027% |    4.973% |   17.011% | 
                      | -257.143  |  257.143  |           | 
----------------------|-----------|-----------|-----------|
                    7 |     3617  |       41  |     3658  | 
                      | 3568.718  |   89.282  |           | 
                      |   98.879% |    1.121% |    6.128% | 
                      |   48.282  |  -48.282  |           | 
----------------------|-----------|-----------|-----------|
                    8 |     4196  |       25  |     4221  | 
                      | 4117.976  |  103.024  |           | 
                      |   99.408% |    0.592% |    7.071% | 
                      |   78.024  |  -78.024  |           | 
----------------------|-----------|-----------|-----------|
                    9 |     5230  |       14  |     5244  | 
                      | 5116.008  |  127.992  |           | 
                      |   99.733% |    0.267% |    8.785% | 
                      |  113.992  | -113.992  |           | 
----------------------|-----------|-----------|-----------|
         Column Total |    58238  |     1457  |    59695  | 
----------------------|-----------|-----------|-----------|

 
Statistics for All Table Factors


Pearson's Chi-squared test 
------------------------------------------------------------
Chi^2 =  1499.458     d.f. =  8     p =  1.758147e-318 


 
       Minimum expected frequency: 65.41184 


step 4 damaged segments for the anova
damaged segments 1457 
mean ln(cost), all damaged 8.902 

step 5 ln(cost) by parent material

descriptives
 group    N   Mean Std_Deviation Std_Error Lower_95 Upper_95 Minimum Maximum
     A  115  9.251         1.934     0.180    8.894    9.608   3.332  14.243
    DT  682  8.738         1.891     0.072    8.595    8.880   3.466  14.404
 DT/GT    8  8.905         2.119     0.749    7.134   10.677   7.249  13.688
    GF  213  9.305         1.817     0.125    9.060    9.551   4.844  14.314
    GL   43  9.803         1.685     0.257    9.285   10.322   7.238  14.567
    GT  386  8.736         1.794     0.091    8.556    8.915   3.871  13.591
 GT/DT    9  9.821         1.841     0.614    8.406   11.236   6.410  11.885
     M    1 11.521            NA        NA       NA       NA  11.521  11.521
 Total 1457  8.902         1.872     0.049    8.805    8.998   3.332  14.567

test of homogeneity of variances
Levene's Test for Homogeneity of Variance (center = mean)
        Df F value Pr(>F)
group    7   0.529 0.8131
      1449               

anova
         Source Sum_of_Squares   df Mean_Square     F   Sig
 Between Groups        127.123    7      18.160 5.287 0.000
  Within Groups       4977.000 1449       3.435            
          Total       5104.123 1456                        

multiple comparisons, bonferroni
     I     J Mean_Difference Std_Error   Sig Lower_95 Upper_95 flag
     A    DT           0.513     0.187 0.170   -0.071    1.098     
     A DT/GT           0.346     0.678 1.000   -1.775    2.467     
     A    GF          -0.054     0.214 1.000   -0.726    0.617     
     A    GL          -0.552     0.331 1.000   -1.589    0.485     
     A    GT           0.515     0.197 0.252   -0.101    1.131     
     A GT/DT          -0.570     0.641 1.000   -2.578    1.437     
     A     M          -2.270     1.861 1.000   -8.095    3.555     
    DT DT/GT          -0.168     0.659 1.000   -2.230    1.895     
    DT    GF          -0.568     0.145 0.003   -1.023   -0.113    *
    DT    GL          -1.066     0.291 0.007   -1.978   -0.154    *
    DT    GT           0.002     0.118 1.000   -0.368    0.371     
    DT GT/DT          -1.084     0.622 1.000   -3.030    0.862     
    DT     M          -2.783     1.855 1.000   -8.587    3.021     
 DT/GT    GF          -0.400     0.667 1.000   -2.489    1.689     
 DT/GT    GL          -0.898     0.714 1.000   -3.131    1.335     
 DT/GT    GT           0.169     0.662 1.000   -1.903    2.241     
 DT/GT GT/DT          -0.916     0.901 1.000   -3.734    1.902     
 DT/GT     M          -2.615     1.966 1.000   -8.767    3.536     
    GF    GL          -0.498     0.310 1.000   -1.467    0.472     
    GF    GT           0.569     0.158 0.009    0.074    1.064    *
    GF GT/DT          -0.516     0.631 1.000   -2.490    1.458     
    GF     M          -2.215     1.858 1.000   -8.029    3.598     
    GL    GT           1.067     0.298 0.010    0.135    2.000    *
    GL GT/DT          -0.018     0.679 1.000   -2.144    2.108     
    GL     M          -1.718     1.875 1.000   -7.585    4.150     
    GT GT/DT          -1.085     0.625 1.000   -3.041    0.870     
    GT     M          -2.785     1.856 1.000   -8.592    3.023     
 GT/DT     M          -1.699     1.954 1.000   -7.813    4.414     

step 6 ln(cost) by hydrologic group

descriptives
     group    N  Mean Std_Deviation Std_Error Lower_95 Upper_95 Minimum Maximum
         A  206 9.325         1.757     0.122    9.084    9.567   4.605  14.243
       A/D    8 9.415         2.146     0.759    7.621   11.209   5.247  12.551
         B  254 8.641         1.906     0.120    8.406    8.877   4.554  14.314
       B/D   75 9.206         1.883     0.217    8.773    9.639   3.332  14.129
         C  341 8.729         1.901     0.103    8.526    8.931   3.871  13.591
       C/D  193 8.790         1.915     0.138    8.518    9.062   4.143  14.567
         D  370 8.964         1.803     0.094    8.780    9.148   3.466  14.404
 not rated   10 9.824         2.098     0.663    8.323   11.325   5.565  11.948
     Total 1457 8.902         1.872     0.049    8.805    8.998   3.332  14.567

test of homogeneity of variances
Levene's Test for Homogeneity of Variance (center = mean)
        Df F value Pr(>F)
group    7  0.8333 0.5595
      1449               

anova
         Source Sum_of_Squares   df Mean_Square     F   Sig
 Between Groups         85.782    7      12.255 3.538 0.001
  Within Groups       5018.341 1449       3.463            
          Total       5104.123 1456                        

multiple comparisons, bonferroni
   I         J Mean_Difference Std_Error   Sig Lower_95 Upper_95 flag
   A       A/D          -0.090     0.671 1.000   -2.189    2.009     
   A         B           0.684     0.174 0.003    0.138    1.230    *
   A       B/D           0.119     0.251 1.000   -0.666    0.905     
   A         C           0.596     0.164 0.008    0.082    1.110    *
   A       C/D           0.535     0.186 0.116   -0.048    1.119     
   A         D           0.361     0.162 0.719   -0.145    0.868     
   A not rated          -0.499     0.603 1.000   -2.385    1.387     
 A/D         B           0.774     0.668 1.000   -1.317    2.865     
 A/D       B/D           0.209     0.692 1.000   -1.957    2.375     
 A/D         C           0.686     0.666 1.000   -1.397    2.769     
 A/D       C/D           0.625     0.671 1.000   -1.476    2.726     
 A/D         D           0.451     0.665 1.000   -1.630    2.532     
 A/D not rated          -0.409     0.883 1.000   -3.172    2.354     
   B       B/D          -0.565     0.245 0.590   -1.330    0.201     
   B         C          -0.088     0.154 1.000   -0.570    0.395     
   B       C/D          -0.149     0.178 1.000   -0.705    0.407     
   B         D          -0.323     0.152 0.936   -0.797    0.152     
   B not rated          -1.183     0.600 1.000   -3.061    0.695     
 B/D         C           0.477     0.237 1.000   -0.266    1.220     
 B/D       C/D           0.416     0.253 1.000   -0.377    1.208     
 B/D         D           0.242     0.236 1.000   -0.496    0.979     
 B/D not rated          -0.618     0.627 1.000   -2.579    1.343     
   C       C/D          -0.061     0.168 1.000   -0.586    0.464     
   C         D          -0.235     0.140 1.000   -0.672    0.202     
   C not rated          -1.095     0.597 1.000   -2.964    0.773     
 C/D         D          -0.174     0.165 1.000   -0.691    0.343     
 C/D not rated          -1.034     0.604 1.000   -2.923    0.855     
   D not rated          -0.860     0.596 1.000   -2.727    1.006     

step 7 ln(cost) by vtrans district

descriptives
 group    N   Mean Std_Deviation Std_Error Lower_95 Upper_95 Minimum Maximum
     1   33  9.184         1.265     0.220    8.735    9.633   6.807  12.082
     2  520  7.921         1.762     0.077    7.770    8.073   3.871  13.112
     3  147 10.700         1.091     0.090   10.522   10.878   7.129  13.372
     4  137  8.798         1.867     0.159    8.482    9.113   4.605  13.528
     5   35  9.167         1.706     0.288    8.582    9.753   4.344  11.778
     6  505  9.315         1.651     0.073    9.171    9.460   4.001  14.567
     7   41  9.540         1.559     0.244    9.048   10.032   3.466  12.785
     8   25  8.909         1.767     0.353    8.179    9.638   6.428  12.605
     9   14  9.313         2.177     0.582    8.056   10.569   3.332  11.589
 Total 1457  8.902         1.872     0.049    8.805    8.998   3.332  14.567

test of homogeneity of variances
Levene's Test for Homogeneity of Variance (center = mean)
        Df F value    Pr(>F)    
group    8  6.8069 8.376e-09 ***
      1448                      
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

anova
         Source Sum_of_Squares   df Mean_Square      F   Sig
 Between Groups       1087.143    8     135.893 48.985 0.000
  Within Groups       4016.980 1448       2.774             
          Total       5104.123 1456                         

multiple comparisons, bonferroni
 I J Mean_Difference Std_Error   Sig Lower_95 Upper_95 flag
 1 2           1.263     0.299 0.001    0.305    2.220    *
 1 3          -1.516     0.321 0.000   -2.544   -0.488    *
 1 4           0.387     0.323 1.000   -0.648    1.421     
 1 5           0.017     0.404 1.000   -1.278    1.311     
 1 6          -0.131     0.299 1.000   -1.090    0.827     
 1 7          -0.356     0.390 1.000   -1.603    0.892     
 1 8           0.275     0.442 1.000   -1.139    1.690     
 1 9          -0.129     0.531 1.000   -1.830    1.573     
 2 3          -2.779     0.156 0.000   -3.277   -2.280    *
 2 4          -0.876     0.160 0.000   -1.388   -0.364    *
 2 5          -1.246     0.291 0.001   -2.178   -0.314    *
 2 6          -1.394     0.104 0.000   -1.727   -1.061    *
 2 7          -1.618     0.270 0.000   -2.484   -0.753    *
 2 8          -0.987     0.341 0.138   -2.080    0.105     
 2 9          -1.391     0.451 0.075   -2.836    0.054     
 3 4           1.903     0.198 0.000    1.269    2.536    *
 3 5           1.533     0.313 0.000    0.529    2.536    *
 3 6           1.385     0.156 0.000    0.885    1.885    *
 3 7           1.160     0.294 0.003    0.218    2.103    *
 3 8           1.791     0.360 0.000    0.637    2.945    *
 3 9           1.387     0.466 0.106   -0.105    2.880     
 4 5          -0.370     0.315 1.000   -1.380    0.640     
 4 6          -0.518     0.160 0.046   -1.032   -0.004    *
 4 7          -0.742     0.296 0.447   -1.692    0.207     
 4 8          -0.111     0.362 1.000   -1.272    1.049     
 4 9          -0.515     0.467 1.000   -2.012    0.982     
 5 6          -0.148     0.291 1.000   -1.080    0.785     
 5 7          -0.372     0.383 1.000   -1.600    0.855     
 5 8           0.259     0.436 1.000   -1.138    1.656     
 5 9          -0.145     0.527 1.000   -1.832    1.542     
 6 7          -0.225     0.270 1.000   -1.091    0.642     
 6 8           0.406     0.341 1.000   -0.687    1.500     
 6 9           0.003     0.451 1.000   -1.443    1.448     
 7 8           0.631     0.423 1.000   -0.723    1.985     
 7 9           0.227     0.516 1.000   -1.424    1.879     
 8 9          -0.404     0.556 1.000   -2.185    1.377     
```

### Figures 2023: repair cost box-and-whisker plots

![Repair cost by parent material, 2023](figures/crosstab_anova_2023_fig1.png)

![Repair cost by hydrologic group, 2023](figures/crosstab_anova_2023_fig2.png)

![Repair cost by VTrans district, 2023](figures/crosstab_anova_2023_fig3.png)

## Output 2024

```
storm 2024 
segments 22571 
damaged 1018 
overall percent damaged 4.51 

step 1 parent material by damage

   Cell Contents
|-------------------------|
|                   Count |
|         Expected Values |
|             Row Percent |
|                Residual |
|-------------------------|

Total Observations in Table:  22571 

             | roads$damaged 
roads$PARENT |        0  |        1  | Row Total | 
-------------|-----------|-----------|-----------|
           A |      930  |       75  |     1005  | 
             |  959.672  |   45.328  |           | 
             |   92.537% |    7.463% |    4.453% | 
             |  -29.672  |   29.672  |           | 
-------------|-----------|-----------|-----------|
          DT |     8117  |      418  |     8535  | 
             | 8150.053  |  384.947  |           | 
             |   95.103% |    4.897% |   37.814% | 
             |  -33.053  |   33.053  |           | 
-------------|-----------|-----------|-----------|
       DT/GT |      123  |       11  |      134  | 
             |  127.956  |    6.044  |           | 
             |   91.791% |    8.209% |    0.594% | 
             |   -4.956  |    4.956  |           | 
-------------|-----------|-----------|-----------|
          GF |     2754  |      157  |     2911  | 
             | 2779.708  |  131.292  |           | 
             |   94.607% |    5.393% |   12.897% | 
             |  -25.708  |   25.708  |           | 
-------------|-----------|-----------|-----------|
       GF/GL |      101  |        3  |      104  | 
             |   99.309  |    4.691  |           | 
             |   97.115% |    2.885% |    0.461% | 
             |    1.691  |   -1.691  |           | 
-------------|-----------|-----------|-----------|
          GL |     3273  |      100  |     3373  | 
             | 3220.871  |  152.129  |           | 
             |   97.035% |    2.965% |   14.944% | 
             |   52.129  |  -52.129  |           | 
-------------|-----------|-----------|-----------|
          GT |     5437  |      237  |     5674  | 
             | 5418.091  |  255.909  |           | 
             |   95.823% |    4.177% |   25.138% | 
             |   18.909  |  -18.909  |           | 
-------------|-----------|-----------|-----------|
       GT/DT |      704  |       11  |      715  | 
             |  682.752  |   32.248  |           | 
             |   98.462% |    1.538% |    3.168% | 
             |   21.248  |  -21.248  |           | 
-------------|-----------|-----------|-----------|
     GT/GT_R |       89  |        2  |       91  | 
             |   86.896  |    4.104  |           | 
             |   97.802% |    2.198% |    0.403% | 
             |    2.104  |   -2.104  |           | 
-------------|-----------|-----------|-----------|
           M |       25  |        4  |       29  | 
             |   27.692  |    1.308  |           | 
             |   86.207% |   13.793% |    0.128% | 
             |   -2.692  |    2.692  |           | 
-------------|-----------|-----------|-----------|
Column Total |    21553  |     1018  |    22571  | 
-------------|-----------|-----------|-----------|

 
Statistics for All Table Factors


Pearson's Chi-squared test 
------------------------------------------------------------
Chi^2 =  75.24351     d.f. =  9     p =  1.414658e-12 


 
       Minimum expected frequency: 1.307962 
Cells with Expected Frequency < 5: 3 of 20 (15%)

Warning message:
In chisq.test(t, correct = FALSE, ...) :
  Chi-squared approximation may be incorrect

step 2 hydrologic group by damage

   Cell Contents
|-------------------------|
|                   Count |
|         Expected Values |
|             Row Percent |
|                Residual |
|-------------------------|

Total Observations in Table:  22571 

                 | roads$damaged 
roads$HYDROGROUP |        0  |        1  | Row Total | 
-----------------|-----------|-----------|-----------|
               A |     2267  |      137  |     2404  | 
                 | 2295.574  |  108.426  |           | 
                 |   94.301% |    5.699% |   10.651% | 
                 |  -28.574  |   28.574  |           | 
-----------------|-----------|-----------|-----------|
             A/D |      346  |       16  |      362  | 
                 |  345.673  |   16.327  |           | 
                 |   95.580% |    4.420% |    1.604% | 
                 |    0.327  |   -0.327  |           | 
-----------------|-----------|-----------|-----------|
               B |     2877  |      166  |     3043  | 
                 | 2905.754  |  137.246  |           | 
                 |   94.545% |    5.455% |   13.482% | 
                 |  -28.754  |   28.754  |           | 
-----------------|-----------|-----------|-----------|
             B/D |      776  |       48  |      824  | 
                 |  786.836  |   37.164  |           | 
                 |   94.175% |    5.825% |    3.651% | 
                 |  -10.836  |   10.836  |           | 
-----------------|-----------|-----------|-----------|
               C |     3901  |      146  |     4047  | 
                 | 3864.472  |  182.528  |           | 
                 |   96.392% |    3.608% |   17.930% | 
                 |   36.528  |  -36.528  |           | 
-----------------|-----------|-----------|-----------|
             C/D |     3330  |      190  |     3520  | 
                 | 3361.241  |  158.759  |           | 
                 |   94.602% |    5.398% |   15.595% | 
                 |  -31.241  |   31.241  |           | 
-----------------|-----------|-----------|-----------|
               D |     7850  |      306  |     8156  | 
                 | 7788.147  |  367.853  |           | 
                 |   96.248% |    3.752% |   36.135% | 
                 |   61.853  |  -61.853  |           | 
-----------------|-----------|-----------|-----------|
       not rated |      206  |        9  |      215  | 
                 |  205.303  |    9.697  |           | 
                 |   95.814% |    4.186% |    0.953% | 
                 |    0.697  |   -0.697  |           | 
-----------------|-----------|-----------|-----------|
    Column Total |    21553  |     1018  |    22571  | 
-----------------|-----------|-----------|-----------|

 
Statistics for All Table Factors


Pearson's Chi-squared test 
------------------------------------------------------------
Chi^2 =  42.5478     d.f. =  7     p =  4.077925e-07 


 
       Minimum expected frequency: 9.696956 


step 3 vtrans district by damage

   Cell Contents
|-------------------------|
|                   Count |
|         Expected Values |
|             Row Percent |
|                Residual |
|-------------------------|

Total Observations in Table:  22571 

                      | roads$damaged 
roads$vtrans_district |        0  |        1  | Row Total | 
----------------------|-----------|-----------|-----------|
                    3 |      137  |        0  |      137  | 
                      |  130.821  |    6.179  |           | 
                      |  100.000% |    0.000% |    0.607% | 
                      |    6.179  |   -6.179  |           | 
----------------------|-----------|-----------|-----------|
                    5 |     5267  |      231  |     5498  | 
                      | 5250.029  |  247.971  |           | 
                      |   95.798% |    4.202% |   24.359% | 
                      |   16.971  |  -16.971  |           | 
----------------------|-----------|-----------|-----------|
                    6 |     5924  |      539  |     6463  | 
                      | 6171.505  |  291.495  |           | 
                      |   91.660% |    8.340% |   28.634% | 
                      | -247.505  |  247.505  |           | 
----------------------|-----------|-----------|-----------|
                    7 |     7034  |      209  |     7243  | 
                      | 6916.325  |  326.675  |           | 
                      |   97.114% |    2.886% |   32.090% | 
                      |  117.675  | -117.675  |           | 
----------------------|-----------|-----------|-----------|
                    8 |      174  |        0  |      174  | 
                      |  166.152  |    7.848  |           | 
                      |  100.000% |    0.000% |    0.771% | 
                      |    7.848  |   -7.848  |           | 
----------------------|-----------|-----------|-----------|
                    9 |     3017  |       39  |     3056  | 
                      | 2918.168  |  137.832  |           | 
                      |   98.724% |    1.276% |   13.539% | 
                      |   98.832  |  -98.832  |           | 
----------------------|-----------|-----------|-----------|
         Column Total |    21553  |     1018  |    22571  | 
----------------------|-----------|-----------|-----------|

 
Statistics for All Table Factors


Pearson's Chi-squared test 
------------------------------------------------------------
Chi^2 =  354.5907     d.f. =  5     p =  1.79761e-74 


 
       Minimum expected frequency: 6.178991 


step 4 damaged segments for the anova
damaged segments 1018 
mean ln(cost), all damaged 9.518 

step 5 ln(cost) by parent material

descriptives
   group    N   Mean Std_Deviation Std_Error Lower_95 Upper_95 Minimum Maximum
       A   75  9.917         1.778     0.205    9.508   10.326   5.434  16.037
      DT  418  9.480         1.620     0.079    9.324    9.636   5.247  14.331
   DT/GT   11  9.099         1.316     0.397    8.215    9.982   6.752  11.082
      GF  157  9.614         1.640     0.131    9.355    9.872   4.984  14.037
   GF/GL    3 10.443         1.176     0.679    7.521   13.365   9.764  11.801
      GL  100  9.123         1.520     0.152    8.822    9.425   4.605  13.109
      GT  237  9.558         1.590     0.103    9.355    9.762   5.808  14.825
   GT/DT   11  9.747         2.111     0.637    8.328   11.165   6.260  13.535
 GT/GT_R    2  9.245         1.360     0.962   -2.978   21.469   8.283  10.207
       M    4  9.599         1.212     0.606    7.670   11.529   8.454  10.963
   Total 1018  9.518         1.623     0.051    9.418    9.617   4.605  16.037

test of homogeneity of variances
Levene's Test for Homogeneity of Variance (center = mean)
        Df F value Pr(>F)
group    9  0.6873  0.721
      1008               

anova
         Source Sum_of_Squares   df Mean_Square     F   Sig
 Between Groups         35.197    9       3.911 1.491 0.146
  Within Groups       2644.051 1008       2.623            
          Total       2679.248 1017                        

multiple comparisons, bonferroni
       I       J Mean_Difference Std_Error   Sig Lower_95 Upper_95 flag
       A      DT           0.437     0.203 1.000   -0.227    1.101     
       A   DT/GT           0.819     0.523 1.000   -0.891    2.529     
       A      GF           0.304     0.227 1.000   -0.440    1.047     
       A   GF/GL          -0.526     0.954 1.000   -3.644    2.593     
       A      GL           0.794     0.247 0.062   -0.015    1.603     
       A      GT           0.359     0.215 1.000   -0.343    1.061     
       A   GT/DT           0.171     0.523 1.000   -1.539    1.881     
       A GT/GT_R           0.672     1.160 1.000   -3.123    4.466     
       A       M           0.318     0.831 1.000   -2.400    3.036     
      DT   DT/GT           0.382     0.495 1.000   -1.236    1.999     
      DT      GF          -0.134     0.152 1.000   -0.629    0.362     
      DT   GF/GL          -0.963     0.938 1.000   -4.032    2.106     
      DT      GL           0.357     0.180 1.000   -0.233    0.946     
      DT      GT          -0.078     0.132 1.000   -0.509    0.353     
      DT   GT/DT          -0.267     0.495 1.000   -1.884    1.351     
      DT GT/GT_R           0.235     1.148 1.000   -3.519    3.989     
      DT       M          -0.119     0.814 1.000   -2.780    2.542     
   DT/GT      GF          -0.515     0.505 1.000   -2.167    1.137     
   DT/GT   GF/GL          -1.344     1.055 1.000   -4.794    2.105     
   DT/GT      GL          -0.025     0.514 1.000   -1.707    1.658     
   DT/GT      GT          -0.460     0.500 1.000   -2.093    1.174     
   DT/GT   GT/DT          -0.648     0.691 1.000   -2.906    1.610     
   DT/GT GT/GT_R          -0.147     1.245 1.000   -4.218    3.924     
   DT/GT       M          -0.501     0.946 1.000   -3.593    2.592     
      GF   GF/GL          -0.829     0.944 1.000   -3.916    2.258     
      GF      GL           0.490     0.207 0.817   -0.187    1.168     
      GF      GT           0.056     0.167 1.000   -0.489    0.601     
      GF   GT/DT          -0.133     0.505 1.000   -1.785    1.519     
      GF GT/GT_R           0.368     1.152 1.000   -3.401    4.137     
      GF       M           0.014     0.820 1.000   -2.667    2.696     
   GF/GL      GL           1.320     0.949 1.000   -1.784    4.423     
   GF/GL      GT           0.885     0.941 1.000   -2.192    3.962     
   GF/GL   GT/DT           0.696     1.055 1.000   -2.753    4.146     
   GF/GL GT/GT_R           1.198     1.478 1.000   -3.637    6.033     
   GF/GL       M           0.844     1.237 1.000   -3.202    4.889     
      GL      GT          -0.435     0.193 1.000   -1.066    0.197     
      GL   GT/DT          -0.623     0.514 1.000   -2.306    1.059     
      GL GT/GT_R          -0.122     1.157 1.000   -3.904    3.660     
      GL       M          -0.476     0.826 1.000   -3.177    2.225     
      GT   GT/DT          -0.189     0.500 1.000   -1.822    1.445     
      GT GT/GT_R           0.313     1.150 1.000   -3.448    4.074     
      GT       M          -0.041     0.817 1.000   -2.712    2.629     
   GT/DT GT/GT_R           0.501     1.245 1.000   -3.570    4.573     
   GT/DT       M           0.147     0.946 1.000   -2.945    3.240     
 GT/GT_R       M          -0.354     1.403 1.000   -4.941    4.233     

step 6 ln(cost) by hydrologic group

descriptives
     group    N   Mean Std_Deviation Std_Error Lower_95 Upper_95 Minimum
         A  137  9.535         1.592     0.136    9.266    9.804   4.984
       A/D   16  9.352         1.785     0.446    8.401   10.304   7.043
         B  166  9.695         1.693     0.131    9.435    9.954   5.808
       B/D   48 10.042         1.931     0.279    9.481   10.603   5.434
         C  146  9.504         1.471     0.122    9.263    9.744   5.247
       C/D  190  9.541         1.576     0.114    9.316    9.767   6.337
         D  306  9.324         1.608     0.092    9.143    9.505   4.605
 not rated    9  9.799         2.100     0.700    8.185   11.413   6.620
     Total 1018  9.518         1.623     0.051    9.418    9.617   4.605
 Maximum
  13.436
  13.091
  14.825
  16.037
  13.230
  14.078
  14.331
  14.037
  16.037

test of homogeneity of variances
Levene's Test for Homogeneity of Variance (center = mean)
        Df F value Pr(>F)
group    7   0.921 0.4892
      1010               

anova
         Source Sum_of_Squares   df Mean_Square     F   Sig
 Between Groups         31.227    7       4.461 1.701 0.105
  Within Groups       2648.022 1010       2.622            
          Total       2679.248 1017                        

multiple comparisons, bonferroni
   I         J Mean_Difference Std_Error   Sig Lower_95 Upper_95 flag
   A       A/D           0.183     0.428 1.000   -1.157    1.523     
   A         B          -0.160     0.187 1.000   -0.745    0.426     
   A       B/D          -0.507     0.272 1.000   -1.357    0.344     
   A         C           0.032     0.193 1.000   -0.572    0.635     
   A       C/D          -0.006     0.181 1.000   -0.574    0.563     
   A         D           0.211     0.166 1.000   -0.310    0.733     
   A not rated          -0.264     0.557 1.000   -2.009    1.481     
 A/D         B          -0.342     0.424 1.000   -1.670    0.985     
 A/D       B/D          -0.690     0.467 1.000   -2.154    0.774     
 A/D         C          -0.151     0.426 1.000   -1.487    1.184     
 A/D       C/D          -0.189     0.421 1.000   -1.509    1.131     
 A/D         D           0.029     0.415 1.000   -1.272    1.329     
 A/D not rated          -0.447     0.675 1.000   -2.560    1.666     
   B       B/D          -0.347     0.265 1.000   -1.178    0.484     
   B         C           0.191     0.184 1.000   -0.384    0.767     
   B       C/D           0.154     0.172 1.000   -0.385    0.692     
   B         D           0.371     0.156 0.495   -0.118    0.860     
   B not rated          -0.105     0.554 1.000   -1.840    1.631     
 B/D         C           0.538     0.269 1.000   -0.305    1.382     
 B/D       C/D           0.501     0.262 1.000   -0.318    1.320     
 B/D         D           0.718     0.251 0.122   -0.069    1.506     
 B/D not rated           0.243     0.588 1.000   -1.600    2.085     
   C       C/D          -0.038     0.178 1.000   -0.596    0.521     
   C         D           0.180     0.163 1.000   -0.330    0.690     
   C not rated          -0.296     0.556 1.000   -2.038    1.446     
 C/D         D           0.217     0.150 1.000   -0.251    0.686     
 C/D not rated          -0.258     0.552 1.000   -1.988    1.472     
   D not rated          -0.476     0.548 1.000   -2.191    1.240     

step 7 ln(cost) by vtrans district

descriptives
 group    N  Mean Std_Deviation Std_Error Lower_95 Upper_95 Minimum Maximum
     5  231 8.916         1.486     0.098    8.723    9.108   4.605  13.535
     6  539 9.633         1.554     0.067    9.502    9.765   5.247  14.331
     7  209 9.903         1.809     0.125    9.656   10.149   5.247  16.037
     9   39 9.420         1.351     0.216    8.982    9.858   5.759  11.269
 Total 1018 9.518         1.623     0.051    9.418    9.617   4.605  16.037

test of homogeneity of variances
Levene's Test for Homogeneity of Variance (center = mean)
        Df F value  Pr(>F)  
group    3  3.0052 0.02956 *
      1014                  
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

anova
         Source Sum_of_Squares   df Mean_Square      F   Sig
 Between Groups        122.215    3      40.738 16.155 0.000
  Within Groups       2557.034 1014       2.522             
          Total       2679.248 1017                         

multiple comparisons, bonferroni
 I J Mean_Difference Std_Error   Sig Lower_95 Upper_95 flag
 5 6          -0.717     0.125 0.000   -1.047   -0.387    *
 5 7          -0.987     0.152 0.000   -1.388   -0.586    *
 5 9          -0.504     0.275 0.403   -1.231    0.223     
 6 7          -0.270     0.129 0.225   -0.612    0.072     
 6 9           0.213     0.263 1.000   -0.483    0.909     
 7 9           0.483     0.277 0.489   -0.249    1.215     
```

### Figures 2024: repair cost box-and-whisker plots

![Repair cost by parent material, 2024](figures/crosstab_anova_2024_fig1.png)

![Repair cost by hydrologic group, 2024](figures/crosstab_anova_2024_fig2.png)

![Repair cost by VTrans district, 2024](figures/crosstab_anova_2024_fig3.png)
