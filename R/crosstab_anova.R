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
