# Data:  ../mrgp-roads-core/twopart_simple_2026-09-17/data/analysis_<STORM>_2026-09-27.csv
#        one row per road segment; compliance-graded segments with at least 3.51 in of rain,
#        Averill excluded; damaged = 1 if the de-duplicated repair cost (cost_dedup) is above zero
# Hydrogroup, Parent, VTRANS District


# step 0  read the data
STORM <- 2023

if (!dir.exists("data") && dir.exists("../data")) setwd("..")
roads <- read.csv(paste0("../mrgp-roads-core/twopart_simple_2026-09-17/data/analysis_", STORM, "_2026-09-27.csv"))
roads <- roads[roads$Precip >= 3.51 & roads$Town != "Averill", ]
roads$damaged <- as.integer(roads$cost_dedup > 0)

# the three categoricals as factors (district arrives as a number in the file)
roads$PARENT          <- factor(roads$PARENT)
roads$HYDROGROUP      <- factor(roads$HYDROGROUP)
roads$vtrans_district <- factor(roads$vtrans_district)

cat("segments\n");                 print(nrow(roads))
cat("damaged\n");                  print(sum(roads$damaged))
cat("overall percent damaged\n");  print(round(100 * mean(roads$damaged), 2))


# step 1  


# step 2  parent material
cat("step 2 parent material by damage\n\n")

tab <- table(Parent = roads$PARENT, Damaged = roads$damaged)

cat("observed counts with totals\n")
print(addmargins(tab))

cat("\nrow percentages\n")
print(round(100 * prop.table(tab, margin = 1), 2))

chi <- chisq.test(tab)

cat("\nexpected counts\n")
print(round(chi$expected, 1))

cat("\npearson chi-square test\n")
print(chi)

cat("\nadjusted residuals\n")
print(round(chi$stdres, 2))

cat("\nthe rule applied\n")
decision <- data.frame(category  = rownames(tab),
                       segments  = as.integer(rowSums(tab)),
                       damaged   = as.integer(tab[, "1"]),
                       expected  = round(chi$expected[, "1"], 1),
                       adj_resid = round(chi$stdres[, "1"], 2),
                       p         = round(2 * pnorm(-abs(chi$stdres[, "1"])), 4))
decision$indicator <- abs(decision$adj_resid) >= 2 & decision$expected >= 5 & decision$damaged >= 1
print(decision, row.names = FALSE)

parent_in <- decision$category[decision$indicator]
cat("\nparent material categories that become indicators:", paste(parent_in, collapse = ", "), "\n")
cat("left off (baseline):", paste(decision$category[!decision$indicator], collapse = ", "), "\n\n")


# step 3  hydrologic group
cat("step 3 hydrologic group by damage\n\n")

tab <- table(Hydrogroup = roads$HYDROGROUP, Damaged = roads$damaged)

cat("observed counts with totals\n")
print(addmargins(tab))

cat("\nrow percentages\n")
print(round(100 * prop.table(tab, margin = 1), 2))

chi <- chisq.test(tab)

cat("\nexpected counts\n")
print(round(chi$expected, 1))

cat("\npearson chi-square test\n")
print(chi)

cat("\nadjusted residuals\n")
print(round(chi$stdres, 2))

cat("\nthe rule applied\n")
decision <- data.frame(category  = rownames(tab),
                       segments  = as.integer(rowSums(tab)),
                       damaged   = as.integer(tab[, "1"]),
                       expected  = round(chi$expected[, "1"], 1),
                       adj_resid = round(chi$stdres[, "1"], 2),
                       p         = round(2 * pnorm(-abs(chi$stdres[, "1"])), 4))
decision$indicator <- abs(decision$adj_resid) >= 2 & decision$expected >= 5 & decision$damaged >= 1
print(decision, row.names = FALSE)

hydro_in <- decision$category[decision$indicator]
cat("\nhydrologic groups that become indicators:", paste(hydro_in, collapse = ", "), "\n")
cat("left off (baseline):", paste(decision$category[!decision$indicator], collapse = ", "), "\n\n")


# step 4  vtrans district
cat("step 4 vtrans district by damage\n\n")

tab <- table(District = roads$vtrans_district, Damaged = roads$damaged)

cat("observed counts with totals\n")
print(addmargins(tab))

cat("\nrow percentages\n")
print(round(100 * prop.table(tab, margin = 1), 2))

chi <- chisq.test(tab)

cat("\nexpected counts\n")
print(round(chi$expected, 1))

cat("\npearson chi-square test\n")
print(chi)

cat("\nadjusted residuals\n")
print(round(chi$stdres, 2))

cat("\nthe rule applied\n")
decision <- data.frame(category  = rownames(tab),
                       segments  = as.integer(rowSums(tab)),
                       damaged   = as.integer(tab[, "1"]),
                       expected  = round(chi$expected[, "1"], 1),
                       adj_resid = round(chi$stdres[, "1"], 2),
                       p         = round(2 * pnorm(-abs(chi$stdres[, "1"])), 4))
decision$indicator <- abs(decision$adj_resid) >= 2 & decision$expected >= 5 & decision$damaged >= 1
print(decision, row.names = FALSE)

district_in <- decision$category[decision$indicator]
cat("\ndistricts that become indicators:", paste(district_in, collapse = ", "), "\n")
cat("left off (baseline):", paste(decision$category[!decision$indicator], collapse = ", "), "\n\n")
# if every district passes, nothing is left off and the indicators rebuild the full factor


# step 5  the cost side: one-way anova of ln(cost) by category
# damaged segments only; cost on the natural-log scale
# the cost-side rule: a category becomes a 0/1 indicator when it is in at least one bonferroni pair with p < .05
dmg <- roads[roads$damaged == 1, ]
dmg$log_cost <- log(dmg$cost_dedup)
cat("step 5 damaged segments for the anova\n")
cat("damaged segments\n");          print(nrow(dmg))
cat("mean ln(cost), all damaged\n"); print(round(mean(dmg$log_cost), 3))


# step 6  parent material anova
cat("\nstep 6 ln(cost) by parent material\n\n")
dmg$PARENT <- factor(dmg$PARENT)

cat("damaged segments per group\n");   print(table(dmg$PARENT))
cat("\nmean ln(cost) per group\n");    print(round(tapply(dmg$log_cost, dmg$PARENT, mean), 3))

# error bar plot: mean +/- 2 standard errors per group
mean_g <- tapply(dmg$log_cost, dmg$PARENT, mean)
se_g   <- tapply(dmg$log_cost, dmg$PARENT, sd) / sqrt(table(dmg$PARENT))
plot(1:nlevels(dmg$PARENT), mean_g, ylim = range(c(mean_g, mean_g - 2 * se_g, mean_g + 2 * se_g), na.rm = TRUE),
     xlim = c(0.5, nlevels(dmg$PARENT) + 0.5), pch = 16, cex = 1.3, xaxt = "n",
     xlab = "Parent material", ylab = "Mean ln(cost)", main = paste("ln(cost) by parent material,", STORM, "storm"))
axis(1, at = 1:nlevels(dmg$PARENT), labels = paste0(levels(dmg$PARENT), "\nn=", table(dmg$PARENT)), cex.axis = 0.8, padj = 0.5)
arrows(1:nlevels(dmg$PARENT), mean_g - 2 * se_g, 1:nlevels(dmg$PARENT), mean_g + 2 * se_g, angle = 90, code = 3, length = 0.05)

cat("\nanova table\n")
print(summary(aov(log_cost ~ PARENT, data = dmg)))

cat("\nlevene test (large p = equal variances)\n")
spread <- abs(dmg$log_cost - ave(dmg$log_cost, dmg$PARENT))
print(anova(lm(spread ~ dmg$PARENT)))

cat("\nbonferroni pairwise comparisons\n")
pw <- pairwise.t.test(dmg$log_cost, dmg$PARENT, p.adjust.method = "bonferroni")
print(pw)

# a category is in a significant pair if any p in its row or its column of the grid is below .05
cat("the rule applied\n")
row_hit <- rowSums(pw$p.value < 0.05, na.rm = TRUE) > 0
col_hit <- colSums(pw$p.value < 0.05, na.rm = TRUE) > 0
decision <- data.frame(category     = levels(dmg$PARENT),
                       n            = as.integer(table(dmg$PARENT)),
                       mean_ln_cost = round(as.numeric(tapply(dmg$log_cost, dmg$PARENT, mean)), 3))
decision$indicator <- decision$category %in% c(names(row_hit)[row_hit], names(col_hit)[col_hit])
print(decision, row.names = FALSE)
cat("\nparent material categories that become cost indicators:",
    paste(decision$category[decision$indicator], collapse = ", "), "\n")
cat("left off (baseline):", paste(decision$category[!decision$indicator], collapse = ", "), "\n\n")


# step 7  hydrologic group anova
cat("step 7 ln(cost) by hydrologic group\n\n")
dmg$HYDROGROUP <- factor(dmg$HYDROGROUP)

cat("damaged segments per group\n");   print(table(dmg$HYDROGROUP))
cat("\nmean ln(cost) per group\n");    print(round(tapply(dmg$log_cost, dmg$HYDROGROUP, mean), 3))

mean_g <- tapply(dmg$log_cost, dmg$HYDROGROUP, mean)
se_g   <- tapply(dmg$log_cost, dmg$HYDROGROUP, sd) / sqrt(table(dmg$HYDROGROUP))
plot(1:nlevels(dmg$HYDROGROUP), mean_g, ylim = range(c(mean_g, mean_g - 2 * se_g, mean_g + 2 * se_g), na.rm = TRUE),
     xlim = c(0.5, nlevels(dmg$HYDROGROUP) + 0.5), pch = 16, cex = 1.3, xaxt = "n",
     xlab = "Hydrologic group", ylab = "Mean ln(cost)", main = paste("ln(cost) by hydrologic group,", STORM, "storm"))
axis(1, at = 1:nlevels(dmg$HYDROGROUP), labels = paste0(levels(dmg$HYDROGROUP), "\nn=", table(dmg$HYDROGROUP)), cex.axis = 0.8, padj = 0.5)
arrows(1:nlevels(dmg$HYDROGROUP), mean_g - 2 * se_g, 1:nlevels(dmg$HYDROGROUP), mean_g + 2 * se_g, angle = 90, code = 3, length = 0.05)

cat("\nanova table\n")
print(summary(aov(log_cost ~ HYDROGROUP, data = dmg)))

cat("\nlevene test (large p = equal variances)\n")
spread <- abs(dmg$log_cost - ave(dmg$log_cost, dmg$HYDROGROUP))
print(anova(lm(spread ~ dmg$HYDROGROUP)))

cat("\nbonferroni pairwise comparisons\n")
pw <- pairwise.t.test(dmg$log_cost, dmg$HYDROGROUP, p.adjust.method = "bonferroni")
print(pw)

cat("the rule applied\n")
row_hit <- rowSums(pw$p.value < 0.05, na.rm = TRUE) > 0
col_hit <- colSums(pw$p.value < 0.05, na.rm = TRUE) > 0
decision <- data.frame(category     = levels(dmg$HYDROGROUP),
                       n            = as.integer(table(dmg$HYDROGROUP)),
                       mean_ln_cost = round(as.numeric(tapply(dmg$log_cost, dmg$HYDROGROUP, mean)), 3))
decision$indicator <- decision$category %in% c(names(row_hit)[row_hit], names(col_hit)[col_hit])
print(decision, row.names = FALSE)
cat("\nhydrologic groups that become cost indicators:",
    paste(decision$category[decision$indicator], collapse = ", "), "\n")
cat("left off (baseline):", paste(decision$category[!decision$indicator], collapse = ", "), "\n\n")


# step 8  vtrans district anova
cat("step 8 ln(cost) by vtrans district\n\n")
dmg$vtrans_district <- factor(dmg$vtrans_district)

cat("damaged segments per group\n");   print(table(dmg$vtrans_district))
cat("\nmean ln(cost) per group\n");    print(round(tapply(dmg$log_cost, dmg$vtrans_district, mean), 3))

mean_g <- tapply(dmg$log_cost, dmg$vtrans_district, mean)
se_g   <- tapply(dmg$log_cost, dmg$vtrans_district, sd) / sqrt(table(dmg$vtrans_district))
plot(1:nlevels(dmg$vtrans_district), mean_g, ylim = range(c(mean_g, mean_g - 2 * se_g, mean_g + 2 * se_g), na.rm = TRUE),
     xlim = c(0.5, nlevels(dmg$vtrans_district) + 0.5), pch = 16, cex = 1.3, xaxt = "n",
     xlab = "VTrans district", ylab = "Mean ln(cost)", main = paste("ln(cost) by VTrans district,", STORM, "storm"))
axis(1, at = 1:nlevels(dmg$vtrans_district), labels = paste0(levels(dmg$vtrans_district), "\nn=", table(dmg$vtrans_district)), cex.axis = 0.8, padj = 0.5)
arrows(1:nlevels(dmg$vtrans_district), mean_g - 2 * se_g, 1:nlevels(dmg$vtrans_district), mean_g + 2 * se_g, angle = 90, code = 3, length = 0.05)

cat("\nanova table\n")
print(summary(aov(log_cost ~ vtrans_district, data = dmg)))

cat("\nlevene test (large p = equal variances)\n")
spread <- abs(dmg$log_cost - ave(dmg$log_cost, dmg$vtrans_district))
print(anova(lm(spread ~ dmg$vtrans_district)))

cat("\nbonferroni pairwise comparisons\n")
pw <- pairwise.t.test(dmg$log_cost, dmg$vtrans_district, p.adjust.method = "bonferroni")
print(pw)

cat("the rule applied\n")
row_hit <- rowSums(pw$p.value < 0.05, na.rm = TRUE) > 0
col_hit <- colSums(pw$p.value < 0.05, na.rm = TRUE) > 0
decision <- data.frame(category     = levels(dmg$vtrans_district),
                       n            = as.integer(table(dmg$vtrans_district)),
                       mean_ln_cost = round(as.numeric(tapply(dmg$log_cost, dmg$vtrans_district, mean)), 3))
decision$indicator <- decision$category %in% c(names(row_hit)[row_hit], names(col_hit)[col_hit])
print(decision, row.names = FALSE)
cat("\ndistricts that become cost indicators:",
    paste(decision$category[decision$indicator], collapse = ", "), "\n")
cat("left off (baseline):", paste(decision$category[!decision$indicator], collapse = ", "), "\n\n")
