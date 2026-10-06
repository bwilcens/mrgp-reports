# =============================================================================
# Crosstabs and one-way ANOVA on the categorical predictors: which categories can become 0/1 indicators?
#
# Data:  ../mrgp-roads-core/twopart_simple_2026-09-17/data/analysis_<STORM>_2026-09-27.csv
#        one row per road segment; compliance-graded segments with at least 3.51 in of rain,
#        Averill excluded; damaged = 1 if the de-duplicated repair cost (cost_dedup) is above zero
#        PARENT           parent material of the soil (A alluvial, DT dense till, GT ablation till,
#                         GF glacio-fluvial, GL glacio-lacustrine; the rest are mixed or rare codes)
#        HYDROGROUP       hydrologic soil group (A well drained ... D poorly drained; A/D etc. dual)
#        vtrans_district  VTrans maintenance district, 1 to 9
#
# The task: cross-tabulate each categorical against damage, then decide which categories carry
# their own signal and can be entered in the model as a 0/1 indicator, leaving the rest off.
#
# Run one line at a time with Ctrl+Enter and read the Console after each line. Set STORM to 2023
# or 2024. Nothing is written to disk.
# =============================================================================


# -----------------------------------------------------------------------------
# STEP 0.  Point R at the data and read it in
# -----------------------------------------------------------------------------
STORM <- 2023

if (!dir.exists("data") && dir.exists("../data")) setwd("..")
roads <- read.csv(paste0("../mrgp-roads-core/twopart_simple_2026-09-17/data/analysis_", STORM, "_2026-09-27.csv"))
roads <- roads[roads$Precip >= 3.51 & roads$Town != "Averill", ]
roads$damaged <- as.integer(roads$cost_dedup > 0)

# the three categoricals as factors (district arrives as a number in the file; a number would be fitted as
# a straight line, not as categories)
roads$PARENT          <- factor(roads$PARENT)
roads$HYDROGROUP      <- factor(roads$HYDROGROUP)
roads$vtrans_district <- factor(roads$vtrans_district)

cat("--- Number of segments ---\n");        print(nrow(roads))
cat("--- Number damaged ---\n");            print(sum(roads$damaged))
cat("--- Overall percent damaged ---\n");   print(round(100 * mean(roads$damaged), 2))


# -----------------------------------------------------------------------------
# STEP 1.  The decision rule, stated before looking at any table
# -----------------------------------------------------------------------------
# A category becomes a 0/1 indicator when all three hold:
#   (a) its adjusted residual in the damaged column is at least 2 in size -- the category has
#       clearly more (positive) or fewer (negative) damaged segments than independence predicts;
#   (b) its expected count in the damaged column is at least 5 -- the usual minimum for the
#       chi-square approximation, and SPSS's footnote threshold;
#   (c) it has at least one damaged segment -- a category with none cannot have a coefficient.
# Every category that fails is left off and becomes part of the baseline.
#
# The adjusted residual is (observed - expected) divided by the standard deviation that difference
# would have if category and damage were unrelated, so it reads like a z-score; its two-sided p-value
# is 2 * P(Z > |residual|), and a residual of 2 is p = .046, which is why 2 is the threshold. R returns it as
# chisq.test(...)$stdres; SPSS prints it when "Adjusted standardized" is ticked under Cells.


# -----------------------------------------------------------------------------
# STEP 2.  Parent material
# -----------------------------------------------------------------------------
# This replaces  Analyze > Descriptive Statistics > Crosstabs  with PARENT as the row variable
# and damaged as the column variable, with Expected and Adjusted standardized ticked under Cells
# and Chi-square ticked under Statistics.
cat("=== STEP 2: parent material by damage ===\n\n")

tab <- table(Parent = roads$PARENT, Damaged = roads$damaged)

cat("--- Observed counts with totals ---\n")
print(addmargins(tab))

cat("\n--- Row percentages (within each parent material) ---\n")
print(round(100 * prop.table(tab, margin = 1), 2))

chi <- chisq.test(tab)

cat("\n--- Expected counts (if parent material and damage were unrelated) ---\n")
print(round(chi$expected, 1))

cat("\n--- Pearson chi-square test ---\n")
print(chi)

cat("\n--- Adjusted residuals ---\n")
print(round(chi$stdres, 2))

# The decision, category by category. Read the damaged column ("1") of the residuals and the
# expected counts against the rule in STEP 1.
cat("\n--- The rule applied ---\n")
decision <- data.frame(category  = rownames(tab),
                       segments  = as.integer(rowSums(tab)),
                       damaged   = as.integer(tab[, "1"]),
                       expected  = round(chi$expected[, "1"], 1),
                       adj_resid = round(chi$stdres[, "1"], 2),
                       p         = round(2 * pnorm(-abs(chi$stdres[, "1"])), 4))
decision$indicator <- abs(decision$adj_resid) >= 2 & decision$expected >= 5 & decision$damaged >= 1
print(decision, row.names = FALSE)

parent_in <- decision$category[decision$indicator]
cat("\nParent material categories that become indicators:", paste(parent_in, collapse = ", "), "\n")
cat("Left off (baseline):", paste(decision$category[!decision$indicator], collapse = ", "), "\n\n")


# -----------------------------------------------------------------------------
# STEP 3.  Hydrologic group
# -----------------------------------------------------------------------------
cat("=== STEP 3: hydrologic group by damage ===\n\n")

tab <- table(Hydrogroup = roads$HYDROGROUP, Damaged = roads$damaged)

cat("--- Observed counts with totals ---\n")
print(addmargins(tab))

cat("\n--- Row percentages (within each hydrologic group) ---\n")
print(round(100 * prop.table(tab, margin = 1), 2))

chi <- chisq.test(tab)

cat("\n--- Expected counts ---\n")
print(round(chi$expected, 1))

cat("\n--- Pearson chi-square test ---\n")
print(chi)

cat("\n--- Adjusted residuals ---\n")
print(round(chi$stdres, 2))

cat("\n--- The rule applied ---\n")
decision <- data.frame(category  = rownames(tab),
                       segments  = as.integer(rowSums(tab)),
                       damaged   = as.integer(tab[, "1"]),
                       expected  = round(chi$expected[, "1"], 1),
                       adj_resid = round(chi$stdres[, "1"], 2),
                       p         = round(2 * pnorm(-abs(chi$stdres[, "1"])), 4))
decision$indicator <- abs(decision$adj_resid) >= 2 & decision$expected >= 5 & decision$damaged >= 1
print(decision, row.names = FALSE)

hydro_in <- decision$category[decision$indicator]
cat("\nHydrologic groups that become indicators:", paste(hydro_in, collapse = ", "), "\n")
cat("Left off (baseline):", paste(decision$category[!decision$indicator], collapse = ", "), "\n\n")


# -----------------------------------------------------------------------------
# STEP 4.  VTrans district
# -----------------------------------------------------------------------------
cat("=== STEP 4: VTrans district by damage ===\n\n")

tab <- table(District = roads$vtrans_district, Damaged = roads$damaged)

cat("--- Observed counts with totals ---\n")
print(addmargins(tab))

cat("\n--- Row percentages (within each district) ---\n")
print(round(100 * prop.table(tab, margin = 1), 2))

chi <- chisq.test(tab)

cat("\n--- Expected counts ---\n")
print(round(chi$expected, 1))

cat("\n--- Pearson chi-square test ---\n")
print(chi)

cat("\n--- Adjusted residuals ---\n")
print(round(chi$stdres, 2))

cat("\n--- The rule applied ---\n")
decision <- data.frame(category  = rownames(tab),
                       segments  = as.integer(rowSums(tab)),
                       damaged   = as.integer(tab[, "1"]),
                       expected  = round(chi$expected[, "1"], 1),
                       adj_resid = round(chi$stdres[, "1"], 2),
                       p         = round(2 * pnorm(-abs(chi$stdres[, "1"])), 4))
decision$indicator <- abs(decision$adj_resid) >= 2 & decision$expected >= 5 & decision$damaged >= 1
print(decision, row.names = FALSE)

district_in <- decision$category[decision$indicator]
cat("\nDistricts that become indicators:", paste(district_in, collapse = ", "), "\n")
cat("Left off (baseline):", paste(decision$category[!decision$indicator], collapse = ", "), "\n\n")
# If every district passes, nothing is left off and the indicators would just rebuild the full
# factor. The district is the storm's footprint, so that is the expected result in 2023.


# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# STEP 5.  The cost side: one-way ANOVA of ln(cost) by category (Exercise 3)
# -----------------------------------------------------------------------------
# The crosstab asks whether a category changes the chance of damage. The ANOVA asks, among the
# segments that were damaged, whether a category changes the repair cost. Cost is taken on the
# natural-log scale because a few very large repairs would otherwise dominate the means.
#
# The decision rule for the cost side: a category becomes a 0/1 indicator when it appears in at
# least one Bonferroni pairwise comparison with p below .05 -- that is, its mean cost differs from
# at least one other category's after correcting for the number of pairs. The rest are left off.
dmg <- roads[roads$damaged == 1, ]
dmg$log_cost <- log(dmg$cost_dedup)
cat("=== STEP 5: damaged segments for the ANOVA ===\n")
cat("--- Number of damaged segments ---\n");  print(nrow(dmg))
cat("--- Mean ln(cost), all damaged ---\n");  print(round(mean(dmg$log_cost), 3))


# -----------------------------------------------------------------------------
# STEP 6.  Parent material: ANOVA
# -----------------------------------------------------------------------------
# This replaces  Analyze > Compare Means > One-Way ANOVA  with log_cost as the dependent variable,
# PARENT as the factor, Bonferroni under Post Hoc, and Homogeneity of variance test under Options.
cat("\n=== STEP 6: ln(cost) by parent material ===\n\n")
dmg$PARENT <- factor(dmg$PARENT)

cat("--- Damaged segments per group ---\n");   print(table(dmg$PARENT))
cat("\n--- Mean ln(cost) per group ---\n");    print(round(tapply(dmg$log_cost, dmg$PARENT, mean), 3))

cat("\n--- ANOVA table ---\n")
print(summary(aov(log_cost ~ PARENT, data = dmg)))

cat("\n--- Levene test (large p = equal variances) ---\n")
spread <- abs(dmg$log_cost - ave(dmg$log_cost, dmg$PARENT))
print(anova(lm(spread ~ dmg$PARENT)))

cat("\n--- Bonferroni pairwise comparisons ---\n")
pw <- pairwise.t.test(dmg$log_cost, dmg$PARENT, p.adjust.method = "bonferroni")
print(pw)

# The decision: a category is in a significant pair if any p in its row or its column of the grid
# is below .05 (rows hold the later categories, columns the earlier ones).
cat("--- The rule applied ---\n")
row_hit <- rowSums(pw$p.value < 0.05, na.rm = TRUE) > 0
col_hit <- colSums(pw$p.value < 0.05, na.rm = TRUE) > 0
decision <- data.frame(category     = levels(dmg$PARENT),
                       n            = as.integer(table(dmg$PARENT)),
                       mean_ln_cost = round(as.numeric(tapply(dmg$log_cost, dmg$PARENT, mean)), 3))
decision$indicator <- decision$category %in% c(names(row_hit)[row_hit], names(col_hit)[col_hit])
print(decision, row.names = FALSE)
cat("\nParent material categories that become cost indicators:",
    paste(decision$category[decision$indicator], collapse = ", "), "\n")
cat("Left off (baseline):", paste(decision$category[!decision$indicator], collapse = ", "), "\n\n")


# -----------------------------------------------------------------------------
# STEP 7.  Hydrologic group: ANOVA
# -----------------------------------------------------------------------------
cat("=== STEP 7: ln(cost) by hydrologic group ===\n\n")
dmg$HYDROGROUP <- factor(dmg$HYDROGROUP)

cat("--- Damaged segments per group ---\n");   print(table(dmg$HYDROGROUP))
cat("\n--- Mean ln(cost) per group ---\n");    print(round(tapply(dmg$log_cost, dmg$HYDROGROUP, mean), 3))

cat("\n--- ANOVA table ---\n")
print(summary(aov(log_cost ~ HYDROGROUP, data = dmg)))

cat("\n--- Levene test (large p = equal variances) ---\n")
spread <- abs(dmg$log_cost - ave(dmg$log_cost, dmg$HYDROGROUP))
print(anova(lm(spread ~ dmg$HYDROGROUP)))

cat("\n--- Bonferroni pairwise comparisons ---\n")
pw <- pairwise.t.test(dmg$log_cost, dmg$HYDROGROUP, p.adjust.method = "bonferroni")
print(pw)

cat("--- The rule applied ---\n")
row_hit <- rowSums(pw$p.value < 0.05, na.rm = TRUE) > 0
col_hit <- colSums(pw$p.value < 0.05, na.rm = TRUE) > 0
decision <- data.frame(category     = levels(dmg$HYDROGROUP),
                       n            = as.integer(table(dmg$HYDROGROUP)),
                       mean_ln_cost = round(as.numeric(tapply(dmg$log_cost, dmg$HYDROGROUP, mean)), 3))
decision$indicator <- decision$category %in% c(names(row_hit)[row_hit], names(col_hit)[col_hit])
print(decision, row.names = FALSE)
cat("\nHydrologic groups that become cost indicators:",
    paste(decision$category[decision$indicator], collapse = ", "), "\n")
cat("Left off (baseline):", paste(decision$category[!decision$indicator], collapse = ", "), "\n\n")


# -----------------------------------------------------------------------------
# STEP 8.  VTrans district: ANOVA
# -----------------------------------------------------------------------------
cat("=== STEP 8: ln(cost) by VTrans district ===\n\n")
dmg$vtrans_district <- factor(dmg$vtrans_district)

cat("--- Damaged segments per group ---\n");   print(table(dmg$vtrans_district))
cat("\n--- Mean ln(cost) per group ---\n");    print(round(tapply(dmg$log_cost, dmg$vtrans_district, mean), 3))

cat("\n--- ANOVA table ---\n")
print(summary(aov(log_cost ~ vtrans_district, data = dmg)))

cat("\n--- Levene test (large p = equal variances) ---\n")
spread <- abs(dmg$log_cost - ave(dmg$log_cost, dmg$vtrans_district))
print(anova(lm(spread ~ dmg$vtrans_district)))

cat("\n--- Bonferroni pairwise comparisons ---\n")
pw <- pairwise.t.test(dmg$log_cost, dmg$vtrans_district, p.adjust.method = "bonferroni")
print(pw)

cat("--- The rule applied ---\n")
row_hit <- rowSums(pw$p.value < 0.05, na.rm = TRUE) > 0
col_hit <- colSums(pw$p.value < 0.05, na.rm = TRUE) > 0
decision <- data.frame(category     = levels(dmg$vtrans_district),
                       n            = as.integer(table(dmg$vtrans_district)),
                       mean_ln_cost = round(as.numeric(tapply(dmg$log_cost, dmg$vtrans_district, mean)), 3))
decision$indicator <- decision$category %in% c(names(row_hit)[row_hit], names(col_hit)[col_hit])
print(decision, row.names = FALSE)
cat("\nDistricts that become cost indicators:",
    paste(decision$category[decision$indicator], collapse = ", "), "\n")
cat("Left off (baseline):", paste(decision$category[!decision$indicator], collapse = ", "), "\n\n")


# =============================================================================
# END
#
# To report:
#   Occurrence side -- for each categorical, the chi-square (value, df, p), the categories chosen as
#   indicators with their adjusted residuals, and the categories left off and why.
#   Cost side -- for each categorical, F (df, p), the Levene test, the Bonferroni pairs below .05,
#   and the categories chosen as cost indicators.
# =============================================================================
