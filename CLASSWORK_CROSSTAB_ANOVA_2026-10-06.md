# Crosstabs and ANOVA, Classwork Style

The task: cross-tabulate each categorical predictor (parent material, hydrologic group, VTrans district) against damage, then run the one-way ANOVA of ln(cost) by category on the damaged segments, and decide from each which categories can enter the two-part model as 0/1 indicators with the rest left off. Written in the step-by-step form of the GEOG 3505 exercises (Exercise 6 for the crosstabs, Exercise 3 for the ANOVA) so it can be stepped through line by line. Graded segments, Precip >= 3.51 in, Averill excluded; damage and cost from cost_dedup. Script run once with STORM = 2023 and once with STORM = 2024.

## Script: crosstab_binaries.R

```r
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
```

## Output 2023

```
--- Number of segments ---
[1] 59695
--- Number damaged ---
[1] 1457
--- Overall percent damaged ---
[1] 2.44
=== STEP 2: parent material by damage ===

--- Observed counts with totals ---
       Damaged
Parent      0     1   Sum
  A      2389   115  2504
  DT    23838   682 24520
  DT/GT   169     8   177
  GF     9947   213 10160
  GL     2456    43  2499
  GT    17767   386 18153
  GT/DT  1509     9  1518
  M       163     1   164
  Sum   58238  1457 59695

--- Row percentages (within each parent material) ---
       Damaged
Parent      0     1
  A     95.41  4.59
  DT    97.22  2.78
  DT/GT 95.48  4.52
  GF    97.90  2.10
  GL    98.28  1.72
  GT    97.87  2.13
  GT/DT 99.41  0.59
  M     99.39  0.61

--- Expected counts (if parent material and damage were unrelated) ---
       Damaged
Parent        0     1
  A      2442.9  61.1
  DT    23921.5 598.5
  DT/GT   172.7   4.3
  GF     9912.0 248.0
  GL     2438.0  61.0
  GT    17709.9 443.1
  GT/DT  1480.9  37.1
  M       160.0   4.0

--- Pearson chi-square test ---

	Pearson's Chi-squared test

data:  tab
X-squared = 105.97, df = 7, p-value < 2.2e-16

--- Adjusted residuals ---
       Damaged
Parent      0     1
  A     -7.13  7.13
  DT    -4.50  4.50
  DT/GT -1.80  1.80
  GF     2.47 -2.47
  GL     2.38 -2.38
  GT     3.29 -3.29
  GT/DT  4.73 -4.73
  M      1.52 -1.52

--- The rule applied ---
 category segments damaged expected adj_resid      p indicator
        A     2504     115     61.1      7.13 0.0000      TRUE
       DT    24520     682    598.5      4.50 0.0000      TRUE
    DT/GT      177       8      4.3      1.80 0.0726     FALSE
       GF    10160     213    248.0     -2.47 0.0136      TRUE
       GL     2499      43     61.0     -2.38 0.0172      TRUE
       GT    18153     386    443.1     -3.29 0.0010      TRUE
    GT/DT     1518       9     37.1     -4.73 0.0000      TRUE
        M      164       1      4.0     -1.52 0.1281     FALSE

Parent material categories that become indicators: A, DT, GF, GL, GT, GT/DT 
Left off (baseline): DT/GT, M 

=== STEP 3: hydrologic group by damage ===

--- Observed counts with totals ---
           Damaged
Hydrogroup      0     1   Sum
  A          9350   206  9556
  A/D         745     8   753
  B          8310   254  8564
  B/D        1976    75  2051
  C         14074   341 14415
  C/D        9331   193  9524
  D         14173   370 14543
  not rated   279    10   289
  Sum       58238  1457 59695

--- Row percentages (within each hydrologic group) ---
           Damaged
Hydrogroup      0     1
  A         97.84  2.16
  A/D       98.94  1.06
  B         97.03  2.97
  B/D       96.34  3.66
  C         97.63  2.37
  C/D       97.97  2.03
  D         97.46  2.54
  not rated 96.54  3.46

--- Expected counts ---
           Damaged
Hydrogroup        0     1
  A          9322.8 233.2
  A/D         734.6  18.4
  B          8355.0 209.0
  B/D        2000.9  50.1
  C         14063.2 351.8
  C/D        9291.5 232.5
  D         14188.0 355.0
  not rated   281.9   7.1

--- Pearson chi-square test ---

	Pearson's Chi-squared test

data:  tab
X-squared = 41.045, df = 7, p-value = 7.936e-07

--- Adjusted residuals ---
           Damaged
Hydrogroup      0     1
  A          1.97 -1.97
  A/D        2.47 -2.47
  B         -3.40  3.40
  B/D       -3.63  3.63
  C          0.67 -0.67
  C/D        2.86 -2.86
  D         -0.93  0.93
  not rated -1.13  1.13

--- The rule applied ---
  category segments damaged expected adj_resid      p indicator
         A     9556     206    233.2     -1.97 0.0488     FALSE
       A/D      753       8     18.4     -2.47 0.0136      TRUE
         B     8564     254    209.0      3.40 0.0007      TRUE
       B/D     2051      75     50.1      3.63 0.0003      TRUE
         C    14415     341    351.8     -0.67 0.5020     FALSE
       C/D     9524     193    232.5     -2.86 0.0043      TRUE
         D    14543     370    355.0      0.93 0.3526     FALSE
 not rated      289      10      7.1      1.13 0.2602     FALSE

Hydrologic groups that become indicators: A/D, B, B/D, C/D 
Left off (baseline): A, C, D, not rated 

=== STEP 4: VTrans district by damage ===

--- Observed counts with totals ---
        Damaged
District     0     1   Sum
     1    6046    33  6079
     2    6069   520  6589
     3   10360   147 10507
     4   10425   137 10562
     5    2645    35  2680
     6    9650   505 10155
     7    3617    41  3658
     8    4196    25  4221
     9    5230    14  5244
     Sum 58238  1457 59695

--- Row percentages (within each district) ---
        Damaged
District     0     1
       1 99.46  0.54
       2 92.11  7.89
       3 98.60  1.40
       4 98.70  1.30
       5 98.69  1.31
       6 95.03  4.97
       7 98.88  1.12
       8 99.41  0.59
       9 99.73  0.27

--- Expected counts ---
        Damaged
District       0     1
       1  5930.6 148.4
       2  6428.2 160.8
       3 10250.6 256.4
       4 10304.2 257.8
       5  2614.6  65.4
       6  9907.1 247.9
       7  3568.7  89.3
       8  4118.0 103.0
       9  5116.0 128.0

--- Pearson chi-square test ---

	Pearson's Chi-squared test

data:  tab
X-squared = 1499.5, df = 8, p-value < 2.2e-16

--- Adjusted residuals ---
        Damaged
District      0      1
       1  10.12 -10.12
       2 -30.40  30.40
       3   7.62  -7.62
       4   8.40  -8.40
       5   3.90  -3.90
       6 -18.15  18.15
       7   5.34  -5.34
       8   8.07  -8.07
       9  10.68 -10.68

--- The rule applied ---
 category segments damaged expected adj_resid     p indicator
        1     6079      33    148.4    -10.12 0e+00      TRUE
        2     6589     520    160.8     30.40 0e+00      TRUE
        3    10507     147    256.4     -7.62 0e+00      TRUE
        4    10562     137    257.8     -8.40 0e+00      TRUE
        5     2680      35     65.4     -3.90 1e-04      TRUE
        6    10155     505    247.9     18.15 0e+00      TRUE
        7     3658      41     89.3     -5.34 0e+00      TRUE
        8     4221      25    103.0     -8.07 0e+00      TRUE
        9     5244      14    128.0    -10.68 0e+00      TRUE

Districts that become indicators: 1, 2, 3, 4, 5, 6, 7, 8, 9 
Left off (baseline):  

=== STEP 5: damaged segments for the ANOVA ===
--- Number of damaged segments ---
[1] 1457
--- Mean ln(cost), all damaged ---
[1] 8.902

=== STEP 6: ln(cost) by parent material ===

--- Damaged segments per group ---

    A    DT DT/GT    GF    GL    GT GT/DT     M 
  115   682     8   213    43   386     9     1 

--- Mean ln(cost) per group ---
     A     DT  DT/GT     GF     GL     GT  GT/DT      M 
 9.251  8.738  8.905  9.305  9.803  8.736  9.821 11.521 

--- ANOVA table ---
              Df Sum Sq Mean Sq F value   Pr(>F)    
PARENT         7    127  18.160   5.287 5.56e-06 ***
Residuals   1449   4977   3.435                     
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

--- Levene test (large p = equal variances) ---
Analysis of Variance Table

Response: spread
             Df  Sum Sq Mean Sq F value Pr(>F)
dmg$PARENT    7    4.73 0.67555   0.529 0.8131
Residuals  1449 1850.40 1.27702               

--- Bonferroni pairwise comparisons ---

	Pairwise comparisons using t tests with pooled SD 

data:  dmg$log_cost and dmg$PARENT 

      A      DT     DT/GT  GF     GL     GT     GT/DT 
DT    0.1700 -      -      -      -      -      -     
DT/GT 1.0000 1.0000 -      -      -      -      -     
GF    1.0000 0.0028 1.0000 -      -      -      -     
GL    1.0000 0.0074 1.0000 1.0000 -      -      -     
GT    0.2521 1.0000 1.0000 0.0092 0.0099 -      -     
GT/DT 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 -     
M     1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000

P value adjustment method: bonferroni 
--- The rule applied ---
 category   n mean_ln_cost indicator
        A 115        9.251     FALSE
       DT 682        8.738      TRUE
    DT/GT   8        8.905     FALSE
       GF 213        9.305      TRUE
       GL  43        9.803      TRUE
       GT 386        8.736      TRUE
    GT/DT   9        9.821     FALSE
        M   1       11.521     FALSE

Parent material categories that become cost indicators: DT, GF, GL, GT 
Left off (baseline): A, DT/GT, GT/DT, M 

=== STEP 7: ln(cost) by hydrologic group ===

--- Damaged segments per group ---

        A       A/D         B       B/D         C       C/D         D not rated 
      206         8       254        75       341       193       370        10 

--- Mean ln(cost) per group ---
        A       A/D         B       B/D         C       C/D         D not rated 
    9.325     9.415     8.641     9.206     8.729     8.790     8.964     9.824 

--- ANOVA table ---
              Df Sum Sq Mean Sq F value   Pr(>F)    
HYDROGROUP     7     86  12.255   3.538 0.000892 ***
Residuals   1449   5018   3.463                     
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

--- Levene test (large p = equal variances) ---
Analysis of Variance Table

Response: spread
                 Df  Sum Sq Mean Sq F value Pr(>F)
dmg$HYDROGROUP    7    7.49  1.0694  0.8333 0.5595
Residuals      1449 1859.49  1.2833               

--- Bonferroni pairwise comparisons ---

	Pairwise comparisons using t tests with pooled SD 

data:  dmg$log_cost and dmg$HYDROGROUP 

          A      A/D    B      B/D    C      C/D    D     
A/D       1.0000 -      -      -      -      -      -     
B         0.0026 1.0000 -      -      -      -      -     
B/D       1.0000 1.0000 0.5898 -      -      -      -     
C         0.0082 1.0000 1.0000 1.0000 -      -      -     
C/D       0.1162 1.0000 1.0000 1.0000 1.0000 -      -     
D         0.7189 1.0000 0.9363 1.0000 1.0000 1.0000 -     
not rated 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000

P value adjustment method: bonferroni 
--- The rule applied ---
  category   n mean_ln_cost indicator
         A 206        9.325      TRUE
       A/D   8        9.415     FALSE
         B 254        8.641      TRUE
       B/D  75        9.206     FALSE
         C 341        8.729      TRUE
       C/D 193        8.790     FALSE
         D 370        8.964     FALSE
 not rated  10        9.824     FALSE

Hydrologic groups that become cost indicators: A, B, C 
Left off (baseline): A/D, B/D, C/D, D, not rated 

=== STEP 8: ln(cost) by VTrans district ===

--- Damaged segments per group ---

  1   2   3   4   5   6   7   8   9 
 33 520 147 137  35 505  41  25  14 

--- Mean ln(cost) per group ---
     1      2      3      4      5      6      7      8      9 
 9.184  7.921 10.700  8.798  9.167  9.315  9.540  8.909  9.313 

--- ANOVA table ---
                  Df Sum Sq Mean Sq F value Pr(>F)    
vtrans_district    8   1087  135.89   48.98 <2e-16 ***
Residuals       1448   4017    2.77                   
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

--- Levene test (large p = equal variances) ---
Analysis of Variance Table

Response: spread
                      Df  Sum Sq Mean Sq F value    Pr(>F)    
dmg$vtrans_district    8   55.92  6.9903  6.8069 8.376e-09 ***
Residuals           1448 1487.01  1.0269                      
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

--- Bonferroni pairwise comparisons ---

	Pairwise comparisons using t tests with pooled SD 

data:  dmg$log_cost and dmg$vtrans_district 

  1       2       3       4       5       6       7       8      
2 0.00092 -       -       -       -       -       -       -      
3 9.1e-05 < 2e-16 -       -       -       -       -       -      
4 1.00000 1.8e-06 < 2e-16 -       -       -       -       -      
5 1.00000 0.00070 4.0e-05 1.00000 -       -       -       -      
6 1.00000 < 2e-16 < 2e-16 0.04608 1.00000 -       -       -      
7 1.00000 9.5e-08 0.00302 0.44661 1.00000 1.00000 -       -      
8 1.00000 0.13842 2.7e-05 1.00000 1.00000 1.00000 1.00000 -      
9 1.00000 0.07492 0.10607 1.00000 1.00000 1.00000 1.00000 1.00000

P value adjustment method: bonferroni 
--- The rule applied ---
 category   n mean_ln_cost indicator
        1  33        9.184      TRUE
        2 520        7.921      TRUE
        3 147       10.700      TRUE
        4 137        8.798      TRUE
        5  35        9.167      TRUE
        6 505        9.315      TRUE
        7  41        9.540      TRUE
        8  25        8.909      TRUE
        9  14        9.313     FALSE

Districts that become cost indicators: 1, 2, 3, 4, 5, 6, 7, 8 
Left off (baseline): 9 

```

## Output 2024

```
--- Number of segments ---
[1] 22571
--- Number damaged ---
[1] 1018
--- Overall percent damaged ---
[1] 4.51
=== STEP 2: parent material by damage ===

--- Observed counts with totals ---
         Damaged
Parent        0     1   Sum
  A         930    75  1005
  DT       8117   418  8535
  DT/GT     123    11   134
  GF       2754   157  2911
  GF/GL     101     3   104
  GL       3273   100  3373
  GT       5437   237  5674
  GT/DT     704    11   715
  GT/GT_R    89     2    91
  M          25     4    29
  Sum     21553  1018 22571

--- Row percentages (within each parent material) ---
         Damaged
Parent        0     1
  A       92.54  7.46
  DT      95.10  4.90
  DT/GT   91.79  8.21
  GF      94.61  5.39
  GF/GL   97.12  2.88
  GL      97.04  2.96
  GT      95.82  4.18
  GT/DT   98.46  1.54
  GT/GT_R 97.80  2.20
  M       86.21 13.79

--- Expected counts (if parent material and damage were unrelated) ---
         Damaged
Parent         0     1
  A        959.7  45.3
  DT      8150.1 384.9
  DT/GT    128.0   6.0
  GF      2779.7 131.3
  GF/GL     99.3   4.7
  GL      3220.9 152.1
  GT      5418.1 255.9
  GT/DT    682.8  32.2
  GT/GT_R   86.9   4.1
  M         27.7   1.3

--- Pearson chi-square test ---

	Pearson's Chi-squared test

data:  tab
X-squared = 75.244, df = 9, p-value = 1.415e-12

--- Adjusted residuals ---
         Damaged
Parent        0     1
  A       -4.61  4.61
  DT      -2.19  2.19
  DT/GT   -2.07  2.07
  GF      -2.46  2.46
  GF/GL    0.80 -0.80
  GL       4.69 -4.69
  GT       1.40 -1.40
  GT/DT    3.89 -3.89
  GT/GT_R  1.07 -1.07
  M       -2.41  2.41

--- The rule applied ---
 category segments damaged expected adj_resid      p indicator
        A     1005      75     45.3      4.61 0.0000      TRUE
       DT     8535     418    384.9      2.19 0.0288      TRUE
    DT/GT      134      11      6.0      2.07 0.0385      TRUE
       GF     2911     157    131.3      2.46 0.0139      TRUE
    GF/GL      104       3      4.7     -0.80 0.4233     FALSE
       GL     3373     100    152.1     -4.69 0.0000      TRUE
       GT     5674     237    255.9     -1.40 0.1621     FALSE
    GT/DT      715      11     32.2     -3.89 0.0001      TRUE
  GT/GT_R       91       2      4.1     -1.07 0.2868     FALSE
        M       29       4      1.3      2.41 0.0159     FALSE

Parent material categories that become indicators: A, DT, DT/GT, GF, GL, GT/DT 
Left off (baseline): GF/GL, GT, GT/GT_R, M 

=== STEP 3: hydrologic group by damage ===

--- Observed counts with totals ---
           Damaged
Hydrogroup      0     1   Sum
  A          2267   137  2404
  A/D         346    16   362
  B          2877   166  3043
  B/D         776    48   824
  C          3901   146  4047
  C/D        3330   190  3520
  D          7850   306  8156
  not rated   206     9   215
  Sum       21553  1018 22571

--- Row percentages (within each hydrologic group) ---
           Damaged
Hydrogroup      0     1
  A         94.30  5.70
  A/D       95.58  4.42
  B         94.54  5.46
  B/D       94.17  5.83
  C         96.39  3.61
  C/D       94.60  5.40
  D         96.25  3.75
  not rated 95.81  4.19

--- Expected counts ---
           Damaged
Hydrogroup       0     1
  A         2295.6 108.4
  A/D        345.7  16.3
  B         2905.8 137.2
  B/D        786.8  37.2
  C         3864.5 182.5
  C/D       3361.2 158.8
  D         7788.1 367.9
  not rated  205.3   9.7

--- Pearson chi-square test ---

	Pearson's Chi-squared test

data:  tab
X-squared = 42.548, df = 7, p-value = 4.078e-07

--- Adjusted residuals ---
           Damaged
Hydrogroup      0     1
  A         -2.97  2.97
  A/D        0.08 -0.08
  B         -2.70  2.70
  B/D       -1.85  1.85
  C          3.05 -3.05
  C/D       -2.76  2.76
  D          4.13 -4.13
  not rated  0.23 -0.23

--- The rule applied ---
  category segments damaged expected adj_resid      p indicator
         A     2404     137    108.4      2.97 0.0030      TRUE
       A/D      362      16     16.3     -0.08 0.9335     FALSE
         B     3043     166    137.2      2.70 0.0069      TRUE
       B/D      824      48     37.2      1.85 0.0639     FALSE
         C     4047     146    182.5     -3.05 0.0023      TRUE
       C/D     3520     190    158.8      2.76 0.0057      TRUE
         D     8156     306    367.9     -4.13 0.0000      TRUE
 not rated      215       9      9.7     -0.23 0.8180     FALSE

Hydrologic groups that become indicators: A, B, C, C/D, D 
Left off (baseline): A/D, B/D, not rated 

=== STEP 4: VTrans district by damage ===

--- Observed counts with totals ---
        Damaged
District     0     1   Sum
     3     137     0   137
     5    5267   231  5498
     6    5924   539  6463
     7    7034   209  7243
     8     174     0   174
     9    3017    39  3056
     Sum 21553  1018 22571

--- Row percentages (within each district) ---
        Damaged
District      0      1
       3 100.00   0.00
       5  95.80   4.20
       6  91.66   8.34
       7  97.11   2.89
       8 100.00   0.00
       9  98.72   1.28

--- Expected counts ---
        Damaged
District      0     1
       3  130.8   6.2
       5 5250.0 248.0
       6 6171.5 291.5
       7 6916.3 326.7
       8  166.2   7.8
       9 2918.2 137.8

--- Pearson chi-square test ---

	Pearson's Chi-squared test

data:  tab
X-squared = 354.59, df = 5, p-value < 2.2e-16

--- Adjusted residuals ---
        Damaged
District      0      1
       3   2.55  -2.55
       5   1.27  -1.27
       6 -17.56  17.56
       7   8.08  -8.08
       8   2.88  -2.88
       9   9.26  -9.26

--- The rule applied ---
 category segments damaged expected adj_resid      p indicator
        3      137       0      6.2     -2.55 0.0107     FALSE
        5     5498     231    248.0     -1.27 0.2048     FALSE
        6     6463     539    291.5     17.56 0.0000      TRUE
        7     7243     209    326.7     -8.08 0.0000      TRUE
        8      174       0      7.8     -2.88 0.0040     FALSE
        9     3056      39    137.8     -9.26 0.0000      TRUE

Districts that become indicators: 6, 7, 9 
Left off (baseline): 3, 5, 8 

=== STEP 5: damaged segments for the ANOVA ===
--- Number of damaged segments ---
[1] 1018
--- Mean ln(cost), all damaged ---
[1] 9.518

=== STEP 6: ln(cost) by parent material ===

--- Damaged segments per group ---

      A      DT   DT/GT      GF   GF/GL      GL      GT   GT/DT GT/GT_R       M 
     75     418      11     157       3     100     237      11       2       4 

--- Mean ln(cost) per group ---
      A      DT   DT/GT      GF   GF/GL      GL      GT   GT/DT GT/GT_R       M 
  9.917   9.480   9.099   9.614  10.443   9.123   9.558   9.747   9.245   9.599 

--- ANOVA table ---
              Df Sum Sq Mean Sq F value Pr(>F)
PARENT         9   35.2   3.911   1.491  0.146
Residuals   1008 2644.1   2.623               

--- Levene test (large p = equal variances) ---
Analysis of Variance Table

Response: spread
             Df Sum Sq Mean Sq F value Pr(>F)
dmg$PARENT    9    5.9 0.65535  0.6873  0.721
Residuals  1008  961.2 0.95357               

--- Bonferroni pairwise comparisons ---

	Pairwise comparisons using t tests with pooled SD 

data:  dmg$log_cost and dmg$PARENT 

        A     DT    DT/GT GF    GF/GL GL    GT    GT/DT GT/GT_R
DT      1.000 -     -     -     -     -     -     -     -      
DT/GT   1.000 1.000 -     -     -     -     -     -     -      
GF      1.000 1.000 1.000 -     -     -     -     -     -      
GF/GL   1.000 1.000 1.000 1.000 -     -     -     -     -      
GL      0.062 1.000 1.000 0.817 1.000 -     -     -     -      
GT      1.000 1.000 1.000 1.000 1.000 1.000 -     -     -      
GT/DT   1.000 1.000 1.000 1.000 1.000 1.000 1.000 -     -      
GT/GT_R 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 -      
M       1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000  

P value adjustment method: bonferroni 
--- The rule applied ---
 category   n mean_ln_cost indicator
        A  75        9.917     FALSE
       DT 418        9.480     FALSE
    DT/GT  11        9.099     FALSE
       GF 157        9.614     FALSE
    GF/GL   3       10.443     FALSE
       GL 100        9.123     FALSE
       GT 237        9.558     FALSE
    GT/DT  11        9.747     FALSE
  GT/GT_R   2        9.245     FALSE
        M   4        9.599     FALSE

Parent material categories that become cost indicators:  
Left off (baseline): A, DT, DT/GT, GF, GF/GL, GL, GT, GT/DT, GT/GT_R, M 

=== STEP 7: ln(cost) by hydrologic group ===

--- Damaged segments per group ---

        A       A/D         B       B/D         C       C/D         D not rated 
      137        16       166        48       146       190       306         9 

--- Mean ln(cost) per group ---
        A       A/D         B       B/D         C       C/D         D not rated 
    9.535     9.352     9.695    10.042     9.504     9.541     9.324     9.799 

--- ANOVA table ---
              Df Sum Sq Mean Sq F value Pr(>F)
HYDROGROUP     7   31.2   4.461   1.701  0.105
Residuals   1010 2648.0   2.622               

--- Levene test (large p = equal variances) ---
Analysis of Variance Table

Response: spread
                 Df Sum Sq Mean Sq F value Pr(>F)
dmg$HYDROGROUP    7   6.22 0.88787   0.921 0.4892
Residuals      1010 973.72 0.96408               

--- Bonferroni pairwise comparisons ---

	Pairwise comparisons using t tests with pooled SD 

data:  dmg$log_cost and dmg$HYDROGROUP 

          A    A/D  B    B/D  C    C/D  D   
A/D       1.00 -    -    -    -    -    -   
B         1.00 1.00 -    -    -    -    -   
B/D       1.00 1.00 1.00 -    -    -    -   
C         1.00 1.00 1.00 1.00 -    -    -   
C/D       1.00 1.00 1.00 1.00 1.00 -    -   
D         1.00 1.00 0.49 0.12 1.00 1.00 -   
not rated 1.00 1.00 1.00 1.00 1.00 1.00 1.00

P value adjustment method: bonferroni 
--- The rule applied ---
  category   n mean_ln_cost indicator
         A 137        9.535     FALSE
       A/D  16        9.352     FALSE
         B 166        9.695     FALSE
       B/D  48       10.042     FALSE
         C 146        9.504     FALSE
       C/D 190        9.541     FALSE
         D 306        9.324     FALSE
 not rated   9        9.799     FALSE

Hydrologic groups that become cost indicators:  
Left off (baseline): A, A/D, B, B/D, C, C/D, D, not rated 

=== STEP 8: ln(cost) by VTrans district ===

--- Damaged segments per group ---

  5   6   7   9 
231 539 209  39 

--- Mean ln(cost) per group ---
    5     6     7     9 
8.916 9.633 9.903 9.420 

--- ANOVA table ---
                  Df Sum Sq Mean Sq F value   Pr(>F)    
vtrans_district    3  122.2   40.74   16.16 2.91e-10 ***
Residuals       1014 2557.0    2.52                     
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

--- Levene test (large p = equal variances) ---
Analysis of Variance Table

Response: spread
                      Df Sum Sq Mean Sq F value  Pr(>F)  
dmg$vtrans_district    3   8.24 2.74569  3.0052 0.02956 *
Residuals           1014 926.44 0.91365                  
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

--- Bonferroni pairwise comparisons ---

	Pairwise comparisons using t tests with pooled SD 

data:  dmg$log_cost and dmg$vtrans_district 

  5       6    7   
6 7.3e-08 -    -   
7 7.1e-10 0.22 -   
9 0.40    1.00 0.49

P value adjustment method: bonferroni 
--- The rule applied ---
 category   n mean_ln_cost indicator
        5 231        8.916      TRUE
        6 539        9.633      TRUE
        7 209        9.903      TRUE
        9  39        9.420     FALSE

Districts that become cost indicators: 5, 6, 7 
Left off (baseline): 9 

```
