# Binary Selection from the Crosstab and ANOVA

The advisor's rule: use the crosstab and ANOVA results to choose which categories of parent material, hydrologic group and VTrans district enter the two-part model as 0/1 indicators, and leave every other category out. Part 1 indicators come from the crosstab (adjusted residual of 2 or more in size, expected count at least 5, at least one damaged segment); Part 2 indicators come from the ANOVA (a category in at least one Bonferroni pair with Sig. below .05). The script prints each category with its numbers and the decision, then refits the two-part model with the indicators and compares it with the full-factor model. Graded segments, Precip >= 3.51 in, Averill excluded; cost_dedup, uncapped; SEs clustered by town. Script run once with STORM = 2023 and once with STORM = 2024.

## Script: binary_select.R

```r
# binary_select.R -- the categorical predictors (parent material, hydrologic group, VTrans district) replaced by
#   0/1 indicators chosen from the crosstab and ANOVA results; every other category is left out (baseline).
#   Part 1 indicators: categories whose adjusted residual in the damaged column is >= 2 in size, with expected
#   count >= 5 and at least one damaged segment. Part 2 indicators: categories in at least one Bonferroni pair
#   with Sig. < .05. Then the two-part model is refitted with the indicators and compared with the full-factor model.
# Set STORM to 2023 or 2024 and run from mrgp-roads-core:
#   "C:/Program Files/R/R-4.5.3/bin/Rscript.exe" twopart_simple_2026-09-17/R/binary_select.R
# Graded segments, Precip >= 3.51 in, Averill excluded; cost_dedup, uncapped; SEs clustered by town.

STORM <- 2023

library(glm2)
library(sandwich)

OTHER1 <- "MeanSlope90m + StreamOrder + PercentImpervious_BaseLC_90m + Precip + RoadGrade_mean_deg + Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m + road_miles_municipal + grand_list_equalized_muni_100M"
OTHER2 <- "MeanSlope90m + StreamCrossing + StreamOrder + K_final + PercentImpervious_BaseLC_90m + Precip + Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m + road_miles_municipal + grand_list_equalized_muni_100M"

d <- read.csv(paste0("twopart_simple_2026-09-17/data/analysis_", STORM, "_2026-09-27.csv"))
d <- d[d$Precip >= 3.51 & d$Town != "Averill", ]
d$damaged   <- as.integer(d$cost_dedup > 0)
d$compliant <- relevel(factor(d$compliant), ref = "Does Not Meet")
k <- d[d$damaged == 1, ]
k$ln_cost <- log(k$cost_dedup)

cat("PART 1 SELECTION", STORM, "-- from the crosstab: |adjusted residual| >= 2, expected count >= 5, at least one damaged\n")
part1_terms <- NULL
for (v in c("PARENT", "HYDROGROUP", "vtrans_district")) {
  tab <- table(d[[v]], d$damaged)
  chi <- chisq.test(tab)
  pick <- abs(chi$stdres[, "1"]) >= 2 & chi$expected[, "1"] >= 5 & tab[, "1"] >= 1
  if (all(pick)) pick[which.min(abs(chi$stdres[, "1"]))] <- FALSE    # every level cannot be in: the weakest becomes the baseline
  cat("\n", v, "\n")
  print(data.frame(category = rownames(tab), segments = as.integer(rowSums(tab)), damaged = as.integer(tab[, "1"]),
                   expected = round(chi$expected[, "1"], 1), adj_resid = round(chi$stdres[, "1"], 2), indicator = pick), row.names = FALSE)
  for (r in rownames(tab)[pick]) {
    nm <- paste0(v, "_", gsub("[^A-Za-z0-9]", "_", r))
    d[[nm]] <- as.integer(d[[v]] == r)
    part1_terms <- c(part1_terms, nm)
  }
}
cat("\nPart 1 indicators:", paste(part1_terms, collapse = ", "), "\n")

cat("\nPART 2 SELECTION", STORM, "-- from the ANOVA: categories in at least one Bonferroni pair with Sig. < .05\n")
part2_terms <- NULL
for (v in c("PARENT", "HYDROGROUP", "vtrans_district")) {
  g <- factor(k[[v]])
  a <- summary(aov(k$ln_cost ~ g))[[1]]
  n_g <- as.numeric(table(g)); mean_g <- as.numeric(tapply(k$ln_cost, g, mean))
  pairs <- nlevels(g) * (nlevels(g) - 1) / 2
  in_pair <- rep(FALSE, nlevels(g))
  for (i in 1:(nlevels(g) - 1)) for (j in (i + 1):nlevels(g)) {
    se <- sqrt(a$`Mean Sq`[2] * (1 / n_g[i] + 1 / n_g[j]))
    p  <- min(1, 2 * pt(-abs((mean_g[i] - mean_g[j]) / se), a$Df[2]) * pairs)
    if (p < 0.05) in_pair[c(i, j)] <- TRUE
  }
  if (all(in_pair)) in_pair[which.min(n_g)] <- FALSE                    # every level cannot be in: the smallest becomes the baseline
  cat("\n", v, "  F =", round(a$`F value`[1], 3), " p =", signif(a$`Pr(>F)`[1], 3), "\n")
  print(data.frame(category = levels(g), n = n_g, mean_ln_cost = round(mean_g, 3), indicator = in_pair), row.names = FALSE)
  for (r in levels(g)[in_pair]) {
    nm <- paste0(v, "_", gsub("[^A-Za-z0-9]", "_", r))
    d[[nm]] <- as.integer(d[[v]] == r)
    k[[nm]] <- as.integer(k[[v]] == r)
    part2_terms <- c(part2_terms, nm)
  }
}
cat("\nPart 2 indicators:", paste(part2_terms, collapse = ", "), "\n")

# ---- the full-factor model for comparison (zero-damage districts pooled, as in the deck script) ----
n <- tapply(d$damaged, d$vtrans_district, sum)
d$district_f <- factor(ifelse(d$vtrans_district %in% names(n)[n == 0], names(which.max(n)), d$vtrans_district))
d$district_f <- relevel(d$district_f, ref = names(which.max(n)))
k$district_f <- d$district_f[d$damaged == 1]
full1 <- glm(as.formula(paste("damaged ~ compliant + PARENT + HYDROGROUP + district_f +", OTHER1)), family = binomial, data = d)
full2 <- glm2(as.formula(paste("cost_dedup ~ PARENT + HYDROGROUP + district_f +", OTHER2)), family = Gamma(link = "log"), data = k)

cat("\nPART 1 WITH INDICATORS", STORM, "\n")
f1 <- paste("damaged ~ compliant +", paste(c(part1_terms, OTHER1), collapse = " + "))
p1 <- glm(as.formula(f1), family = binomial, data = d)
print(summary(p1))
se <- sqrt(diag(vcovCL(p1, cluster = d$Town, type = "HC1")))
print(round(exp(cbind(odds_ratio = coef(p1), lo_95 = coef(p1) - 1.96 * se, hi_95 = coef(p1) + 1.96 * se)), 3))
n1 <- sum(d$damaged); n0 <- sum(d$damaged == 0)
cat("AUC, indicators\n");    print(round((sum(rank(fitted(p1))[d$damaged == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0), 4))
cat("AUC, full factors\n");  print(round((sum(rank(fitted(full1))[d$damaged == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0), 4))
cat("AIC, indicators / full factors\n"); print(round(c(AIC(p1), AIC(full1)), 1))
cat("drop in deviance, indicators vs full factors\n"); print(anova(p1, full1, test = "Chisq"))

cat("\nPART 2 WITH INDICATORS", STORM, "\n")
f2 <- paste("cost_dedup ~", paste(c(part2_terms, OTHER2), collapse = " + "))
p2 <- glm2(as.formula(f2), family = Gamma(link = "log"), data = k)
print(summary(p2))
se <- sqrt(diag(vcovCL(p2, cluster = k$Town, type = "HC1")))
print(round(exp(cbind(cost_ratio = coef(p2), lo_95 = coef(p2) - 1.96 * se, hi_95 = coef(p2) + 1.96 * se)), 3))
cat("deviance R2, indicators / full factors\n"); print(round(c(1 - p2$deviance / p2$null.deviance, 1 - full2$deviance / full2$null.deviance), 4))
cat("AIC, indicators / full factors\n"); print(round(c(AIC(p2), AIC(full2)), 1))
cat("F test, indicators vs full factors\n"); print(anova(p2, full2, test = "F"))

cat("\nAVOIDED COST", STORM, "\n")
d0 <- d
d0$compliant <- factor("Does Not Meet", levels = levels(d$compliant))
p_observed <- predict(p1, newdata = d,  type = "response")
p_recoded  <- predict(p1, newdata = d0, type = "response")
cost_hat   <- predict(p2, newdata = d,  type = "response")
avoided    <- (p_recoded - p_observed) * cost_hat
compliant  <- d$compliant == "Compliant"
cat("statewide avoided cost, indicators\n");   print(round(sum(avoided[compliant])))
cat("statewide avoided cost, full factors\n")
print(round(sum(((predict(full1, newdata = d0, type = "response") - predict(full1, newdata = d, type = "response")) * predict(full2, newdata = d, type = "response"))[compliant])))

# town-clustered interval for the indicator model's avoided cost (delta method, joint HC1 sandwich as in R/105)
k1 <- names(coef(p1)); k2 <- names(coef(p2))
X1 <- model.matrix(delete.response(terms(p1)), d,  xlev = p1$xlevels)[, k1]
X0 <- model.matrix(delete.response(terms(p1)), d0, xlev = p1$xlevels)[, k1]
Z  <- model.matrix(delete.response(terms(p2)), d,  xlev = p2$xlevels)[, k2]
h  <- c(colSums(((cost_hat * p_recoded * (1 - p_recoded)) * X0 - (cost_hat * p_observed * (1 - p_observed)) * X1)[compliant, ]),
        colSums(((p_recoded - p_observed) * cost_hat * Z)[compliant, ]))
B1 <- summary(p1)$cov.unscaled
w  <- weights(p2, "working"); rr <- residuals(p2, "working") * w
B2 <- summary(p2)$cov.unscaled * sum(rr^2) / sum(w)
towns <- sort(unique(d$Town))
U1 <- matrix(0, length(towns), length(k1), dimnames = list(towns, k1)); s1 <- rowsum(estfun(p1), d$Town); U1[rownames(s1), ] <- s1
U2 <- matrix(0, length(towns), length(k2), dimnames = list(towns, k2)); s2 <- rowsum(estfun(p2), k$Town); U2[rownames(s2), ] <- s2
G1 <- length(towns); G2 <- length(unique(k$Town))
a1 <- G1 / (G1 - 1) * (nrow(d) - 1) / (nrow(d) - length(k1)); a2 <- G2 / (G2 - 1) * (nrow(k) - 1) / (nrow(k) - length(k2))
V  <- crossprod(cbind(sqrt(a1) * U1 %*% B1, sqrt(a2) * U2 %*% B2))
se_avoided <- sqrt(drop(t(h) %*% V %*% h))
cat("clustered 95% interval, indicators\n"); print(round(sum(avoided[compliant]) + c(-1.96, 1.96) * se_avoided))
cat("ten sample towns, indicators\n")
print(round(tapply(avoided[compliant], d$Town[compliant], sum)[c("Brattleboro", "Corinth", "Wallingford", "Hardwick",
      "West Windsor", "Wolcott", "Washington", "Richmond", "Stamford", "Starksboro")]))
```

## Output 2023

```
PART 1 SELECTION 2023 -- from the crosstab: |adjusted residual| >= 2, expected count >= 5, at least one damaged

 PARENT 
 category segments damaged expected adj_resid indicator
        A     2504     115     61.1      7.13      TRUE
       DT    24520     682    598.5      4.50      TRUE
    DT/GT      177       8      4.3      1.80     FALSE
       GF    10160     213    248.0     -2.47      TRUE
       GL     2499      43     61.0     -2.38      TRUE
       GT    18153     386    443.1     -3.29      TRUE
    GT/DT     1518       9     37.1     -4.73      TRUE
        M      164       1      4.0     -1.52     FALSE

 HYDROGROUP 
  category segments damaged expected adj_resid indicator
         A     9556     206    233.2     -1.97     FALSE
       A/D      753       8     18.4     -2.47      TRUE
         B     8564     254    209.0      3.40      TRUE
       B/D     2051      75     50.1      3.63      TRUE
         C    14415     341    351.8     -0.67     FALSE
       C/D     9524     193    232.5     -2.86      TRUE
         D    14543     370    355.0      0.93     FALSE
 not rated      289      10      7.1      1.13     FALSE

 vtrans_district 
 category segments damaged expected adj_resid indicator
        1     6079      33    148.4    -10.12      TRUE
        2     6589     520    160.8     30.40      TRUE
        3    10507     147    256.4     -7.62      TRUE
        4    10562     137    257.8     -8.40      TRUE
        5     2680      35     65.4     -3.90     FALSE
        6    10155     505    247.9     18.15      TRUE
        7     3658      41     89.3     -5.34      TRUE
        8     4221      25    103.0     -8.07      TRUE
        9     5244      14    128.0    -10.68      TRUE

Part 1 indicators: PARENT_A, PARENT_DT, PARENT_GF, PARENT_GL, PARENT_GT, PARENT_GT_DT, HYDROGROUP_A_D, HYDROGROUP_B, HYDROGROUP_B_D, HYDROGROUP_C_D, vtrans_district_1, vtrans_district_2, vtrans_district_3, vtrans_district_4, vtrans_district_6, vtrans_district_7, vtrans_district_8, vtrans_district_9 

PART 2 SELECTION 2023 -- from the ANOVA: categories in at least one Bonferroni pair with Sig. < .05

 PARENT   F = 5.287  p = 5.56e-06 
 category   n mean_ln_cost indicator
        A 115        9.251     FALSE
       DT 682        8.738      TRUE
    DT/GT   8        8.905     FALSE
       GF 213        9.305      TRUE
       GL  43        9.803      TRUE
       GT 386        8.736      TRUE
    GT/DT   9        9.821     FALSE
        M   1       11.521     FALSE

 HYDROGROUP   F = 3.538  p = 0.000892 
  category   n mean_ln_cost indicator
         A 206        9.325      TRUE
       A/D   8        9.415     FALSE
         B 254        8.641      TRUE
       B/D  75        9.206     FALSE
         C 341        8.729      TRUE
       C/D 193        8.790     FALSE
         D 370        8.964     FALSE
 not rated  10        9.824     FALSE

 vtrans_district   F = 48.985  p = 3.05e-70 
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

Part 2 indicators: PARENT_DT, PARENT_GF, PARENT_GL, PARENT_GT, HYDROGROUP_A, HYDROGROUP_B, HYDROGROUP_C, vtrans_district_1, vtrans_district_2, vtrans_district_3, vtrans_district_4, vtrans_district_5, vtrans_district_6, vtrans_district_7, vtrans_district_8 

PART 1 WITH INDICATORS 2023 

Call:
glm(formula = as.formula(f1), family = binomial, data = d)

Coefficients:
                                Estimate Std. Error z value Pr(>|z|)    
(Intercept)                    -9.064810   0.482553 -18.785  < 2e-16 ***
compliantCompliant             -0.224982   0.066028  -3.407 0.000656 ***
PARENT_A                        0.077213   0.377760   0.204 0.838044    
PARENT_DT                      -0.032458   0.351919  -0.092 0.926515    
PARENT_GF                      -0.331783   0.354917  -0.935 0.349883    
PARENT_GL                      -0.113242   0.382634  -0.296 0.767266    
PARENT_GT                      -0.513881   0.355225  -1.447 0.147999    
PARENT_GT_DT                   -1.119658   0.493939  -2.267 0.023403 *  
HYDROGROUP_A_D                 -0.094033   0.371895  -0.253 0.800386    
HYDROGROUP_B                    0.178927   0.090280   1.982 0.047489 *  
HYDROGROUP_B_D                 -0.065974   0.179492  -0.368 0.713203    
HYDROGROUP_C_D                 -0.199535   0.088945  -2.243 0.024875 *  
vtrans_district_1              -0.956361   0.248265  -3.852 0.000117 ***
vtrans_district_2               1.348164   0.185355   7.273 3.51e-13 ***
vtrans_district_3              -0.144056   0.196521  -0.733 0.463540    
vtrans_district_4              -0.674493   0.200996  -3.356 0.000791 ***
vtrans_district_6               0.546291   0.188567   2.897 0.003767 ** 
vtrans_district_7              -0.313914   0.237973  -1.319 0.187130    
vtrans_district_8              -1.005367   0.267729  -3.755 0.000173 ***
vtrans_district_9              -2.196452   0.323005  -6.800 1.05e-11 ***
MeanSlope90m                    0.064328   0.007531   8.542  < 2e-16 ***
StreamOrder                     0.286928   0.026845  10.688  < 2e-16 ***
PercentImpervious_BaseLC_90m    1.483485   0.456924   3.247 0.001168 ** 
Precip                          0.384849   0.024983  15.404  < 2e-16 ***
RoadGrade_mean_deg              0.049478   0.015045   3.289 0.001007 ** 
SurfaceUnpaved                  0.387465   0.080448   4.816 1.46e-06 ***
Culvert_Any                     0.509731   0.059421   8.578  < 2e-16 ***
Driveway_Count                  0.008350   0.037718   0.221 0.824796    
Mean_TopoConvergence_30m        0.259465   0.025946  10.000  < 2e-16 ***
road_miles_municipal            0.001553   0.001518   1.023 0.306438    
grand_list_equalized_muni_100M -0.033701   0.006762  -4.984 6.23e-07 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 13697  on 59694  degrees of freedom
Residual deviance: 11499  on 59664  degrees of freedom
AIC: 11561

Number of Fisher Scoring iterations: 8

                               odds_ratio lo_95  hi_95
(Intercept)                         0.000 0.000  0.001
compliantCompliant                  0.799 0.674  0.946
PARENT_A                            1.080 0.432  2.701
PARENT_DT                           0.968 0.385  2.432
PARENT_GF                           0.718 0.293  1.758
PARENT_GL                           0.893 0.329  2.422
PARENT_GT                           0.598 0.240  1.489
PARENT_GT_DT                        0.326 0.088  1.214
HYDROGROUP_A_D                      0.910 0.445  1.860
HYDROGROUP_B                        1.196 0.958  1.493
HYDROGROUP_B_D                      0.936 0.643  1.364
HYDROGROUP_C_D                      0.819 0.645  1.040
vtrans_district_1                   0.384 0.095  1.548
vtrans_district_2                   3.850 1.177 12.597
vtrans_district_3                   0.866 0.249  3.015
vtrans_district_4                   0.509 0.140  1.850
vtrans_district_6                   1.727 0.530  5.631
vtrans_district_7                   0.731 0.202  2.639
vtrans_district_8                   0.366 0.105  1.270
vtrans_district_9                   0.111 0.031  0.395
MeanSlope90m                        1.066 1.042  1.092
StreamOrder                         1.332 1.250  1.420
PercentImpervious_BaseLC_90m        4.408 1.628 11.936
Precip                              1.469 1.256  1.719
RoadGrade_mean_deg                  1.051 1.014  1.088
SurfaceUnpaved                      1.473 1.132  1.917
Culvert_Any                         1.665 1.447  1.916
Driveway_Count                      1.008 0.900  1.130
Mean_TopoConvergence_30m            1.296 1.198  1.402
road_miles_municipal                1.002 0.992  1.011
grand_list_equalized_muni_100M      0.967 0.916  1.020
AUC, indicators
[1] 0.8229
AUC, full factors
[1] 0.8244
AIC, indicators / full factors
[1] 11561.3 11556.3
drop in deviance, indicators vs full factors
Analysis of Deviance Table

Model 1: damaged ~ compliant + PARENT_A + PARENT_DT + PARENT_GF + PARENT_GL + 
    PARENT_GT + PARENT_GT_DT + HYDROGROUP_A_D + HYDROGROUP_B + 
    HYDROGROUP_B_D + HYDROGROUP_C_D + vtrans_district_1 + vtrans_district_2 + 
    vtrans_district_3 + vtrans_district_4 + vtrans_district_6 + 
    vtrans_district_7 + vtrans_district_8 + vtrans_district_9 + 
    MeanSlope90m + StreamOrder + PercentImpervious_BaseLC_90m + 
    Precip + RoadGrade_mean_deg + Surface + Culvert_Any + Driveway_Count + 
    Mean_TopoConvergence_30m + road_miles_municipal + grand_list_equalized_muni_100M
Model 2: damaged ~ compliant + PARENT + HYDROGROUP + district_f + MeanSlope90m + 
    StreamOrder + PercentImpervious_BaseLC_90m + Precip + RoadGrade_mean_deg + 
    Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m + 
    road_miles_municipal + grand_list_equalized_muni_100M
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)  
1     59664      11499                       
2     59660      11486  4   13.067  0.01095 *
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

PART 2 WITH INDICATORS 2023 

Call:
glm2(formula = as.formula(f2), family = Gamma(link = "log"), 
    data = k)

Coefficients:
                                Estimate Std. Error t value Pr(>|t|)    
(Intercept)                     6.892481   0.927515   7.431 1.85e-13 ***
PARENT_DT                      -0.257744   0.243751  -1.057 0.290505    
PARENT_GF                      -0.155527   0.316727  -0.491 0.623470    
PARENT_GL                      -0.183410   0.413068  -0.444 0.657096    
PARENT_GT                      -0.972020   0.269414  -3.608 0.000319 ***
HYDROGROUP_A                    0.350879   0.308123   1.139 0.254994    
HYDROGROUP_B                    0.429700   0.238637   1.801 0.071970 .  
HYDROGROUP_C                    0.266466   0.178967   1.489 0.136732    
vtrans_district_1              -0.311720   0.739212  -0.422 0.673313    
vtrans_district_2              -1.211042   0.618879  -1.957 0.050562 .  
vtrans_district_3               0.835206   0.643958   1.297 0.194845    
vtrans_district_4              -0.655653   0.637585  -1.028 0.303965    
vtrans_district_5               0.366209   0.730441   0.501 0.616200    
vtrans_district_6               0.110851   0.613144   0.181 0.856557    
vtrans_district_7               0.166796   0.701788   0.238 0.812168    
vtrans_district_8               0.204565   0.755834   0.271 0.786701    
MeanSlope90m                    0.072319   0.016010   4.517 6.79e-06 ***
StreamCrossing                  0.090595   0.256598   0.353 0.724094    
StreamOrder                     0.219385   0.103655   2.116 0.034475 *  
K_final                         0.937063   0.581094   1.613 0.107056    
PercentImpervious_BaseLC_90m    2.373979   1.111955   2.135 0.032934 *  
Precip                          0.226436   0.057932   3.909 9.72e-05 ***
SurfaceUnpaved                 -0.366277   0.173967  -2.105 0.035427 *  
Culvert_Any                    -0.135023   0.132383  -1.020 0.307925    
Driveway_Count                 -0.262524   0.085405  -3.074 0.002153 ** 
Mean_TopoConvergence_30m        0.155423   0.059959   2.592 0.009635 ** 
road_miles_municipal            0.007507   0.003191   2.352 0.018795 *  
grand_list_equalized_muni_100M -0.002626   0.016701  -0.157 0.875100    
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

(Dispersion parameter for Gamma family taken to be 5.039538)

    Null deviance: 4936.4  on 1456  degrees of freedom
Residual deviance: 3630.1  on 1429  degrees of freedom
AIC: 32003

Number of Fisher Scoring iterations: 16

                               cost_ratio   lo_95    hi_95
(Intercept)                       984.842 280.519 3457.570
PARENT_DT                           0.773   0.435    1.371
PARENT_GF                           0.856   0.463    1.584
PARENT_GL                           0.832   0.362    1.913
PARENT_GT                           0.378   0.234    0.612
HYDROGROUP_A                        1.420   0.931    2.167
HYDROGROUP_B                        1.537   0.957    2.467
HYDROGROUP_C                        1.305   0.903    1.887
vtrans_district_1                   0.732   0.366    1.463
vtrans_district_2                   0.298   0.184    0.481
vtrans_district_3                   2.305   1.375    3.865
vtrans_district_4                   0.519   0.317    0.850
vtrans_district_5                   1.442   0.726    2.865
vtrans_district_6                   1.117   0.653    1.911
vtrans_district_7                   1.182   0.681    2.051
vtrans_district_8                   1.227   0.342    4.405
MeanSlope90m                        1.075   1.039    1.112
StreamCrossing                      1.095   0.448    2.673
StreamOrder                         1.245   0.820    1.892
K_final                             2.552   0.694    9.385
PercentImpervious_BaseLC_90m       10.740   0.277  417.093
Precip                              1.254   1.089    1.444
SurfaceUnpaved                      0.693   0.411    1.171
Culvert_Any                         0.874   0.673    1.134
Driveway_Count                      0.769   0.650    0.910
Mean_TopoConvergence_30m            1.168   1.042    1.310
road_miles_municipal                1.008   0.999    1.016
grand_list_equalized_muni_100M      0.997   0.965    1.031
deviance R2, indicators / full factors
[1] 0.2646 0.2666
AIC, indicators / full factors
[1] 32003.0 32011.7
F test, indicators vs full factors
Analysis of Deviance Table

Model 1: cost_dedup ~ PARENT_DT + PARENT_GF + PARENT_GL + PARENT_GT + 
    HYDROGROUP_A + HYDROGROUP_B + HYDROGROUP_C + vtrans_district_1 + 
    vtrans_district_2 + vtrans_district_3 + vtrans_district_4 + 
    vtrans_district_5 + vtrans_district_6 + vtrans_district_7 + 
    vtrans_district_8 + MeanSlope90m + StreamCrossing + StreamOrder + 
    K_final + PercentImpervious_BaseLC_90m + Precip + Surface + 
    Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m + 
    road_miles_municipal + grand_list_equalized_muni_100M
Model 2: cost_dedup ~ PARENT + HYDROGROUP + district_f + MeanSlope90m + 
    StreamCrossing + StreamOrder + K_final + PercentImpervious_BaseLC_90m + 
    Precip + Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m + 
    road_miles_municipal + grand_list_equalized_muni_100M
  Resid. Df Resid. Dev Df Deviance      F Pr(>F)
1      1429     3630.1                          
2      1422     3620.2  7   9.9296 0.2816 0.9613

AVOIDED COST 2023 
statewide avoided cost, indicators
[1] 11437970
statewide avoided cost, full factors
[1] 11492123
clustered 95% interval, indicators
[1]   664120 22211819
ten sample towns, indicators
 Brattleboro      Corinth  Wallingford     Hardwick West Windsor      Wolcott 
       41819        59907       114741        67863          830       158196 
  Washington     Richmond     Stamford   Starksboro 
      189822        14261         7962        43047 
```

## Output 2024

```
PART 1 SELECTION 2024 -- from the crosstab: |adjusted residual| >= 2, expected count >= 5, at least one damaged

 PARENT 
 category segments damaged expected adj_resid indicator
        A     1005      75     45.3      4.61      TRUE
       DT     8535     418    384.9      2.19      TRUE
    DT/GT      134      11      6.0      2.07      TRUE
       GF     2911     157    131.3      2.46      TRUE
    GF/GL      104       3      4.7     -0.80     FALSE
       GL     3373     100    152.1     -4.69      TRUE
       GT     5674     237    255.9     -1.40     FALSE
    GT/DT      715      11     32.2     -3.89      TRUE
  GT/GT_R       91       2      4.1     -1.07     FALSE
        M       29       4      1.3      2.41     FALSE

 HYDROGROUP 
  category segments damaged expected adj_resid indicator
         A     2404     137    108.4      2.97      TRUE
       A/D      362      16     16.3     -0.08     FALSE
         B     3043     166    137.2      2.70      TRUE
       B/D      824      48     37.2      1.85     FALSE
         C     4047     146    182.5     -3.05      TRUE
       C/D     3520     190    158.8      2.76      TRUE
         D     8156     306    367.9     -4.13      TRUE
 not rated      215       9      9.7     -0.23     FALSE

 vtrans_district 
 category segments damaged expected adj_resid indicator
        3      137       0      6.2     -2.55     FALSE
        5     5498     231    248.0     -1.27     FALSE
        6     6463     539    291.5     17.56      TRUE
        7     7243     209    326.7     -8.08      TRUE
        8      174       0      7.8     -2.88     FALSE
        9     3056      39    137.8     -9.26      TRUE

Part 1 indicators: PARENT_A, PARENT_DT, PARENT_DT_GT, PARENT_GF, PARENT_GL, PARENT_GT_DT, HYDROGROUP_A, HYDROGROUP_B, HYDROGROUP_C, HYDROGROUP_C_D, HYDROGROUP_D, vtrans_district_6, vtrans_district_7, vtrans_district_9 

PART 2 SELECTION 2024 -- from the ANOVA: categories in at least one Bonferroni pair with Sig. < .05

 PARENT   F = 1.491  p = 0.146 
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

 HYDROGROUP   F = 1.701  p = 0.105 
  category   n mean_ln_cost indicator
         A 137        9.535     FALSE
       A/D  16        9.352     FALSE
         B 166        9.695     FALSE
       B/D  48       10.042     FALSE
         C 146        9.504     FALSE
       C/D 190        9.541     FALSE
         D 306        9.324     FALSE
 not rated   9        9.799     FALSE

 vtrans_district   F = 16.155  p = 2.91e-10 
 category   n mean_ln_cost indicator
        5 231        8.916      TRUE
        6 539        9.633      TRUE
        7 209        9.903      TRUE
        9  39        9.420     FALSE

Part 2 indicators: vtrans_district_5, vtrans_district_6, vtrans_district_7 

PART 1 WITH INDICATORS 2024 

Call:
glm(formula = as.formula(f1), family = binomial, data = d)

Coefficients:
                                Estimate Std. Error z value Pr(>|z|)    
(Intercept)                    -9.517773   0.508039 -18.734  < 2e-16 ***
compliantCompliant             -0.238145   0.077115  -3.088  0.00201 ** 
PARENT_A                        0.092025   0.187518   0.491  0.62360    
PARENT_DT                       0.539394   0.136586   3.949 7.84e-05 ***
PARENT_DT_GT                    0.755905   0.351927   2.148  0.03172 *  
PARENT_GF                      -0.098645   0.211795  -0.466  0.64139    
PARENT_GL                      -0.086536   0.149193  -0.580  0.56190    
PARENT_GT_DT                   -1.257483   0.328695  -3.826  0.00013 ***
HYDROGROUP_A                    0.256561   0.199742   1.284  0.19898    
HYDROGROUP_B                    0.160613   0.187417   0.857  0.39145    
HYDROGROUP_C                   -0.364531   0.200321  -1.820  0.06880 .  
HYDROGROUP_C_D                 -0.416218   0.227246  -1.832  0.06701 .  
HYDROGROUP_D                   -0.547710   0.216050  -2.535  0.01124 *  
vtrans_district_6               0.580989   0.098285   5.911 3.39e-09 ***
vtrans_district_7              -0.486490   0.116871  -4.163 3.15e-05 ***
vtrans_district_9              -1.026873   0.189639  -5.415 6.13e-08 ***
MeanSlope90m                    0.078224   0.009849   7.942 1.99e-15 ***
StreamOrder                     0.233770   0.032269   7.244 4.34e-13 ***
PercentImpervious_BaseLC_90m    1.221736   0.541792   2.255  0.02413 *  
Precip                          0.443932   0.063248   7.019 2.24e-12 ***
RoadGrade_mean_deg              0.080814   0.018957   4.263 2.02e-05 ***
SurfaceUnpaved                  0.560769   0.095735   5.858 4.70e-09 ***
Culvert_Any                     0.616133   0.073011   8.439  < 2e-16 ***
Driveway_Count                  0.089044   0.042229   2.109  0.03498 *  
Mean_TopoConvergence_30m        0.383597   0.034984  10.965  < 2e-16 ***
road_miles_municipal           -0.005076   0.001767  -2.873  0.00407 ** 
grand_list_equalized_muni_100M  0.005401   0.003652   1.479  0.13910    
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 8298.6  on 22570  degrees of freedom
Residual deviance: 7289.8  on 22544  degrees of freedom
AIC: 7343.8

Number of Fisher Scoring iterations: 7

                               odds_ratio lo_95  hi_95
(Intercept)                         0.000 0.000  0.001
compliantCompliant                  0.788 0.625  0.994
PARENT_A                            1.096 0.699  1.719
PARENT_DT                           1.715 1.198  2.456
PARENT_DT_GT                        2.130 1.237  3.665
PARENT_GF                           0.906 0.524  1.565
PARENT_GL                           0.917 0.594  1.416
PARENT_GT_DT                        0.284 0.113  0.715
HYDROGROUP_A                        1.292 0.740  2.256
HYDROGROUP_B                        1.174 0.737  1.871
HYDROGROUP_C                        0.695 0.411  1.174
HYDROGROUP_C_D                      0.660 0.420  1.036
HYDROGROUP_D                        0.578 0.346  0.968
vtrans_district_6                   1.788 1.070  2.988
vtrans_district_7                   0.615 0.304  1.242
vtrans_district_9                   0.358 0.135  0.947
MeanSlope90m                        1.081 1.048  1.116
StreamOrder                         1.263 1.179  1.354
PercentImpervious_BaseLC_90m        3.393 1.068 10.783
Precip                              1.559 1.063  2.285
RoadGrade_mean_deg                  1.084 1.034  1.137
SurfaceUnpaved                      1.752 1.257  2.442
Culvert_Any                         1.852 1.508  2.274
Driveway_Count                      1.093 0.963  1.241
Mean_TopoConvergence_30m            1.468 1.340  1.607
road_miles_municipal                0.995 0.983  1.007
grand_list_equalized_muni_100M      1.005 0.994  1.017
AUC, indicators
[1] 0.7797
AUC, full factors
[1] 0.7782
AIC, indicators / full factors
[1] 7343.8 7361.6
drop in deviance, indicators vs full factors
Analysis of Deviance Table

Model 1: damaged ~ compliant + PARENT_A + PARENT_DT + PARENT_DT_GT + PARENT_GF + 
    PARENT_GL + PARENT_GT_DT + HYDROGROUP_A + HYDROGROUP_B + 
    HYDROGROUP_C + HYDROGROUP_C_D + HYDROGROUP_D + vtrans_district_6 + 
    vtrans_district_7 + vtrans_district_9 + MeanSlope90m + StreamOrder + 
    PercentImpervious_BaseLC_90m + Precip + RoadGrade_mean_deg + 
    Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m + 
    road_miles_municipal + grand_list_equalized_muni_100M
Model 2: damaged ~ compliant + PARENT + HYDROGROUP + district_f + MeanSlope90m + 
    StreamOrder + PercentImpervious_BaseLC_90m + Precip + RoadGrade_mean_deg + 
    Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m + 
    road_miles_municipal + grand_list_equalized_muni_100M
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)
1     22544     7289.8                     
2     22539     7297.6  5   -7.752         

PART 2 WITH INDICATORS 2024 

Call:
glm2(formula = as.formula(f2), family = Gamma(link = "log"), 
    data = k)

Coefficients:
                                Estimate Std. Error t value Pr(>|t|)    
(Intercept)                    10.178221   1.014552  10.032  < 2e-16 ***
vtrans_district_5               0.015853   0.437085   0.036  0.97107    
vtrans_district_6               0.717658   0.406952   1.763  0.07812 .  
vtrans_district_7               1.963994   0.452614   4.339 1.57e-05 ***
MeanSlope90m                    0.089585   0.021781   4.113 4.23e-05 ***
StreamCrossing                 -0.645655   0.316656  -2.039  0.04171 *  
StreamOrder                     0.656045   0.131543   4.987 7.21e-07 ***
K_final                        -0.316371   0.527525  -0.600  0.54882    
PercentImpervious_BaseLC_90m    1.671158   1.320721   1.265  0.20604    
Precip                         -0.184231   0.150615  -1.223  0.22155    
SurfaceUnpaved                 -0.304925   0.209671  -1.454  0.14618    
Culvert_Any                     0.227386   0.169005   1.345  0.17879    
Driveway_Count                 -0.113912   0.102681  -1.109  0.26753    
Mean_TopoConvergence_30m        0.048835   0.074718   0.654  0.51353    
road_miles_municipal           -0.013866   0.004529  -3.061  0.00226 ** 
grand_list_equalized_muni_100M  0.014429   0.008779   1.643  0.10060    
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

(Dispersion parameter for Gamma family taken to be 5.595851)

    Null deviance: 3157.1  on 1017  degrees of freedom
Residual deviance: 2405.4  on 1002  degrees of freedom
AIC: 23537

Number of Fisher Scoring iterations: 14

                               cost_ratio    lo_95      hi_95
(Intercept)                     26323.608 5410.284 128076.885
vtrans_district_5                   1.016    0.587      1.759
vtrans_district_6                   2.050    1.427      2.945
vtrans_district_7                   7.128    4.077     12.461
MeanSlope90m                        1.094    1.050      1.139
StreamCrossing                      0.524    0.269      1.023
StreamOrder                         1.927    1.439      2.582
K_final                             0.729    0.307      1.729
PercentImpervious_BaseLC_90m        5.318    0.390     72.529
Precip                              0.832    0.606      1.141
SurfaceUnpaved                      0.737    0.553      0.982
Culvert_Any                         1.255    0.920      1.712
Driveway_Count                      0.892    0.752      1.059
Mean_TopoConvergence_30m            1.050    0.940      1.173
road_miles_municipal                0.986    0.975      0.998
grand_list_equalized_muni_100M      1.015    1.001      1.028
deviance R2, indicators / full factors
[1] 0.2381 0.2662
AIC, indicators / full factors
[1] 23537.2 23519.5
F test, indicators vs full factors
Analysis of Deviance Table

Model 1: cost_dedup ~ vtrans_district_5 + vtrans_district_6 + vtrans_district_7 + 
    MeanSlope90m + StreamCrossing + StreamOrder + K_final + PercentImpervious_BaseLC_90m + 
    Precip + Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m + 
    road_miles_municipal + grand_list_equalized_muni_100M
Model 2: cost_dedup ~ PARENT + HYDROGROUP + district_f + MeanSlope90m + 
    StreamCrossing + StreamOrder + K_final + PercentImpervious_BaseLC_90m + 
    Precip + Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m + 
    road_miles_municipal + grand_list_equalized_muni_100M
  Resid. Df Resid. Dev Df Deviance      F Pr(>F)
1      1002     2405.4                          
2       986     2316.7 16    88.71 1.1607 0.2939

AVOIDED COST 2024 
statewide avoided cost, indicators
[1] 12265670
statewide avoided cost, full factors
[1] 12714930
clustered 95% interval, indicators
[1] -2215460 26746800
ten sample towns, indicators
      <NA>       <NA>       <NA>   Hardwick       <NA>    Wolcott       <NA> 
        NA         NA         NA     226097         NA     170552         NA 
  Richmond       <NA> Starksboro 
    204163         NA     149107 
```
