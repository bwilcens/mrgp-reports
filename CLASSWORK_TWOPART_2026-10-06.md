# Two-Part Model, Classwork Style

The selected two-part model written as a step-by-step exercise script: Part 1 logistic regression for whether a segment was damaged (Exercise 7), town-clustered standard errors, Part 2 Gamma regression with a log link for the repair cost, the avoided cost from recoding every compliant segment to Does Not Meet with its clustered interval, and the same model on the towns holding both FEMA and VTrans records run on the pooled de-duplicated cost and on the VTrans records alone. Graded segments, Precip >= 3.51 in, Averill excluded; cost_dedup, uncapped. Script run once with STORM = 2023 and once with STORM = 2024.

## Script: twopart_classwork.R

```r
# =============================================================================
# The two-part model of road damage, step by step
#
# Data:  ../mrgp-roads-core/twopart_simple_2026-09-17/data/analysis_<STORM>_2026-09-27.csv
#        one row per road segment; compliance-graded segments with at least 3.51 in of rain,
#        Averill excluded (it has no equalized grand list); damaged = 1 if cost_dedup > 0
#        compliant        "Compliant" (Fully or Partially Meets) or "Does Not Meet"
#        cost_dedup       repair cost in dollars, FEMA and VTrans records de-duplicated, uncapped
#        vt_cost          the VTrans records alone (STEP 5)
#        the covariates   slope, stream order, parent material, hydrologic group, imperviousness,
#                         rainfall, road grade, surface, culvert, driveways, topographic wetness,
#                         VTrans district, town road miles, town grand list
#
# Part 1 is a logistic regression: did the segment flood? (Exercise 7)
# Part 2 is a Gamma regression with a log link: given that it flooded, what did it cost?
# The avoided cost recodes every compliant segment to Does Not Meet and adds up the change in
# (probability of damage) x (predicted cost). Compliance is in Part 1 only.
#
# Run one line at a time with Ctrl+Enter. Set STORM to 2023 or 2024. Nothing is written to disk.
# =============================================================================


# -----------------------------------------------------------------------------
# STEP 0.  Point R at the data and read it in
# -----------------------------------------------------------------------------
STORM <- 2023

library(glm2)        # a steadier fitting routine for the Gamma model (plain glm can stall on it)
library(sandwich)    # town-clustered standard errors

if (!dir.exists("data") && dir.exists("../data")) setwd("..")
roads <- read.csv(paste0("../mrgp-roads-core/twopart_simple_2026-09-17/data/analysis_", STORM, "_2026-09-27.csv"))
roads <- roads[roads$Precip >= 3.51 & roads$Town != "Averill", ]
roads$damaged   <- as.integer(roads$cost_dedup > 0)
roads$compliant <- relevel(factor(roads$compliant), ref = "Does Not Meet")   # Does Not Meet is the reference

# a district with no damaged segment cannot have a coefficient: pool it into the busiest district
damaged_by_district <- tapply(roads$damaged, roads$vtrans_district, sum)
busiest <- names(which.max(damaged_by_district))
roads$vtrans_district[roads$vtrans_district %in% names(damaged_by_district)[damaged_by_district == 0]] <- busiest
roads$vtrans_district <- relevel(factor(roads$vtrans_district), ref = busiest)

cat("--- Number of segments ---\n");        print(nrow(roads))
cat("--- Number damaged ---\n");            print(sum(roads$damaged))
cat("--- Number of towns ---\n");           print(length(unique(roads$Town)))
cat("--- Compliance status ---\n");         print(table(roads$compliant))


# -----------------------------------------------------------------------------
# STEP 1.  Part 1: logistic regression for whether the segment was damaged (Exercise 7)
# -----------------------------------------------------------------------------
# This replaces  Analyze > Regression > Binary Logistic  with damaged as the dependent variable
# and compliant entered as a categorical covariate (Indicator contrast, reference = Does Not Meet).
cat("=== STEP 1: Part 1, did the segment flood? ===\n\n")

part1 <- glm(damaged ~ compliant + MeanSlope90m + StreamOrder + PARENT + HYDROGROUP +
               PercentImpervious_BaseLC_90m + Precip + RoadGrade_mean_deg + Surface + Culvert_Any +
               Driveway_Count + Mean_TopoConvergence_30m + vtrans_district +
               road_miles_municipal + grand_list_equalized_muni_100M,
             family = binomial, data = roads)

cat("--- The fitted model (coefficients are log-odds) ---\n")
print(summary(part1))

cat("\n--- Odds ratios: exp(coefficient) ---\n")
print(round(exp(coef(part1)), 3))

# Drop in deviance for the compliance term: the same model without it, then the change in deviance
# is a chi-square on 1 df (Exercise 7, Question on nested models).
cat("\n--- Does compliance belong in the model? Drop in deviance ---\n")
part1_no_compliance <- update(part1, . ~ . - compliant)
print(anova(part1_no_compliance, part1, test = "Chisq"))

# AUC: the chance that a randomly chosen damaged segment gets a higher fitted probability than a
# randomly chosen undamaged one; 0.5 is a coin flip. Percent correctly classified is not used
# because only 2 to 4 segments in 100 are damaged, so "predict nobody" scores 97 percent.
cat("\n--- AUC ---\n")
n1 <- sum(roads$damaged); n0 <- sum(roads$damaged == 0)
print(round((sum(rank(fitted(part1))[roads$damaged == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0), 4))


# -----------------------------------------------------------------------------
# STEP 2.  Standard errors clustered by town
# -----------------------------------------------------------------------------
# The regular standard errors treat every segment as an independent observation. Segments in one
# town share a road crew, a budget, soils and a position under the storm, so they are not. The
# clustered standard error adds up each town's contribution first and treats the towns as the
# independent units. The coefficient does not change; its standard error gets larger.
# (SPSS has no checkbox for this in Binary Logistic; it lives in Complex Samples or GEE.)
cat("=== STEP 2: town-clustered standard errors ===\n\n")

se_regular   <- sqrt(diag(vcov(part1)))
se_clustered <- sqrt(diag(vcovCL(part1, cluster = roads$Town, type = "HC1")))

cat("--- Compliance coefficient: regular vs clustered ---\n")
b <- coef(part1)[["compliantCompliant"]]
cat("  coefficient          :", round(b, 4), "\n")
cat("  regular SE           :", round(se_regular[["compliantCompliant"]], 4), "\n")
cat("  clustered SE         :", round(se_clustered[["compliantCompliant"]], 4), "\n")
cat("  odds ratio           :", round(exp(b), 3), "\n")
cat("  95% interval, regular:", round(exp(b + c(-1.96, 1.96) * se_regular[["compliantCompliant"]]), 3), "\n")
cat("  95% interval, clustered:", round(exp(b + c(-1.96, 1.96) * se_clustered[["compliantCompliant"]]), 3), "\n")

cat("\n--- Every coefficient as an odds ratio with its clustered 95% interval ---\n")
print(round(exp(cbind(odds_ratio = coef(part1),
                      lo_95 = coef(part1) - 1.96 * se_clustered,
                      hi_95 = coef(part1) + 1.96 * se_clustered)), 3))


# -----------------------------------------------------------------------------
# STEP 3.  Part 2: Gamma regression for the repair cost of the damaged segments
# -----------------------------------------------------------------------------
# This replaces  Analyze > Generalized Linear Models  with Type of Model = Gamma with log link,
# cost_dedup as the dependent variable, fitted on the damaged segments only. Compliance is not a
# predictor here; the cost model describes what a repair costs once it is needed.
cat("=== STEP 3: Part 2, what did it cost? ===\n\n")

damaged <- roads[roads$damaged == 1, ]
cat("--- Damaged segments ---\n");          print(nrow(damaged))
cat("--- Mean and median cost ---\n");      print(round(c(mean = mean(damaged$cost_dedup), median = median(damaged$cost_dedup))))

part2 <- glm2(cost_dedup ~ MeanSlope90m + StreamCrossing + StreamOrder + PARENT + K_final + HYDROGROUP +
                PercentImpervious_BaseLC_90m + Precip + Surface + Culvert_Any + Driveway_Count +
                Mean_TopoConvergence_30m + vtrans_district + road_miles_municipal + grand_list_equalized_muni_100M,
              family = Gamma(link = "log"), data = damaged)

cat("\n--- The fitted model (coefficients are on the log-cost scale) ---\n")
print(summary(part2))

# With a log link, exp(coefficient) is a cost ratio: 1.20 means 20 percent higher repair cost.
cat("\n--- Cost ratios with clustered 95% intervals ---\n")
se_clustered2 <- sqrt(diag(vcovCL(part2, cluster = damaged$Town, type = "HC1")))
print(round(exp(cbind(cost_ratio = coef(part2),
                      lo_95 = coef(part2) - 1.96 * se_clustered2,
                      hi_95 = coef(part2) + 1.96 * se_clustered2)), 3))

# Deviance R-squared: the share of the cost variation the model accounts for.
cat("\n--- Deviance R-squared ---\n")
print(round(1 - part2$deviance / part2$null.deviance, 4))
cat("--- Mean predicted / mean observed cost (1 = the model reproduces the total) ---\n")
print(round(mean(fitted(part2)) / mean(damaged$cost_dedup), 4))


# -----------------------------------------------------------------------------
# STEP 4.  The avoided cost: recode every compliant segment to Does Not Meet
# -----------------------------------------------------------------------------
# For every segment: probability of damage as rated, probability if it were rated Does Not Meet,
# and its predicted repair cost. The avoided cost of a compliant segment is
#   (probability if Does Not Meet - probability as rated) x predicted cost
# and the statewide figure is the sum over the compliant segments.
cat("=== STEP 4: the avoided cost ===\n\n")

recoded <- roads
recoded$compliant <- factor("Does Not Meet", levels = levels(roads$compliant))

p_as_rated   <- predict(part1, newdata = roads,   type = "response")
p_if_dnm     <- predict(part1, newdata = recoded, type = "response")
cost_hat     <- predict(part2, newdata = roads,   type = "response")
avoided      <- (p_if_dnm - p_as_rated) * cost_hat
is_compliant <- roads$compliant == "Compliant"

cat("--- Compliant segments recoded ---\n");          print(sum(is_compliant))
cat("--- Mean probability as rated, compliant ---\n"); print(round(mean(p_as_rated[is_compliant]), 4))
cat("--- Mean probability if Does Not Meet ---\n");    print(round(mean(p_if_dnm[is_compliant]), 4))
cat("--- Statewide avoided cost, dollars ---\n");      print(round(sum(avoided[is_compliant])))

cat("\n--- The ten sample towns (NA = not in this storm's exposed population) ---\n")
print(round(tapply(avoided[is_compliant], roads$Town[is_compliant], sum)[c("Brattleboro", "Corinth", "Wallingford",
      "Hardwick", "West Windsor", "Wolcott", "Washington", "Richmond", "Stamford", "Starksboro")]))

# The town-clustered interval for the statewide total. The total is a function of both parts'
# coefficients, so its standard error comes from the delta method: the gradient of the total with
# respect to every coefficient, and the joint clustered covariance of both parts, with each town's
# score contributions added up before squaring. This is the block the slide deck's intervals use.
cat("\n--- Town-clustered 95% interval for the statewide total ---\n")
k1 <- names(coef(part1)); k2 <- names(coef(part2))
X1 <- model.matrix(delete.response(terms(part1)), roads,   xlev = part1$xlevels)[, k1]
X0 <- model.matrix(delete.response(terms(part1)), recoded, xlev = part1$xlevels)[, k1]
Z  <- model.matrix(delete.response(terms(part2)), roads,   xlev = part2$xlevels)[, k2]
gradient <- c(colSums(((cost_hat * p_if_dnm * (1 - p_if_dnm)) * X0 - (cost_hat * p_as_rated * (1 - p_as_rated)) * X1)[is_compliant, ]),
              colSums(((p_if_dnm - p_as_rated) * cost_hat * Z)[is_compliant, ]))
B1 <- summary(part1)$cov.unscaled
w  <- weights(part2, "working"); rr <- residuals(part2, "working") * w
B2 <- summary(part2)$cov.unscaled * sum(rr^2) / sum(w)
towns <- sort(unique(roads$Town))
U1 <- matrix(0, length(towns), length(k1), dimnames = list(towns, k1)); s1 <- rowsum(estfun(part1), roads$Town);   U1[rownames(s1), ] <- s1
U2 <- matrix(0, length(towns), length(k2), dimnames = list(towns, k2)); s2 <- rowsum(estfun(part2), damaged$Town); U2[rownames(s2), ] <- s2
G1 <- length(towns); G2 <- length(unique(damaged$Town))
a1 <- G1 / (G1 - 1) * (nrow(roads) - 1) / (nrow(roads) - length(k1))
a2 <- G2 / (G2 - 1) * (nrow(damaged) - 1) / (nrow(damaged) - length(k2))
V_joint <- crossprod(cbind(sqrt(a1) * U1 %*% B1, sqrt(a2) * U2 %*% B2))
se_total <- sqrt(drop(t(gradient) %*% V_joint %*% gradient))
cat("  clustered SE :", round(se_total), "\n")
cat("  95% interval :", round(sum(avoided[is_compliant]) + c(-1.96, 1.96) * se_total), "\n")


# -----------------------------------------------------------------------------
# STEP 5.  The same model on the towns with both FEMA and VTrans records: pooled vs VTrans only
# -----------------------------------------------------------------------------
# Keep only the towns that hold at least one FEMA-recorded and one VTrans-recorded damage, so that
# the two runs below use identical segments. Run A counts a segment as damaged by the pooled,
# de-duplicated record (cost_dedup); run B by the VTrans records alone (vt_cost). Everything else is
# the same, so any difference is the record source. Intervals on the dollar totals for these runs
# are in both_source_towns.R; here the compliance odds ratio carries the clustered interval.
cat("=== STEP 5: both-source towns ===\n\n")

both <- roads[roads$Town %in% intersect(unique(roads$Town[roads$fema_cost > 0]), unique(roads$Town[roads$vt_cost > 0])), ]
cat("--- Towns with both record sources ---\n");  print(length(unique(both$Town)))
cat("--- Segments ---\n");                        print(nrow(both))

# ---- Run A: pooled, de-duplicated ----
cat("\n--- Run A: pooled de-duplicated cost ---\n")
both$damaged <- as.integer(both$cost_dedup > 0)
both$cost    <- both$cost_dedup
cat("damaged segments:", sum(both$damaged), "\n")
# categories with no damaged segment in this smaller population are pooled into the busiest one
n <- tapply(both$damaged, both$vtrans_district, sum); both$vtrans_district <- as.character(both$vtrans_district)
both$vtrans_district[both$vtrans_district %in% names(n)[n == 0]] <- names(which.max(n)); both$vtrans_district <- relevel(factor(both$vtrans_district), ref = names(which.max(n)))
n <- tapply(both$damaged, both$PARENT, sum);     both$PARENT[both$PARENT %in% names(n)[n == 0]] <- names(which.max(n));         both$PARENT <- factor(both$PARENT)
n <- tapply(both$damaged, both$HYDROGROUP, sum); both$HYDROGROUP[both$HYDROGROUP %in% names(n)[n == 0]] <- names(which.max(n)); both$HYDROGROUP <- factor(both$HYDROGROUP)

p1A <- glm(damaged ~ compliant + MeanSlope90m + StreamOrder + PARENT + HYDROGROUP + PercentImpervious_BaseLC_90m + Precip +
             RoadGrade_mean_deg + Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m + vtrans_district +
             road_miles_municipal + grand_list_equalized_muni_100M, family = binomial, data = both)
dA  <- both[both$damaged == 1, ]
p2A <- glm2(cost ~ MeanSlope90m + StreamCrossing + StreamOrder + PARENT + K_final + HYDROGROUP + PercentImpervious_BaseLC_90m +
              Precip + Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m + vtrans_district +
              road_miles_municipal + grand_list_equalized_muni_100M, family = Gamma(link = "log"), data = dA)
bA  <- coef(p1A)[["compliantCompliant"]]; seA <- sqrt(vcovCL(p1A, cluster = both$Town, type = "HC1")["compliantCompliant", "compliantCompliant"])
n1 <- sum(both$damaged); n0 <- sum(both$damaged == 0)
aucA <- (sum(rank(fitted(p1A))[both$damaged == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
recodedA <- both; recodedA$compliant <- factor("Does Not Meet", levels = levels(both$compliant))
avoidedA <- sum(((predict(p1A, newdata = recodedA, type = "response") - predict(p1A, newdata = both, type = "response")) *
                 predict(p2A, newdata = both, type = "response"))[both$compliant == "Compliant"])
cat("compliance odds ratio (clustered 95%):", round(exp(bA), 3), round(exp(bA + c(-1.96, 1.96) * seA), 3), "\n")
cat("AUC:", round(aucA, 4), "  deviance R2:", round(1 - p2A$deviance / p2A$null.deviance, 4), "\n")
cat("avoided cost, these towns:", round(avoidedA), "\n")

# ---- Run B: VTrans records only ----
cat("\n--- Run B: VTrans records only ---\n")
both <- roads[roads$Town %in% intersect(unique(roads$Town[roads$fema_cost > 0]), unique(roads$Town[roads$vt_cost > 0])), ]
both$damaged <- as.integer(both$vt_cost > 0)
both$cost    <- both$vt_cost
cat("damaged segments:", sum(both$damaged), "\n")
n <- tapply(both$damaged, both$vtrans_district, sum); both$vtrans_district <- as.character(both$vtrans_district)
both$vtrans_district[both$vtrans_district %in% names(n)[n == 0]] <- names(which.max(n)); both$vtrans_district <- relevel(factor(both$vtrans_district), ref = names(which.max(n)))
n <- tapply(both$damaged, both$PARENT, sum);     both$PARENT[both$PARENT %in% names(n)[n == 0]] <- names(which.max(n));         both$PARENT <- factor(both$PARENT)
n <- tapply(both$damaged, both$HYDROGROUP, sum); both$HYDROGROUP[both$HYDROGROUP %in% names(n)[n == 0]] <- names(which.max(n)); both$HYDROGROUP <- factor(both$HYDROGROUP)

p1B <- glm(damaged ~ compliant + MeanSlope90m + StreamOrder + PARENT + HYDROGROUP + PercentImpervious_BaseLC_90m + Precip +
             RoadGrade_mean_deg + Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m + vtrans_district +
             road_miles_municipal + grand_list_equalized_muni_100M, family = binomial, data = both)
dB  <- both[both$damaged == 1, ]
p2B <- glm2(cost ~ MeanSlope90m + StreamCrossing + StreamOrder + PARENT + K_final + HYDROGROUP + PercentImpervious_BaseLC_90m +
              Precip + Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m + vtrans_district +
              road_miles_municipal + grand_list_equalized_muni_100M, family = Gamma(link = "log"), data = dB)
bB  <- coef(p1B)[["compliantCompliant"]]; seB <- sqrt(vcovCL(p1B, cluster = both$Town, type = "HC1")["compliantCompliant", "compliantCompliant"])
n1 <- sum(both$damaged); n0 <- sum(both$damaged == 0)
aucB <- (sum(rank(fitted(p1B))[both$damaged == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
recodedB <- both; recodedB$compliant <- factor("Does Not Meet", levels = levels(both$compliant))
avoidedB <- sum(((predict(p1B, newdata = recodedB, type = "response") - predict(p1B, newdata = both, type = "response")) *
                 predict(p2B, newdata = both, type = "response"))[both$compliant == "Compliant"])
cat("compliance odds ratio (clustered 95%):", round(exp(bB), 3), round(exp(bB + c(-1.96, 1.96) * seB), 3), "\n")
cat("AUC:", round(aucB, 4), "  deviance R2:", round(1 - p2B$deviance / p2B$null.deviance, 4), "\n")
cat("avoided cost, these towns:", round(avoidedB), "\n")

# ---- Side by side ----
cat("\n--- Side by side, same towns and segments ---\n")
print(data.frame(run = c("A pooled de-duplicated", "B VTrans only"),
                 damaged = c(sum(dA$damaged), sum(dB$damaged)),
                 odds_ratio = round(c(exp(bA), exp(bB)), 3),
                 lo_95 = round(c(exp(bA - 1.96 * seA), exp(bB - 1.96 * seB)), 3),
                 hi_95 = round(c(exp(bA + 1.96 * seA), exp(bB + 1.96 * seB)), 3),
                 AUC = round(c(aucA, aucB), 3),
                 dev_R2 = round(c(1 - p2A$deviance / p2A$null.deviance, 1 - p2B$deviance / p2B$null.deviance), 3),
                 avoided = round(c(avoidedA, avoidedB))), row.names = FALSE)


# =============================================================================
# END
#
# To report:
#   Part 1 -- the compliance odds ratio with its clustered interval, the drop-in-deviance test for
#   the compliance term, and AUC. Part 2 -- deviance R-squared and the predicted/observed ratio.
#   The avoided cost statewide with its clustered interval, and for the ten sample towns.
#   The both-source comparison: how the odds ratio and avoided cost change when only VTrans
#   records count as damage, on identical segments.
# =============================================================================
```

## Output 2023

```
--- Number of segments ---
[1] 59695
--- Number damaged ---
[1] 1457
--- Number of towns ---
[1] 204
--- Compliance status ---

Does Not Meet     Compliant 
         9520         50175 
=== STEP 1: Part 1, did the segment flood? ===

--- The fitted model (coefficients are log-odds) ---

Call:
glm(formula = damaged ~ compliant + MeanSlope90m + StreamOrder + 
    PARENT + HYDROGROUP + PercentImpervious_BaseLC_90m + Precip + 
    RoadGrade_mean_deg + Surface + Culvert_Any + Driveway_Count + 
    Mean_TopoConvergence_30m + vtrans_district + road_miles_municipal + 
    grand_list_equalized_muni_100M, family = binomial, data = roads)

Coefficients:
                                Estimate Std. Error z value Pr(>|z|)    
(Intercept)                    -7.466475   0.353294 -21.134  < 2e-16 ***
compliantCompliant             -0.224594   0.066154  -3.395 0.000686 ***
MeanSlope90m                    0.064833   0.007613   8.516  < 2e-16 ***
StreamOrder                     0.287403   0.026874  10.694  < 2e-16 ***
PARENTDT                        0.032226   0.184585   0.175 0.861407    
PARENTDT/GT                     0.378696   0.416743   0.909 0.363507    
PARENTGF                       -0.557078   0.171563  -3.247 0.001166 ** 
PARENTGL                       -0.080598   0.220206  -0.366 0.714356    
PARENTGT                       -0.481129   0.165242  -2.912 0.003595 ** 
PARENTGT/DT                    -1.082519   0.378519  -2.860 0.004238 ** 
PARENTM                        -2.004703   1.064968  -1.882 0.059781 .  
HYDROGROUPA/D                  -0.146623   0.374130  -0.392 0.695130    
HYDROGROUPB                    -0.068516   0.147348  -0.465 0.641935    
HYDROGROUPB/D                  -0.225437   0.196462  -1.147 0.251182    
HYDROGROUPC                    -0.364976   0.160233  -2.278 0.022739 *  
HYDROGROUPC/D                  -0.509089   0.187974  -2.708 0.006763 ** 
HYDROGROUPD                    -0.279821   0.177497  -1.576 0.114914    
HYDROGROUPnot rated             0.761443   0.368852   2.064 0.038984 *  
PercentImpervious_BaseLC_90m    1.514083   0.466348   3.247 0.001168 ** 
Precip                          0.384915   0.025079  15.348  < 2e-16 ***
RoadGrade_mean_deg              0.049024   0.015068   3.254 0.001140 ** 
SurfaceUnpaved                  0.392556   0.080622   4.869 1.12e-06 ***
Culvert_Any                     0.510180   0.059454   8.581  < 2e-16 ***
Driveway_Count                  0.005928   0.037783   0.157 0.875322    
Mean_TopoConvergence_30m        0.258015   0.026016   9.918  < 2e-16 ***
vtrans_district1               -2.304813   0.184375 -12.501  < 2e-16 ***
vtrans_district3               -1.477504   0.102912 -14.357  < 2e-16 ***
vtrans_district4               -2.031914   0.104894 -19.371  < 2e-16 ***
vtrans_district5               -1.371802   0.186370  -7.361 1.83e-13 ***
vtrans_district6               -0.803284   0.079783 -10.068  < 2e-16 ***
vtrans_district7               -1.682827   0.170807  -9.852  < 2e-16 ***
vtrans_district8               -2.386332   0.210167 -11.354  < 2e-16 ***
vtrans_district9               -3.553555   0.275906 -12.880  < 2e-16 ***
road_miles_municipal            0.001585   0.001521   1.042 0.297406    
grand_list_equalized_muni_100M -0.033924   0.006725  -5.045 4.54e-07 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 13697  on 59694  degrees of freedom
Residual deviance: 11486  on 59660  degrees of freedom
AIC: 11556

Number of Fisher Scoring iterations: 8

--- Odds ratios: exp(coefficient) ---
                   (Intercept)             compliantCompliant 
                         0.001                          0.799 
                  MeanSlope90m                    StreamOrder 
                         1.067                          1.333 
                      PARENTDT                    PARENTDT/GT 
                         1.033                          1.460 
                      PARENTGF                       PARENTGL 
                         0.573                          0.923 
                      PARENTGT                    PARENTGT/DT 
                         0.618                          0.339 
                       PARENTM                  HYDROGROUPA/D 
                         0.135                          0.864 
                   HYDROGROUPB                  HYDROGROUPB/D 
                         0.934                          0.798 
                   HYDROGROUPC                  HYDROGROUPC/D 
                         0.694                          0.601 
                   HYDROGROUPD            HYDROGROUPnot rated 
                         0.756                          2.141 
  PercentImpervious_BaseLC_90m                         Precip 
                         4.545                          1.469 
            RoadGrade_mean_deg                 SurfaceUnpaved 
                         1.050                          1.481 
                   Culvert_Any                 Driveway_Count 
                         1.666                          1.006 
      Mean_TopoConvergence_30m               vtrans_district1 
                         1.294                          0.100 
              vtrans_district3               vtrans_district4 
                         0.228                          0.131 
              vtrans_district5               vtrans_district6 
                         0.254                          0.448 
              vtrans_district7               vtrans_district8 
                         0.186                          0.092 
              vtrans_district9           road_miles_municipal 
                         0.029                          1.002 
grand_list_equalized_muni_100M 
                         0.967 

--- Does compliance belong in the model? Drop in deviance ---
Analysis of Deviance Table

Model 1: damaged ~ MeanSlope90m + StreamOrder + PARENT + HYDROGROUP + 
    PercentImpervious_BaseLC_90m + Precip + RoadGrade_mean_deg + 
    Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m + 
    vtrans_district + road_miles_municipal + grand_list_equalized_muni_100M
Model 2: damaged ~ compliant + MeanSlope90m + StreamOrder + PARENT + HYDROGROUP + 
    PercentImpervious_BaseLC_90m + Precip + RoadGrade_mean_deg + 
    Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m + 
    vtrans_district + road_miles_municipal + grand_list_equalized_muni_100M
  Resid. Df Resid. Dev Df Deviance  Pr(>Chi)    
1     59661      11498                          
2     59660      11486  1   11.228 0.0008058 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

--- AUC ---
[1] 0.8244
=== STEP 2: town-clustered standard errors ===

--- Compliance coefficient: regular vs clustered ---
  coefficient          : -0.2246 
  regular SE           : 0.0662 
  clustered SE         : 0.0852 
  odds ratio           : 0.799 
  95% interval, regular: 0.702 0.909 
  95% interval, clustered: 0.676 0.944 

--- Every coefficient as an odds ratio with its clustered 95% interval ---
                               odds_ratio lo_95  hi_95
(Intercept)                         0.001 0.000  0.002
compliantCompliant                  0.799 0.676  0.944
MeanSlope90m                        1.067 1.043  1.091
StreamOrder                         1.333 1.250  1.421
PARENTDT                            1.033 0.695  1.535
PARENTDT/GT                         1.460 0.501  4.256
PARENTGF                            0.573 0.400  0.821
PARENTGL                            0.923 0.590  1.442
PARENTGT                            0.618 0.424  0.901
PARENTGT/DT                         0.339 0.120  0.960
PARENTM                             0.135 0.011  1.650
HYDROGROUPA/D                       0.864 0.420  1.776
HYDROGROUPB                         0.934 0.626  1.394
HYDROGROUPB/D                       0.798 0.503  1.268
HYDROGROUPC                         0.694 0.463  1.042
HYDROGROUPC/D                       0.601 0.374  0.967
HYDROGROUPD                         0.756 0.475  1.204
HYDROGROUPnot rated                 2.141 0.748  6.130
PercentImpervious_BaseLC_90m        4.545 1.705 12.116
Precip                              1.469 1.257  1.718
RoadGrade_mean_deg                  1.050 1.014  1.088
SurfaceUnpaved                      1.481 1.139  1.926
Culvert_Any                         1.666 1.449  1.915
Driveway_Count                      1.006 0.900  1.125
Mean_TopoConvergence_30m            1.294 1.197  1.400
vtrans_district1                    0.100 0.040  0.250
vtrans_district3                    0.228 0.118  0.440
vtrans_district4                    0.131 0.062  0.277
vtrans_district5                    0.254 0.077  0.831
vtrans_district6                    0.448 0.286  0.702
vtrans_district7                    0.186 0.085  0.406
vtrans_district8                    0.092 0.047  0.181
vtrans_district9                    0.029 0.015  0.056
road_miles_municipal                1.002 0.992  1.011
grand_list_equalized_muni_100M      0.967 0.917  1.019
=== STEP 3: Part 2, what did it cost? ===

--- Damaged segments ---
[1] 1457
--- Mean and median cost ---
  mean median 
 39960   7795 

--- The fitted model (coefficients are on the log-cost scale) ---

Call:
glm2(formula = cost_dedup ~ MeanSlope90m + StreamCrossing + StreamOrder + 
    PARENT + K_final + HYDROGROUP + PercentImpervious_BaseLC_90m + 
    Precip + Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m + 
    vtrans_district + road_miles_municipal + grand_list_equalized_muni_100M, 
    family = Gamma(link = "log"), data = damaged)

Coefficients:
                                Estimate Std. Error t value Pr(>|t|)    
(Intercept)                     5.707725   0.793585   7.192 1.03e-12 ***
MeanSlope90m                    0.072232   0.016182   4.464 8.69e-06 ***
StreamCrossing                  0.059876   0.258136   0.232 0.816606    
StreamOrder                     0.227261   0.104513   2.174 0.029834 *  
PARENTDT                        0.170580   0.414485   0.412 0.680733    
PARENTDT/GT                     0.994072   0.932999   1.065 0.286849    
PARENTGF                        0.003897   0.384789   0.010 0.991920    
PARENTGL                        0.132984   0.495213   0.269 0.788323    
PARENTGT                       -0.662830   0.373139  -1.776 0.075887 .  
PARENTGT/DT                     0.110840   0.854064   0.130 0.896759    
PARENTM                        -0.390565   2.293107  -0.170 0.864782    
K_final                         0.951998   0.611062   1.558 0.119470    
HYDROGROUPA/D                  -0.143556   0.837695  -0.171 0.863957    
HYDROGROUPB                    -0.004477   0.312535  -0.014 0.988573    
HYDROGROUPB/D                  -0.135864   0.430539  -0.316 0.752377    
HYDROGROUPC                    -0.271026   0.342537  -0.791 0.428942    
HYDROGROUPC/D                  -0.567432   0.408276  -1.390 0.164800    
HYDROGROUPD                    -0.673869   0.400711  -1.682 0.092849 .  
HYDROGROUPnot rated            -0.244306   0.793885  -0.308 0.758329    
PercentImpervious_BaseLC_90m    2.389935   1.136576   2.103 0.035663 *  
Precip                          0.227222   0.058973   3.853 0.000122 ***
SurfaceUnpaved                 -0.337897   0.175204  -1.929 0.053980 .  
Culvert_Any                    -0.133432   0.132689  -1.006 0.314780    
Driveway_Count                 -0.260546   0.085929  -3.032 0.002473 ** 
Mean_TopoConvergence_30m        0.166419   0.060433   2.754 0.005966 ** 
vtrans_district1                0.930248   0.433477   2.146 0.032041 *  
vtrans_district3                2.106675   0.252597   8.340  < 2e-16 ***
vtrans_district4                0.616431   0.230988   2.669 0.007702 ** 
vtrans_district5                1.624267   0.428341   3.792 0.000156 ***
vtrans_district6                1.384926   0.180058   7.692 2.70e-14 ***
vtrans_district7                1.427096   0.381584   3.740 0.000191 ***
vtrans_district8                1.415119   0.492936   2.871 0.004155 ** 
vtrans_district9                1.258631   0.622373   2.022 0.043331 *  
road_miles_municipal            0.007656   0.003223   2.375 0.017663 *  
grand_list_equalized_muni_100M -0.005518   0.016880  -0.327 0.743802    
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

(Dispersion parameter for Gamma family taken to be 5.037085)

    Null deviance: 4936.4  on 1456  degrees of freedom
Residual deviance: 3620.2  on 1422  degrees of freedom
AIC: 32012

Number of Fisher Scoring iterations: 19

--- Cost ratios with clustered 95% intervals ---
                               cost_ratio  lo_95    hi_95
(Intercept)                       301.185 78.197 1160.058
MeanSlope90m                        1.075  1.037    1.115
StreamCrossing                      1.062  0.424    2.660
StreamOrder                         1.255  0.816    1.931
PARENTDT                            1.186  0.590    2.385
PARENTDT/GT                         2.702  0.427   17.103
PARENTGF                            1.004  0.510    1.976
PARENTGL                            1.142  0.518    2.517
PARENTGT                            0.515  0.276    0.963
PARENTGT/DT                         1.117  0.408    3.063
PARENTM                             0.677  0.268    1.711
K_final                             2.591  0.682    9.843
HYDROGROUPA/D                       0.866  0.373    2.013
HYDROGROUPB                         0.996  0.530    1.870
HYDROGROUPB/D                       0.873  0.406    1.875
HYDROGROUPC                         0.763  0.431    1.349
HYDROGROUPC/D                       0.567  0.330    0.974
HYDROGROUPD                         0.510  0.269    0.966
HYDROGROUPnot rated                 0.783  0.381    1.609
PercentImpervious_BaseLC_90m       10.913  0.263  452.680
Precip                              1.255  1.088    1.448
SurfaceUnpaved                      0.713  0.422    1.207
Culvert_Any                         0.875  0.674    1.136
Driveway_Count                      0.771  0.655    0.907
Mean_TopoConvergence_30m            1.181  1.051    1.327
vtrans_district1                    2.535  1.262    5.091
vtrans_district3                    8.221  4.850   13.935
vtrans_district4                    1.852  1.141    3.007
vtrans_district5                    5.075  2.561   10.056
vtrans_district6                    3.995  2.481    6.433
vtrans_district7                    4.167  2.419    7.177
vtrans_district8                    4.117  0.903   18.774
vtrans_district9                    3.521  2.100    5.903
road_miles_municipal                1.008  0.999    1.016
grand_list_equalized_muni_100M      0.994  0.962    1.028

--- Deviance R-squared ---
[1] 0.2666
--- Mean predicted / mean observed cost (1 = the model reproduces the total) ---
[1] 1.1006
=== STEP 4: the avoided cost ===

--- Compliant segments recoded ---
[1] 50175
--- Mean probability as rated, compliant ---
[1] 0.0216
--- Mean probability if Does Not Meet ---
[1] 0.0266
--- Statewide avoided cost, dollars ---
[1] 11492123

--- The ten sample towns (NA = not in this storm's exposed population) ---
 Brattleboro      Corinth  Wallingford     Hardwick West Windsor      Wolcott 
       40441        59487       114854        69499          929       175214 
  Washington     Richmond     Stamford   Starksboro 
      196540        14877         8352        42818 

--- Town-clustered 95% interval for the statewide total ---
  clustered SE : 5428829 
  95% interval : 851619 22132628 
=== STEP 5: both-source towns ===

--- Towns with both record sources ---
[1] 24
--- Segments ---
[1] 8641

--- Run A: pooled de-duplicated cost ---
damaged segments: 667 
compliance odds ratio (clustered 95%): 0.889 0.713 1.11 
AUC: 0.7491   deviance R2: 0.3079 
avoided cost, these towns: 2258820 

--- Run B: VTrans records only ---
damaged segments: 261 
compliance odds ratio (clustered 95%): 1.117 0.822 1.519 
AUC: 0.7925   deviance R2: 0.3893 
avoided cost, these towns: -790326 

--- Side by side, same towns and segments ---
                    run damaged odds_ratio lo_95 hi_95   AUC dev_R2 avoided
 A pooled de-duplicated     667      0.889 0.713 1.110 0.749  0.308 2258820
          B VTrans only     261      1.117 0.822 1.519 0.793  0.389 -790326
```

## Output 2024

```
--- Number of segments ---
[1] 22571
--- Number damaged ---
[1] 1018
--- Number of towns ---
[1] 94
--- Compliance status ---

Does Not Meet     Compliant 
         4448         18123 
=== STEP 1: Part 1, did the segment flood? ===

--- The fitted model (coefficients are log-odds) ---

Call:
glm(formula = damaged ~ compliant + MeanSlope90m + StreamOrder + 
    PARENT + HYDROGROUP + PercentImpervious_BaseLC_90m + Precip + 
    RoadGrade_mean_deg + Surface + Culvert_Any + Driveway_Count + 
    Mean_TopoConvergence_30m + vtrans_district + road_miles_municipal + 
    grand_list_equalized_muni_100M, family = binomial, data = roads)

Coefficients:
                                Estimate Std. Error z value Pr(>|z|)    
(Intercept)                    -8.670501   0.512459 -16.919  < 2e-16 ***
compliantCompliant             -0.247991   0.077253  -3.210  0.00133 ** 
MeanSlope90m                    0.080641   0.009854   8.183 2.76e-16 ***
StreamOrder                     0.234168   0.032226   7.266 3.69e-13 ***
PARENTDT                        0.422059   0.228742   1.845  0.06502 .  
PARENTDT/GT                     0.631193   0.395277   1.597  0.11030    
PARENTGF                       -0.189904   0.237092  -0.801  0.42315    
PARENTGF/GL                     0.555528   0.802132   0.693  0.48858    
PARENTGL                       -0.240439   0.218113  -1.102  0.27031    
PARENTGT                       -0.126837   0.203175  -0.624  0.53245    
PARENTGT/DT                    -1.445285   0.367454  -3.933 8.38e-05 ***
PARENTGT/GT_R                  -0.474274   0.757607  -0.626  0.53130    
PARENTM                         1.187802   0.650971   1.825  0.06805 .  
HYDROGROUPA/D                  -0.152670   0.284666  -0.536  0.59174    
HYDROGROUPB                    -0.034265   0.232244  -0.148  0.88271    
HYDROGROUPB/D                  -0.241252   0.261172  -0.924  0.35563    
HYDROGROUPC                    -0.549553   0.250132  -2.197  0.02802 *  
HYDROGROUPC/D                  -0.623265   0.267331  -2.331  0.01973 *  
HYDROGROUPD                    -0.745123   0.260271  -2.863  0.00420 ** 
HYDROGROUPnot rated            -0.738621   0.516483  -1.430  0.15269    
PercentImpervious_BaseLC_90m    1.243215   0.546534   2.275  0.02292 *  
Precip                          0.443346   0.064256   6.900 5.21e-12 ***
RoadGrade_mean_deg              0.080014   0.018948   4.223 2.41e-05 ***
SurfaceUnpaved                  0.584786   0.095998   6.092 1.12e-09 ***
Culvert_Any                     0.619183   0.073124   8.468  < 2e-16 ***
Driveway_Count                  0.087219   0.042311   2.061  0.03927 *  
Mean_TopoConvergence_30m        0.380955   0.035011  10.881  < 2e-16 ***
vtrans_district5               -0.482150   0.100780  -4.784 1.72e-06 ***
vtrans_district7               -1.030241   0.096111 -10.719  < 2e-16 ***
vtrans_district9               -1.564506   0.175244  -8.928  < 2e-16 ***
road_miles_municipal           -0.004963   0.001775  -2.797  0.00516 ** 
grand_list_equalized_muni_100M  0.006021   0.003663   1.644  0.10020    
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 8298.6  on 22570  degrees of freedom
Residual deviance: 7297.6  on 22539  degrees of freedom
AIC: 7361.6

Number of Fisher Scoring iterations: 7

--- Odds ratios: exp(coefficient) ---
                   (Intercept)             compliantCompliant 
                         0.000                          0.780 
                  MeanSlope90m                    StreamOrder 
                         1.084                          1.264 
                      PARENTDT                    PARENTDT/GT 
                         1.525                          1.880 
                      PARENTGF                    PARENTGF/GL 
                         0.827                          1.743 
                      PARENTGL                       PARENTGT 
                         0.786                          0.881 
                   PARENTGT/DT                  PARENTGT/GT_R 
                         0.236                          0.622 
                       PARENTM                  HYDROGROUPA/D 
                         3.280                          0.858 
                   HYDROGROUPB                  HYDROGROUPB/D 
                         0.966                          0.786 
                   HYDROGROUPC                  HYDROGROUPC/D 
                         0.577                          0.536 
                   HYDROGROUPD            HYDROGROUPnot rated 
                         0.475                          0.478 
  PercentImpervious_BaseLC_90m                         Precip 
                         3.467                          1.558 
            RoadGrade_mean_deg                 SurfaceUnpaved 
                         1.083                          1.795 
                   Culvert_Any                 Driveway_Count 
                         1.857                          1.091 
      Mean_TopoConvergence_30m               vtrans_district5 
                         1.464                          0.617 
              vtrans_district7               vtrans_district9 
                         0.357                          0.209 
          road_miles_municipal grand_list_equalized_muni_100M 
                         0.995                          1.006 

--- Does compliance belong in the model? Drop in deviance ---
Analysis of Deviance Table

Model 1: damaged ~ MeanSlope90m + StreamOrder + PARENT + HYDROGROUP + 
    PercentImpervious_BaseLC_90m + Precip + RoadGrade_mean_deg + 
    Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m + 
    vtrans_district + road_miles_municipal + grand_list_equalized_muni_100M
Model 2: damaged ~ compliant + MeanSlope90m + StreamOrder + PARENT + HYDROGROUP + 
    PercentImpervious_BaseLC_90m + Precip + RoadGrade_mean_deg + 
    Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m + 
    vtrans_district + road_miles_municipal + grand_list_equalized_muni_100M
  Resid. Df Resid. Dev Df Deviance Pr(>Chi)   
1     22540     7307.6                        
2     22539     7297.6  1   10.068 0.001509 **
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

--- AUC ---
[1] 0.7782
=== STEP 2: town-clustered standard errors ===

--- Compliance coefficient: regular vs clustered ---
  coefficient          : -0.248 
  regular SE           : 0.0773 
  clustered SE         : 0.1179 
  odds ratio           : 0.78 
  95% interval, regular: 0.671 0.908 
  95% interval, clustered: 0.619 0.983 

--- Every coefficient as an odds ratio with its clustered 95% interval ---
                               odds_ratio lo_95  hi_95
(Intercept)                         0.000 0.000  0.002
compliantCompliant                  0.780 0.619  0.983
MeanSlope90m                        1.084 1.050  1.119
StreamOrder                         1.264 1.180  1.354
PARENTDT                            1.525 0.896  2.596
PARENTDT/GT                         1.880 0.939  3.761
PARENTGF                            0.827 0.353  1.936
PARENTGF/GL                         1.743 0.378  8.028
PARENTGL                            0.786 0.490  1.262
PARENTGT                            0.881 0.540  1.437
PARENTGT/DT                         0.236 0.093  0.597
PARENTGT/GT_R                       0.622 0.341  1.136
PARENTM                             3.280 0.617 17.422
HYDROGROUPA/D                       0.858 0.521  1.414
HYDROGROUPB                         0.966 0.456  2.049
HYDROGROUPB/D                       0.786 0.317  1.946
HYDROGROUPC                         0.577 0.244  1.368
HYDROGROUPC/D                       0.536 0.254  1.130
HYDROGROUPD                         0.475 0.208  1.083
HYDROGROUPnot rated                 0.478 0.109  2.098
PercentImpervious_BaseLC_90m        3.467 1.086 11.066
Precip                              1.558 1.064  2.282
RoadGrade_mean_deg                  1.083 1.033  1.136
SurfaceUnpaved                      1.795 1.287  2.503
Culvert_Any                         1.857 1.510  2.284
Driveway_Count                      1.091 0.962  1.238
Mean_TopoConvergence_30m            1.464 1.336  1.604
vtrans_district5                    0.617 0.372  1.025
vtrans_district7                    0.357 0.200  0.637
vtrans_district9                    0.209 0.088  0.496
road_miles_municipal                0.995 0.983  1.007
grand_list_equalized_muni_100M      1.006 0.995  1.018
=== STEP 3: Part 2, what did it cost? ===

--- Damaged segments ---
[1] 1018
--- Mean and median cost ---
  mean median 
 64100  11922 

--- The fitted model (coefficients are on the log-cost scale) ---

Call:
glm2(formula = cost_dedup ~ MeanSlope90m + StreamCrossing + StreamOrder + 
    PARENT + K_final + HYDROGROUP + PercentImpervious_BaseLC_90m + 
    Precip + Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m + 
    vtrans_district + road_miles_municipal + grand_list_equalized_muni_100M, 
    family = Gamma(link = "log"), data = damaged)

Coefficients:
                                Estimate Std. Error t value Pr(>|t|)    
(Intercept)                    10.722225   1.017810  10.535  < 2e-16 ***
MeanSlope90m                    0.075822   0.021686   3.496 0.000493 ***
StreamCrossing                 -0.708385   0.297651  -2.380 0.017506 *  
StreamOrder                     0.643796   0.124201   5.184 2.64e-07 ***
PARENTDT                       -0.067974   0.492960  -0.138 0.890356    
PARENTDT/GT                    -0.763631   0.838437  -0.911 0.362635    
PARENTGF                       -0.038010   0.425319  -0.089 0.928808    
PARENTGF/GL                    -1.347058   1.656224  -0.813 0.416225    
PARENTGL                       -0.599336   0.460664  -1.301 0.193553    
PARENTGT                       -0.226473   0.434726  -0.521 0.602515    
PARENTGT/DT                     0.742920   0.776483   0.957 0.338915    
PARENTGT/GT_R                  -1.528109   1.631363  -0.937 0.349139    
PARENTM                        -0.571195   1.248553  -0.457 0.647423    
K_final                         0.850870   0.898843   0.947 0.344060    
HYDROGROUPA/D                  -0.175704   0.603361  -0.291 0.770953    
HYDROGROUPB                     0.299445   0.440709   0.679 0.497005    
HYDROGROUPB/D                   0.176840   0.510081   0.347 0.728899    
HYDROGROUPC                    -0.521416   0.500838  -1.041 0.298090    
HYDROGROUPC/D                  -0.032566   0.553832  -0.059 0.953122    
HYDROGROUPD                    -0.275236   0.557586  -0.494 0.621684    
HYDROGROUPnot rated             0.698255   1.044927   0.668 0.504141    
PercentImpervious_BaseLC_90m    1.471794   1.233753   1.193 0.233180    
Precip                         -0.137135   0.143743  -0.954 0.340303    
SurfaceUnpaved                 -0.372482   0.200623  -1.857 0.063662 .  
Culvert_Any                     0.250207   0.157692   1.587 0.112904    
Driveway_Count                 -0.080170   0.095928  -0.836 0.403511    
Mean_TopoConvergence_30m        0.051087   0.071697   0.713 0.476295    
vtrans_district5               -0.822385   0.215583  -3.815 0.000145 ***
vtrans_district7                1.113055   0.224199   4.965 8.11e-07 ***
vtrans_district9               -0.892969   0.381751  -2.339 0.019527 *  
road_miles_municipal           -0.014900   0.004244  -3.511 0.000467 ***
grand_list_equalized_muni_100M  0.014098   0.008339   1.691 0.091219 .  
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

(Dispersion parameter for Gamma family taken to be 4.776752)

    Null deviance: 3157.1  on 1017  degrees of freedom
Residual deviance: 2316.7  on  986  degrees of freedom
AIC: 23519

Number of Fisher Scoring iterations: 16

--- Cost ratios with clustered 95% intervals ---
                               cost_ratio    lo_95      hi_95
(Intercept)                     45352.717 8407.748 244639.704
MeanSlope90m                        1.079    1.034      1.125
StreamCrossing                      0.492    0.265      0.916
StreamOrder                         1.904    1.426      2.542
PARENTDT                            0.934    0.417      2.095
PARENTDT/GT                         0.466    0.136      1.593
PARENTGF                            0.963    0.364      2.549
PARENTGF/GL                         0.260    0.035      1.910
PARENTGL                            0.549    0.319      0.947
PARENTGT                            0.797    0.393      1.619
PARENTGT/DT                         2.102    0.620      7.122
PARENTGT/GT_R                       0.217    0.064      0.739
PARENTM                             0.565    0.112      2.858
K_final                             2.342    0.456     12.026
HYDROGROUPA/D                       0.839    0.325      2.163
HYDROGROUPB                         1.349    0.504      3.611
HYDROGROUPB/D                       1.193    0.444      3.211
HYDROGROUPC                         0.594    0.192      1.836
HYDROGROUPC/D                       0.968    0.343      2.732
HYDROGROUPD                         0.759    0.283      2.040
HYDROGROUPnot rated                 2.010    0.221     18.285
PercentImpervious_BaseLC_90m        4.357    0.387     49.009
Precip                              0.872    0.649      1.171
SurfaceUnpaved                      0.689    0.518      0.917
Culvert_Any                         1.284    0.940      1.755
Driveway_Count                      0.923    0.781      1.090
Mean_TopoConvergence_30m            1.052    0.955      1.159
vtrans_district5                    0.439    0.265      0.728
vtrans_district7                    3.044    1.912      4.845
vtrans_district9                    0.409    0.295      0.569
road_miles_municipal                0.985    0.974      0.997
grand_list_equalized_muni_100M      1.014    1.000      1.029

--- Deviance R-squared ---
[1] 0.2662
--- Mean predicted / mean observed cost (1 = the model reproduces the total) ---
[1] 1.1088
=== STEP 4: the avoided cost ===

--- Compliant segments recoded ---
[1] 18123
--- Mean probability as rated, compliant ---
[1] 0.0398
--- Mean probability if Does Not Meet ---
[1] 0.0498
--- Statewide avoided cost, dollars ---
[1] 12714930

--- The ten sample towns (NA = not in this storm's exposed population) ---
      <NA>       <NA>       <NA>   Hardwick       <NA>    Wolcott       <NA> 
        NA         NA         NA     222483         NA     180992         NA 
  Richmond       <NA> Starksboro 
    206667         NA     238569 

--- Town-clustered 95% interval for the statewide total ---
  clustered SE : 7423237 
  95% interval : -1834614 27264475 
=== STEP 5: both-source towns ===

--- Towns with both record sources ---
[1] 26
--- Segments ---
[1] 9473

--- Run A: pooled de-duplicated cost ---
damaged segments: 747 
compliance odds ratio (clustered 95%): 0.705 0.552 0.901 
AUC: 0.7438   deviance R2: 0.1938 
avoided cost, these towns: 10490094 

--- Run B: VTrans records only ---
damaged segments: 262 
compliance odds ratio (clustered 95%): 0.618 0.474 0.807 
AUC: 0.7988   deviance R2: 0.3235 
avoided cost, these towns: 5926812 

--- Side by side, same towns and segments ---
                    run damaged odds_ratio lo_95 hi_95   AUC dev_R2  avoided
 A pooled de-duplicated     747      0.705 0.552 0.901 0.744  0.194 10490094
          B VTrans only     262      0.618 0.474 0.807 0.799  0.324  5926812
```
