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
