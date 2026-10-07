# Two-Part Model with Indicators, by Damage Record Source

Parent material and hydrologic group enter as the 0/1 indicators chosen from the crosstab and ANOVA; VTrans district stays as a factor. The model is fitted statewide, then on the towns that hold both FEMA and VTrans records with damage counted from the pooled de-duplicated records (run A) and from the VTrans records alone (run B), and the three runs are put side by side. Every coefficient is shown as B with its clustered SE and as Exp(B) with its interval. Script run once with STORM = 2023 and once with STORM = 2024.

## Script: twopart_sources.R

```r
# =============================================================================
# The two-part model with indicator variables, and what changes with the damage record source
#
# Data:  ../mrgp-roads-core/twopart_simple_2026-09-17/data/analysis_<STORM>_2026-09-27.csv
#        one row per road segment; compliance-graded segments with at least 3.51 in of rain,
#        Averill excluded; damage and cost from cost_dedup (FEMA and VTrans, de-duplicated)
#        or from vt_cost (VTrans records alone)
#
# Parent material and hydrologic group enter as 0/1 indicators chosen from the crosstab and ANOVA:
#   occurrence model: parent_A, parent_DT, parent_GL, parent_GT_DT, hydro_B
#   cost model:       parent_till (dense till or ablation till)
# Every other category is the baseline. VTrans district stays as a factor: it marks where the storm
# went, and without it the compliance effect absorbs the storm's path.
#
# STEP 1 fits the model statewide. STEPs 2 and 3 fit it on the towns that hold both a FEMA and a
# VTrans record, once counting damage from the pooled records and once from VTrans alone -- same
# segments, same terms, so the difference is the record source. STEP 4 puts the runs side by side.
#
# Run one line at a time with Ctrl+Enter. Set STORM to 2023 or 2024. Nothing is written to disk.
# =============================================================================


# -----------------------------------------------------------------------------
# STEP 0.  Read the data
# -----------------------------------------------------------------------------
STORM <- 2023

library(glm2)        # a steadier fitting routine for the Gamma model
library(sandwich)    # town-clustered standard errors

if (!dir.exists("data") && dir.exists("../data")) setwd("..")
roads <- read.csv(paste0("../mrgp-roads-core/twopart_simple_2026-09-17/data/analysis_", STORM, "_2026-09-27.csv"))
roads <- roads[roads$Precip >= 3.51 & roads$Town != "Averill", ]
roads$compliant <- relevel(factor(roads$compliant), ref = "Does Not Meet")

# the model, written once
PART1 <- "compliant + parent_A + parent_DT + parent_GL + parent_GT_DT + hydro_B + MeanSlope90m + StreamOrder + PercentImpervious_BaseLC_90m + Precip + RoadGrade_mean_deg + Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m + vtrans_district + road_miles_municipal + grand_list_equalized_muni_100M"
PART2 <- "parent_till + MeanSlope90m + StreamCrossing + StreamOrder + K_final + PercentImpervious_BaseLC_90m + Precip + Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m + vtrans_district + road_miles_municipal + grand_list_equalized_muni_100M"

cat("--- Segments ---\n");            print(nrow(roads))
cat("--- Towns ---\n");               print(length(unique(roads$Town)))
cat("--- Compliance status ---\n");   print(table(roads$compliant))


# -----------------------------------------------------------------------------
# STEP 1.  Statewide, pooled de-duplicated records
# -----------------------------------------------------------------------------
cat("=== STEP 1: statewide ===\n\n")
d <- roads
d$damaged <- as.integer(d$cost_dedup > 0)
d$cost    <- d$cost_dedup
n <- tapply(d$damaged, d$vtrans_district, sum)                    # a district with no damage joins the busiest
d$vtrans_district[d$vtrans_district %in% names(n)[n == 0]] <- names(which.max(n))
d$vtrans_district <- relevel(factor(d$vtrans_district), ref = names(which.max(n)))
k <- d[d$damaged == 1, ]
cat("damaged segments:", sum(d$damaged), "  total cost:", format(sum(k$cost), big.mark = ","), "\n")

# Part 1 -- did the segment flood?
p1 <- glm(as.formula(paste("damaged ~", PART1)), family = binomial, data = d)
se <- sqrt(diag(vcovCL(p1, cluster = d$Town, type = "HC1")))
# Every term: B is the coefficient on the log-odds scale (SPSS's B column), Exp(B) is the odds
# ratio (SPSS's Exp(B)); the interval is built on B with the clustered SE and then exponentiated.
cat("\nPart 1: every term -- B, clustered SE, Exp(B) = odds ratio, clustered 95% interval\n")
print(round(cbind(B = coef(p1), SE_clustered = se, Exp_B = exp(coef(p1)),
                  lo_95 = exp(coef(p1) - 1.96 * se), hi_95 = exp(coef(p1) + 1.96 * se)), 3))
n1 <- sum(d$damaged); n0 <- sum(d$damaged == 0)
cat("AUC:", round((sum(rank(fitted(p1))[d$damaged == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0), 4), "\n")

# Part 2 -- what did it cost?
p2 <- glm2(as.formula(paste("cost ~", PART2)), family = Gamma(link = "log"), data = k)
se2 <- sqrt(diag(vcovCL(p2, cluster = k$Town, type = "HC1")))
# Every term: B is the coefficient on the log-cost scale, Exp(B) is the cost ratio (1.20 = 20% higher cost)
cat("\nPart 2: every term -- B, clustered SE, Exp(B) = cost ratio, clustered 95% interval\n")
print(round(cbind(B = coef(p2), SE_clustered = se2, Exp_B = exp(coef(p2)),
                  lo_95 = exp(coef(p2) - 1.96 * se2), hi_95 = exp(coef(p2) + 1.96 * se2)), 3))
cat("deviance R2:", round(1 - p2$deviance / p2$null.deviance, 4), "\n")

# Avoided cost -- recode every compliant segment to Does Not Meet
d0 <- d; d0$compliant <- factor("Does Not Meet", levels = levels(d$compliant))
p_rated <- predict(p1, newdata = d,  type = "response")
p_dnm   <- predict(p1, newdata = d0, type = "response")
cost_hat <- predict(p2, newdata = d, type = "response")
avoided <- (p_dnm - p_rated) * cost_hat
cc <- d$compliant == "Compliant"
# its town-clustered standard error: delta method with the joint sandwich of both parts
k1 <- names(coef(p1)); k2 <- names(coef(p2))
X1 <- model.matrix(delete.response(terms(p1)), d, xlev = p1$xlevels)[, k1]; X0 <- model.matrix(delete.response(terms(p1)), d0, xlev = p1$xlevels)[, k1]
Z  <- model.matrix(delete.response(terms(p2)), d, xlev = p2$xlevels)[, k2]
h  <- c(colSums(((cost_hat * p_dnm * (1 - p_dnm)) * X0 - (cost_hat * p_rated * (1 - p_rated)) * X1)[cc, ]), colSums(((p_dnm - p_rated) * cost_hat * Z)[cc, ]))
B1 <- summary(p1)$cov.unscaled; w <- weights(p2, "working"); rr <- residuals(p2, "working") * w; B2 <- summary(p2)$cov.unscaled * sum(rr^2) / sum(w)
tw <- sort(unique(d$Town))
U1 <- matrix(0, length(tw), length(k1), dimnames = list(tw, k1)); s1 <- rowsum(estfun(p1), d$Town); U1[rownames(s1), ] <- s1
U2 <- matrix(0, length(tw), length(k2), dimnames = list(tw, k2)); s2 <- rowsum(estfun(p2), k$Town); U2[rownames(s2), ] <- s2
G1 <- length(tw); G2 <- length(unique(k$Town)); a1 <- G1 / (G1 - 1) * (nrow(d) - 1) / (nrow(d) - length(k1)); a2 <- G2 / (G2 - 1) * (nrow(k) - 1) / (nrow(k) - length(k2))
V  <- crossprod(cbind(sqrt(a1) * U1 %*% B1, sqrt(a2) * U2 %*% B2)); se_av <- sqrt(drop(t(h) %*% V %*% h))
cat("\nAvoided cost statewide:", format(round(sum(avoided[cc])), big.mark = ","),
    " clustered 95%:", format(round(sum(avoided[cc]) - 1.96 * se_av), big.mark = ","), "to", format(round(sum(avoided[cc]) + 1.96 * se_av), big.mark = ","), "\n")
cat("as a share of the observed damage cost:", round(100 * sum(avoided[cc]) / sum(k$cost), 1), "%\n")
cat("ten sample towns (NA = not in this storm's exposed population):\n")
print(round(tapply(avoided[cc], d$Town[cc], sum)[c("Brattleboro", "Corinth", "Wallingford", "Hardwick", "West Windsor", "Wolcott", "Washington", "Richmond", "Stamford", "Starksboro")]))
state <- c(segments = nrow(d), damaged = sum(d$damaged), odds_ratio = exp(coef(p1)[["compliantCompliant"]]),
           lo = exp(coef(p1)[["compliantCompliant"]] - 1.96 * se[["compliantCompliant"]]), hi = exp(coef(p1)[["compliantCompliant"]] + 1.96 * se[["compliantCompliant"]]),
           avoided = sum(avoided[cc]), avoided_lo = sum(avoided[cc]) - 1.96 * se_av, avoided_hi = sum(avoided[cc]) + 1.96 * se_av, pct_of_cost = 100 * sum(avoided[cc]) / sum(k$cost))


# -----------------------------------------------------------------------------
# STEP 2.  Both-source towns, run A: damage from the pooled de-duplicated records
# -----------------------------------------------------------------------------
# Only the towns with at least one FEMA-recorded and one VTrans-recorded damage, so that run B
# below uses exactly the same segments.
cat("\n=== STEP 2: both-source towns, pooled records ===\n\n")
both <- roads[roads$Town %in% intersect(unique(roads$Town[roads$fema_cost > 0]), unique(roads$Town[roads$vt_cost > 0])), ]
cat("towns:", length(unique(both$Town)), "  segments:", nrow(both), "\n")

d <- both
d$damaged <- as.integer(d$cost_dedup > 0)
d$cost    <- d$cost_dedup
n <- tapply(d$damaged, d$vtrans_district, sum)
d$vtrans_district[d$vtrans_district %in% names(n)[n == 0]] <- names(which.max(n))
d$vtrans_district <- relevel(factor(d$vtrans_district), ref = names(which.max(n)))
k <- d[d$damaged == 1, ]
cat("damaged segments:", sum(d$damaged), "  total cost:", format(sum(k$cost), big.mark = ","), "\n")
# an indicator with no damaged segment in these towns cannot be estimated and is left out of this run
ind <- c("parent_A", "parent_DT", "parent_GL", "parent_GT_DT", "hydro_B")
ind <- ind[colSums(k[, ind]) > 0]
cat("indicators estimable here:", paste(ind, collapse = ", "), "\n")
f1 <- paste("damaged ~ compliant +", paste(ind, collapse = " + "), "+ MeanSlope90m + StreamOrder + PercentImpervious_BaseLC_90m + Precip + RoadGrade_mean_deg + Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m + vtrans_district + road_miles_municipal + grand_list_equalized_muni_100M")

p1 <- glm(as.formula(f1), family = binomial, data = d)
se <- sqrt(diag(vcovCL(p1, cluster = d$Town, type = "HC1")))
cat("\nPart 1: every term -- B, clustered SE, Exp(B) = odds ratio, clustered 95% interval\n")
print(round(cbind(B = coef(p1), SE_clustered = se, Exp_B = exp(coef(p1)),
                  lo_95 = exp(coef(p1) - 1.96 * se), hi_95 = exp(coef(p1) + 1.96 * se)), 3))
n1 <- sum(d$damaged); n0 <- sum(d$damaged == 0)
cat("AUC:", round((sum(rank(fitted(p1))[d$damaged == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0), 4), "\n")
p2 <- glm2(as.formula(paste("cost ~", PART2)), family = Gamma(link = "log"), data = k)
se2 <- sqrt(diag(vcovCL(p2, cluster = k$Town, type = "HC1")))
cat("\nPart 2: every term -- B, clustered SE, Exp(B) = cost ratio, clustered 95% interval\n")
print(round(cbind(B = coef(p2), SE_clustered = se2, Exp_B = exp(coef(p2)),
                  lo_95 = exp(coef(p2) - 1.96 * se2), hi_95 = exp(coef(p2) + 1.96 * se2)), 3))
cat("Part 2 deviance R2:", round(1 - p2$deviance / p2$null.deviance, 4), "\n")

d0 <- d; d0$compliant <- factor("Does Not Meet", levels = levels(d$compliant))
p_rated <- predict(p1, newdata = d, type = "response"); p_dnm <- predict(p1, newdata = d0, type = "response"); cost_hat <- predict(p2, newdata = d, type = "response")
avoided <- (p_dnm - p_rated) * cost_hat; cc <- d$compliant == "Compliant"
k1 <- names(coef(p1)); k2 <- names(coef(p2))
X1 <- model.matrix(delete.response(terms(p1)), d, xlev = p1$xlevels)[, k1]; X0 <- model.matrix(delete.response(terms(p1)), d0, xlev = p1$xlevels)[, k1]; Z <- model.matrix(delete.response(terms(p2)), d, xlev = p2$xlevels)[, k2]
h  <- c(colSums(((cost_hat * p_dnm * (1 - p_dnm)) * X0 - (cost_hat * p_rated * (1 - p_rated)) * X1)[cc, ]), colSums(((p_dnm - p_rated) * cost_hat * Z)[cc, ]))
B1 <- summary(p1)$cov.unscaled; w <- weights(p2, "working"); rr <- residuals(p2, "working") * w; B2 <- summary(p2)$cov.unscaled * sum(rr^2) / sum(w)
tw <- sort(unique(d$Town)); U1 <- matrix(0, length(tw), length(k1), dimnames = list(tw, k1)); s1 <- rowsum(estfun(p1), d$Town); U1[rownames(s1), ] <- s1
U2 <- matrix(0, length(tw), length(k2), dimnames = list(tw, k2)); s2 <- rowsum(estfun(p2), k$Town); U2[rownames(s2), ] <- s2
G1 <- length(tw); G2 <- length(unique(k$Town)); a1 <- G1 / (G1 - 1) * (nrow(d) - 1) / (nrow(d) - length(k1)); a2 <- G2 / (G2 - 1) * (nrow(k) - 1) / (nrow(k) - length(k2))
V  <- crossprod(cbind(sqrt(a1) * U1 %*% B1, sqrt(a2) * U2 %*% B2)); se_av <- sqrt(drop(t(h) %*% V %*% h))
cat("\nAvoided cost, these towns:", format(round(sum(avoided[cc])), big.mark = ","), " clustered 95%:", format(round(sum(avoided[cc]) - 1.96 * se_av), big.mark = ","), "to", format(round(sum(avoided[cc]) + 1.96 * se_av), big.mark = ","), "\n")
cat("as a share of the observed damage cost:", round(100 * sum(avoided[cc]) / sum(k$cost), 1), "%\n")
runA <- c(segments = nrow(d), damaged = sum(d$damaged), odds_ratio = exp(coef(p1)[["compliantCompliant"]]),
          lo = exp(coef(p1)[["compliantCompliant"]] - 1.96 * se[["compliantCompliant"]]), hi = exp(coef(p1)[["compliantCompliant"]] + 1.96 * se[["compliantCompliant"]]),
          avoided = sum(avoided[cc]), avoided_lo = sum(avoided[cc]) - 1.96 * se_av, avoided_hi = sum(avoided[cc]) + 1.96 * se_av, pct_of_cost = 100 * sum(avoided[cc]) / sum(k$cost))


# -----------------------------------------------------------------------------
# STEP 3.  Both-source towns, run B: damage from the VTrans records alone
# -----------------------------------------------------------------------------
# The same lines as STEP 2 with vt_cost in place of cost_dedup. A segment that only FEMA recorded
# counts as undamaged here.
cat("\n=== STEP 3: both-source towns, VTrans records only ===\n\n")
d <- both
d$damaged <- as.integer(d$vt_cost > 0)
d$cost    <- d$vt_cost
n <- tapply(d$damaged, d$vtrans_district, sum)
d$vtrans_district[d$vtrans_district %in% names(n)[n == 0]] <- names(which.max(n))
d$vtrans_district <- relevel(factor(d$vtrans_district), ref = names(which.max(n)))
k <- d[d$damaged == 1, ]
cat("damaged segments:", sum(d$damaged), "  total cost:", format(sum(k$cost), big.mark = ","), "\n")
ind <- c("parent_A", "parent_DT", "parent_GL", "parent_GT_DT", "hydro_B")
ind <- ind[colSums(k[, ind]) > 0]
cat("indicators estimable here:", paste(ind, collapse = ", "), "\n")
f1 <- paste("damaged ~ compliant +", paste(ind, collapse = " + "), "+ MeanSlope90m + StreamOrder + PercentImpervious_BaseLC_90m + Precip + RoadGrade_mean_deg + Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m + vtrans_district + road_miles_municipal + grand_list_equalized_muni_100M")

p1 <- glm(as.formula(f1), family = binomial, data = d)
se <- sqrt(diag(vcovCL(p1, cluster = d$Town, type = "HC1")))
cat("\nPart 1: every term -- B, clustered SE, Exp(B) = odds ratio, clustered 95% interval\n")
print(round(cbind(B = coef(p1), SE_clustered = se, Exp_B = exp(coef(p1)),
                  lo_95 = exp(coef(p1) - 1.96 * se), hi_95 = exp(coef(p1) + 1.96 * se)), 3))
n1 <- sum(d$damaged); n0 <- sum(d$damaged == 0)
cat("AUC:", round((sum(rank(fitted(p1))[d$damaged == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0), 4), "\n")
p2 <- glm2(as.formula(paste("cost ~", PART2)), family = Gamma(link = "log"), data = k)
se2 <- sqrt(diag(vcovCL(p2, cluster = k$Town, type = "HC1")))
cat("\nPart 2: every term -- B, clustered SE, Exp(B) = cost ratio, clustered 95% interval\n")
print(round(cbind(B = coef(p2), SE_clustered = se2, Exp_B = exp(coef(p2)),
                  lo_95 = exp(coef(p2) - 1.96 * se2), hi_95 = exp(coef(p2) + 1.96 * se2)), 3))
cat("Part 2 deviance R2:", round(1 - p2$deviance / p2$null.deviance, 4), "\n")

d0 <- d; d0$compliant <- factor("Does Not Meet", levels = levels(d$compliant))
p_rated <- predict(p1, newdata = d, type = "response"); p_dnm <- predict(p1, newdata = d0, type = "response"); cost_hat <- predict(p2, newdata = d, type = "response")
avoided <- (p_dnm - p_rated) * cost_hat; cc <- d$compliant == "Compliant"
k1 <- names(coef(p1)); k2 <- names(coef(p2))
X1 <- model.matrix(delete.response(terms(p1)), d, xlev = p1$xlevels)[, k1]; X0 <- model.matrix(delete.response(terms(p1)), d0, xlev = p1$xlevels)[, k1]; Z <- model.matrix(delete.response(terms(p2)), d, xlev = p2$xlevels)[, k2]
h  <- c(colSums(((cost_hat * p_dnm * (1 - p_dnm)) * X0 - (cost_hat * p_rated * (1 - p_rated)) * X1)[cc, ]), colSums(((p_dnm - p_rated) * cost_hat * Z)[cc, ]))
B1 <- summary(p1)$cov.unscaled; w <- weights(p2, "working"); rr <- residuals(p2, "working") * w; B2 <- summary(p2)$cov.unscaled * sum(rr^2) / sum(w)
tw <- sort(unique(d$Town)); U1 <- matrix(0, length(tw), length(k1), dimnames = list(tw, k1)); s1 <- rowsum(estfun(p1), d$Town); U1[rownames(s1), ] <- s1
U2 <- matrix(0, length(tw), length(k2), dimnames = list(tw, k2)); s2 <- rowsum(estfun(p2), k$Town); U2[rownames(s2), ] <- s2
G1 <- length(tw); G2 <- length(unique(k$Town)); a1 <- G1 / (G1 - 1) * (nrow(d) - 1) / (nrow(d) - length(k1)); a2 <- G2 / (G2 - 1) * (nrow(k) - 1) / (nrow(k) - length(k2))
V  <- crossprod(cbind(sqrt(a1) * U1 %*% B1, sqrt(a2) * U2 %*% B2)); se_av <- sqrt(drop(t(h) %*% V %*% h))
cat("\nAvoided cost, these towns:", format(round(sum(avoided[cc])), big.mark = ","), " clustered 95%:", format(round(sum(avoided[cc]) - 1.96 * se_av), big.mark = ","), "to", format(round(sum(avoided[cc]) + 1.96 * se_av), big.mark = ","), "\n")
cat("as a share of the observed damage cost:", round(100 * sum(avoided[cc]) / sum(k$cost), 1), "%\n")
runB <- c(segments = nrow(d), damaged = sum(d$damaged), odds_ratio = exp(coef(p1)[["compliantCompliant"]]),
          lo = exp(coef(p1)[["compliantCompliant"]] - 1.96 * se[["compliantCompliant"]]), hi = exp(coef(p1)[["compliantCompliant"]] + 1.96 * se[["compliantCompliant"]]),
          avoided = sum(avoided[cc]), avoided_lo = sum(avoided[cc]) - 1.96 * se_av, avoided_hi = sum(avoided[cc]) + 1.96 * se_av, pct_of_cost = 100 * sum(avoided[cc]) / sum(k$cost))


# -----------------------------------------------------------------------------
# STEP 4.  Side by side
# -----------------------------------------------------------------------------
# Runs A and B share towns, segments and terms; only the record source differs. Compare the
# compliance odds ratio (scale-free) and the avoided cost as a share of each run's own damage
# cost; the raw dollars are not comparable because VTrans dollars are a subset of the pooled ones.
cat("\n=== STEP 4: side by side,", STORM, "storm ===\n\n")
side <- rbind(state, runA, runB)
rownames(side) <- c("statewide, pooled", "both-source towns, pooled (A)", "both-source towns, VTrans only (B)")
print(data.frame(segments   = side[, "segments"], damaged = side[, "damaged"],
                 odds_ratio = round(side[, "odds_ratio"], 3),
                 clustered_95 = paste0(round(side[, "lo"], 3), " to ", round(side[, "hi"], 3)),
                 avoided    = format(round(side[, "avoided"]), big.mark = ","),
                 avoided_95 = paste0(format(round(side[, "avoided_lo"]), big.mark = ","), " to ", format(round(side[, "avoided_hi"]), big.mark = ",")),
                 pct_of_cost = round(side[, "pct_of_cost"], 1)))


# =============================================================================
# END
#
# To report: the compliance odds ratio and the avoided cost (with its clustered interval and as a
# share of observed cost) for the statewide run, and for runs A and B on the both-source towns.
# Where A and B disagree, the difference is what the two record systems captured, not the roads.
# =============================================================================
```

## Output 2023

```
--- Segments ---
[1] 59695
--- Towns ---
[1] 204
--- Compliance status ---

Does Not Meet     Compliant 
         9520         50175 
=== STEP 1: statewide ===

damaged segments: 1457   total cost: 58,222,056 

Part 1: every term -- B, clustered SE, Exp(B) = odds ratio, clustered 95% interval
                                    B SE_clustered Exp_B lo_95  hi_95
(Intercept)                    -8.159        0.602 0.000 0.000  0.001
compliantCompliant             -0.231        0.086 0.794 0.670  0.940
parent_A                        0.462        0.111 1.587 1.275  1.974
parent_DT                       0.345        0.102 1.412 1.156  1.724
parent_GL                       0.293        0.225 1.340 0.862  2.083
parent_GT_DT                   -0.642        0.517 0.526 0.191  1.449
hydro_B                         0.132        0.100 1.141 0.938  1.388
MeanSlope90m                    0.064        0.012 1.066 1.041  1.091
StreamOrder                     0.289        0.032 1.335 1.253  1.423
PercentImpervious_BaseLC_90m    1.597        0.524 4.939 1.770 13.781
Precip                          0.380        0.079 1.462 1.252  1.709
RoadGrade_mean_deg              0.047        0.018 1.048 1.012  1.086
SurfaceUnpaved                  0.377        0.133 1.458 1.124  1.892
Culvert_Any                     0.502        0.072 1.652 1.435  1.902
Driveway_Count                  0.005        0.058 1.005 0.897  1.127
Mean_TopoConvergence_30m        0.268        0.041 1.308 1.207  1.417
vtrans_district1               -2.323        0.468 0.098 0.039  0.245
vtrans_district3               -1.519        0.329 0.219 0.115  0.417
vtrans_district4               -2.037        0.385 0.130 0.061  0.277
vtrans_district5               -1.371        0.603 0.254 0.078  0.828
vtrans_district6               -0.820        0.231 0.440 0.280  0.692
vtrans_district7               -1.689        0.394 0.185 0.085  0.400
vtrans_district8               -2.380        0.351 0.093 0.047  0.184
vtrans_district9               -3.570        0.338 0.028 0.015  0.055
road_miles_municipal            0.002        0.005 1.002 0.992  1.011
grand_list_equalized_muni_100M -0.033        0.027 0.967 0.917  1.021
AUC: 0.8223 

Part 2: every term -- B, clustered SE, Exp(B) = cost ratio, clustered 95% interval
                                    B SE_clustered   Exp_B   lo_95    hi_95
(Intercept)                     6.050        0.695 424.307 108.711 1656.099
parent_till                    -0.453        0.183   0.636   0.444    0.910
MeanSlope90m                    0.064        0.016   1.066   1.033    1.100
StreamCrossing                  0.081        0.457   1.085   0.443    2.658
StreamOrder                     0.216        0.210   1.242   0.823    1.873
K_final                         0.779        0.625   2.179   0.640    7.418
PercentImpervious_BaseLC_90m    1.938        1.856   6.947   0.183  263.945
Precip                          0.239        0.072   1.269   1.103    1.461
SurfaceUnpaved                 -0.433        0.301   0.649   0.360    1.169
Culvert_Any                    -0.120        0.139   0.887   0.675    1.164
Driveway_Count                 -0.268        0.081   0.765   0.653    0.896
Mean_TopoConvergence_30m        0.164        0.065   1.178   1.037    1.338
vtrans_district1                0.751        0.365   2.119   1.035    4.337
vtrans_district3                1.913        0.274   6.776   3.957   11.604
vtrans_district4                0.487        0.251   1.628   0.995    2.663
vtrans_district5                1.450        0.331   4.262   2.228    8.154
vtrans_district6                1.174        0.219   3.234   2.106    4.966
vtrans_district7                1.253        0.253   3.503   2.135    5.746
vtrans_district8                1.202        0.635   3.327   0.958   11.552
vtrans_district9                1.007        0.268   2.737   1.620    4.624
road_miles_municipal            0.007        0.004   1.007   0.998    1.015
grand_list_equalized_muni_100M -0.004        0.018   0.996   0.962    1.031
deviance R2: 0.2515 

Avoided cost statewide: 11,368,894  clustered 95%: 899,243 to 21,838,545 
as a share of the observed damage cost: 19.5 %
ten sample towns (NA = not in this storm's exposed population):
 Brattleboro      Corinth  Wallingford     Hardwick West Windsor      Wolcott 
       42251        59902       113580        71885         1009       159527 
  Washington     Richmond     Stamford   Starksboro 
      173891        14667         8320        37902 

=== STEP 2: both-source towns, pooled records ===

towns: 24   segments: 8641 
damaged segments: 667   total cost: 29,577,443 
indicators estimable here: parent_A, parent_DT, parent_GL, hydro_B 

Part 1: every term -- B, clustered SE, Exp(B) = odds ratio, clustered 95% interval
                                    B SE_clustered Exp_B lo_95  hi_95
(Intercept)                    -6.430        0.663 0.002 0.000  0.006
compliantCompliant             -0.133        0.115 0.876 0.699  1.097
parent_A                        0.464        0.161 1.590 1.159  2.181
parent_DT                       0.478        0.127 1.613 1.258  2.069
parent_GL                       0.096        0.314 1.101 0.594  2.038
hydro_B                         0.002        0.144 1.002 0.756  1.329
MeanSlope90m                    0.047        0.016 1.048 1.015  1.081
StreamOrder                     0.310        0.049 1.364 1.238  1.502
PercentImpervious_BaseLC_90m    0.636        0.982 1.888 0.275 12.943
Precip                          0.025        0.072 1.025 0.891  1.180
RoadGrade_mean_deg              0.080        0.023 1.083 1.036  1.132
SurfaceUnpaved                  0.334        0.267 1.396 0.828  2.356
Culvert_Any                     0.559        0.090 1.749 1.468  2.085
Driveway_Count                  0.029        0.087 1.030 0.868  1.222
Mean_TopoConvergence_30m        0.285        0.045 1.330 1.218  1.452
vtrans_district1               -0.351        0.291 0.704 0.398  1.244
vtrans_district2                0.777        0.174 2.174 1.547  3.055
vtrans_district3               -0.766        0.195 0.465 0.317  0.681
vtrans_district4               -0.806        0.241 0.447 0.279  0.717
vtrans_district5                0.083        0.540 1.087 0.377  3.134
vtrans_district7               -0.493        0.145 0.611 0.460  0.812
vtrans_district8               -1.402        0.297 0.246 0.137  0.441
road_miles_municipal            0.002        0.005 1.002 0.992  1.012
grand_list_equalized_muni_100M -0.081        0.057 0.922 0.824  1.032
AUC: 0.7436 

Part 2: every term -- B, clustered SE, Exp(B) = cost ratio, clustered 95% interval
                                    B SE_clustered     Exp_B    lo_95     hi_95
(Intercept)                     9.414        0.737 12255.535 2891.854 51938.354
parent_till                    -0.468        0.364     0.626    0.307     1.278
MeanSlope90m                    0.062        0.016     1.063    1.030     1.098
StreamCrossing                 -0.459        0.443     0.632    0.265     1.508
StreamOrder                     0.557        0.190     1.745    1.203     2.530
K_final                         0.741        0.844     2.098    0.401    10.969
PercentImpervious_BaseLC_90m   -2.152        1.179     0.116    0.012     1.171
Precip                          0.035        0.096     1.036    0.858     1.250
SurfaceUnpaved                 -0.575        0.325     0.563    0.298     1.064
Culvert_Any                    -0.104        0.191     0.901    0.619     1.311
Driveway_Count                 -0.178        0.104     0.837    0.683     1.026
Mean_TopoConvergence_30m        0.125        0.069     1.133    0.990     1.298
vtrans_district1               -0.475        0.323     0.622    0.330     1.170
vtrans_district2               -0.685        0.157     0.504    0.370     0.686
vtrans_district3               -0.568        0.225     0.566    0.365     0.880
vtrans_district4               -1.639        0.186     0.194    0.135     0.279
vtrans_district5                0.488        0.601     1.629    0.502     5.292
vtrans_district7               -0.272        0.126     0.762    0.596     0.975
vtrans_district8                1.157        0.521     3.181    1.145     8.839
road_miles_municipal            0.005        0.006     1.005    0.993     1.016
grand_list_equalized_muni_100M -0.112        0.069     0.894    0.781     1.024
Part 2 deviance R2: 0.2682 

Avoided cost, these towns: 2,669,712  clustered 95%: -2,160,509 to 7,499,932 
as a share of the observed damage cost: 9 %

=== STEP 3: both-source towns, VTrans records only ===

damaged segments: 261   total cost: 9,223,297 
indicators estimable here: parent_A, parent_DT, parent_GL, hydro_B 

Part 1: every term -- B, clustered SE, Exp(B) = odds ratio, clustered 95% interval
                                    B SE_clustered Exp_B lo_95  hi_95
(Intercept)                    -7.413        0.936 0.001 0.000  0.004
compliantCompliant              0.083        0.164 1.087 0.788  1.499
parent_A                       -0.638        0.376 0.529 0.253  1.105
parent_DT                       0.349        0.188 1.417 0.981  2.048
parent_GL                       0.520        0.449 1.683 0.698  4.060
hydro_B                        -0.282        0.252 0.754 0.460  1.237
MeanSlope90m                    0.086        0.029 1.090 1.029  1.154
StreamOrder                     0.298        0.111 1.347 1.084  1.673
PercentImpervious_BaseLC_90m   -0.268        1.889 0.765 0.019 30.982
Precip                          0.004        0.093 1.004 0.836  1.206
RoadGrade_mean_deg              0.082        0.047 1.086 0.990  1.190
SurfaceUnpaved                 -0.305        0.323 0.737 0.391  1.389
Culvert_Any                     0.340        0.159 1.405 1.029  1.918
Driveway_Count                  0.105        0.084 1.110 0.941  1.310
Mean_TopoConvergence_30m        0.342        0.066 1.408 1.238  1.603
vtrans_district1                0.243        0.365 1.276 0.624  2.609
vtrans_district2                1.203        0.297 3.329 1.859  5.960
vtrans_district3               -0.332        0.273 0.718 0.420  1.226
vtrans_district4               -0.810        0.258 0.445 0.268  0.738
vtrans_district5                0.336        0.617 1.400 0.418  4.689
vtrans_district7               -2.059        0.226 0.128 0.082  0.199
vtrans_district8               -0.887        0.324 0.412 0.218  0.777
road_miles_municipal            0.001        0.006 1.001 0.989  1.013
grand_list_equalized_muni_100M -0.146        0.066 0.865 0.759  0.985
AUC: 0.7807 

Part 2: every term -- B, clustered SE, Exp(B) = cost ratio, clustered 95% interval
                                    B SE_clustered     Exp_B    lo_95
(Intercept)                    10.737        1.390 46038.820 3018.333
parent_till                    -0.532        0.315     0.587    0.317
MeanSlope90m                    0.002        0.028     1.002    0.948
StreamCrossing                 -0.381        0.388     0.683    0.319
StreamOrder                     0.156        0.175     1.168    0.829
K_final                         2.166        0.951     8.727    1.354
PercentImpervious_BaseLC_90m   -1.615        1.765     0.199    0.006
Precip                         -0.048        0.124     0.953    0.747
SurfaceUnpaved                 -0.291        0.227     0.748    0.479
Culvert_Any                     0.024        0.269     1.024    0.604
Driveway_Count                 -0.581        0.213     0.559    0.368
Mean_TopoConvergence_30m        0.066        0.119     1.068    0.846
vtrans_district1               -0.950        0.388     0.387    0.181
vtrans_district2               -1.946        0.298     0.143    0.080
vtrans_district3               -0.043        0.266     0.958    0.569
vtrans_district4               -1.297        0.292     0.273    0.154
vtrans_district5               -0.111        0.782     0.895    0.193
vtrans_district7                0.030        0.486     1.031    0.398
vtrans_district8               -0.239        0.444     0.788    0.330
road_miles_municipal            0.005        0.007     1.005    0.992
grand_list_equalized_muni_100M -0.150        0.075     0.861    0.742
                                    hi_95
(Intercept)                    702233.019
parent_till                         1.089
MeanSlope90m                        1.058
StreamCrossing                      1.462
StreamOrder                         1.648
K_final                            56.236
PercentImpervious_BaseLC_90m        6.319
Precip                              1.216
SurfaceUnpaved                      1.167
Culvert_Any                         1.737
Driveway_Count                      0.849
Mean_TopoConvergence_30m            1.347
vtrans_district1                    0.828
vtrans_district2                    0.256
vtrans_district3                    1.613
vtrans_district4                    0.485
vtrans_district5                    4.145
vtrans_district7                    2.671
vtrans_district8                    1.880
road_miles_municipal                1.018
grand_list_equalized_muni_100M      0.998
Part 2 deviance R2: 0.3078 

Avoided cost, these towns: -522,629  clustered 95%: -2,456,174 to 1,410,916 
as a share of the observed damage cost: -5.7 %

=== STEP 4: side by side, 2023 storm ===

                                   segments damaged odds_ratio   clustered_95
statewide, pooled                     59695    1457      0.794   0.67 to 0.94
both-source towns, pooled (A)          8641     667      0.876 0.699 to 1.097
both-source towns, VTrans only (B)     8641     261      1.087 0.788 to 1.499
                                      avoided               avoided_95
statewide, pooled                  11,368,894    899,243 to 21,838,545
both-source towns, pooled (A)       2,669,712 -2,160,509 to  7,499,932
both-source towns, VTrans only (B)   -522,629 -2,456,174 to  1,410,916
                                   pct_of_cost
statewide, pooled                         19.5
both-source towns, pooled (A)              9.0
both-source towns, VTrans only (B)        -5.7
```

## Output 2024

```
--- Segments ---
[1] 22571
--- Towns ---
[1] 94
--- Compliance status ---

Does Not Meet     Compliant 
         4448         18123 
=== STEP 1: statewide ===

damaged segments: 1018   total cost: 65,254,042 

Part 1: every term -- B, clustered SE, Exp(B) = odds ratio, clustered 95% interval
                                    B SE_clustered Exp_B lo_95  hi_95
(Intercept)                    -9.086        1.208 0.000 0.000  0.001
compliantCompliant             -0.263        0.118 0.769 0.610  0.970
parent_A                        0.147        0.205 1.158 0.775  1.730
parent_DT                       0.225        0.123 1.252 0.984  1.593
parent_GL                      -0.376        0.175 0.686 0.487  0.968
parent_GT_DT                   -1.389        0.455 0.249 0.102  0.609
hydro_B                         0.345        0.121 1.412 1.114  1.790
MeanSlope90m                    0.083        0.016 1.087 1.054  1.121
StreamOrder                     0.236        0.036 1.266 1.179  1.359
PercentImpervious_BaseLC_90m    1.525        0.580 4.593 1.475 14.309
Precip                          0.421        0.200 1.524 1.030  2.255
RoadGrade_mean_deg              0.077        0.025 1.080 1.029  1.133
SurfaceUnpaved                  0.562        0.166 1.755 1.266  2.432
Culvert_Any                     0.599        0.105 1.821 1.483  2.236
Driveway_Count                  0.084        0.066 1.088 0.957  1.237
Mean_TopoConvergence_30m        0.388        0.046 1.475 1.348  1.614
vtrans_district5               -0.437        0.265 0.646 0.384  1.085
vtrans_district7               -1.026        0.301 0.358 0.199  0.647
vtrans_district9               -1.554        0.439 0.211 0.089  0.500
road_miles_municipal           -0.006        0.006 0.994 0.982  1.007
grand_list_equalized_muni_100M  0.009        0.006 1.009 0.996  1.021
AUC: 0.7747 

Part 2: every term -- B, clustered SE, Exp(B) = cost ratio, clustered 95% interval
                                    B SE_clustered     Exp_B    lo_95
(Intercept)                    10.871        0.912 52635.351 8804.891
parent_till                     0.062        0.205     1.064    0.712
MeanSlope90m                    0.090        0.021     1.094    1.051
StreamCrossing                 -0.651        0.345     0.521    0.265
StreamOrder                     0.659        0.152     1.933    1.435
K_final                        -0.423        0.496     0.655    0.248
PercentImpervious_BaseLC_90m    1.716        1.327     5.564    0.413
Precip                         -0.183        0.162     0.833    0.606
SurfaceUnpaved                 -0.310        0.151     0.734    0.546
Culvert_Any                     0.224        0.158     1.251    0.918
Driveway_Count                 -0.116        0.087     0.890    0.751
Mean_TopoConvergence_30m        0.052        0.060     1.053    0.936
vtrans_district5               -0.693        0.235     0.500    0.316
vtrans_district7                1.237        0.265     3.445    2.048
vtrans_district9               -0.722        0.185     0.486    0.338
road_miles_municipal           -0.014        0.006     0.986    0.975
grand_list_equalized_muni_100M  0.014        0.007     1.015    1.001
                                    hi_95
(Intercept)                    314652.404
parent_till                         1.590
MeanSlope90m                        1.139
StreamCrossing                      1.026
StreamOrder                         2.605
K_final                             1.733
PercentImpervious_BaseLC_90m       74.953
Precip                              1.144
SurfaceUnpaved                      0.986
Culvert_Any                         1.706
Driveway_Count                      1.056
Mean_TopoConvergence_30m            1.185
vtrans_district5                    0.793
vtrans_district7                    5.795
vtrans_district9                    0.698
road_miles_municipal                0.998
grand_list_equalized_muni_100M      1.028
deviance R2: 0.2382 

Avoided cost statewide: 13,603,182  clustered 95%: -1,240,559 to 28,446,922 
as a share of the observed damage cost: 20.8 %
ten sample towns (NA = not in this storm's exposed population):
      <NA>       <NA>       <NA>   Hardwick       <NA>    Wolcott       <NA> 
        NA         NA         NA     236877         NA     160096         NA 
  Richmond       <NA> Starksboro 
    222257         NA     170113 

=== STEP 2: both-source towns, pooled records ===

towns: 26   segments: 9473 
damaged segments: 747   total cost: 42,092,428 
indicators estimable here: parent_A, parent_DT, parent_GL, parent_GT_DT, hydro_B 

Part 1: every term -- B, clustered SE, Exp(B) = odds ratio, clustered 95% interval
                                    B SE_clustered Exp_B lo_95  hi_95
(Intercept)                    -6.111        1.185 0.002 0.000  0.023
compliantCompliant             -0.355        0.125 0.701 0.549  0.895
parent_A                       -0.044        0.270 0.957 0.564  1.623
parent_DT                       0.048        0.143 1.049 0.793  1.389
parent_GL                      -0.103        0.207 0.902 0.601  1.353
parent_GT_DT                   -1.938        0.292 0.144 0.081  0.255
hydro_B                         0.187        0.138 1.205 0.919  1.580
MeanSlope90m                    0.072        0.019 1.075 1.035  1.116
StreamOrder                     0.285        0.036 1.330 1.240  1.427
PercentImpervious_BaseLC_90m    1.704        0.797 5.495 1.153 26.189
Precip                         -0.065        0.236 0.937 0.591  1.487
RoadGrade_mean_deg              0.067        0.030 1.070 1.008  1.135
SurfaceUnpaved                  0.701        0.205 2.016 1.349  3.014
Culvert_Any                     0.595        0.122 1.813 1.427  2.305
Driveway_Count                  0.093        0.069 1.097 0.958  1.256
Mean_TopoConvergence_30m        0.370        0.054 1.447 1.301  1.609
vtrans_district5               -0.023        0.253 0.977 0.595  1.604
vtrans_district7                0.393        0.390 1.482 0.690  3.185
road_miles_municipal           -0.015        0.008 0.985 0.969  1.000
grand_list_equalized_muni_100M  0.012        0.009 1.012 0.995  1.030
AUC: 0.7411 

Part 2: every term -- B, clustered SE, Exp(B) = cost ratio, clustered 95% interval
                                    B SE_clustered     Exp_B    lo_95
(Intercept)                    11.055        1.135 63248.542 6832.032
parent_till                     0.079        0.247     1.083    0.668
MeanSlope90m                    0.064        0.029     1.066    1.008
StreamCrossing                 -0.115        0.425     0.891    0.388
StreamOrder                     0.495        0.158     1.640    1.203
K_final                        -0.454        0.670     0.635    0.171
PercentImpervious_BaseLC_90m    0.042        1.262     1.043    0.088
Precip                          0.012        0.206     1.012    0.675
SurfaceUnpaved                 -0.363        0.158     0.696    0.510
Culvert_Any                     0.089        0.200     1.093    0.738
Driveway_Count                 -0.054        0.074     0.947    0.819
Mean_TopoConvergence_30m        0.019        0.069     1.020    0.890
vtrans_district5               -0.626        0.205     0.535    0.358
vtrans_district7                1.339        0.373     3.816    1.838
road_miles_municipal           -0.022        0.009     0.978    0.961
grand_list_equalized_muni_100M  0.021        0.010     1.021    1.001
                                    hi_95
(Intercept)                    585532.726
parent_till                         1.756
MeanSlope90m                        1.127
StreamCrossing                      2.048
StreamOrder                         2.236
K_final                             2.359
PercentImpervious_BaseLC_90m       12.383
Precip                              1.517
SurfaceUnpaved                      0.949
Culvert_Any                         1.619
Driveway_Count                      1.096
Mean_TopoConvergence_30m            1.168
vtrans_district5                    0.800
vtrans_district7                    7.923
road_miles_municipal                0.996
grand_list_equalized_muni_100M      1.041
Part 2 deviance R2: 0.1589 

Avoided cost, these towns: 10,716,687  clustered 95%: 1,272,943 to 20,160,431 
as a share of the observed damage cost: 25.5 %

=== STEP 3: both-source towns, VTrans records only ===

damaged segments: 262   total cost: 15,319,997 
indicators estimable here: parent_A, parent_DT, parent_GL, parent_GT_DT, hydro_B 

Part 1: every term -- B, clustered SE, Exp(B) = odds ratio, clustered 95% interval
                                    B SE_clustered Exp_B lo_95 hi_95
(Intercept)                    -5.416        1.364 0.004 0.000 0.064
compliantCompliant             -0.485        0.142 0.616 0.466 0.813
parent_A                        0.429        0.394 1.536 0.709 3.328
parent_DT                       0.163        0.167 1.177 0.849 1.633
parent_GL                      -0.036        0.287 0.965 0.550 1.694
parent_GT_DT                   -1.466        0.430 0.231 0.099 0.536
hydro_B                        -0.095        0.238 0.909 0.570 1.449
MeanSlope90m                    0.093        0.022 1.098 1.052 1.146
StreamOrder                     0.390        0.052 1.477 1.334 1.636
PercentImpervious_BaseLC_90m    0.117        1.036 1.124 0.148 8.564
Precip                         -0.414        0.243 0.661 0.410 1.065
RoadGrade_mean_deg             -0.034        0.035 0.966 0.902 1.035
SurfaceUnpaved                  0.066        0.235 1.068 0.674 1.692
Culvert_Any                     0.694        0.191 2.002 1.375 2.913
Driveway_Count                  0.098        0.071 1.103 0.960 1.267
Mean_TopoConvergence_30m        0.409        0.065 1.505 1.325 1.711
vtrans_district5                0.734        0.346 2.084 1.058 4.103
vtrans_district7                0.218        0.372 1.244 0.600 2.577
road_miles_municipal           -0.020        0.008 0.980 0.966 0.995
grand_list_equalized_muni_100M  0.024        0.009 1.024 1.007 1.042
AUC: 0.7944 

Part 2: every term -- B, clustered SE, Exp(B) = cost ratio, clustered 95% interval
                                    B SE_clustered    Exp_B   lo_95     hi_95
(Intercept)                     8.277        1.283 3931.780 318.121 48594.401
parent_till                     0.053        0.363    1.055   0.518     2.147
MeanSlope90m                    0.110        0.025    1.117   1.063     1.172
StreamCrossing                 -0.271        0.309    0.762   0.416     1.397
StreamOrder                     0.235        0.146    1.265   0.951     1.684
K_final                         0.117        0.866    1.124   0.206     6.133
PercentImpervious_BaseLC_90m    5.614        2.613  274.372   1.638 45962.834
Precip                          0.109        0.248    1.115   0.685     1.814
SurfaceUnpaved                  0.304        0.209    1.355   0.899     2.042
Culvert_Any                     0.188        0.208    1.206   0.802     1.815
Driveway_Count                 -0.199        0.126    0.820   0.641     1.049
Mean_TopoConvergence_30m        0.082        0.077    1.085   0.933     1.262
vtrans_district5               -0.804        0.271    0.448   0.263     0.761
vtrans_district7                1.675        0.700    5.338   1.355    21.033
road_miles_municipal           -0.009        0.009    0.991   0.973     1.009
grand_list_equalized_muni_100M  0.000        0.011    1.000   0.980     1.022
Part 2 deviance R2: 0.2922 

Avoided cost, these towns: 5,941,244  clustered 95%: 1,363,655 to 10,518,832 
as a share of the observed damage cost: 38.8 %

=== STEP 4: side by side, 2024 storm ===

                                   segments damaged odds_ratio   clustered_95
statewide, pooled                     22571    1018      0.769   0.61 to 0.97
both-source towns, pooled (A)          9473     747      0.701 0.549 to 0.895
both-source towns, VTrans only (B)     9473     262      0.616 0.466 to 0.813
                                      avoided               avoided_95
statewide, pooled                  13,603,182 -1,240,559 to 28,446,922
both-source towns, pooled (A)      10,716,687  1,272,943 to 20,160,431
both-source towns, VTrans only (B)  5,941,244  1,363,655 to 10,518,832
                                   pct_of_cost
statewide, pooled                         20.8
both-source towns, pooled (A)             25.5
both-source towns, VTrans only (B)        38.8
```
