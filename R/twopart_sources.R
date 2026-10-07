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
SHOW  <- c("compliantCompliant", "parent_A", "parent_DT", "parent_GL", "parent_GT_DT", "hydro_B")   # the rows to print

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
cat("\nPart 1: odds ratios with town-clustered 95% intervals (compliance and the indicators)\n")
print(round(exp(cbind(odds_ratio = coef(p1), lo_95 = coef(p1) - 1.96 * se, hi_95 = coef(p1) + 1.96 * se))[SHOW, ], 3))
n1 <- sum(d$damaged); n0 <- sum(d$damaged == 0)
cat("AUC:", round((sum(rank(fitted(p1))[d$damaged == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0), 4), "\n")

# Part 2 -- what did it cost?
p2 <- glm2(as.formula(paste("cost ~", PART2)), family = Gamma(link = "log"), data = k)
se2 <- sqrt(diag(vcovCL(p2, cluster = k$Town, type = "HC1")))
cat("\nPart 2: cost ratio for the till indicator, clustered 95%\n")
print(round(exp(c(cost_ratio = coef(p2)[["parent_till"]], lo_95 = coef(p2)[["parent_till"]] - 1.96 * se2[["parent_till"]], hi_95 = coef(p2)[["parent_till"]] + 1.96 * se2[["parent_till"]])), 3))
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
cat("\nPart 1: odds ratios with clustered 95% intervals\n")
print(round(exp(cbind(odds_ratio = coef(p1), lo_95 = coef(p1) - 1.96 * se, hi_95 = coef(p1) + 1.96 * se))[c("compliantCompliant", ind), ], 3))
n1 <- sum(d$damaged); n0 <- sum(d$damaged == 0)
cat("AUC:", round((sum(rank(fitted(p1))[d$damaged == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0), 4), "\n")
p2 <- glm2(as.formula(paste("cost ~", PART2)), family = Gamma(link = "log"), data = k)
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
cat("\nPart 1: odds ratios with clustered 95% intervals\n")
print(round(exp(cbind(odds_ratio = coef(p1), lo_95 = coef(p1) - 1.96 * se, hi_95 = coef(p1) + 1.96 * se))[c("compliantCompliant", ind), ], 3))
n1 <- sum(d$damaged); n0 <- sum(d$damaged == 0)
cat("AUC:", round((sum(rank(fitted(p1))[d$damaged == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0), 4), "\n")
p2 <- glm2(as.formula(paste("cost ~", PART2)), family = Gamma(link = "log"), data = k)
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
