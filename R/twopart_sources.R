# Data:  ../mrgp-roads-core/twopart_simple_2026-09-17/data/analysis_<STORM>_2026-09-27.csv
#        one row per road segment; compliance-graded segments with at least 3.51 in of rain,
#        Averill excluded; damage and cost from cost_dedup (FEMA and VTrans, de-duplicated)
#        or from vt_cost (VTrans records alone)
#        indicators from the crosstab and anova: parent_A, parent_DT, parent_GL, parent_GT_DT,
#        hydro_B (occurrence side); parent_till (cost side)


# step 0  read the data
STORM <- 2023

library(glm2)        # glm package failed, glm2 did not
library(sandwich)    # town-clustered standard errors

if (!dir.exists("data") && dir.exists("../data")) setwd("..")
roads <- read.csv(paste0("../mrgp-roads-core/twopart_simple_2026-09-17/data/analysis_", STORM, "_2026-09-27.csv"))
roads <- roads[roads$Precip >= 3.51 & roads$Town != "Averill", ]
roads$compliant <- relevel(factor(roads$compliant), ref = "Does Not Meet")

# the model, written once; district stays as a factor
PART1 <- "compliant + parent_A + parent_DT + parent_GL + parent_GT_DT + hydro_B + MeanSlope90m + StreamOrder + PercentImpervious_BaseLC_90m + Precip + RoadGrade_mean_deg + Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m + vtrans_district + road_miles_municipal + grand_list_equalized_muni_100M"
PART2 <- "parent_till + MeanSlope90m + StreamCrossing + StreamOrder + K_final + PercentImpervious_BaseLC_90m + Precip + Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m + vtrans_district + road_miles_municipal + grand_list_equalized_muni_100M"

cat("segments\n");           print(nrow(roads))
cat("towns\n");              print(length(unique(roads$Town)))
cat("compliance status\n");  print(table(roads$compliant))


# step 1  statewide, pooled de-duplicated records
cat("step 1 statewide\n\n")
d <- roads
d$damaged <- as.integer(d$cost_dedup > 0)
d$cost    <- d$cost_dedup
n <- tapply(d$damaged, d$vtrans_district, sum)                    # a district with no damage joins the busiest
d$vtrans_district[d$vtrans_district %in% names(n)[n == 0]] <- names(which.max(n))
d$vtrans_district <- relevel(factor(d$vtrans_district), ref = names(which.max(n)))
k <- d[d$damaged == 1, ]
cat("damaged segments:", sum(d$damaged), "  total cost:", format(sum(k$cost), big.mark = ","), "\n")

# part 1  did the segment flood
p1 <- glm(as.formula(paste("damaged ~", PART1)), family = binomial, data = d)
se <- sqrt(diag(vcovCL(p1, cluster = d$Town, type = "HC1")))
# B is the log-odds coefficient, Exp(B) the odds ratio; the interval is built on B then exponentiated
cat("\npart 1, every term: B, clustered SE, Exp(B) = odds ratio, clustered 95% interval\n")
print(round(cbind(B = coef(p1), SE_clustered = se, Exp_B = exp(coef(p1)),
                  lo_95 = exp(coef(p1) - 1.96 * se), hi_95 = exp(coef(p1) + 1.96 * se)), 3))
n1 <- sum(d$damaged); n0 <- sum(d$damaged == 0)
cat("AUC:", round((sum(rank(fitted(p1))[d$damaged == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0), 4), "\n")

# part 2  what did it cost
p2 <- glm2(as.formula(paste("cost ~", PART2)), family = Gamma(link = "log"), data = k)
se2 <- sqrt(diag(vcovCL(p2, cluster = k$Town, type = "HC1")))
# B is the log-cost coefficient, Exp(B) the cost ratio (1.20 = 20% higher cost)
cat("\npart 2, every term: B, clustered SE, Exp(B) = cost ratio, clustered 95% interval\n")
print(round(cbind(B = coef(p2), SE_clustered = se2, Exp_B = exp(coef(p2)),
                  lo_95 = exp(coef(p2) - 1.96 * se2), hi_95 = exp(coef(p2) + 1.96 * se2)), 3))
cat("deviance R2:", round(1 - p2$deviance / p2$null.deviance, 4), "\n")

# avoided cost  recode every compliant segment to does not meet
d0 <- d; d0$compliant <- factor("Does Not Meet", levels = levels(d$compliant))
p_rated <- predict(p1, newdata = d,  type = "response")
p_dnm   <- predict(p1, newdata = d0, type = "response")
cost_hat <- predict(p2, newdata = d, type = "response")
avoided <- (p_dnm - p_rated) * cost_hat
cc <- d$compliant == "Compliant"
# town-clustered standard error: delta method with the joint sandwich of both parts
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
cat("\navoided cost statewide:", format(round(sum(avoided[cc])), big.mark = ","),
    " clustered 95%:", format(round(sum(avoided[cc]) - 1.96 * se_av), big.mark = ","), "to", format(round(sum(avoided[cc]) + 1.96 * se_av), big.mark = ","), "\n")
cat("as a share of the observed damage cost:", round(100 * sum(avoided[cc]) / sum(k$cost), 1), "%\n")
cat("ten sample towns (NA = not in this storm's exposed population):\n")
print(round(tapply(avoided[cc], d$Town[cc], sum)[c("Brattleboro", "Corinth", "Wallingford", "Hardwick", "West Windsor", "Wolcott", "Washington", "Richmond", "Stamford", "Starksboro")]))
state <- c(segments = nrow(d), damaged = sum(d$damaged), odds_ratio = exp(coef(p1)[["compliantCompliant"]]),
           lo = exp(coef(p1)[["compliantCompliant"]] - 1.96 * se[["compliantCompliant"]]), hi = exp(coef(p1)[["compliantCompliant"]] + 1.96 * se[["compliantCompliant"]]),
           avoided = sum(avoided[cc]), avoided_lo = sum(avoided[cc]) - 1.96 * se_av, avoided_hi = sum(avoided[cc]) + 1.96 * se_av, pct_of_cost = 100 * sum(avoided[cc]) / sum(k$cost))


# step 2  both-source towns, run A: damage from the pooled de-duplicated records
# only the towns with at least one FEMA-recorded and one VTrans-recorded damage, so run B uses the same segments
cat("\nstep 2 both-source towns, pooled records\n\n")
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
cat("\npart 1, every term: B, clustered SE, Exp(B) = odds ratio, clustered 95% interval\n")
print(round(cbind(B = coef(p1), SE_clustered = se, Exp_B = exp(coef(p1)),
                  lo_95 = exp(coef(p1) - 1.96 * se), hi_95 = exp(coef(p1) + 1.96 * se)), 3))
n1 <- sum(d$damaged); n0 <- sum(d$damaged == 0)
cat("AUC:", round((sum(rank(fitted(p1))[d$damaged == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0), 4), "\n")
p2 <- glm2(as.formula(paste("cost ~", PART2)), family = Gamma(link = "log"), data = k)
se2 <- sqrt(diag(vcovCL(p2, cluster = k$Town, type = "HC1")))
cat("\npart 2, every term: B, clustered SE, Exp(B) = cost ratio, clustered 95% interval\n")
print(round(cbind(B = coef(p2), SE_clustered = se2, Exp_B = exp(coef(p2)),
                  lo_95 = exp(coef(p2) - 1.96 * se2), hi_95 = exp(coef(p2) + 1.96 * se2)), 3))
cat("deviance R2:", round(1 - p2$deviance / p2$null.deviance, 4), "\n")

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
cat("\navoided cost, these towns:", format(round(sum(avoided[cc])), big.mark = ","), " clustered 95%:", format(round(sum(avoided[cc]) - 1.96 * se_av), big.mark = ","), "to", format(round(sum(avoided[cc]) + 1.96 * se_av), big.mark = ","), "\n")
cat("as a share of the observed damage cost:", round(100 * sum(avoided[cc]) / sum(k$cost), 1), "%\n")
runA <- c(segments = nrow(d), damaged = sum(d$damaged), odds_ratio = exp(coef(p1)[["compliantCompliant"]]),
          lo = exp(coef(p1)[["compliantCompliant"]] - 1.96 * se[["compliantCompliant"]]), hi = exp(coef(p1)[["compliantCompliant"]] + 1.96 * se[["compliantCompliant"]]),
          avoided = sum(avoided[cc]), avoided_lo = sum(avoided[cc]) - 1.96 * se_av, avoided_hi = sum(avoided[cc]) + 1.96 * se_av, pct_of_cost = 100 * sum(avoided[cc]) / sum(k$cost))


# step 3  both-source towns, run B: damage from the VTrans records alone
# the same lines as step 2 with vt_cost in place of cost_dedup; a segment only FEMA recorded counts as undamaged
cat("\nstep 3 both-source towns, VTrans records only\n\n")
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
cat("\npart 1, every term: B, clustered SE, Exp(B) = odds ratio, clustered 95% interval\n")
print(round(cbind(B = coef(p1), SE_clustered = se, Exp_B = exp(coef(p1)),
                  lo_95 = exp(coef(p1) - 1.96 * se), hi_95 = exp(coef(p1) + 1.96 * se)), 3))
n1 <- sum(d$damaged); n0 <- sum(d$damaged == 0)
cat("AUC:", round((sum(rank(fitted(p1))[d$damaged == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0), 4), "\n")
p2 <- glm2(as.formula(paste("cost ~", PART2)), family = Gamma(link = "log"), data = k)
se2 <- sqrt(diag(vcovCL(p2, cluster = k$Town, type = "HC1")))
cat("\npart 2, every term: B, clustered SE, Exp(B) = cost ratio, clustered 95% interval\n")
print(round(cbind(B = coef(p2), SE_clustered = se2, Exp_B = exp(coef(p2)),
                  lo_95 = exp(coef(p2) - 1.96 * se2), hi_95 = exp(coef(p2) + 1.96 * se2)), 3))
cat("deviance R2:", round(1 - p2$deviance / p2$null.deviance, 4), "\n")

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
cat("\navoided cost, these towns:", format(round(sum(avoided[cc])), big.mark = ","), " clustered 95%:", format(round(sum(avoided[cc]) - 1.96 * se_av), big.mark = ","), "to", format(round(sum(avoided[cc]) + 1.96 * se_av), big.mark = ","), "\n")
cat("as a share of the observed damage cost:", round(100 * sum(avoided[cc]) / sum(k$cost), 1), "%\n")
runB <- c(segments = nrow(d), damaged = sum(d$damaged), odds_ratio = exp(coef(p1)[["compliantCompliant"]]),
          lo = exp(coef(p1)[["compliantCompliant"]] - 1.96 * se[["compliantCompliant"]]), hi = exp(coef(p1)[["compliantCompliant"]] + 1.96 * se[["compliantCompliant"]]),
          avoided = sum(avoided[cc]), avoided_lo = sum(avoided[cc]) - 1.96 * se_av, avoided_hi = sum(avoided[cc]) + 1.96 * se_av, pct_of_cost = 100 * sum(avoided[cc]) / sum(k$cost))


# step 4  side by side
# runs A and B share towns, segments and terms; compare the odds ratio and the avoided cost as a share of each
# run's own damage cost; the raw dollars are not comparable because VTrans dollars are a subset of the pooled ones
cat("\nstep 4 side by side,", STORM, "storm\n\n")
side <- rbind(state, runA, runB)
rownames(side) <- c("statewide, pooled", "both-source towns, pooled (A)", "both-source towns, VTrans only (B)")
print(data.frame(segments   = side[, "segments"], damaged = side[, "damaged"],
                 odds_ratio = round(side[, "odds_ratio"], 3),
                 clustered_95 = paste0(round(side[, "lo"], 3), " to ", round(side[, "hi"], 3)),
                 avoided    = format(round(side[, "avoided"]), big.mark = ","),
                 avoided_95 = paste0(format(round(side[, "avoided_lo"]), big.mark = ","), " to ", format(round(side[, "avoided_hi"]), big.mark = ",")),
                 pct_of_cost = round(side[, "pct_of_cost"], 1)))
