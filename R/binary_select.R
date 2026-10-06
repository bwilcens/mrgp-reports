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
