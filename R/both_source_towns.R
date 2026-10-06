# both_source_towns.R -- the selected two-part model on the towns that hold BOTH a FEMA and a VTrans damage record,
#   fitted on the pooled de-duplicated cost (OUTCOME = "cost_dedup") or on the VTrans records alone ("vt_cost").
#   Same towns, same segments, same terms either way, so the two runs are comparable.
# Set STORM and OUTCOME and run from mrgp-roads-core:
#   "C:/Program Files/R/R-4.5.3/bin/Rscript.exe" twopart_simple_2026-09-17/R/both_source_towns.R
# Graded segments, Precip >= 3.51 in, Averill excluded; uncapped cost; SEs clustered by town.

STORM   <- 2023
OUTCOME <- "cost_dedup"      # "cost_dedup" = FEMA + VTrans pooled and de-duplicated; "vt_cost" = VTrans only

library(glm2)
library(sandwich)

d <- read.csv(paste0("twopart_simple_2026-09-17/data/analysis_", STORM, "_2026-09-27.csv"))
d <- d[d$Precip >= 3.51 & d$Town != "Averill", ]

fema_towns <- unique(d$Town[d$fema_cost > 0])
vt_towns   <- unique(d$Town[d$vt_cost > 0])
d <- d[d$Town %in% intersect(fema_towns, vt_towns), ]

d$damaged   <- as.integer(d[[OUTCOME]] > 0)
d$cost      <- d[[OUTCOME]]
d$compliant <- relevel(factor(d$compliant), ref = "Does Not Meet")

# a category with no damaged segment cannot be estimated: pool it into the category with the most
n <- tapply(d$damaged, d$vtrans_district, sum)
d$vtrans_district[d$vtrans_district %in% names(n)[n == 0]] <- names(which.max(n))
d$vtrans_district <- relevel(factor(d$vtrans_district), ref = names(which.max(n)))
n <- tapply(d$damaged, d$PARENT, sum)
d$PARENT[d$PARENT %in% names(n)[n == 0]] <- names(which.max(n))
d$PARENT <- factor(d$PARENT)
n <- tapply(d$damaged, d$HYDROGROUP, sum)
d$HYDROGROUP[d$HYDROGROUP %in% names(n)[n == 0]] <- names(which.max(n))
d$HYDROGROUP <- factor(d$HYDROGROUP)

cat("BOTH-SOURCE TOWNS", STORM, OUTCOME, "\n")
cat("towns\n");                      print(length(unique(d$Town)))
cat("segments\n");                   print(nrow(d))
cat("damaged\n");                    print(sum(d$damaged))
cat("total cost\n");                 print(sum(d$cost))

cat("PART 1\n")
p1 <- glm(damaged ~ compliant + MeanSlope90m + StreamOrder + PARENT + HYDROGROUP + PercentImpervious_BaseLC_90m +
            Precip + RoadGrade_mean_deg + Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m +
            vtrans_district + road_miles_municipal + grand_list_equalized_muni_100M,
          family = binomial, data = d)
print(summary(p1))
se <- sqrt(diag(vcovCL(p1, cluster = d$Town, type = "HC1")))[["compliantCompliant"]]
b  <- coef(p1)[["compliantCompliant"]]
cat("compliance odds ratio, clustered 95%\n"); print(round(exp(c(odds_ratio = b, lo = b - 1.96 * se, hi = b + 1.96 * se)), 3))
n1 <- sum(d$damaged); n0 <- sum(d$damaged == 0)
cat("AUC\n"); print(round((sum(rank(fitted(p1))[d$damaged == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0), 4))

cat("PART 2\n")
k <- d[d$damaged == 1, ]
p2 <- glm2(cost ~ MeanSlope90m + StreamCrossing + StreamOrder + PARENT + K_final + HYDROGROUP +
             PercentImpervious_BaseLC_90m + Precip + Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m +
             vtrans_district + road_miles_municipal + grand_list_equalized_muni_100M,
           family = Gamma(link = "log"), data = k)
print(summary(p2))
cat("deviance R2\n"); print(round(1 - p2$deviance / p2$null.deviance, 4))

cat("AVOIDED COST\n")
d0 <- d
d0$compliant <- factor("Does Not Meet", levels = levels(d$compliant))
p_observed <- predict(p1, newdata = d,  type = "response")   # probability of damage as rated
p_recoded  <- predict(p1, newdata = d0, type = "response")   # probability if rated Does Not Meet
cost_hat   <- predict(p2, newdata = d,  type = "response")   # predicted repair cost
avoided    <- (p_recoded - p_observed) * cost_hat
compliant  <- d$compliant == "Compliant"
cat("compliant segments recoded\n");  print(sum(compliant))
cat("avoided cost, these towns\n");   print(round(sum(avoided[compliant])))

# town-clustered standard error of the avoided cost, as in the deck script (R/105): the delta method with the
# joint HC1 sandwich of both parts, scores added within each town so the cross-part block is kept
k1 <- names(coef(p1))[!is.na(coef(p1))]
k2 <- names(coef(p2))[!is.na(coef(p2))]
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
a1 <- G1 / (G1 - 1) * (nrow(d) - 1) / (nrow(d) - length(k1))
a2 <- G2 / (G2 - 1) * (nrow(k) - 1) / (nrow(k) - length(k2))
V  <- crossprod(cbind(sqrt(a1) * U1 %*% B1, sqrt(a2) * U2 %*% B2))
se_avoided <- sqrt(drop(t(h) %*% V %*% h))
cat("clustered SE\n");                print(round(se_avoided))
cat("clustered 95% interval\n");      print(round(sum(avoided[compliant]) + c(-1.96, 1.96) * se_avoided))

cat("study towns\n")
print(round(tapply(avoided[compliant], d$Town[compliant], sum)[c("Brattleboro", "Corinth", "Wallingford", "Hardwick",
      "West Windsor", "Wolcott", "Washington", "Richmond", "Stamford", "Starksboro")]))

cat("SUMMARY", STORM, OUTCOME, "\n")
print(c(towns = length(unique(d$Town)), damaged = sum(d$damaged), odds_ratio = round(exp(b), 3),
        lo = round(exp(b - 1.96 * se), 3), hi = round(exp(b + 1.96 * se), 3),
        AUC = round((sum(rank(fitted(p1))[d$damaged == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0), 3),
        dev_r2 = round(1 - p2$deviance / p2$null.deviance, 3), avoided = round(sum(avoided[compliant])),
        avoided_lo = round(sum(avoided[compliant]) - 1.96 * se_avoided), avoided_hi = round(sum(avoided[compliant]) + 1.96 * se_avoided)))
