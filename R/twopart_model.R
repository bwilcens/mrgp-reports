# twopart_model.R -- the selected two-part model: Part 1 logit (did the segment flood), Part 2 Gamma log-link
#   (what it cost), then the avoided cost from recoding every compliant segment to Does Not Meet.
# Set STORM to 2023 or 2024 and run from mrgp-roads-core:
#   "C:/Program Files/R/R-4.5.3/bin/Rscript.exe" twopart_simple_2026-09-17/R/twopart_model.R
# Graded segments, Precip >= 3.51 in, Averill excluded; damage and cost from cost_dedup, uncapped; SEs clustered by town.

STORM <- 2023

library(glm2)
library(sandwich)

d <- read.csv(paste0("twopart_simple_2026-09-17/data/analysis_", STORM, "_2026-09-27.csv"))
d <- d[d$Precip >= 3.51 & d$Town != "Averill", ]
d$damaged   <- as.integer(d$cost_dedup > 0)
d$compliant <- relevel(factor(d$compliant), ref = "Does Not Meet")

# a district with no damaged segment cannot be estimated: pool it into the district with the most
damaged_by_district <- tapply(d$damaged, d$vtrans_district, sum)
busiest <- names(which.max(damaged_by_district))
d$vtrans_district[d$vtrans_district %in% names(damaged_by_district)[damaged_by_district == 0]] <- busiest
d$vtrans_district <- relevel(factor(d$vtrans_district), ref = busiest)

cat("PART 1", STORM, " segments", nrow(d), " damaged", sum(d$damaged), " towns", length(unique(d$Town)), "\n")

p1 <- glm(damaged ~ compliant + MeanSlope90m + StreamOrder + PARENT + HYDROGROUP + PercentImpervious_BaseLC_90m +
            Precip + RoadGrade_mean_deg + Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m +
            vtrans_district + road_miles_municipal + grand_list_equalized_muni_100M,
          family = binomial, data = d)
print(summary(p1))

se <- sqrt(diag(vcovCL(p1, cluster = d$Town, type = "HC1")))
print(round(exp(cbind(odds_ratio = coef(p1), lo_95 = coef(p1) - 1.96 * se, hi_95 = coef(p1) + 1.96 * se)), 3))

n1 <- sum(d$damaged); n0 <- sum(d$damaged == 0)
cat("AUC\n"); print(round((sum(rank(fitted(p1))[d$damaged == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0), 4))

cat("PART 2", STORM, "\n")

k <- d[d$damaged == 1, ]
p2 <- glm2(cost_dedup ~ MeanSlope90m + StreamCrossing + StreamOrder + PARENT + K_final + HYDROGROUP +
             PercentImpervious_BaseLC_90m + Precip + Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m +
             vtrans_district + road_miles_municipal + grand_list_equalized_muni_100M,
           family = Gamma(link = "log"), data = k)
print(summary(p2))

se <- sqrt(diag(vcovCL(p2, cluster = k$Town, type = "HC1")))
print(round(exp(cbind(cost_ratio = coef(p2), lo_95 = coef(p2) - 1.96 * se, hi_95 = coef(p2) + 1.96 * se)), 3))

cat("deviance R2\n"); print(round(1 - p2$deviance / p2$null.deviance, 4))

cat("AVOIDED COST", STORM, "\n")

d0 <- d
d0$compliant <- factor("Does Not Meet", levels = levels(d$compliant))
p_observed <- predict(p1, newdata = d,  type = "response")   # probability of damage as rated
p_recoded  <- predict(p1, newdata = d0, type = "response")   # probability if rated Does Not Meet
cost       <- predict(p2, newdata = d,  type = "response")   # predicted repair cost
avoided    <- (p_recoded - p_observed) * cost
compliant  <- d$compliant == "Compliant"

cat("compliant segments recoded\n");  print(sum(compliant))
cat("statewide avoided cost\n");      print(round(sum(avoided[compliant])))
cat("ten sample towns\n")
print(round(tapply(avoided[compliant], d$Town[compliant], sum)[c("Brattleboro", "Corinth", "Wallingford", "Hardwick",
      "West Windsor", "Wolcott", "Washington", "Richmond", "Stamford", "Starksboro")]))
