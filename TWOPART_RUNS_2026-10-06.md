# Two-Part Model Runs

The selected model. Graded segments, Precip >= 3.51 in, Averill excluded. Damage and cost = cost_dedup from analysis_{storm}_2026-09-27.csv, uncapped. Standard errors clustered by town. Script run once with STORM = 2023 and once with STORM = 2024.

## Script: twopart_model.R

```r
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
```

## Output 2023

```
PART 1 2023  segments 59695  damaged 1457  towns 204 

Call:
glm(formula = damaged ~ compliant + MeanSlope90m + StreamOrder + 
    PARENT + HYDROGROUP + PercentImpervious_BaseLC_90m + Precip + 
    RoadGrade_mean_deg + Surface + Culvert_Any + Driveway_Count + 
    Mean_TopoConvergence_30m + vtrans_district + road_miles_municipal + 
    grand_list_equalized_muni_100M, family = binomial, data = d)

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
AUC
[1] 0.8244
PART 2 2023 

Call:
glm2(formula = cost_dedup ~ MeanSlope90m + StreamCrossing + StreamOrder + 
    PARENT + K_final + HYDROGROUP + PercentImpervious_BaseLC_90m + 
    Precip + Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m + 
    vtrans_district + road_miles_municipal + grand_list_equalized_muni_100M, 
    family = Gamma(link = "log"), data = k)

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
deviance R2
[1] 0.2666
AVOIDED COST 2023 
compliant segments recoded
[1] 50175
statewide avoided cost
[1] 11492123
ten sample towns
 Brattleboro      Corinth  Wallingford     Hardwick West Windsor      Wolcott 
       40441        59487       114854        69499          929       175214 
  Washington     Richmond     Stamford   Starksboro 
      196540        14877         8352        42818 
```

## Output 2024

```
PART 1 2024  segments 22571  damaged 1018  towns 94 

Call:
glm(formula = damaged ~ compliant + MeanSlope90m + StreamOrder + 
    PARENT + HYDROGROUP + PercentImpervious_BaseLC_90m + Precip + 
    RoadGrade_mean_deg + Surface + Culvert_Any + Driveway_Count + 
    Mean_TopoConvergence_30m + vtrans_district + road_miles_municipal + 
    grand_list_equalized_muni_100M, family = binomial, data = d)

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
AUC
[1] 0.7782
PART 2 2024 
Warning message:
step size truncated due to increasing deviance 

Call:
glm2(formula = cost_dedup ~ MeanSlope90m + StreamCrossing + StreamOrder + 
    PARENT + K_final + HYDROGROUP + PercentImpervious_BaseLC_90m + 
    Precip + Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m + 
    vtrans_district + road_miles_municipal + grand_list_equalized_muni_100M, 
    family = Gamma(link = "log"), data = k)

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
deviance R2
[1] 0.2662
AVOIDED COST 2024 
compliant segments recoded
[1] 18123
statewide avoided cost
[1] 12714930
ten sample towns
      <NA>       <NA>       <NA>   Hardwick       <NA>    Wolcott       <NA> 
        NA         NA         NA     222483         NA     180992         NA 
  Richmond       <NA> Starksboro 
    206667         NA     238569 
```

# Both-source towns: pooled vs VTrans only

The same model on only the towns that hold both a FEMA and a VTrans damage record in that storm, fitted twice on the same segments: once on the pooled de-duplicated cost (OUTCOME = cost_dedup) and once on the VTrans records alone (OUTCOME = vt_cost). The only thing that changes between the two runs is which records count as damage. The avoided cost carries the same town-clustered delta-method interval as the deck. Script run four times: each storm with each OUTCOME.

## Script: both_source_towns.R

```r
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
```

## Summary of the four runs

```
SUMMARY 2023 cost_dedup 
       towns      damaged   odds_ratio           lo           hi          AUC 
      24.000      667.000        0.889        0.713        1.110        0.749 
      dev_r2      avoided   avoided_lo   avoided_hi 
       0.308  2258820.000 -2200659.000  6718298.000 
```

```
SUMMARY 2023 vt_cost 
       towns      damaged   odds_ratio           lo           hi          AUC 
      24.000      261.000        1.117        0.822        1.519        0.793 
      dev_r2      avoided   avoided_lo   avoided_hi 
       0.389  -790326.000 -2947316.000  1366663.000 
```

```
SUMMARY 2024 cost_dedup 
       towns      damaged   odds_ratio           lo           hi          AUC 
      26.000      747.000        0.705        0.552        0.901        0.744 
      dev_r2      avoided   avoided_lo   avoided_hi 
       0.194 10490094.000   911257.000 20068931.000 
```

```
SUMMARY 2024 vt_cost 
       towns      damaged   odds_ratio           lo           hi          AUC 
      26.000      262.000        0.618        0.474        0.807        0.799 
      dev_r2      avoided   avoided_lo   avoided_hi 
       0.324  5926812.000  1612337.000 10241287.000 
```

## Output 2023, pooled de-duplicated cost

```
BOTH-SOURCE TOWNS 2023 cost_dedup 
towns
[1] 24
segments
[1] 8641
damaged
[1] 667
total cost
[1] 29577443
PART 1

Call:
glm(formula = damaged ~ compliant + MeanSlope90m + StreamOrder + 
    PARENT + HYDROGROUP + PercentImpervious_BaseLC_90m + Precip + 
    RoadGrade_mean_deg + Surface + Culvert_Any + Driveway_Count + 
    Mean_TopoConvergence_30m + vtrans_district + road_miles_municipal + 
    grand_list_equalized_muni_100M, family = binomial, data = d)

Coefficients:
                                 Estimate Std. Error z value Pr(>|z|)    
(Intercept)                    -5.7635061  0.6017820  -9.577  < 2e-16 ***
compliantCompliant             -0.1173741  0.0995292  -1.179 0.238282    
MeanSlope90m                    0.0473108  0.0123264   3.838 0.000124 ***
StreamOrder                     0.3103634  0.0421535   7.363 1.80e-13 ***
PARENTDT                        0.2361871  0.2941640   0.803 0.422027    
PARENTDT/GT                     0.6104717  0.5765367   1.059 0.289664    
PARENTGF                       -0.4998556  0.2734275  -1.828 0.067533 .  
PARENTGL                       -0.2362696  0.3505452  -0.674 0.500307    
PARENTGT                       -0.4631921  0.2653646  -1.745 0.080899 .  
HYDROGROUPA/D                  -0.5207669  0.6222171  -0.837 0.402619    
HYDROGROUPB                    -0.2370073  0.2141567  -1.107 0.268423    
HYDROGROUPB/D                  -0.1761800  0.3012001  -0.585 0.558597    
HYDROGROUPC                    -0.4445239  0.2254312  -1.972 0.048623 *  
HYDROGROUPC/D                  -0.7417464  0.2755274  -2.692 0.007100 ** 
HYDROGROUPD                    -0.3349818  0.2559923  -1.309 0.190683    
PercentImpervious_BaseLC_90m    0.5882629  0.7642427   0.770 0.441458    
Precip                          0.0466962  0.0454210   1.028 0.303914    
RoadGrade_mean_deg              0.0850051  0.0239152   3.554 0.000379 ***
SurfaceUnpaved                  0.3604836  0.1288777   2.797 0.005156 ** 
Culvert_Any                     0.5709837  0.0905150   6.308 2.82e-10 ***
Driveway_Count                  0.0333406  0.0657207   0.507 0.611939    
Mean_TopoConvergence_30m        0.2682358  0.0442156   6.067 1.31e-09 ***
vtrans_district1               -0.3245800  0.3635464  -0.893 0.371956    
vtrans_district2                0.7982949  0.1238809   6.444 1.16e-10 ***
vtrans_district3               -0.7585243  0.2890102  -2.625 0.008676 ** 
vtrans_district4               -0.7406319  0.1995773  -3.711 0.000206 ***
vtrans_district5                0.1288874  0.2579293   0.500 0.617286    
vtrans_district7               -0.4514389  0.2811330  -1.606 0.108321    
vtrans_district8               -1.4178808  0.3065400  -4.625 3.74e-06 ***
road_miles_municipal            0.0009416  0.0027054   0.348 0.727797    
grand_list_equalized_muni_100M -0.0665756  0.0268685  -2.478 0.013218 *  
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 4698.2  on 8640  degrees of freedom
Residual deviance: 4180.3  on 8610  degrees of freedom
AIC: 4242.3

Number of Fisher Scoring iterations: 6

compliance odds ratio, clustered 95%
odds_ratio         lo         hi 
     0.889      0.713      1.110 
AUC
[1] 0.7491
PART 2

Call:
glm2(formula = cost ~ MeanSlope90m + StreamCrossing + StreamOrder + 
    PARENT + K_final + HYDROGROUP + PercentImpervious_BaseLC_90m + 
    Precip + Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m + 
    vtrans_district + road_miles_municipal + grand_list_equalized_muni_100M, 
    family = Gamma(link = "log"), data = k)

Coefficients:
                                Estimate Std. Error t value Pr(>|t|)    
(Intercept)                     9.122501   1.177190   7.749 3.67e-14 ***
MeanSlope90m                    0.060740   0.023071   2.633  0.00867 ** 
StreamCrossing                 -0.466074   0.348814  -1.336  0.18197    
StreamOrder                     0.560680   0.140780   3.983 7.60e-05 ***
PARENTDT                        0.700326   0.576433   1.215  0.22484    
PARENTDT/GT                     2.115465   1.094304   1.933  0.05366 .  
PARENTGF                        0.232856   0.553256   0.421  0.67398    
PARENTGL                        0.610237   0.689833   0.885  0.37670    
PARENTGT                       -0.376259   0.521238  -0.722  0.47065    
K_final                         0.370659   0.940385   0.394  0.69360    
HYDROGROUPA/D                  -0.395455   1.240542  -0.319  0.75000    
HYDROGROUPB                     0.353181   0.438734   0.805  0.42112    
HYDROGROUPB/D                  -0.308940   0.607423  -0.509  0.61120    
HYDROGROUPC                    -0.373439   0.460809  -0.810  0.41801    
HYDROGROUPC/D                  -0.484635   0.548444  -0.884  0.37722    
HYDROGROUPD                    -0.935536   0.542436  -1.725  0.08507 .  
PercentImpervious_BaseLC_90m   -1.939952   1.665015  -1.165  0.24440    
Precip                          0.015932   0.096539   0.165  0.86897    
SurfaceUnpaved                 -0.512585   0.248551  -2.062  0.03959 *  
Culvert_Any                    -0.035854   0.173057  -0.207  0.83594    
Driveway_Count                 -0.157643   0.137778  -1.144  0.25298    
Mean_TopoConvergence_30m        0.140326   0.083074   1.689  0.09168 .  
vtrans_district1               -0.694910   0.718036  -0.968  0.33352    
vtrans_district2               -1.048682   0.238516  -4.397 1.29e-05 ***
vtrans_district3               -0.869043   0.575845  -1.509  0.13175    
vtrans_district4               -1.821015   0.387590  -4.698 3.22e-06 ***
vtrans_district5                0.666863   0.552930   1.206  0.22825    
vtrans_district7               -0.188323   0.543473  -0.347  0.72907    
vtrans_district8                0.958414   0.626513   1.530  0.12657    
road_miles_municipal            0.007666   0.005589   1.372  0.17065    
grand_list_equalized_muni_100M -0.142469   0.054513  -2.614  0.00917 ** 
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

(Dispersion parameter for Gamma family taken to be 4.033662)

    Null deviance: 2016.9  on 666  degrees of freedom
Residual deviance: 1395.8  on 636  degrees of freedom
AIC: 14919

Number of Fisher Scoring iterations: 18

deviance R2
[1] 0.3079
AVOIDED COST
compliant segments recoded
[1] 6926
avoided cost, these towns
[1] 2258820
clustered SE
[1] 2275244
clustered 95% interval
[1] -2200659  6718298
study towns
       <NA>        <NA> Wallingford    Hardwick        <NA>     Wolcott 
         NA          NA       40173       77545          NA      121610 
 Washington    Richmond    Stamford        <NA> 
     252246       44755       49184          NA 
SUMMARY 2023 cost_dedup 
       towns      damaged   odds_ratio           lo           hi          AUC 
      24.000      667.000        0.889        0.713        1.110        0.749 
      dev_r2      avoided   avoided_lo   avoided_hi 
       0.308  2258820.000 -2200659.000  6718298.000 
```

## Output 2023, VTrans only

```
BOTH-SOURCE TOWNS 2023 vt_cost 
towns
[1] 24
segments
[1] 8641
damaged
[1] 261
total cost
[1] 9223297
PART 1

Call:
glm(formula = damaged ~ compliant + MeanSlope90m + StreamOrder + 
    PARENT + HYDROGROUP + PercentImpervious_BaseLC_90m + Precip + 
    RoadGrade_mean_deg + Surface + Culvert_Any + Driveway_Count + 
    Mean_TopoConvergence_30m + vtrans_district + road_miles_municipal + 
    grand_list_equalized_muni_100M, family = binomial, data = d)

Coefficients:
                                Estimate Std. Error z value Pr(>|z|)    
(Intercept)                    -7.008916   0.929283  -7.542 4.62e-14 ***
compliantCompliant              0.110748   0.161424   0.686 0.492668    
MeanSlope90m                    0.088515   0.017714   4.997 5.83e-07 ***
StreamOrder                     0.294380   0.063217   4.657 3.21e-06 ***
PARENTDT                        0.924831   0.502541   1.840 0.065723 .  
PARENTDT/GT                     1.486395   0.892127   1.666 0.095688 .  
PARENTGF                       -0.438971   0.482109  -0.911 0.362547    
PARENTGL                        0.988011   0.569382   1.735 0.082700 .  
PARENTGT                        0.335945   0.475996   0.706 0.480329    
HYDROGROUPB                    -1.183474   0.321661  -3.679 0.000234 ***
HYDROGROUPB/D                  -1.737017   0.652479  -2.662 0.007764 ** 
HYDROGROUPC                    -1.297569   0.323458  -4.012 6.03e-05 ***
HYDROGROUPC/D                  -1.389748   0.383166  -3.627 0.000287 ***
HYDROGROUPD                    -0.934766   0.342228  -2.731 0.006306 ** 
PercentImpervious_BaseLC_90m   -0.162324   1.302899  -0.125 0.900850    
Precip                          0.031653   0.064394   0.492 0.623035    
RoadGrade_mean_deg              0.083968   0.035361   2.375 0.017569 *  
SurfaceUnpaved                 -0.268911   0.175683  -1.531 0.125852    
Culvert_Any                     0.357121   0.138787   2.573 0.010077 *  
Driveway_Count                  0.119620   0.098286   1.217 0.223582    
Mean_TopoConvergence_30m        0.331922   0.067408   4.924 8.48e-07 ***
vtrans_district1                0.241803   0.434975   0.556 0.578280    
vtrans_district2                1.360119   0.185388   7.337 2.19e-13 ***
vtrans_district3               -0.255235   0.413550  -0.617 0.537116    
vtrans_district4               -0.782046   0.327846  -2.385 0.017060 *  
vtrans_district5                0.324414   0.368905   0.879 0.379185    
vtrans_district7               -1.917400   1.010715  -1.897 0.057818 .  
vtrans_district8               -0.797424   0.440245  -1.811 0.070092 .  
road_miles_municipal           -0.001875   0.004185  -0.448 0.654061    
grand_list_equalized_muni_100M -0.112589   0.047668  -2.362 0.018178 *  
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 2340.9  on 8640  degrees of freedom
Residual deviance: 2037.5  on 8611  degrees of freedom
AIC: 2097.5

Number of Fisher Scoring iterations: 8

compliance odds ratio, clustered 95%
odds_ratio         lo         hi 
     1.117      0.822      1.519 
AUC
[1] 0.7925
PART 2
Warning messages:
1: step size truncated due to increasing deviance 
2: step size truncated due to increasing deviance 
3: step size truncated due to increasing deviance 

Call:
glm2(formula = cost ~ MeanSlope90m + StreamCrossing + StreamOrder + 
    PARENT + K_final + HYDROGROUP + PercentImpervious_BaseLC_90m + 
    Precip + Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m + 
    vtrans_district + road_miles_municipal + grand_list_equalized_muni_100M, 
    family = Gamma(link = "log"), data = k)

Coefficients:
                                Estimate Std. Error t value Pr(>|t|)    
(Intercept)                     9.001988   1.725200   5.218 4.02e-07 ***
MeanSlope90m                    0.010976   0.030058   0.365  0.71533    
StreamCrossing                 -0.182536   0.499896  -0.365  0.71534    
StreamOrder                     0.130051   0.191353   0.680  0.49741    
PARENTDT                        1.096335   0.955108   1.148  0.25221    
PARENTDT/GT                     1.595954   1.589666   1.004  0.31645    
PARENTGF                        1.951232   0.947253   2.060  0.04053 *  
PARENTGL                        0.597287   1.038767   0.575  0.56585    
PARENTGT                       -0.568927   0.872089  -0.652  0.51481    
K_final                         0.626557   1.335951   0.469  0.63951    
HYDROGROUPB                     2.013688   0.733026   2.747  0.00649 ** 
HYDROGROUPB/D                  -0.026914   1.322491  -0.020  0.98378    
HYDROGROUPC                     1.011973   0.634184   1.596  0.11192    
HYDROGROUPC/D                   1.307375   0.740717   1.765  0.07888 .  
HYDROGROUPD                    -0.029148   0.700751  -0.042  0.96686    
PercentImpervious_BaseLC_90m   -3.159320   2.354528  -1.342  0.18098    
Precip                         -0.040406   0.132311  -0.305  0.76035    
SurfaceUnpaved                 -0.215209   0.301756  -0.713  0.47645    
Culvert_Any                     0.252203   0.239654   1.052  0.29373    
Driveway_Count                 -0.588517   0.195665  -3.008  0.00292 ** 
Mean_TopoConvergence_30m        0.075214   0.113278   0.664  0.50736    
vtrans_district1               -1.724548   0.790371  -2.182  0.03012 *  
vtrans_district2               -2.821714   0.367299  -7.682 4.42e-13 ***
vtrans_district3               -0.816872   0.725462  -1.126  0.26133    
vtrans_district4               -1.698382   0.566620  -2.997  0.00302 ** 
vtrans_district5               -0.553038   0.743260  -0.744  0.45759    
vtrans_district7                1.168720   2.020267   0.578  0.56349    
vtrans_district8               -0.592688   0.813762  -0.728  0.46715    
road_miles_municipal            0.013506   0.008041   1.680  0.09438 .  
grand_list_equalized_muni_100M -0.192617   0.089032  -2.163  0.03153 *  
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

(Dispersion parameter for Gamma family taken to be 2.834837)

    Null deviance: 935.24  on 260  degrees of freedom
Residual deviance: 571.17  on 231  degrees of freedom
AIC: 5625.3

Number of Fisher Scoring iterations: 12

deviance R2
[1] 0.3893
AVOIDED COST
compliant segments recoded
[1] 6926
avoided cost, these towns
[1] -790326
clustered SE
[1] 1100505
clustered 95% interval
[1] -2947316  1366663
study towns
       <NA>        <NA> Wallingford    Hardwick        <NA>     Wolcott 
         NA          NA      -20879      -26204          NA      -54291 
 Washington    Richmond    Stamford        <NA> 
    -127674       -5306      -11387          NA 
SUMMARY 2023 vt_cost 
       towns      damaged   odds_ratio           lo           hi          AUC 
      24.000      261.000        1.117        0.822        1.519        0.793 
      dev_r2      avoided   avoided_lo   avoided_hi 
       0.389  -790326.000 -2947316.000  1366663.000 
```

## Output 2024, pooled de-duplicated cost

```
BOTH-SOURCE TOWNS 2024 cost_dedup 
towns
[1] 26
segments
[1] 9473
damaged
[1] 747
total cost
[1] 42092428
PART 1

Call:
glm(formula = damaged ~ compliant + MeanSlope90m + StreamOrder + 
    PARENT + HYDROGROUP + PercentImpervious_BaseLC_90m + Precip + 
    RoadGrade_mean_deg + Surface + Culvert_Any + Driveway_Count + 
    Mean_TopoConvergence_30m + vtrans_district + road_miles_municipal + 
    grand_list_equalized_muni_100M, family = binomial, data = d)

Coefficients:
                                Estimate Std. Error z value Pr(>|z|)    
(Intercept)                    -5.941380   0.633277  -9.382  < 2e-16 ***
compliantCompliant             -0.349819   0.089867  -3.893 9.92e-05 ***
MeanSlope90m                    0.071721   0.011651   6.156 7.48e-10 ***
StreamOrder                     0.283632   0.039737   7.138 9.48e-13 ***
PARENTDT                        0.251904   0.274097   0.919 0.358079    
PARENTDT/GT                     0.388240   0.424744   0.914 0.360687    
PARENTGF                        0.030092   0.276721   0.109 0.913406    
PARENTGF/GL                     0.802140   1.198426   0.669 0.503287    
PARENTGL                        0.062495   0.259707   0.241 0.809836    
PARENTGT                        0.050909   0.246105   0.207 0.836118    
PARENTGT/DT                    -1.847844   0.564264  -3.275 0.001057 ** 
PARENTGT/GT_R                   0.028458   1.077737   0.026 0.978934    
PARENTM                         1.954810   0.961454   2.033 0.042034 *  
HYDROGROUPA/D                  -0.301157   0.355106  -0.848 0.396396    
HYDROGROUPB                    -0.091362   0.254945  -0.358 0.720074    
HYDROGROUPB/D                  -0.129642   0.317607  -0.408 0.683139    
HYDROGROUPC                    -0.496965   0.282205  -1.761 0.078236 .  
HYDROGROUPC/D                  -0.403981   0.303201  -1.332 0.182733    
HYDROGROUPD                    -0.434404   0.294268  -1.476 0.139885    
HYDROGROUPnot rated            -1.186214   0.990486  -1.198 0.231070    
PercentImpervious_BaseLC_90m    1.566727   0.655339   2.391 0.016816 *  
Precip                         -0.052566   0.087419  -0.601 0.547634    
RoadGrade_mean_deg              0.073328   0.022063   3.324 0.000889 ***
SurfaceUnpaved                  0.708616   0.116464   6.084 1.17e-09 ***
Culvert_Any                     0.608915   0.086287   7.057 1.70e-12 ***
Driveway_Count                  0.099800   0.048745   2.047 0.040621 *  
Mean_TopoConvergence_30m        0.363108   0.042177   8.609  < 2e-16 ***
vtrans_district5               -0.067667   0.115760  -0.585 0.558852    
vtrans_district7                0.396191   0.167301   2.368 0.017878 *  
road_miles_municipal           -0.014768   0.002668  -5.535 3.12e-08 ***
grand_list_equalized_muni_100M  0.010059   0.004265   2.358 0.018366 *  
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 5228.4  on 9472  degrees of freedom
Residual deviance: 4711.2  on 9442  degrees of freedom
AIC: 4773.2

Number of Fisher Scoring iterations: 6

compliance odds ratio, clustered 95%
odds_ratio         lo         hi 
     0.705      0.552      0.901 
AUC
[1] 0.7438
PART 2
Warning messages:
1: step size truncated due to increasing deviance 
2: step size truncated due to increasing deviance 
3: step size truncated due to increasing deviance 

Call:
glm2(formula = cost ~ MeanSlope90m + StreamCrossing + StreamOrder + 
    PARENT + K_final + HYDROGROUP + PercentImpervious_BaseLC_90m + 
    Precip + Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m + 
    vtrans_district + road_miles_municipal + grand_list_equalized_muni_100M, 
    family = Gamma(link = "log"), data = k)

Coefficients:
                                Estimate Std. Error t value Pr(>|t|)    
(Intercept)                    10.748312   1.219501   8.814  < 2e-16 ***
MeanSlope90m                    0.050750   0.024799   2.046 0.041076 *  
StreamCrossing                 -0.351821   0.349229  -1.007 0.314072    
StreamOrder                     0.529825   0.143009   3.705 0.000228 ***
PARENTDT                       -0.288908   0.561097  -0.515 0.606783    
PARENTDT/GT                    -1.155433   0.873695  -1.322 0.186435    
PARENTGF                        0.032308   0.464062   0.070 0.944516    
PARENTGF/GL                    -1.959306   2.126754  -0.921 0.357222    
PARENTGL                       -0.632956   0.524480  -1.207 0.227898    
PARENTGT                       -0.437028   0.503724  -0.868 0.385907    
PARENTGT/DT                     1.184707   1.198099   0.989 0.323084    
PARENTGT/GT_R                  -2.041367   2.249824  -0.907 0.364529    
PARENTM                        -0.467626   1.396806  -0.335 0.737887    
K_final                         1.147690   1.084065   1.059 0.290098    
HYDROGROUPA/D                   0.056385   0.706044   0.080 0.936371    
HYDROGROUPB                     0.310304   0.460706   0.674 0.500822    
HYDROGROUPB/D                  -0.495546   0.599749  -0.826 0.408934    
HYDROGROUPC                    -0.372036   0.556277  -0.669 0.503841    
HYDROGROUPC/D                   0.230887   0.614311   0.376 0.707142    
HYDROGROUPD                    -0.072268   0.623309  -0.116 0.907730    
HYDROGROUPnot rated             1.608972   1.551061   1.037 0.299929    
PercentImpervious_BaseLC_90m   -0.368237   1.384572  -0.266 0.790349    
Precip                          0.001267   0.173086   0.007 0.994161    
SurfaceUnpaved                 -0.406511   0.240113  -1.693 0.090890 .  
Culvert_Any                     0.220269   0.179472   1.227 0.220107    
Driveway_Count                 -0.044020   0.106898  -0.412 0.680616    
Mean_TopoConvergence_30m        0.041226   0.084768   0.486 0.626879    
vtrans_district5               -0.795779   0.240389  -3.310 0.000978 ***
vtrans_district7                1.240777   0.328496   3.777 0.000172 ***
road_miles_municipal           -0.021602   0.005318  -4.062 5.41e-05 ***
grand_list_equalized_muni_100M  0.021747   0.009073   2.397 0.016788 *  
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

(Dispersion parameter for Gamma family taken to be 4.608719)

    Null deviance: 2059.3  on 746  degrees of freedom
Residual deviance: 1660.1  on 716  degrees of freedom
AIC: 17316

Number of Fisher Scoring iterations: 18

deviance R2
[1] 0.1938
AVOIDED COST
compliant segments recoded
[1] 7435
avoided cost, these towns
[1] 10490094
clustered SE
[1] 4887162
clustered 95% interval
[1]   911257 20068931
study towns
      <NA>       <NA>       <NA>   Hardwick       <NA>    Wolcott       <NA> 
        NA         NA         NA     705430         NA     233715         NA 
  Richmond       <NA> Starksboro 
    400347         NA     403974 
SUMMARY 2024 cost_dedup 
       towns      damaged   odds_ratio           lo           hi          AUC 
      26.000      747.000        0.705        0.552        0.901        0.744 
      dev_r2      avoided   avoided_lo   avoided_hi 
       0.194 10490094.000   911257.000 20068931.000 
```

## Output 2024, VTrans only

```
BOTH-SOURCE TOWNS 2024 vt_cost 
towns
[1] 26
segments
[1] 9473
damaged
[1] 262
total cost
[1] 15319997
PART 1

Call:
glm(formula = damaged ~ compliant + MeanSlope90m + StreamOrder + 
    PARENT + HYDROGROUP + PercentImpervious_BaseLC_90m + Precip + 
    RoadGrade_mean_deg + Surface + Culvert_Any + Driveway_Count + 
    Mean_TopoConvergence_30m + vtrans_district + road_miles_municipal + 
    grand_list_equalized_muni_100M, family = binomial, data = d)

Coefficients:
                                Estimate Std. Error z value Pr(>|z|)    
(Intercept)                    -4.584912   1.033167  -4.438 9.09e-06 ***
compliantCompliant             -0.480540   0.151813  -3.165 0.001549 ** 
MeanSlope90m                    0.100628   0.018280   5.505 3.70e-08 ***
StreamOrder                     0.384866   0.055790   6.898 5.26e-12 ***
PARENTDT                       -0.211702   0.389427  -0.544 0.586700    
PARENTDT/GT                    -0.425974   0.648531  -0.657 0.511291    
PARENTGF                       -0.691175   0.395663  -1.747 0.080659 .  
PARENTGF/GL                    -0.132932   1.428370  -0.093 0.925852    
PARENTGL                       -0.424946   0.363588  -1.169 0.242501    
PARENTGT                       -0.472640   0.346566  -1.364 0.172637    
PARENTGT/DT                    -1.953231   0.784642  -2.489 0.012799 *  
PARENTM                         0.132317   1.146309   0.115 0.908106    
HYDROGROUPA/D                  -0.368221   0.480009  -0.767 0.443015    
HYDROGROUPB                    -0.564899   0.384427  -1.469 0.141708    
HYDROGROUPB/D                  -0.395553   0.457669  -0.864 0.387435    
HYDROGROUPC                    -1.079006   0.440383  -2.450 0.014279 *  
HYDROGROUPC/D                  -0.768845   0.456174  -1.685 0.091907 .  
HYDROGROUPD                    -0.357419   0.435839  -0.820 0.412175    
HYDROGROUPnot rated            -0.204033   0.957284  -0.213 0.831219    
PercentImpervious_BaseLC_90m    0.029509   1.115595   0.026 0.978897    
Precip                         -0.409657   0.156910  -2.611 0.009034 ** 
RoadGrade_mean_deg             -0.023780   0.038292  -0.621 0.534589    
SurfaceUnpaved                  0.083646   0.162578   0.515 0.606902    
Culvert_Any                     0.699727   0.145762   4.800 1.58e-06 ***
Driveway_Count                  0.107623   0.078949   1.363 0.172816    
Mean_TopoConvergence_30m        0.411314   0.066819   6.156 7.48e-10 ***
vtrans_district5                0.705762   0.173132   4.076 4.57e-05 ***
vtrans_district7                0.212943   0.323256   0.659 0.510060    
road_miles_municipal           -0.020645   0.004738  -4.357 1.32e-05 ***
grand_list_equalized_muni_100M  0.024045   0.007004   3.433 0.000597 ***
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 2396.7  on 9472  degrees of freedom
Residual deviance: 2085.3  on 9443  degrees of freedom
AIC: 2145.3

Number of Fisher Scoring iterations: 7

compliance odds ratio, clustered 95%
odds_ratio         lo         hi 
     0.618      0.474      0.807 
AUC
[1] 0.7988
PART 2
Warning message:
step size truncated due to increasing deviance 

Call:
glm2(formula = cost ~ MeanSlope90m + StreamCrossing + StreamOrder + 
    PARENT + K_final + HYDROGROUP + PercentImpervious_BaseLC_90m + 
    Precip + Surface + Culvert_Any + Driveway_Count + Mean_TopoConvergence_30m + 
    vtrans_district + road_miles_municipal + grand_list_equalized_muni_100M, 
    family = Gamma(link = "log"), data = k)

Coefficients:
                                Estimate Std. Error t value Pr(>|t|)    
(Intercept)                     8.335150   1.574710   5.293 2.78e-07 ***
MeanSlope90m                    0.096267   0.030551   3.151  0.00184 ** 
StreamCrossing                 -0.476636   0.396492  -1.202  0.23054    
StreamOrder                     0.327674   0.156202   2.098  0.03701 *  
PARENTDT                       -0.210676   0.584224  -0.361  0.71872    
PARENTDT/GT                     0.484069   1.037157   0.467  0.64113    
PARENTGF                        0.110313   0.446440   0.247  0.80505    
PARENTGF/GL                    -2.932508   2.456585  -1.194  0.23380    
PARENTGL                       -0.414496   0.525733  -0.788  0.43126    
PARENTGT                        0.074707   0.533868   0.140  0.88883    
PARENTGT/DT                    -2.675828   1.269476  -2.108  0.03612 *  
PARENTM                        -1.872780   2.145322  -0.873  0.38359    
K_final                         1.820946   1.417339   1.285  0.20016    
HYDROGROUPA/D                  -0.205890   0.749759  -0.275  0.78386    
HYDROGROUPB                    -0.116571   0.535709  -0.218  0.82793    
HYDROGROUPB/D                  -1.015999   0.624243  -1.628  0.10497    
HYDROGROUPC                    -0.876670   0.694023  -1.263  0.20780    
HYDROGROUPC/D                  -0.170979   0.736165  -0.232  0.81654    
HYDROGROUPD                    -0.469762   0.731723  -0.642  0.52151    
HYDROGROUPnot rated             0.886362   1.761025   0.503  0.61522    
PercentImpervious_BaseLC_90m    5.425926   2.219432   2.445  0.01524 *  
Precip                          0.105597   0.236760   0.446  0.65601    
SurfaceUnpaved                  0.141167   0.255458   0.553  0.58106    
Culvert_Any                     0.255243   0.242519   1.052  0.29368    
Driveway_Count                 -0.158423   0.143156  -1.107  0.26959    
Mean_TopoConvergence_30m        0.077072   0.109135   0.706  0.48077    
vtrans_district5               -0.929228   0.279187  -3.328  0.00102 ** 
vtrans_district7                1.465181   0.492015   2.978  0.00321 ** 
road_miles_municipal           -0.006416   0.006656  -0.964  0.33609    
grand_list_equalized_muni_100M -0.009623   0.010750  -0.895  0.37164    
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

(Dispersion parameter for Gamma family taken to be 2.46314)

    Null deviance: 720.54  on 261  degrees of freedom
Residual deviance: 487.42  on 232  degrees of freedom
AIC: 6074.6

Number of Fisher Scoring iterations: 17

deviance R2
[1] 0.3235
AVOIDED COST
compliant segments recoded
[1] 7435
avoided cost, these towns
[1] 5926812
clustered SE
[1] 2201263
clustered 95% interval
[1]  1612337 10241287
study towns
      <NA>       <NA>       <NA>   Hardwick       <NA>    Wolcott       <NA> 
        NA         NA         NA     399424         NA     152250         NA 
  Richmond       <NA> Starksboro 
    310168         NA     152020 
SUMMARY 2024 vt_cost 
       towns      damaged   odds_ratio           lo           hi          AUC 
      26.000      262.000        0.618        0.474        0.807        0.799 
      dev_r2      avoided   avoided_lo   avoided_hi 
       0.324  5926812.000  1612337.000 10241287.000 
```

