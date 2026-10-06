# Road Damage Category Screen

Graded segments, Precip >= 3.51 in, Averill excluded. Damage and cost = cost_dedup from analysis_{storm}_2026-09-27.csv, uncapped. Each script run once with STORM = 2023 and once with STORM = 2024.

## Script: R/114_crosstab_anova_simple.R

```r
# 114_crosstab_anova_simple.R -- crosstab of each category against damage, then ANOVA of ln(cost) by category.
# Set STORM to 2023 or 2024 and run from mrgp-roads-core:
#   "C:/Program Files/R/R-4.5.3/bin/Rscript.exe" twopart_simple_2026-09-17/R/114_crosstab_anova_simple.R
# Graded segments, Precip >= 3.51 in, Averill excluded; damage and cost from cost_dedup, uncapped.

STORM <- 2023

d <- read.csv(paste0("twopart_simple_2026-09-17/data/analysis_", STORM, "_2026-09-27.csv"))
d <- d[d$Precip >= 3.51 & d$Town != "Averill", ]
d$damaged <- as.integer(d$cost_dedup > 0)

cat("CROSSTABS", STORM, "\n")

tab <- table(d$PARENT, d$damaged)
print(tab)
print(chisq.test(tab))
print(round(chisq.test(tab)$stdres, 2))

tab <- table(d$HYDROGROUP, d$damaged)
print(tab)
print(chisq.test(tab))
print(round(chisq.test(tab)$stdres, 2))

tab <- table(d$vtrans_district, d$damaged)
print(tab)
print(chisq.test(tab))
print(round(chisq.test(tab)$stdres, 2))

cat("ANOVA", STORM, "\n")

k <- d[d$damaged == 1, ]
k$ln_cost <- log(k$cost_dedup)

print(round(tapply(k$ln_cost, k$PARENT, mean), 2))
print(summary(aov(ln_cost ~ PARENT, data = k)))
print(pairwise.t.test(k$ln_cost, k$PARENT, p.adjust.method = "bonferroni"))

print(round(tapply(k$ln_cost, k$HYDROGROUP, mean), 2))
print(summary(aov(ln_cost ~ HYDROGROUP, data = k)))
print(pairwise.t.test(k$ln_cost, k$HYDROGROUP, p.adjust.method = "bonferroni"))

print(round(tapply(k$ln_cost, k$vtrans_district, mean), 2))
print(summary(aov(ln_cost ~ vtrans_district, data = k)))
print(pairwise.t.test(k$ln_cost, k$vtrans_district, p.adjust.method = "bonferroni"))
```

## Crosstabs 2023

```
CROSSTABS 2023 
       
            0     1
  A      2389   115
  DT    23838   682
  DT/GT   169     8
  GF     9947   213
  GL     2456    43
  GT    17767   386
  GT/DT  1509     9
  M       163     1

	Pearson's Chi-squared test

data:  tab
X-squared = 105.97, df = 7, p-value < 2.2e-16

Warning message:
In chisq.test(tab) : Chi-squared approximation may be incorrect
       
            0     1
  A     -7.13  7.13
  DT    -4.50  4.50
  DT/GT -1.80  1.80
  GF     2.47 -2.47
  GL     2.38 -2.38
  GT     3.29 -3.29
  GT/DT  4.73 -4.73
  M      1.52 -1.52
Warning message:
In chisq.test(tab) : Chi-squared approximation may be incorrect
           
                0     1
  A          9350   206
  A/D         745     8
  B          8310   254
  B/D        1976    75
  C         14074   341
  C/D        9331   193
  D         14173   370
  not rated   279    10

	Pearson's Chi-squared test

data:  tab
X-squared = 41.045, df = 7, p-value = 7.936e-07

           
                0     1
  A          1.97 -1.97
  A/D        2.47 -2.47
  B         -3.40  3.40
  B/D       -3.63  3.63
  C          0.67 -0.67
  C/D        2.86 -2.86
  D         -0.93  0.93
  not rated -1.13  1.13
   
        0     1
  1  6046    33
  2  6069   520
  3 10360   147
  4 10425   137
  5  2645    35
  6  9650   505
  7  3617    41
  8  4196    25
  9  5230    14

	Pearson's Chi-squared test

data:  tab
X-squared = 1499.5, df = 8, p-value < 2.2e-16

   
         0      1
  1  10.12 -10.12
  2 -30.40  30.40
  3   7.62  -7.62
  4   8.40  -8.40
  5   3.90  -3.90
  6 -18.15  18.15
  7   5.34  -5.34
  8   8.07  -8.07
  9  10.68 -10.68
```

## Crosstabs 2024

```
CROSSTABS 2024 
         
             0    1
  A        930   75
  DT      8117  418
  DT/GT    123   11
  GF      2754  157
  GF/GL    101    3
  GL      3273  100
  GT      5437  237
  GT/DT    704   11
  GT/GT_R   89    2
  M         25    4

	Pearson's Chi-squared test

data:  tab
X-squared = 75.244, df = 9, p-value = 1.415e-12

Warning message:
In chisq.test(tab) : Chi-squared approximation may be incorrect
         
              0     1
  A       -4.61  4.61
  DT      -2.19  2.19
  DT/GT   -2.07  2.07
  GF      -2.46  2.46
  GF/GL    0.80 -0.80
  GL       4.69 -4.69
  GT       1.40 -1.40
  GT/DT    3.89 -3.89
  GT/GT_R  1.07 -1.07
  M       -2.41  2.41
Warning message:
In chisq.test(tab) : Chi-squared approximation may be incorrect
           
               0    1
  A         2267  137
  A/D        346   16
  B         2877  166
  B/D        776   48
  C         3901  146
  C/D       3330  190
  D         7850  306
  not rated  206    9

	Pearson's Chi-squared test

data:  tab
X-squared = 42.548, df = 7, p-value = 4.078e-07

           
                0     1
  A         -2.97  2.97
  A/D        0.08 -0.08
  B         -2.70  2.70
  B/D       -1.85  1.85
  C          3.05 -3.05
  C/D       -2.76  2.76
  D          4.13 -4.13
  not rated  0.23 -0.23
   
       0    1
  3  137    0
  5 5267  231
  6 5924  539
  7 7034  209
  8  174    0
  9 3017   39

	Pearson's Chi-squared test

data:  tab
X-squared = 354.59, df = 5, p-value < 2.2e-16

   
         0      1
  3   2.55  -2.55
  5   1.27  -1.27
  6 -17.56  17.56
  7   8.08  -8.08
  8   2.88  -2.88
  9   9.26  -9.26
```

## ANOVA 2023

```
ANOVA 2023 
    A    DT DT/GT    GF    GL    GT GT/DT     M 
 9.25  8.74  8.91  9.31  9.80  8.74  9.82 11.52 
              Df Sum Sq Mean Sq F value   Pr(>F)    
PARENT         7    127  18.160   5.287 5.56e-06 ***
Residuals   1449   4977   3.435                     
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

	Pairwise comparisons using t tests with pooled SD 

data:  k$ln_cost and k$PARENT 

      A      DT     DT/GT  GF     GL     GT     GT/DT 
DT    0.1700 -      -      -      -      -      -     
DT/GT 1.0000 1.0000 -      -      -      -      -     
GF    1.0000 0.0028 1.0000 -      -      -      -     
GL    1.0000 0.0074 1.0000 1.0000 -      -      -     
GT    0.2521 1.0000 1.0000 0.0092 0.0099 -      -     
GT/DT 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 -     
M     1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000

P value adjustment method: bonferroni 
        A       A/D         B       B/D         C       C/D         D not rated 
     9.33      9.42      8.64      9.21      8.73      8.79      8.96      9.82 
              Df Sum Sq Mean Sq F value   Pr(>F)    
HYDROGROUP     7     86  12.255   3.538 0.000892 ***
Residuals   1449   5018   3.463                     
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

	Pairwise comparisons using t tests with pooled SD 

data:  k$ln_cost and k$HYDROGROUP 

          A      A/D    B      B/D    C      C/D    D     
A/D       1.0000 -      -      -      -      -      -     
B         0.0026 1.0000 -      -      -      -      -     
B/D       1.0000 1.0000 0.5898 -      -      -      -     
C         0.0082 1.0000 1.0000 1.0000 -      -      -     
C/D       0.1162 1.0000 1.0000 1.0000 1.0000 -      -     
D         0.7189 1.0000 0.9363 1.0000 1.0000 1.0000 -     
not rated 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000 1.0000

P value adjustment method: bonferroni 
    1     2     3     4     5     6     7     8     9 
 9.18  7.92 10.70  8.80  9.17  9.32  9.54  8.91  9.31 
                  Df Sum Sq Mean Sq F value Pr(>F)    
vtrans_district    1    280  279.66   84.34 <2e-16 ***
Residuals       1455   4824    3.32                   
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

	Pairwise comparisons using t tests with pooled SD 

data:  k$ln_cost and k$vtrans_district 

  1       2       3       4       5       6       7       8      
2 0.00092 -       -       -       -       -       -       -      
3 9.1e-05 < 2e-16 -       -       -       -       -       -      
4 1.00000 1.8e-06 < 2e-16 -       -       -       -       -      
5 1.00000 0.00070 4.0e-05 1.00000 -       -       -       -      
6 1.00000 < 2e-16 < 2e-16 0.04608 1.00000 -       -       -      
7 1.00000 9.5e-08 0.00302 0.44661 1.00000 1.00000 -       -      
8 1.00000 0.13842 2.7e-05 1.00000 1.00000 1.00000 1.00000 -      
9 1.00000 0.07492 0.10607 1.00000 1.00000 1.00000 1.00000 1.00000

P value adjustment method: bonferroni 
```

## ANOVA 2024

```
ANOVA 2024 
      A      DT   DT/GT      GF   GF/GL      GL      GT   GT/DT GT/GT_R       M 
   9.92    9.48    9.10    9.61   10.44    9.12    9.56    9.75    9.25    9.60 
              Df Sum Sq Mean Sq F value Pr(>F)
PARENT         9   35.2   3.911   1.491  0.146
Residuals   1008 2644.1   2.623               

	Pairwise comparisons using t tests with pooled SD 

data:  k$ln_cost and k$PARENT 

        A     DT    DT/GT GF    GF/GL GL    GT    GT/DT GT/GT_R
DT      1.000 -     -     -     -     -     -     -     -      
DT/GT   1.000 1.000 -     -     -     -     -     -     -      
GF      1.000 1.000 1.000 -     -     -     -     -     -      
GF/GL   1.000 1.000 1.000 1.000 -     -     -     -     -      
GL      0.062 1.000 1.000 0.817 1.000 -     -     -     -      
GT      1.000 1.000 1.000 1.000 1.000 1.000 -     -     -      
GT/DT   1.000 1.000 1.000 1.000 1.000 1.000 1.000 -     -      
GT/GT_R 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 -      
M       1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000 1.000  

P value adjustment method: bonferroni 
        A       A/D         B       B/D         C       C/D         D not rated 
     9.54      9.35      9.69     10.04      9.50      9.54      9.32      9.80 
              Df Sum Sq Mean Sq F value Pr(>F)
HYDROGROUP     7   31.2   4.461   1.701  0.105
Residuals   1010 2648.0   2.622               

	Pairwise comparisons using t tests with pooled SD 

data:  k$ln_cost and k$HYDROGROUP 

          A    A/D  B    B/D  C    C/D  D   
A/D       1.00 -    -    -    -    -    -   
B         1.00 1.00 -    -    -    -    -   
B/D       1.00 1.00 1.00 -    -    -    -   
C         1.00 1.00 1.00 1.00 -    -    -   
C/D       1.00 1.00 1.00 1.00 1.00 -    -   
D         1.00 1.00 0.49 0.12 1.00 1.00 -   
not rated 1.00 1.00 1.00 1.00 1.00 1.00 1.00

P value adjustment method: bonferroni 
   5    6    7    9 
8.92 9.63 9.90 9.42 
                  Df Sum Sq Mean Sq F value   Pr(>F)    
vtrans_district    1   55.3   55.35   21.43 4.14e-06 ***
Residuals       1016 2623.9    2.58                     
---
Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

	Pairwise comparisons using t tests with pooled SD 

data:  k$ln_cost and k$vtrans_district 

  5       6    7   
6 7.3e-08 -    -   
7 7.1e-10 0.22 -   
9 0.40    1.00 0.49

P value adjustment method: bonferroni 
```

## Plots: mean ln(cost) +/- 2 SE by category

### Script: R/116_category_errorbars.R

```r
# 116_category_errorbars.R -- mean ln(cost) +/- 2 standard errors for each category of parent material,
#   hydrologic group and VTrans district, both storms, on the damaged segments the ANOVA uses.
# RUN from mrgp-roads-core: "C:/Program Files/R/R-4.5.3/bin/Rscript.exe" twopart_simple_2026-09-17/R/116_category_errorbars.R
# Graded segments, Precip >= 3.51 in, Averill excluded; cost_dedup, uncapped. Writes figures/category_errorbars_2026-10-06.png

png("twopart_simple_2026-09-17/figures/category_errorbars_2026-10-06.png", width = 2400, height = 2700, res = 200)
par(mfrow = c(3, 2), mar = c(6, 5, 3, 1))

for (v in c("PARENT", "HYDROGROUP", "vtrans_district")) {
  for (STORM in c(2023, 2024)) {
    d <- read.csv(paste0("twopart_simple_2026-09-17/data/analysis_", STORM, "_2026-09-27.csv"))
    k <- d[d$Precip >= 3.51 & d$Town != "Averill" & d$cost_dedup > 0, ]
    k$ln_cost <- log(k$cost_dedup)
    g <- factor(k[[v]])
    n_g    <- as.numeric(table(g))
    mean_g <- as.numeric(tapply(k$ln_cost, g, mean))
    se_g   <- as.numeric(tapply(k$ln_cost, g, sd)) / sqrt(n_g)
    plot(1:nlevels(g), mean_g, ylim = c(6.5, 11.5), xlim = c(0.5, nlevels(g) + 0.5), pch = 16, xaxt = "n",
         xlab = "", ylab = "mean ln(cost)", main = paste(v, STORM))
    axis(1, at = 1:nlevels(g), labels = paste0(levels(g), "\nn=", n_g), cex.axis = 0.8, padj = 0.5)
    arrows(1:nlevels(g), mean_g - 2 * se_g, 1:nlevels(g), mean_g + 2 * se_g, angle = 90, code = 3, length = 0.05)
  }
}

dev.off()
```

![mean ln(cost) with 2 SE bars for each category, both storms](category_errorbars_2026-10-06.png)
