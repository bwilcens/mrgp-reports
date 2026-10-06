# MRGP road-damage model reports

Meeting reports for the two-part road-damage model: the script that produced each result, then the R output unedited.

- [Road Damage Category Screen](CATEGORY_SCREEN_2026-10-06.md) — crosstabs and one-way ANOVA of parent material, hydrologic group and VTrans district, 2023 and 2024, with error-bar plots.
- [Two-Part Model Runs](TWOPART_RUNS_2026-10-06.md) — the selected logit and Gamma models and the avoided cost, 2023 and 2024; then the same model on the both-source towns, pooled de-duplicated cost against VTrans records alone.
- [Binary Selection from the Crosstab and ANOVA](BINARY_SELECTION_2026-10-06.md) — the categories chosen as 0/1 indicators by the crosstab (Part 1) and ANOVA (Part 2) rules, the refitted two-part model, and the comparison with the full-factor model, 2023 and 2024.
- [Crosstabs and ANOVA, Classwork Style](CLASSWORK_CROSSTAB_ANOVA_2026-10-06.md) — the same crosstabs and ANOVAs as a step-by-step exercise script, with the indicator decision printed next to each category, 2023 and 2024.

Data: compliance-graded road segments with at least 3.51 in of rain; damage and cost from the de-duplicated cost field, uncapped. Scripts live in the thesis repository under `twopart_simple_2026-09-17/R/` as `crosstab.R`, `anova_spss_tables.R`, `twopart_model.R`, `both_source_towns.R` and `category_errorbars.R`.

## Sources for the code

- Crosstab workflow (`table`, `margin.table`, `prop.table`, `chisq.test`): Kabacoff, *Quick-R: Frequencies and Crosstabs*, https://www.datacamp.com/doc/r/frequencies
- Expected counts and adjusted (standardized) residuals: R documentation for `chisq.test`, components `expected` and `stdres`, following Agresti (2007) *An Introduction to Categorical Data Analysis*, section 2.4.5, https://stat.ethz.ch/R-manual/R-devel/library/stats/html/chisq.test.html
- One-way ANOVA (`aov`, `summary`) and pairwise comparisons: Datanovia, *One-Way ANOVA in R*, https://www.datanovia.com/learn/biostatistics/anova/anova-in-r ; R documentation for `pairwise.t.test` (pooled SD, `p.adjust.method = "bonferroni"`), https://stat.ethz.ch/R-manual/R-devel/library/stats/html/pairwise.t.test.html
- SPSS ONEWAY post hoc layout (mean difference, standard error, Bonferroni significance and confidence interval): IBM SPSS Statistics Algorithms, chapter *ONEWAY Algorithms*
- Town-clustered standard errors (`sandwich::vcovCL`): Zeileis, Köll and Graham (2020), *Various Versatile Variances: An Object-Oriented Implementation of Clustered Covariances in R*, Journal of Statistical Software 95(1), https://doi.org/10.18637/jss.v095.i01
- Two-part model and its statistics: Deb and Norton (2018), *Modeling Health Care Expenditures and Use*, Annual Review of Public Health 39, 489–505; Belotti, Deb, Manning and Norton (2015), *twopm: Two-part models*, Stata Journal 15(1), 3–20
