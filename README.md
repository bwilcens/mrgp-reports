# MRGP road-damage model reports

- [Crosstab and ANOVA on Categoricals](CROSSTAB_ANOVA_2026-10-06.md) — crosstabs and one-way ANOVAs on parent material, hydrologic group and VTrans district, with the indicator decision beside each category and the mean-and-whisker figures, 2023 and 2024.
- [Two-Part Model with Indicators, by Damage Record Source](TWOPART_SOURCES_2026-10-06.md) — the two-part model with the chosen 0/1 indicators, statewide and on the both-source towns under pooled and VTrans-only records, side by side, 2023 and 2024.

Scripts are in `R/`: `crosstab_anova.R`, `twopart_sources.R`, and `add_indicators.R` (adds the indicator columns to the data files).

Data: compliance-graded road segments with at least 3.51 in of rain; damage and cost from the de-duplicated cost field, uncapped.

## Sources for the code

- Crosstab workflow (`table`, `margin.table`, `prop.table`, `chisq.test`): Kabacoff, *Quick-R: Frequencies and Crosstabs*, https://www.datacamp.com/doc/r/frequencies
- Expected counts and adjusted (standardized) residuals: R documentation for `chisq.test`, components `expected` and `stdres`, following Agresti (2007) *An Introduction to Categorical Data Analysis*, section 2.4.5, https://stat.ethz.ch/R-manual/R-devel/library/stats/html/chisq.test.html
- One-way ANOVA (`aov`, `summary`) and pairwise comparisons: Datanovia, *One-Way ANOVA in R*, https://www.datanovia.com/learn/biostatistics/anova/anova-in-r ; R documentation for `pairwise.t.test` (pooled SD, `p.adjust.method = "bonferroni"`), https://stat.ethz.ch/R-manual/R-devel/library/stats/html/pairwise.t.test.html
- Town-clustered standard errors (`sandwich::vcovCL`): Zeileis, Köll and Graham (2020), *Various Versatile Variances: An Object-Oriented Implementation of Clustered Covariances in R*, Journal of Statistical Software 95(1), https://doi.org/10.18637/jss.v095.i01
- Two-part model and its statistics: Deb and Norton (2018), *Modeling Health Care Expenditures and Use*, Annual Review of Public Health 39, 489–505; Belotti, Deb, Manning and Norton (2015), *twopm: Two-part models*, Stata Journal 15(1), 3–20
