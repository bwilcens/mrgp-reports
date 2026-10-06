# MRGP road-damage model reports

Meeting reports for the two-part road-damage model: the script that produced each result, then the R output unedited.

- [Road Damage Category Screen](CATEGORY_SCREEN_2026-10-06.md) — crosstabs and one-way ANOVA of parent material, hydrologic group and VTrans district, 2023 and 2024, with error-bar plots.
- [Two-Part Model Runs](TWOPART_RUNS_2026-10-06.md) — the selected logit and Gamma models and the avoided cost, 2023 and 2024.

Data: compliance-graded road segments with at least 3.51 in of rain; damage and cost from the de-duplicated cost field, uncapped. Scripts live in the thesis repository under `twopart_simple_2026-09-17/R/` as `crosstab_anova.R`, `twopart_model.R` and `category_errorbars.R`.
