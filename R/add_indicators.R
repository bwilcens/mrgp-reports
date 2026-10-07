# add_indicators.R -- adds the 0/1 indicator columns chosen from the crosstab and ANOVA screens to the two
#   analysis files, so the two-part scripts can use them in place of the parent-material and hydrologic-group
#   factors. Occurrence side: parent_A, parent_DT, parent_GL, parent_GT_DT, hydro_B. Cost side: parent_till.
# RUN from mrgp-roads-core: "C:/Program Files/R/R-4.5.3/bin/Rscript.exe" twopart_simple_2026-09-17/R/add_indicators.R
# The original files are copied to data/_backup_2026-10-06_pre_indicators/ before being overwritten.

dir.create("twopart_simple_2026-09-17/data/_backup_2026-10-06_pre_indicators", showWarnings = FALSE)

for (STORM in c(2023, 2024)) {
  f <- paste0("twopart_simple_2026-09-17/data/analysis_", STORM, "_2026-09-27.csv")
  file.copy(f, "twopart_simple_2026-09-17/data/_backup_2026-10-06_pre_indicators/", overwrite = FALSE)
  d <- read.csv(f)

  d$parent_A     <- as.integer(d$PARENT == "A")                    # alluvium
  d$parent_DT    <- as.integer(d$PARENT == "DT")                   # dense till
  d$parent_GL    <- as.integer(d$PARENT == "GL")                   # glacio-lacustrine
  d$parent_GT_DT <- as.integer(d$PARENT == "GT/DT")                # ablation till over dense till
  d$parent_till  <- as.integer(d$PARENT %in% c("DT", "GT"))        # either till (cost side)
  d$hydro_B      <- as.integer(d$HYDROGROUP == "B")                # hydrologic group B

  write.csv(d, f, row.names = FALSE)
  cat("storm", STORM, "rows", nrow(d), "columns", ncol(d), "\n")
  print(colSums(d[, c("parent_A", "parent_DT", "parent_GL", "parent_GT_DT", "parent_till", "hydro_B")]))
}
