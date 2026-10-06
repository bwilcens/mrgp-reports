# category_errorbars.R -- mean ln(cost) +/- 2 standard errors for each category of parent material,
#   hydrologic group and VTrans district, both storms, on the damaged segments the ANOVA uses.
# RUN from mrgp-roads-core: "C:/Program Files/R/R-4.5.3/bin/Rscript.exe" twopart_simple_2026-09-17/R/category_errorbars.R
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
