# |-----------------------------------------------------------------------------------|
# | Project:  Kidney RO1                                                              |
# | Script:   Sample size and power calculations                                      |
# | Author:   Davit Sargsyan                                                          |
# | Created:  05/27/2018                                                              |
# |-----------------------------------------------------------------------------------|
# sink(file = "tmp/log_kidney_ro1_power_v3.R")
date()

require(data.table)
require(ggplot2)

# Effect size estimates based on Research Strategy Kidney RO1 Ref #188:
# Zhou G, Johansson U, Peng XR, Bamberg K, Huang Y. An additive effect of 
# eplerenone to ACE inhibitor on slowing the progression of diabetic nephropathy
# in the db/db mice. Am J Transl Res. 2016;8(3):1339-54. 
# PubMed PMID: 27186263; PMCID: PMC4859623.
# Table 1: db/db-22 wk and db/db+Epl-22 wk, UAE (ug/24h)
# Power Curves----
dt1 <- data.table(time = rep(c(0, 22), 2),
                  trt = rep(c("dbdb", "dbdb+Epi"), each = 2),
                  mean = c(374.03, 475.3, 341.13, 202.2),
                  sd = c(327.8, 355.5, 176.2, 160.8),
                  n = rep(8:9, each = 2))
dt1
write.csv(dt1,
          file = "tmp/dt1.csv")

pval <- list()
for (i in 1:1000) {
  tmp <- list(x1 = rnorm(n = 8,
                               mean = dt1$mean[2],
                               sd = dt1$sd[2]),
                    x2 = rnorm(n = 9,
                               mean = dt1$mean[4],
                               sd = dt1$sd[4]))
  pval[[i]] <- t.test(x = tmp$x1,
                      y = tmp$x2,
                      alternative = "greater",
                      var.equal = FALSE)$p.value
}
pval <- do.call("c", pval)
hist(pval, 100)
sum(pval < 0.05)/1000

# Different sample sizes----
nn <- 5:20
out <- list()
for (j in nn) {
  pval <- list()
  for (i in 1:1000) {
    tmp <- list(x1 = rnorm(n = j,
                           mean = dt1$mean[2],
                           sd = dt1$sd[2]),
                x2 = rnorm(n = j,
                           mean = dt1$mean[4],
                           sd = dt1$sd[4]))
    pval[[i]] <- t.test(x = tmp$x1,
                        y = tmp$x2,
                        alternative = "greater",
                        var.equal = FALSE)$p.value
  }
  pval <- do.call("c", pval)
  out[[j]] <- sum(pval < 0.05)/1000
}
out <- data.table(n = nn,
                  power = do.call("c", out))
out

# Plot----
p1 <- ggplot(out,
             aes(x = n,
                 y = power)) +
  geom_rect(xmin = 13,
            xmax = Inf,
            ymin = -Inf,
            ymax = Inf,
            fill = "light green",
            alpha = 0.1,
            color = "black",
            linetype = "dashed") +
  geom_smooth(method = "loess",
              size = 1) + 
  geom_point(fill = "red",
             size = 2,
             shape = 21) +
  geom_hline(yintercept = 0.8,
             linetype = "dashed") +
  scale_x_continuous("Number of Animals per Group",
                     breaks = 5:20) + 
  scale_y_continuous("Power",
                     breaks = seq(0, 
                                  1, 
                                  by = 0.1)) + 
  ggtitle("Power Curve For One-Sided t-Test") +
  theme(plot.title = element_text(hjust = 0.5))
p1

tiff(filename = "tmp/power_curve_from_published_kidney_ro1.tiff",
     height = 5,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# Assuming log-normal distribution----
dt1$lmean <- log(dt1$mean)
dt1$lsd <- log(dt1$sd)
dt1$lse <- dt1$lsd/sqrt(dt1$n)
dt1

# sessionInfo()
# sink()