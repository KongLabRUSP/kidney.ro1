# |-----------------------------------------------------------------------------------|
# | Project:  Kidney RO1                                                              |
# | Script:   Sample size and power calculations                                      |
# | Author:   Davit Sargsyan                                                          |
# | Created:  08/23/2017                                                              |
# |-----------------------------------------------------------------------------------|
# sink(file = "tmp/log_kidney_ro1_power_v1.R")
date()

require(data.table)
require(pwr)
require(bit64)
require(ggplot2)

# Power Curves----
# Standard deviations (% mean)
std <- seq(20, 80, by = 20)
std

# Mean differences (%)
delta <- matrix(rep(seq(20, 100, by = 20),
                    length(std)),
                nrow = length(std),
                byrow = TRUE)
delta

# Effect size
h <- delta/std
h
colnames(h) <- paste("Mean Diff=", 
                     delta[1, ],
                     "%",
                     sep = "")
rownames(h) <- paste("SD=", 
                     std,
                     "%",
                     sep = "")
h

# Number of animals
n <- seq(3, 20, by = 1)
n

# Aplha = 0.05
res <- list()
for(i in 1:nrow(h)) {
  out <- list()
  for(j in 1:ncol(h)) {
    out[[j]] <- pwr.2p.test(h = h[i, j],
                            n = n,
                            sig.level = 0.05/4,
                            alternative = "two.sided")$power
  }
  res[[i]] <- data.table(mu = rep(colnames(h), 
                                  each = length(n)),
                         n = rep(n,
                                 ncol(h)), 
                         power = do.call("c", out))
}
dt1 <- data.table(std = rep(rownames(h),
                            each = nrow(res[[1]])),
                  do.call("rbind",
                          res))
dt1$std <- factor(dt1$std,
                  levels = unique(dt1$std))
dt1$mu <- factor(dt1$mu,
                 levels = rev(unique(dt1$mu)))
dt1

tmp <- subset(dt1,
              mu == "Mean Diff=100%" &
              std == "SD=20%")
tmp

p1 <- ggplot(dt1) +
  facet_wrap(~ std,
             nrow = 2) +
  geom_line(aes(x = n,
                y = power,
                group = mu,
                colour = mu),
            size = 1) + 
  geom_hline(yintercept = 0.8,
             linetype = "dashed") +
  scale_x_continuous("Number of Animals per Group",
                     breaks = 3:20) + 
  scale_y_continuous("Power",
                     breaks = seq(0, 
                                  1, 
                                  by = 0.1)) + 
  ggtitle("") +
  guides(colour = guide_legend(title = "Mean Difference (%)")) +
  theme(plot.title = element_text(hjust = 0.5))
p1

tiff(filename = "tmp/power_curves_kidney_ro1.tiff",
     height = 7,
     width = 9,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()


x1 <- 2.39
sd1 <- 0.35
x2 <- 0.054
sd2 <- 0.02

log2(x1*sd2/x2*sd1)

pwr.2p.test(h = log2(x1*sd2/x2*sd1),
            n = 10,
            sig.level = 0.05/5,
            alternative = "two.sided")

# sessionInfo()
# sink()
