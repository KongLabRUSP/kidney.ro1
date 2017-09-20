# Project: GM-CSF
# Description: sample size calculation
# Author: Davit Sargsyan
# Date: 08/23/2017
#****************************************************
require(data.table)
require(pwr)
require(bit64)
require(ggplot2)

# Power Curves----
# Standard deviations (% mean)
std <- seq(10, 80, by = 10)
std

# Mean differences (%): 2-fold change
delta <- 100

# Effect size
h <- delta/std
names(h) <- paste("SD=", 
                     std,
                     "%",
                     sep = "")
h

# Number of animals
n <- seq(3, 20, by = 1)
n

# Aplha = 0.05
# 5 pairwise comparisons: ctrl- vs. ctrl+, each of 4 arms vs. ctrl+
res <- list()
for(i in 1:length(h)) {
    res[[i]] <- pwr.2p.test(h = h[i],
                            n = n,
                            sig.level = 0.05/5,
                            alternative = "g")$power
}

dt1 <- data.table(std = rep(names(h),
                            each = length(res[[1]])),
                  n = rep(n, 
                          length(res)),
                  power = do.call("c",
                          res))
dt1$std <- factor(dt1$std,
                  levels = unique(dt1$std))
dt1

p1 <- ggplot(dt1) +
  geom_rect(xmin = 4,
            xmax = 8,
            ymin = -Inf,
            ymax = Inf,
            fill = "light green",
            alpha = 0.1,
            color = "black",
            linetype = "dashed") +
  geom_line(aes(x = n,
                y = power,
                group = std,
                colour = std),
            size = 1) + 
  geom_hline(yintercept = 0.8,
             linetype = "dashed") +
  scale_x_continuous("Number of Animals per Group",
                     breaks = 3:20) + 
  scale_y_continuous("Power",
                     breaks = seq(0, 
                                  1, 
                                  by = 0.1)) + 
  ggtitle("Power Curves For One-Sided Test\nAt Least Two-Fold Difference Between Group Means") +
  guides(colour = guide_legend(title = "Standard Deviation\n(% Treatment)")) +
  theme(plot.title = element_text(hjust = 0.5))
p1

tiff(filename = "tmp/power_curves_kidney_ro1.tiff",
     height = 5,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()