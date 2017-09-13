# |----------------------------------------------------------------------------------|
# | Project: Study of Diabetes in MES13 cells                                        |
# | Script: Methyl-seq data analysis and visualization                               |
# | Author: Davit Sargsyan                                                           |
# | Created: 09/12/2017                                                              |
# |----------------------------------------------------------------------------------|
# Header----
# Source: file:///C:/R/R-3.3.2/library/ChIPseeker/doc/ChIPseeker.html
# Source: http://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html

# # NOTE: FOR LINUX, RUN THIS IN THE TERMINAL AS SUDO!
# # Source: https://support.bioconductor.org/p/70093/
# sudo R
# source("http://bioconductor.org/biocLite.R")
# biocLite("TxDb.Mmusculus.UCSC.mm9.knownGene",
#          suppressUpdates = TRUE)
# biocLite("org.Mm.eg.db")

# biocLite("ReactomePA")

# Save consol output to a log file
# sink(file = "tmp/log_mes13_methylseq_data_analysis_v1.txt")

require(data.table)
require(ggplot2)
require(ChIPseeker)
require(TxDb.Mmusculus.UCSC.mm9.knownGene)
require(knitr)

# Load data----
peakAnno1 <- annotatePeak(peak = "mes13/data/combined.csv", 
                          tssRegion = c(-3000, 3000), 
                          TxDb = TxDb.Mmusculus.UCSC.mm9.knownGene,
                          annoDb = "org.Mm.eg.db")
dt1 <- data.table(as.data.frame(peakAnno1@anno@elementMetadata@listData))

# dt1[, p.fdr := p.adjust(p = Control..Exptl.pval,
#                         method = "fdr")]
setkey(dt1)
dt1

# Remove unmapped regions
dt1 <- dt1[!is.na(dt1$SYMBOL == "NA"), ]

# Subset data----
dt1 <- data.table(gene = dt1$SYMBOL,
                  anno = dt1$annotation,
                  reg = NA,
                  CpG = dt1$CpG,
                  dt1[, 2:11])

unique(dt1$anno)
kable(data.table(table(substr(dt1$anno, 1, 9))))
  # |V1        |    N|
  # |:---------|----:|
  # |3' UTR    |   38|
  # |5' UTR    |   15|
  # |Distal In |  838|
  # |Downstrea |   22|
  # |Exon (uc0 |  212|
  # |Intron (u |  468|
  # |Promoter  | 2581|

# Separate Promoter, Body and Downstream; remove everything else
# a. Promoter: up to 3kb upstream
dt1$reg[substr(dt1$anno, 
               1,
               8) == "Promoter"] <- "Promoter"

# b. Body: exons and introns
dt1$reg[substr(dt1$anno, 
               1, 
               4) %in% c("Exon",
                         "Intr")] <- "Body"

# c. Downstream: Distal Intergenic and  Downstream
dt1$reg[substr(dt1$anno, 
               1, 
               4) %in% c("Dist",
                         "Down")] <- "Downstream"
dt1$reg <- factor(dt1$reg,
                  levels = c("Promoter",
                             "Body",
                             "Downstream"))

# NOTE: disregarded 5' and 3' (only 53 regions combined)
dt1 <- droplevels(subset(dt1,
                         !is.na(reg)))
dt1[, anno := NULL]
summary(dt1)

# Part I: total methylation----
# Aggregate data by region
out <- list()
for (i in 3:13) {
  out[[i - 2]] <- aggregate(dt1[, i, with = FALSE],
                            by = list(dt1$reg),
                            FUN = sum,
                            na.rm = TRUE)
}
dt.reg <- data.table(Reduce(merge, out))
dt.reg

# Calculate percent methylation in each sample----
dt.reg <- data.table(Region = dt.reg$Group.1,
                     CpG = dt.reg$CpG,
                     pct = round(100*dt.reg[, c(4, 6, 8, 10, 12)]/
                                   dt.reg[, c(3, 5, 7, 9, 11)],
                                 1))
summary(dt.reg)

# Melt data
dt.reg.l <- melt.data.table(data = dt.reg,
                            id.vars = 1:2,
                            measure.vars = 3:7,
                            variable.name = "Treatment",
                            value.name = "Percent Methylation")
dt.reg.l$reg.id <- paste(dt.reg.l$Region,
                         "(",
                         dt.reg.l$CpG,
                         " CpG)",
                         sep = "")
dt.reg.l$reg.id <- factor(dt.reg.l$reg.id,
                          levels = c("Promoter(22856 CpG)",
                                     "Body(5153 CpG)",
                                     "Downstream(6458 CpG)"))

dt.reg.l$Treatment <- factor(dt.reg.l$Treatment)
levels(dt.reg.l$Treatment) <- c("LG",
                                "HG",
                                "MIC 1.5uM",
                                "FX 1uM",
                                "Ber 6uM")

# Plot
p1 <- ggplot(dt.reg.l,
             aes(x = reg.id,
                 y = `Percent Methylation`,
                 group = Treatment,
                 fill = Treatment)) +
  geom_bar(position = position_dodge(),
           stat="identity",
           color = "black") +
  scale_x_discrete("Region") +
  scale_y_continuous(limits = c(0, 100)) +
  ggtitle("Total Methylation (%)") +
  theme(plot.title = element_text(hjust = 0.5))

tiff(filename = "mes13/tmp/total_pct_methyl.tiff",
     height = 5,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# Part II: gene methylation----
# Collapse by gene
out <- list()
for (i in 3:13) {
  out[[i - 2]] <- aggregate(dt1[, i, with = FALSE],
                            by = list(dt1$gene,
                                      dt1$reg),
                            FUN = sum,
                            na.rm = TRUE)
}

dt.gene <- data.table(Reduce(merge, out))
dt.gene

# Calculate percent methylation in each sample----
dt.gene <- data.table(gene = dt.gene$Group.1,
                      reg = dt.gene$Group.2,
                      CpG = dt.gene$CpG,
                      pct = round(100*dt.gene[, c(5, 7, 9, 11, 13)]/
                                    dt.gene[, c(4, 6, 8, 10, 12)],
                                  1))
summary(dt.gene)

# Melt data
dt.gene.l <- melt.data.table(data = dt.gene,
                             id.vars = 1:3,
                             measure.vars = 4:8,
                             variable.name = "Treatment",
                             value.name = "Percent Methylation")
dt.gene.l$gene.id <- paste(dt.gene.l $gene,
                           "(",
                           dt.gene.l $CpG,
                           " CpG)",
                           sep = "")
dt.gene.l$Treatment <- factor(dt.gene.l$Treatment)
levels(dt.gene.l$Treatment) <- c("LG",
                                 "HG",
                                 "MIC 1.5uM",
                                 "FX 1uM",
                                 "Ber 6uM")

# Heatmap----
# Plot top 50 genes
dt.gene.l[, summ := sum(`Percent Methylation`,
                        na.rm = TRUE),
          by = c("gene",
                 "reg")]

dt.gene.l <- dt.gene.l[order(dt.gene.l$summ,
                             decreasing = TRUE)]
dt.gene.l

# Keep top 50 genes from each region
gene.keep <- c(as.character(dt.gene.l$gene[dt.gene.l$reg == "Promoter"][1:50]),
               as.character(dt.gene.l$gene[dt.gene.l$reg == "Body"][1:50]),
               as.character(dt.gene.l$gene[dt.gene.l$reg == "Downstream"][1:50]))

tmp <- subset(dt.gene.l,
              gene %in% gene.keep)

p2 <- ggplot(data = tmp) +
  facet_wrap(~ reg,
             # scales = "free_y",
             nrow = 1) +
  geom_tile(aes(x =  Treatment,
                y = gene.id,
                fill = `Percent Methylation`),
            color = "black") +
  scale_fill_gradient2(high = "red", 
                       limit = c(0, 100), 
                       name = "Percent Methylation") +
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete("Gene and CpG",
                   expand = c(0, 0)) +
  ggtitle("Methylation in MES13 by Region, Top 50 Genes") +
  theme(axis.text.x = element_text(angle = 30,
                                   hjust = 1),
        # legend.position = "top",
        plot.title = element_text(hjust = 0.5))
print(p2)

tiff(filename = "mes13/tmp/top50_pct_methyl_heatmap.tiff",
     height = 8,
     width = 8,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p2)
graphics.off()

# Save the table----
write.csv(dt.gene,
          file = "mes13/tmp/dt.gene.csv")

# Log2 change (2-fold change)
# NOTE: add 1 to all non-NA values so log can be calculated
dt.log2 <- apply(dt.gene[, 4:8],
                 MARGIN = 2,
                 FUN = function(a){
                   a <- log2(a + 1)
                   a[is.na(a)] <- 0
                   return(a)
                 })

# Comute log2 differences with positive control (High Glucose)
log2.diff <- dt.log2[, -3] - dt.log2[, 3]
colnames(log2.diff) <- paste(levels(dt.gene.l$Treatment)[-2],
                             levels(dt.gene.l$Treatment)[2],
                             sep = " - ")
log2.diff <- data.table(dt.gene[, 1:2],
                        log2.diff)
log2.diff 

# Heatmap----
dt.diff.l <- melt.data.table(data = log2.diff,
                             id.vars = 1:2,
                             measure.vars = 3:6,
                             variable.name = "Treatment",
                             value.name = "Log2 Difference")
dt.diff.l$Treatment <- factor(dt.diff.l$Treatment)

# Plot top 50 genes
dt.diff.l[, summ := sum(abs(`Log2 Difference`),
                        na.rm = TRUE),
          by = c("gene",
                 "reg")]

dt.diff.l <- dt.diff.l[order(dt.diff.l$summ,
                             decreasing = TRUE)]
dt.diff.l

gene.keep <- c(as.character(dt.diff.l$gene[dt.diff.l$reg == "Promoter"][1:50]),
               as.character(dt.diff.l$gene[dt.diff.l$reg == "Body"][1:50]),
               as.character(dt.diff.l$gene[dt.diff.l$reg == "Downstream"][1:50]))

tmp <- subset(dt.diff.l,
              gene %in% gene.keep)

p3 <- ggplot(data = tmp) +
  facet_wrap(~ reg,
             # scales = "free_y",
             nrow = 1) +
  geom_tile(aes(x =  Treatment,
                y = gene,
                fill = `Log2 Difference`),
            color = "black") +
  scale_fill_gradient2(high = "green", 
                       mid = "black", 
                       low = "red", 
                       midpoint = 0, 
                       name = "Difference of Log2(% Methyl)") +
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete("Gene",
                   expand = c(0, 0)) +
  ggtitle("Difference in Methylation in MES13, Top 50 Genes") +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        plot.title = element_text(hjust = 0.5))
print(p3)

tiff(filename = "mes13/tmp/top50_diff_log2_pct_methyl_heatmap.tiff",
     height = 8,
     width = 8,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p3)
graphics.off()

# Save the table----
write.csv(dt.diff.l,
          file = "mes13/tmp/dt.diff.l.csv")
