# |----------------------------------------------------------------------------------|
# | Project: Study of Diabetes in MES13 cells                                        |
# | Script: Methyl-seq data analysis and visualization                               |
# | Author: Davit Sargsyan                                                           |
# | Created: 09/11/2017                                                              |
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

length(unique(dt1$SYMBOL))
unique(substr(dt1$annotation, 1, 4))
dt2 <- subset(dt1,
              substr(annotation, 1, 4) == "Exon")
length(unique(dt2$SYMBOL))

# Remove unmapped regions
dt1 <- dt1[!is.na(dt1$SYMBOL == "NA"), ]

# Subset data----
# a. Exons Only----
dt2 <- subset(dt1,
              substr(annotation, 1, 4) == "Exon")
length(unique(dt2$SYMBOL))

# # b. Exons and Introns Only----
# dt2 <- subset(dt1,
#               substr(annotation, 1, 4) %in%"Exon")
# length(unique(dt2$SYMBOL))

#...

# Calculate percent methylation in each sample----
dt2 <- data.table(gene = dt2$SYMBOL,
                  toTSS = dt2$distanceToTSS,
                  CpG = dt2$CpG,
                  pct = round(100*dt2[, c(3, 5, 7, 9, 11)]/
                                dt2[, c(2, 4, 6, 8, 10)],
                              1))
summary(dt2)

# Sort by distance to TSS and gene
dt2 <- dt2[order(dt2$toTSS)]
dt2 <- dt2[order(dt2$gene)]
dt2

# Heatmap
dt.l <- melt.data.table(data = dt2,
                        id.vars = 1:3,
                        measure.vars = 4:8,
                        variable.name = "Treatment",
                        value.name = "Percent Methylation")
dt.l$gene.id <- paste(dt.l$gene,
                      "(dist =",
                      dt.l$toTSS,
                      ", CpG = ",
                      dt.l$CpG,
                      ")",
                      sep = "")
dt.l$Treatment <- factor(dt.l$Treatment)
levels(dt.l$Treatment) <- c("LG",
                            "HG",
                            "MIC 1.5uM",
                            "FX 1uM",
                            "Ber 6uM")

# Plot----
p1 <- ggplot(data = dt.l) +
  geom_tile(aes(x =  Treatment,
                y = gene.id,
                fill = `Percent Methylation`),
            color = "black") +
  # geom_text(aes(x =  Treatment,
  #               y = gene.id,
  #               label = CpG)) +
  scale_fill_gradient2(high = "red", 
                       #mid = "black", 
                       #low = "yellow", 
                       #midpoint = 50, 
                       limit = c(0, 100), 
                       name = "Percent Methylation") +
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete("Gene and Distance from TSS",
                   expand = c(0, 0)) +
  ggtitle("Methylation in MES13 Exons") +
  theme(axis.text.x = element_text(angle = 0,
                                   hjust = 1),
        # legend.position = "top",
        plot.title = element_text(hjust = 0.5))
print(p1)

tiff(filename = "mes13/tmp/pct_methyl_heatmap.tiff",
     height = 25,
     width = 8,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# Save the table----
write.csv(dt2,
          file = "mes13/tmp/dt2.csv")

# Log2 change (2-fold change)
# NOTE: add 1 to all non-NA values so log can be calculated
dt.log2 <- apply(dt2[, 4:8],
                 MARGIN = 2,
                 FUN = function(a){
                   a <- log2(a + 1)
                   a[is.na(a)] <- 0
                   return(a)
                 })
# Comute log2 differences with positive control (High Glucose)
log2.diff <- dt.log2[, -2] - dt.log2[, 2]
colnames(log2.diff) <- paste(levels(dt.l$Treatment)[-2],
                             levels(dt.l$Treatment)[2],
                             sep = " - ")
log2.diff <- data.table(dt2[, 1:3],
                        log2.diff)
log2.diff 

# Plot----
dt.diff.l <- melt.data.table(data = log2.diff,
                             id.vars = 1:3,
                             measure.vars = 4:7,
                             variable.name = "Treatment",
                             value.name = "Log2 Difference")
dt.diff.l$gene.id <- paste(dt.diff.l$gene,
                           "(dist =",
                           dt.diff.l$toTSS,
                           ", CpG = ",
                           dt.diff.l$CpG,
                           ")",
                           sep = "")
dt.diff.l$Treatment <- factor(dt.diff.l$Treatment)

p2 <- ggplot(data = dt.diff.l) +
  geom_tile(aes(x =  Treatment,
                y = gene.id,
                fill = `Log2 Difference`),
            color = "black") +
  scale_fill_gradient2(high = "green", 
                       mid = "black", 
                       low = "red", 
                       midpoint = 0, 
                       # limit = c(0, 100), 
                       name = "Difference of Log2(% Methyl)") +
  scale_x_discrete(expand = c(0, 0)) + 
  scale_y_discrete("Gene, Distance from TSS, and # of CpG",
                   expand = c(0, 0)) +
  ggtitle("Difference in Methylation in MES13 Exons") +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        plot.title = element_text(hjust = 0.5))
print(p2)

tiff(filename = "mes13/tmp/diff_log2_pct_methyl_heatmap.tiff",
     height = 25,
     width = 8,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p2)
graphics.off()

# Save the table----
write.csv(dt.diff.l,
          file = "mes13/tmp/dt.diff.l.csv")