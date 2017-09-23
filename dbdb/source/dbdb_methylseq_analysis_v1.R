# |----------------------------------------------------------------------------------|
# | Project: Study of DiabetesDb/Db                                       |
# | Script: Methyl-seq data analysis and visualization                               |
# | Author: Davit Sargsyan                                                           |
# | Created: 09/23/2017                                                              |
# |----------------------------------------------------------------------------------|
# Header----
# Save consol output to a log file
# sink(file = "tmp/log_mes13_methylseq_data_analysis_v1.txt")

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

require(data.table)
require(ggplot2)
require(ChIPseeker)
require(TxDb.Mmusculus.UCSC.mm9.knownGene)
require(knitr)
require(ReactomePA)

# Treatment legend----
trt.names <- c("db/m,16w",
               "db/m,16w",
               "db/db,16w",
               "db/db,16w",
               "db/m,21w",
               "db/m,21w",
               "db/m,21w",
               "db/db,21w",
               "db/db,21w",
               "db/db,21w")

# Load data----
peakAnno1 <- annotatePeak(peak = "dbdb/data/methyl_seq/combined.dedup.csv", 
                          tssRegion = c(-3000, 3000), 
                          TxDb = TxDb.Mmusculus.UCSC.mm9.knownGene,
                          annoDb = "org.Mm.eg.db")
dt1 <- data.table(as.data.frame(peakAnno1@anno@elementMetadata@listData))

# dt1[, p.fdr := p.adjust(p = Control..Exptl.pval,
#                         method = "fdr")]
dt1

# Remove unmapped regions
dt1 <- dt1[!is.na(dt1$SYMBOL == "NA"), ]

# Subset data----
dt1 <- data.table(gene = dt1$SYMBOL,
                  anno = dt1$annotation,
                  geneId = dt1$geneId,
                  reg = NA,
                  CpG = dt1$CpG,
                  dt1[, 2:21])

unique(dt1$anno)
kable(data.table(table(substr(dt1$anno, 1, 9))))
  # |V1        |    N|
  # |:---------|----:|
  # |3' UTR    |  126|
  # |5' UTR    |   53|
  # |Distal In | 2318|
  # |Downstrea |   83|
  # |Exon (uc0 |  764|
  # |Intron (u | 1400|
  # |Promoter  | 8476|

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
for (i in 4:24) {
  out[[i - 3]] <- aggregate(dt1[, i, with = FALSE],
                            by = list(dt1$reg),
                            FUN = sum,
                            na.rm = TRUE)
}
dt.reg <- data.table(Reduce(merge, out))

# Change row order
dt.reg <- dt.reg[c(3, 1, 2), ]
dt.reg

# Hits per CpG average (i.e. vertical coverage)
t1 <- dt.reg[, c(1, 2, seq(3, 21, 2))]
t1[, 3:12] <- do.call("cbind",
                     lapply(t1[, 3:12],
                   function(a) {
                     return(round(a/t1[, 2],
                                  1))
                   }))
colnames(t1) <- c("Gene Region",
                  "Total CpG Count",
                  trt.names)
kable(t1)
  # |Gene Region | Total CpG Count| db/m,16w| db/m,16w| db/db,16w| db/db,16w| db/m,21w| db/m,21w| db/m,21w| db/db,21w| db/db,21w| db/db,21w|
  # |:-----------|---------------:|--------:|--------:|---------:|---------:|--------:|--------:|--------:|---------:|---------:|---------:|
  # |Promoter    |           89043|      0.7|      0.7|       1.0|       0.8|      0.6|      0.9|      0.9|       0.7|       0.8|       0.7|
  # |Body        |           18649|      0.9|      0.8|       1.3|       1.0|      0.8|      1.1|      1.1|       0.8|       0.9|       0.9|
  # |Downstream  |           20414|      1.5|      1.4|       1.7|       1.6|      1.3|      1.7|      1.6|       1.3|       1.5|       1.5|
write.csv(t1,
          file = "dbdb/tmp/t1.csv",
          row.names = FALSE)

# Calculate percent methylation in each sample----
# collapse samples
dt.reg <- data.table(Region = dt.reg$Group.1,
                     CpG = dt.reg$CpG,
                     `db_m_16w_N` = dt.reg$X2.O_S1_R1_001.N + 
                       dt.reg$X2.R_S2_R1_001.N,
                     `db_m_16w_X` = dt.reg$X2.O_S1_R1_001.X + 
                       dt.reg$X2.R_S2_R1_001.X,
                     `db_db_16w_N` = dt.reg$X4.L_S3_R1_001.N +
                       dt.reg$X4.LL_S4_R1_001.N,
                     `db_db_16w_X` = dt.reg$X4.L_S3_R1_001.X +
                       dt.reg$X4.LL_S4_R1_001.X,
                     `db_m_21w_N` = dt.reg$X6.L_S5_R1_001.N +
                       dt.reg$X6.LL_S6_R1_001.N +
                       dt.reg$X7.L_S1_R1_001.N,
                     `db_m_21w_X` = dt.reg$X6.L_S5_R1_001.X +
                       dt.reg$X6.LL_S6_R1_001.X +
                       dt.reg$X7.L_S1_R1_001.X,
                     `db_db_21w_N` = dt.reg$X5.O_S2_R1_001.N +
                       dt.reg$X5.R_S7_R1_001.N +
                       dt.reg$X5.RR_S8_R1_001.N,
                     `db_db_21w_X` = dt.reg$X5.O_S2_R1_001.X +
                       dt.reg$X5.R_S7_R1_001.X +
                       dt.reg$X5.RR_S8_R1_001.X)

dt.reg <- data.table(Region = dt.reg$Region,
                     CpG = dt.reg$CpG,
                     pct = round(100*dt.reg[, c(4, 6, 8, 10)]/
                                   dt.reg[, c(3, 5, 7, 9)],
                                 1))
dt.reg

# Melt data
dt.reg.l <- melt.data.table(data = dt.reg,
                            id.vars = 1:2,
                            measure.vars = 3:6,
                            variable.name = "Treatment",
                            value.name = "Methylation(%)")

dt.reg.l$Treatment <- factor(dt.reg.l$Treatment)
levels(dt.reg.l$Treatment) <- unique(trt.names)
dt.reg.l

# Plot
p1 <- ggplot(dt.reg.l,
             aes(x = Region,
                 y = `Methylation(%)`,
                 group = Treatment,
                 fill = Treatment)) +
  geom_bar(position = position_dodge(),
           stat="identity",
           color = "black") +
  scale_x_discrete("Region") +
  scale_y_continuous(limits = c(0, 100)) +
  ggtitle("Total Methylation (%)") +
  theme(plot.title = element_text(hjust = 0.5))
p1

tiff(filename = "dbdb/tmp/total_pct_methyl.tiff",
     height = 5,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# Part II: gene methylation----
# NOTE: collapse samples
dt2 <- data.table(dt1[, 1:4],
                  `db_m_16w_N` = dt1$X2.O_S1_R1_001.N + 
                    dt1$X2.R_S2_R1_001.N,
                  `db_m_16w_X` = dt1$X2.O_S1_R1_001.X + 
                    dt1$X2.R_S2_R1_001.X,
                  `db_db_16w_N` = dt1$X4.L_S3_R1_001.N +
                    dt1$X4.LL_S4_R1_001.N,
                  `db_db_16w_X` = dt1$X4.L_S3_R1_001.X +
                    dt1$X4.LL_S4_R1_001.X,
                  `db_m_21w_N` = dt1$X6.L_S5_R1_001.N +
                    dt1$X6.LL_S6_R1_001.N +
                    dt1$X7.L_S1_R1_001.N,
                  `db_m_21w_X` = dt1$X6.L_S5_R1_001.X +
                    dt1$X6.LL_S6_R1_001.X +
                    dt1$X7.L_S1_R1_001.X,
                  `db_db_21w_N` = dt1$X5.O_S2_R1_001.N +
                    dt1$X5.R_S7_R1_001.N +
                    dt1$X5.RR_S8_R1_001.N,
                  `db_db_21w_X` = dt1$X5.O_S2_R1_001.X +
                    dt1$X5.R_S7_R1_001.X +
                    dt1$X5.RR_S8_R1_001.X)

# Collapse by gene 
out <- list()
for (i in 4:12) {
  out[[i - 3]] <- aggregate(dt2[, i, with = FALSE],
                            by = list(dt2$gene,
                                      dt2$reg),
                            FUN = sum,
                            na.rm = TRUE)
}

dt.gene <- data.table(Reduce(merge, out))
dt.gene
summary(dt.gene)

# Calculate percent methylation in each sample----
dt.gene <- data.table(gene = dt.gene$Group.1,
                      reg = dt.gene$Group.2,
                      CpG = dt.gene$CpG,
                      pct = round(100*dt.gene[, c(5, 7, 9, 11)]/
                                    dt.gene[, c(4, 6, 8, 10)],
                                  1))

# Melt data
dt.gene.l <- melt.data.table(data = dt.gene,
                             id.vars = 1:3,
                             measure.vars = 4:7,
                             variable.name = "Treatment",
                             value.name = "Methylation(%)")
dt.gene.l$gene.id <- paste(dt.gene.l $gene,
                           "(",
                           dt.gene.l $CpG,
                           " CpG)",
                           sep = "")
dt.gene.l$Treatment <- factor(dt.gene.l$Treatment)
levels(dt.gene.l$Treatment) <- c("db/m, 16 weeks",
                                 "db/db, 16 weeks",
                                 "db/m, 21 weeks",
                                 "db/db, 21 weeks")
summary(dt.gene.l)
dt.gene.l

# Genes with largest differences at 16 weeks (db/db - db/m methylation)----
dt.gene$dbdb_dbm_16w <- dt.gene$pct.db_db_16w_X - dt.gene$pct.db_m_16w_X
m.diff <- dt.gene[!is.nan(dbdb_dbm_16w), c(1, 2, 8)]
m.diff <- m.diff[order(dbdb_dbm_16w), ]
m.diff

# Genes with largest change in methylation (db/db - db/m)----
gene.sorted <- unique(as.character(m.diff$gene))

# a. Highest Positive Difference at 16 weeks (db/db - db/m)----
gene.up <- gene.sorted[(length(gene.sorted) - 49):length(gene.sorted)]
tmp <- subset(dt.gene.l,
              as.character(gene) %in% gene.up)
tmp

p2a <- ggplot(data = tmp) +
  facet_wrap(~ reg,
             #scales = "free_y",
             nrow = 1) +
  geom_tile(aes(x =  Treatment,
                y = gene,
                fill = `Methylation(%)`),
            color = "black") +
  scale_fill_gradient2(high = "red",
                       limit = c(0, 100),
                       name = "Methylation(%)") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete("Gene",
                   expand = c(0, 0)) +
  ggtitle("Top 50 Genes With Increased Methylation \n Db/Db - Db/M, at 16 Weeks") +
  theme(axis.text.x = element_text(angle = 30,
                                   hjust = 1),
        # legend.position = "top",
        plot.title = element_text(hjust = 0.5))
p2a

tiff(filename = "dbdb/tmp/top50_dbdb_dbm_16w_pos_meth_diff.tiff",
     height = 8,
     width = 7,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p2a)
graphics.off()

# b. Genes with smallest change in methylation (db/db - db/m)----
gene.dn <- gene.sorted[1:50]
tmp <- subset(dt.gene.l,
              as.character(gene) %in% gene.dn)
tmp

p2b <- ggplot(data = tmp) +
  facet_wrap(~ reg,
             #scales = "free_y",
             nrow = 1) +
  geom_tile(aes(x =  Treatment,
                y = gene,
                fill = `Methylation(%)`),
            color = "black") +
  scale_fill_gradient2(high = "red",
                       limit = c(0, 100),
                       name = "Methylation(%)") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete("Gene",
                   expand = c(0, 0)) +
  ggtitle("Top 50 Genes With Decreased Methylation \n Db/Db - Db/M, at 16 Weeks") +
  theme(axis.text.x = element_text(angle = 30,
                                   hjust = 1),
        # legend.position = "top",
        plot.title = element_text(hjust = 0.5))
p2b

tiff(filename = "dbdb/tmp/top50_dbdb_dbm_16w_neg_meth_diff.tiff",
     height = 8,
     width = 7,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p2b)
graphics.off()

# Save the tables----
write.csv(dt.gene,
          file = "dbdb/tmp/dt.gene.csv")
 
# # Log2 change (2-fold change)----
# # NOTE: add 1 to all non-NA values so log can be calculated
# dt.log2 <- apply(dt.gene[, 4:8],
#                             MARGIN = 2,
#                             FUN = function(a){
#                               a <- log2(a + 1)
#                               return(a)
#                             })
# head(dt.log2)
# 
# # Compute log2 differences with positive control (High Glucose)
# log2.diff <- dt.log2[, -2] - dt.log2[, 2]
# colnames(log2.diff) <- paste(levels(dt.gene.l$Treatment)[-2],
#                              levels(dt.gene.l$Treatment)[2],
#                              sep = " - ")
# log2.diff <- data.table(dt.gene[, 1:2],
#                         log2.diff)
# log2.diff 
# 
# # Remove gene if all treatment differences are NaN
# ndx.rm <- apply(X = log2.diff[, 3:6],
#                 MARGIN = 1,
#                 FUN = function(a) {
#                   return(all(is.nan(a)))
#                 })
# sum(ndx.rm)
# dt2 <- log2.diff[!ndx.rm, ]
# dt2
# 
# # Save the table----
# write.csv(dt2,
#           file = "mes13/tmp/dt2.csv")
# 
# # Long data----
# dt2.l <- melt.data.table(data = dt2,
#                          id.vars = 1:2,
#                          measure.vars = 3:6,
#                          variable.name = "Treatment",
#                          value.name = "Log2 Difference")
# dt2.l$Treatment <- factor(dt2.l$Treatment)
# 
# # a. Highest Positive  or Negative Difference (HG - LG)----
# tmp <- subset(dt2.l,
#               as.character(gene) %in% unique(c(gene.up,
#                                                gene.dn)))
# tmp
# 
# p3 <- ggplot(data = tmp) +
#   coord_flip() +
#   facet_wrap(~ reg,
#              # scales = "free_y",
#              ncol = 1) +
#   geom_tile(aes(x =  Treatment,
#                 y = gene,
#                 fill = `Log2 Difference`),
#             color = "black") +
#   scale_fill_gradient2(high = "green",
#                        mid = "black",
#                        low = "red",
#                        midpoint = 0,
#                        name = "Difference of Log2(% Methyl)") +
#   scale_x_discrete(expand = c(0, 0)) +
#   scale_y_discrete("Gene",
#                    expand = c(0, 0)) +
#   ggtitle("MES13: Differences of Log2(% Methylation), Top 50 Largest Posotive and Top50 Largest Negative HG-LG Methylation Genes") +
#   theme(axis.text.x = element_text(angle = 45,
#                                    hjust = 1),
#         legend.position = "top",
#         plot.title = element_text(hjust = 0.5))
# print(p3)
# 
# tiff(filename = "mes13/tmp/top50_lg-hg_hg-lg_log2_pct_methyl.tiff",
#      height = 6,
#      width = 12,
#      units = 'in',
#      res = 300,
#      compression = "lzw+p")
# print(p3)
# graphics.off()
# 
# # b. Genes with defined (HG - LG) and log2 change of at least +/-1----
# tmp <- subset(dt2.l,
#               as.character(gene) %in% gene.sorted)
# tmp
# 
# tmp[, Keep := (max(abs(`Log2 Difference`),
#                    na.rm = TRUE) >= 2),
#     by = "gene"]
# sum(tmp$Keep)
# gene.keep <- unique(as.character(tmp$gene[tmp$Keep]))
# tmp <- subset(tmp,
#               gene %in% gene.keep)
# length(unique(tmp$gene))
# 
# # Plot
# p4 <- ggplot(data = tmp) +
#   coord_flip() +
#   facet_wrap(~ reg,
#              # scales = "free_y",
#              ncol = 1) +
#   geom_tile(aes(x =  Treatment,
#                 y = gene,
#                 fill = `Log2 Difference`),
#             color = "black") +
#   scale_fill_gradient2(high = "green",
#                        mid = "black",
#                        low = "red",
#                        midpoint = 0,
#                        name = "Difference of\nLog2(% Methyl)") +
#   scale_x_discrete(expand = c(0, 0)) +
#   scale_y_discrete("Gene",
#                    expand = c(0, 0)) +
#   ggtitle("MES13: Differences of Log2(% Methylation)\nGenes With At Least One Log2 Diff >= +/-2") +
#   theme(axis.text.x = element_text(angle = 45,
#                                    hjust = 1),
#         legend.position = "top",
#         plot.title = element_text(hjust = 0.5))
# print(p4)
# 
# tiff(filename = "mes13/tmp/log2diff_atleast2.tiff",
#      height = 6,
#      width = 9,
#      units = 'in',
#      res = 300,
#      compression = "lzw+p")
# print(p4)
# graphics.off()
# 
# # Part III: Pathway analysis with Reactome----
# # Get Entrezgene IDs
# geneID <- unique(dt1$geneId[as.character(dt1$gene) %in% gene.sorted])
# 
# # Save Entrezgene IDs and gene names
# out <- data.table(geneID = geneID,
#                   geneName = gene.sorted)
# write.csv(out,
#           file = "mes13/tmp/gene.sorted.csv",
#           row.names = FALSE)
# 
# # Get pathways----
# m1 <- enrichPathway(gene = geneID,
#                     pvalueCutoff = 1, 
#                     readable = T, 
#                     organism = "mouse")
# t2 <- as.data.table(m1)
# t2
# write.csv(t2,
#           file = "mes13/tmp/pathways.csv",
#           row.names = FALSE)
# 
# # Pathway barplot----
# tiff(filename = "mes13/tmp/pathway.tiff",
#      height = 3,
#      width = 10,
#      units = 'in',
#      res = 300,
#      compression = "lzw+p")
# barplot(m1, 
#         showCategory = 10)
# graphics.off()
# 
# # Other plots etc.----
# # Source: https://bioconductor.org/packages/release/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html
# dotplot(m1, showCategory=15)
# 
# enrichMap(m1, 
#           layout = igraph::layout.kamada.kawai, 
#           vertex.label.cex = 1)
# 
# cnetplot(m1, categorySize="pvalue", foldChange = geneID)
# 
# m2 <- gsePathway(geneID,
#                  pvalueCutoff = 1)
# 
# enrichMap(t2)
