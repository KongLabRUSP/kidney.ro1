# |----------------------------------------------------------------------------------|
# | Project: Study of Diabetes in MES13 cells                                        |
# | Script: RNA-seq data analysis and visualization                                  |
# | Author: Davit Sargsyan                                                           |
# | Created: 09/16/2017                                                              |
# |----------------------------------------------------------------------------------|
# Header----
# Save consol output to a log file
# sink(file = "tmp/log_mes13_rnaseq_figures_v1.txt")

require(data.table)
require(ggplot2)
require(org.Mm.eg.db)
require(ReactomePA)

# Load data----
dt1 <- fread("mes13/data/rna_seq/gene_exp.csv")
dt1

# Map  gene names to Entrez Gene ID
gene_names <- unique(dt1$gene)

entraz_id <- list()
for (i in 1:length(gene_names)) {
  out <- try(unlist(mget(x = gene_names[i],
                         envir = org.Mm.egALIAS2EG)),
             silent = TRUE)
  if (class(out)[1] != "try-error") {
    entraz_id[i] <- out
  } else {
    entraz_id[i] <- NA
  }
}

gene_map <- data.table(gene = gene_names,
                       entraz_id = do.call("c", 
                                           entraz_id))

# Merge Entraz IDs with the data
dt1 <- merge(dt1, 
             gene_map,
             by = "gene")

# a. Subset LG vs. HG----
lg_hg <- subset(dt1,
                sample_1 == "LG" & 
                  sample_2 == "HG",
                select = c(1, 15, 5, 6, 8, 9))
lg_hg$value_1 <- log2(lg_hg$value_1 + 1)
lg_hg$value_2 <- log2(lg_hg$value_2 + 1)
lg_hg$log2diff <- lg_hg$value_1 - lg_hg$value_2
lg_hg

# Subset to at least 2-fold change----
dt2 <- subset(lg_hg,
              log2diff >= 1|
                log2diff <= -1)

# Pathways----
m1 <- enrichPathway(gene = dt2$entraz_id,
                    pvalueCutoff = 1, 
                    readable = TRUE, 
                    organism = "mouse")
t1 <- as.data.table(m1)
t1
write.csv(t1,
          file = "mes13/tmp/mes13_rnaseq_lg-hg_pathway.csv",
          row.names = FALSE)

tiff(filename = "mes13/tmp/rnaseq_lg-hg_pathway.tiff",
     height = 3,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
barplot(m1, 
        showCategory = 30)
graphics.off()

# Network plot----
geneList <- dt2$log2diff[!duplicated(dt2$entraz_id)]
names(geneList) <- dt2$entraz_id[!duplicated(dt2$entraz_id)]

tiff(filename = "mes13/tmp/rnaseq_lg-hg_network.tiff",
     height = 10,
     width = 12,
     units = 'in',
     res = 300,
     compression = "lzw+p")
viewPathway(pathName = "RNA Polymerase III Transcription",
            readable = TRUE,
            organism = "mouse",
            foldChange = geneList)
graphics.off()

# Heatmap of genes involved in known pathways-----
genes.keep <- unique(do.call("c",
                             strsplit(m1@result$geneID,
                                      split = "/")))
tmp <- subset(dt1,
              gene %in% genes.keep,
              select = c(1, 5, 6, 8, 9))
tmp$contr <- paste(tmp$sample_2,
                   "-",
                   tmp$sample_1)
tmp <- subset(tmp,
              contr %in% c("HG - LG",
                           "FX_1_uM - HG",
                           "TIIA_5_uM - HG",
                           "Ber_6_uM - HG",
                           "MIC_1_5_uM - HG"))
tmp$contr <- factor(tmp$contr,
                    levels = c("HG - LG",
                               "FX_1_uM - HG",
                               "TIIA_5_uM - HG",
                               "Ber_6_uM - HG",
                               "MIC_1_5_uM - HG"))
tmp$log2diff <- log2(tmp$value_2 + 1) - log2(tmp$value_1 ++ 1)
  
p1a <- ggplot(data = tmp) +
  geom_tile(aes(x =  contr,
                y = gene,
                fill = log2diff),
            color = "black") +
  scale_fill_gradient2(high = "green",
                       mid = "black",
                       low = "red",
                       midpoint = 0,
                       name = "Difference of\nLog2(RNA)") + 
  scale_x_discrete(expand = c(0, 0),
                   "Treatment Difference") +
  scale_y_discrete("Gene",
                   expand = c(0, 0)) +
  ggtitle("MES13 RNA-Seq: Genes With At Least 2-Fold Change\nAnd Mapped To Known Pathways") +
  theme(axis.text.x = element_text(angle = 30,
                                   hjust = 1),
        plot.title = element_text(hjust = 0.5))
p1a

tiff(filename = "mes13/tmp/rnaseq_lg-hg_heatmap.tiff",
     height = 7,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1a)
graphics.off()

# b. Subset FX_1_uM - HG----
fx_hg <- subset(dt1,
                sample_1 == "HG"  & 
                  sample_2 == "FX_1_uM",
                select = c(1, 15, 5, 6, 8, 9))
fx_hg$value_1 <- log2(fx_hg$value_1 + 1)
fx_hg$value_2 <- log2(fx_hg$value_2 + 1)
fx_hg$log2diff <- fx_hg$value_2 - fx_hg$value_1
fx_hg

# Subset to at least 2-fold change
dt2 <- subset(fx_hg,
              log2diff >= 2 |
                log2diff <= -2)

# Pathways----
m1 <- enrichPathway(gene = dt2$entraz_id,
                    pvalueCutoff = 1, 
                    readable = TRUE, 
                    organism = "mouse")
t1 <- as.data.table(m1)
t1
write.csv(t1,
          file = "mes13/tmp/rnaseq_fx-hg_pathway.csv",
          row.names = FALSE)

tiff(filename = "mes13/tmp/rnaseq_fx-hg_pathway.tiff",
     height = 3,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
barplot(m1, 
        showCategory = 10)
graphics.off()

# Network plot----
geneList <- dt2$log2diff[!duplicated(dt2$entraz_id)]
names(geneList) <- dt2$entraz_id[!duplicated(dt2$entraz_id)]
geneList[geneList > 3] <- 3
geneList[geneList < -3] <- -3

tiff(filename = "mes13/tmp/rnaseq_fx-hg_network.tiff",
     height = 30,
     width = 35,
     units = 'in',
     res = 300,
     compression = "lzw+p")
viewPathway(pathName = "Biological oxidations",
            readable = TRUE,
            organism = "mouse",
            foldChange = geneList)
graphics.off()

# Heatmap of genes involved in known pathways-----
genes.keep <- unique(do.call("c",
                             strsplit(m1@result$geneID,
                                      split = "/")))
tmp <- subset(dt1,
              gene %in% genes.keep,
              select = c(1, 5, 6, 8, 9))
tmp$contr <- paste(tmp$sample_2,
                   "-",
                   tmp$sample_1)
tmp <- subset(tmp,
              contr %in% c("HG - LG",
                           "FX_1_uM - HG",
                           "TIIA_5_uM - HG",
                           "Ber_6_uM - HG",
                           "MIC_1_5_uM - HG"))
tmp$contr <- factor(tmp$contr,
                    levels = c("HG - LG",
                               "FX_1_uM - HG",
                               "TIIA_5_uM - HG",
                               "Ber_6_uM - HG",
                               "MIC_1_5_uM - HG"))
tmp$log2diff <- log2(tmp$value_2 + 1) - log2(tmp$value_1 ++ 1)

p1b <- ggplot(data = tmp) +
  geom_tile(aes(x =  contr,
                y = gene,
                fill = log2diff),
            color = "black") +
  scale_fill_gradient2(high = "green",
                       mid = "black",
                       low = "red",
                       midpoint = 0,
                       name = "Difference of\nLog2(RNA)") + 
  scale_x_discrete(expand = c(0, 0),
                   "Treatment Difference") +
  scale_y_discrete("Gene",
                   expand = c(0, 0)) +
  ggtitle("MES13 RNA-Seq: Genes With At Least 8-Fold Change\nAnd Mapped To Known Pathways") +
  theme(axis.text.x = element_text(angle = 30,
                                   hjust = 1),
        plot.title = element_text(hjust = 0.5))
p1b

tiff(filename = "mes13/tmp/rnaseq_fx-hg_heatmap.tiff",
     height = 7,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1b)
graphics.off()

# c. Subset TIIA_5_uM vs. HG----
ti_hg <- subset(dt1,
                sample_1 == "HG" & 
                  sample_2 == "TIIA_5_uM",
                select = c(1, 15, 5, 6, 8, 9))
ti_hg$value_1 <- log2(ti_hg$value_1 + 1)
ti_hg$value_2 <- log2(ti_hg$value_2 + 1)
ti_hg$log2diff <- ti_hg$value_1 - ti_hg$value_2
ti_hg

# Subset to at least 2-fold change----
dt2 <- subset(ti_hg,
              log2diff >= 2.5|
                log2diff <= -2.5)

# Pathways----
m1 <- enrichPathway(gene = dt2$entraz_id,
                    pvalueCutoff = 1, 
                    readable = TRUE, 
                    organism = "mouse")
t1 <- as.data.table(m1)
t1
write.csv(t1,
          file = "mes13/tmp/rnaseq_ti-hg_pathway.csv",
          row.names = FALSE)

tiff(filename = "mes13/tmp/rnaseq_ti-hg_pathway.tiff",
     height = 3,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
barplot(m1, 
        showCategory = 30)
graphics.off()

# Network plot----
geneList <- dt2$log2diff[!duplicated(dt2$entraz_id)]
names(geneList) <- dt2$entraz_id[!duplicated(dt2$entraz_id)]
geneList[geneList > 3] <- 3
geneList[geneList < -3] <- -3

tiff(filename = "mes13/tmp/rnaseq_ti-hg_network.tiff",
     height = 30,
     width = 35,
     units = 'in',
     res = 300,
     compression = "lzw+p")
viewPathway(pathName = "Biological oxidations",
            readable = TRUE,
            organism = "mouse",
            foldChange = geneList)
graphics.off()

# Heatmap of genes involved in known pathways-----
genes.keep <- unique(do.call("c",
                             strsplit(m1@result$geneID,
                                      split = "/")))
tmp <- subset(dt1,
              gene %in% genes.keep,
              select = c(1, 5, 6, 8, 9))
tmp$contr <- paste(tmp$sample_2,
                   "-",
                   tmp$sample_1)
tmp <- subset(tmp,
              contr %in% c("HG - LG",
                           "FX_1_uM - HG",
                           "TIIA_5_uM - HG",
                           "Ber_6_uM - HG",
                           "MIC_1_5_uM - HG"))
tmp$contr <- factor(tmp$contr,
                    levels = c("HG - LG",
                               "FX_1_uM - HG",
                               "TIIA_5_uM - HG",
                               "Ber_6_uM - HG",
                               "MIC_1_5_uM - HG"))
tmp$log2diff <- log2(tmp$value_2 + 1) - log2(tmp$value_1 ++ 1)

p1c <- ggplot(data = tmp) +
  geom_tile(aes(x =  contr,
                y = gene,
                fill = log2diff),
            color = "black") +
  scale_fill_gradient2(high = "green",
                       mid = "black",
                       low = "red",
                       midpoint = 0,
                       name = "Difference of\nLog2(RNA)") + 
  scale_x_discrete(expand = c(0, 0),
                   "Treatment Difference") +
  scale_y_discrete("Gene",
                   expand = c(0, 0)) +
  ggtitle("MES13 RNA-Seq: Genes With At Least 2.5 Log2 Difference\nAnd Mapped To Known Pathways") +
  theme(axis.text.x = element_text(angle = 30,
                                   hjust = 1),
        plot.title = element_text(hjust = 0.5))
p1c

tiff(filename = "mes13/tmp/rnaseq_ti-hg_heatmap.tiff",
     height = 7,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1a)
graphics.off()

# sink()