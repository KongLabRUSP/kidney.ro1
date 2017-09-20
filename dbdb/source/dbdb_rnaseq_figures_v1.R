# |----------------------------------------------------------------------------------|
# | Project: DB/DB                                                                   |
# | Script: RNA-seq data analysis and visualization                                  |
# | Author: Davit Sargsyan                                                           |
# | Created: 09/19/2017                                                              |
# |----------------------------------------------------------------------------------|
# Header----
# Save consol output to a log file
# sink(file = "tmp/log_dbdb_rnaseq_analysis_v1.txt")

require(data.table)
require(ggplot2)
require(ReactomePA)

# Load data----
# Data Script: .../kidney.ro1/dbdb/source/dbdb_rnaseq_deseq2_v1.R
load("dbdb/data/rna_seq/dbdb.rnaseq.deseq2.RData")
dt1

# a. DB vs. Control at 16 Weeks----
dt.16w <- subset(dt1,
                  contr == "DB vs. Control" &
                   time == "16 Weeks")
dt.16w

# Subset to at least log2 +/-2----
dt2 <- subset(dt.16w,
              log2diff >= 1.5|
                log2diff <= -1.5)

# Save gene list for Reactome
# geneList <- dt2$log2diff
# names(geneList) <- dt2$entraz_id
geneList <- as.numeric(dt2$entraz_id)
write.csv(geneList,
          file = "dbdb/tmp/geneList_dbdb_16w.csv",
          row.names = FALSE)

# Pathways----
m1 <- enrichPathway(gene = dt2$entraz_id,
                    pvalueCutoff = 1, 
                    readable = TRUE, 
                    organism = "mouse")
t1 <- as.data.table(m1)
t1
write.csv(t1,
          file = "dbdb/tmp/rnaseq_db_ctrl_16w_pathway.csv",
          row.names = FALSE)

tiff(filename = "dbdb/tmp/rnaseq_db_ctrl_16w_pathway.tiff",
     height = 10,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
barplot(m1, 
        showCategory = 30,
        main = "DB vs. Control at 16 Weeks")
graphics.off()

# Network plot----
geneList <- dt2$log2diff[!duplicated(dt2$entraz_id)]
names(geneList) <- dt2$entraz_id[!duplicated(dt2$entraz_id)]
geneList[geneList > 3] <- 3
geneList[geneList < -3] <- -3

tiff(filename = "dbdb/tmp/rnaseq_db_ctrl_16w_network.tiff",
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
              gene %in% genes.keep &
                contr == "DB vs. Control")
  
p1a <- ggplot(data = tmp) +
  geom_tile(aes(x =  time,
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
  ggtitle("DB/DB RNA-Seq: Genes With Log2 Change Of At Least 1.5 \nAt 16 Weeks And Mapped To Known Pathways") +
  theme(axis.text.x = element_text(angle = 30,
                                   hjust = 1),
        plot.title = element_text(hjust = 0.5))
p1a

tiff(filename = "dbdb/tmp/rnaseq_db_ctrl_16w_heatmap.tiff",
     height = 9,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1a)
graphics.off()

# b. DB vs. Control at 21 Weeks----
dt.21w <- subset(dt1,
                 contr == "DB vs. Control" &
                   time == "21 Weeks")
dt.21w

# Subset to at least log2 +/-1.5----
dt2 <- subset(dt.21w,
              log2diff >= 1.5|
                log2diff <= -1.5)

# Pathways----
m1 <- enrichPathway(gene = dt2$entraz_id,
                    pvalueCutoff = 1, 
                    readable = TRUE, 
                    organism = "mouse")
t1 <- as.data.table(m1)
t1
write.csv(t1,
          file = "dbdb/tmp/rnaseq_db_ctrl_21w_pathway.csv",
          row.names = FALSE)

tiff(filename = "dbdb/tmp/rnaseq_db_ctrl_21w_pathway.tiff",
     height = 10,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
barplot(m1, 
        showCategory = 30,
        main = "DB vs. Control at 21 Weeks")
graphics.off()

# Network plot----
geneList <- dt2$log2diff[!duplicated(dt2$entraz_id)]
names(geneList) <- dt2$entraz_id[!duplicated(dt2$entraz_id)]
geneList[geneList > 3] <- 3
geneList[geneList < -3] <- -3

tiff(filename = "dbdb/tmp/rnaseq_db_ctrl_21w_network.tiff",
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
              gene %in% genes.keep &
                contr == "DB vs. Control")

p1a <- ggplot(data = tmp) +
  geom_tile(aes(x =  time,
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
  ggtitle("DB/DB RNA-Seq: Genes With Log2 Change Of At Least 1.5 \nAt 21 Weeks And Mapped To Known Pathways") +
  theme(axis.text.x = element_text(angle = 30,
                                   hjust = 1),
        plot.title = element_text(hjust = 0.5))
p1a

tiff(filename = "dbdb/tmp/rnaseq_db_ctrl_21w_heatmap.tiff",
     height = 9,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1a)
graphics.off()

# sink()