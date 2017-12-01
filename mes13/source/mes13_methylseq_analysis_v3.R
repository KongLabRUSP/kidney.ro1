# |----------------------------------------------------------------------------------|
# | Project: Study of Diabetes in MES13 cells                                        |
# | Script: Methyl-seq data analysis and visualization                               |
# | Author: Davit Sargsyan                                                           |
# | Created: 09/12/2017                                                              |
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
trt.names <- c("LG",
               "HG",
               "MIC 1.5uM",
               "FX 1uM",
               "Ber 6uM")

# Load data----
peakAnno1 <- annotatePeak(peak = "mes13/data/methyl_seq/combined.csv", 
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
                  dt1[, 2:11], 
                  geneName = dt1$GENENAME)

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
for (i in 4:14) {
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
t1 <- dt.reg[, c(1, 2, 3, 5, 7, 9, 11)]
t1[, 3:7] <- do.call("cbind",
                     lapply(t1[, 3:7],
                   function(a) {
                     return(round(a/t1[, 2],
                                  1))
                   }))
colnames(t1) <- c("Gene Region",
                  "Total CpG Count",
                  trt.names)
kable(t1)
  # |Gene Region | Total CpG Count|  LG|  HG| MIC 1.5uM| FX 1uM| Ber 6uM|
  # |:-----------|---------------:|---:|---:|---------:|------:|-------:|
  # |Promoter    |           22856| 0.6| 1.3|       1.2|    0.9|     1.2|
  # |Body        |            5153| 1.0| 1.9|       1.6|    1.5|     1.6|
  # |Downstream  |            6458| 2.2| 2.9|       3.0|    4.1|     3.2|
write.csv(t1,
          file = "mes13/tmp/t1.csv",
          row.names = FALSE)

# Calculate percent methylation in each sample----
dt.reg <- data.table(Region = dt.reg$Group.1,
                     CpG = dt.reg$CpG,
                     pct = round(100*dt.reg[, c(4, 6, 8, 10, 12)]/
                                   dt.reg[, c(3, 5, 7, 9, 11)],
                                 1))
dt.reg

# Melt data
dt.reg.l <- melt.data.table(data = dt.reg,
                            id.vars = 1:2,
                            measure.vars = 3:7,
                            variable.name = "Treatment",
                            value.name = "Methylation(%)")

dt.reg.l$Treatment <- factor(dt.reg.l$Treatment)
levels(dt.reg.l$Treatment) <- trt.names
dt.reg.l

# NEW (09/27/2017): remove Ber
dt.reg.l <- droplevels(subset(dt.reg.l,
                              Treatment != "Ber 6uM"))

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

# Part II: gene methylation----
# Collapse by gene
out <- list()
for (i in 4:14) {
  out[[i - 3]] <- aggregate(dt1[, i, with = FALSE],
                            by = list(dt1$gene,
                                      dt1$reg),
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
                      pct = round(100*dt.gene[, c(5, 7, 9, 11, 13)]/
                                    dt.gene[, c(4, 6, 8, 10, 12)],
                                  1))
colnames(dt.gene)[4:8] <- trt.names
dt.gene

# NEW (09/27/2017): remove Ber and keep promoter only----
dt.gene <- droplevels(dt.gene[dt.gene$reg == "Promoter", -c("Ber 6uM")])
dt.gene

# SAVE THIS VERSION
dt2 <- dt.gene

# NEW (09/27/2017): keep genes with values for  HG and 
# at least one of the treatments
dt.gene <- subset(dt.gene,
                  !is.nan(HG) &
                    (!is.nan(`MIC 1.5uM`) | 
                       !is.nan(`FX 1uM`)))

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
summary(dt.gene.l)

# Genes with largest window (Positive - negative control methylation)----
dt.gene$hg_lg <- dt.gene$HG - dt.gene$LG
m.diff <- dt.gene[!is.nan(hg_lg), c(1, 2, 8)]
m.diff <- m.diff[order(hg_lg), ]
m.diff

write.csv(m.diff,
          file = "mes13/tmp/mes13_methyl_hg-lg_all_genes.csv")

# Genes with largest change in methylation (HG - LG)----
gene.sorted <- unique(as.character(m.diff$gene))

# a. Highest Positive Difference (HG - LG)----
gene.up <- gene.sorted[(length(gene.sorted) - 19):length(gene.sorted)]
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
  #ggtitle("MES13: Top 30 Largest Positive Diff in Methylation \n High - Low Glucose") +
  theme(axis.text.x = element_text(angle = 30,
                                   hjust = 1),
        # legend.position = "top",
        plot.title = element_text(hjust = 0.5))
p2a

# b. Highest Negative Difference (HG - LG)----
gene.dn <- gene.sorted[1:20]
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
  # ggtitle("MES13: Top 30 Largest Negative Diff in Methylation \n High - Low Glucose") +
  theme(axis.text.x = element_text(angle = 30,
                                   hjust = 1),
        # legend.position = "top",
        plot.title = element_text(hjust = 0.5))
p2b

tiff(filename = "mes13/tmp/top_meth_diff.tiff",
     height = 5,
     width = 12,
     units = 'in',
     res = 300,
     compression = "lzw+p")
gridExtra::grid.arrange(p2a, p2b, p1, nrow = 1)
graphics.off()

# c. Dr. Kong's request (09/27/2017): heatmap of all (HG-LG)----
gene.dn <- gene.sorted
tmp <- subset(dt.gene.l,
              as.character(gene) %in% gene.dn)
tmp

p2c <- ggplot(data = tmp) +
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
  # ggtitle("MES13: Top 30 Largest Negative Diff in Methylation \n High - Low Glucose") +
  theme(axis.text.x = element_text(angle = 30,
                                   hjust = 1),
        # legend.position = "top",
        plot.title = element_text(hjust = 0.5))
p2c

tiff(filename = "mes13/tmp/mes13_meth_hg-lg_all.tiff",
     height = 15,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p2c)
graphics.off()
write.csv(tmp,
          file = "mes13/tmp/mes13_meth_hg-lg_all.csv")

# Save the tables----
gene.meth.save <- dt.gene
dt1[, c("gene",
        "geneName")]
gene.meth.save <- merge(unique(dt1[, c("gene",
                                       "geneId",
                                       "geneName")]),
                        gene.meth.save,
                        by = "gene")
gene.meth.save

write.csv(gene.meth.save ,
          file = "mes13/tmp/gene.meth.save.csv")

# Part III: Pathway analysis with Reactome----
# Get Entrezgene IDs
geneID <- as.character(unique(dt1$geneId[as.character(dt1$gene) %in% gene.sorted]))

# Save Entrezgene IDs and gene names
out <- data.table(geneID = geneID,
                  geneName = gene.sorted)
write.csv(out,
          file = "mes13/tmp/gene.entrezID.csv",
          row.names = FALSE)

# Get pathways----
m1 <- enrichPathway(gene = geneID,
                    pvalueCutoff = 1, 
                    readable = T, 
                    organism = "mouse")
t2 <- as.data.table(m1)
t2
write.csv(t2,
          file = "mes13/tmp/pathways.csv",
          row.names = FALSE)

# Pathway barplot----
tiff(filename = "mes13/tmp/pathway.tiff",
     height = 3,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
barplot(m1)
graphics.off()

# Part IV: get specific pathways----
path.names <- do.call("rbind", 
                      as.list(reactome.db::reactomePATHID2NAME))
path.names <- data.table(reactome_id = rownames(path.names),
                         path = path.names[, 1])
path.names

# Keep mouse pathways only
path.mm <- subset(path.names,
                  substr(path, 1, 12) == "Mus musculus")
path.mm$path <- gsub(pattern = "Mus musculus: ",
                     replacement = "",
                     x = path.mm$path)
path.mm

##  David's selection of pathways, 09/20/2017----
# TNF signaling
# Activation of NF-kappaB in B cells
# Activation of the AP-1 family of trTNF signalinganscription factors
# Biological oxidations
# COX reactions
# Cell-extracellular matrix interactions
# Extracellular matrix organization
# IL-6-type cytokine receptor ligand interactions
# Interleukin-6 family signaling
# Regulation of TNFR1 signaling
# TNFR2 non-canonical NF-kB pathway
# NF-kB is activated and signals survival
# Inflammasomes
# Oxidative Stress Induced Senescence
# TGF-beta receptor signaling activates SMADs
# TGF-beta receptor signaling in EMT (epithelial to mesenchymal transition)

# Separate TNF signaling pathway----
path.list <- c("TNF signaling",
               "Activation of NF-kappaB in B cells",
               "Activation of the AP-1 family of trTNF signalinganscription factors",
               "Biological oxidations",
               "COX reactions",
               "Cell-extracellular matrix interactions",
               "Extracellular matrix organization",
               "IL-6-type cytokine receptor ligand interactions",
               "Interleukin-6 family signaling",
               "Regulation of TNFR1 signaling",
               "TNFR2 non-canonical NF-kB pathway",
               "NF-kB is activated and signals survival",
               "Inflammasomes",
               "Oxidative Stress Induced Senescence",
               "TGF-beta receptor signaling activates SMADs",
               "TGF-beta receptor signaling in EMT (epithelial to mesenchymal transition)")
path.mm[path %in% path.list]

# Map path IDs to Entrez gene IDs----
path2gene <- as.list(reactome.db::reactomePATHID2EXTID)
path2gene <- path2gene[names(path2gene) %in% path.mm$reactome_id]

# Get TNF signaling gene Entrez IDs----
get.genes <- path2gene[names(path2gene) %in% path.mm$reactome_id[path.mm$path %in% path.list]]
get.genes <- unique(do.call("c", get.genes))

dt2 <- merge(dt1[, 1:2],
             dt2, 
             by = "gene")
genes.from.path <- unique(dt2[dt2$geneId %in% get.genes, ])
write.csv(genes.from.path,
          file = "mes13/tmp/mes13_meth_genes.from.path.csv")

# CHECKPOINT: do the genes I found map to pathways correctly?
m1 <- enrichPathway(gene = as.character(genes.from.path$geneId),
                    pvalueCutoff = 0.01, 
                    readable = TRUE, 
                    organism = "mouse")
t1 <- as.data.table(m1)
t1
barplot(m1, 
        showCategory = 10)

tiff(filename = "mes13/tmp/mes13_meth_genes.from.path_pathways.tiff",
     height = 6,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
barplot(m1, 
        showCategory = 20)
graphics.off()

# Heatmap of the selected genes----
dt3 <- melt.data.table(data = genes.from.path,
                       id.vars = c(1, 3),
                       measure.vars = 5:8,
                       variable.name = "Treatment",
                       value.name = "Methylation(%)")

dt3$Treatment <- factor(dt3$Treatment)
dt3

p3 <- ggplot(data = dt3) +
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
  # ggtitle("MES13: Top 30 Largest Negative Diff in Methylation \n High - Low Glucose") +
  theme(axis.text.x = element_text(angle = 30,
                                   hjust = 1),
        # legend.position = "top",
        plot.title = element_text(hjust = 0.5))
p3

tiff(filename = "mes13/tmp/mes13_meth_genes.from.pathways_heatmap.tiff",
     height = 15,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p3)
graphics.off()
write.csv(genes.from.path,
          file = "mes13/tmp/mes13_meth_genes.from.pathways_heatmap.csv")