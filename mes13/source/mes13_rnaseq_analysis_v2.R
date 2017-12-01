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
dt1

# Get specific pathways----
path.names <- do.call("rbind", 
                      as.list(reactome.db::reactomePATHID2NAME))
path.names <- data.table(reactome_id = rownames(path.names),
                         path = path.names[, 1])
path.names

# Keep only the genes from the following pathways----
# Keep mouse pathways only
path.mm <- subset(path.names,
                  substr(path, 1, 12) == "Mus musculus")
path.mm$path <- gsub(pattern = "Mus musculus: ",
                     replacement = "",
                     x = path.mm$path)
path.mm

##  David's selection of pathways, 09/20/2017----
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
genes.from.path <- dt1[dt1$entraz_id %in% get.genes, ]
genes.from.path <- subset(genes.from.path,
                          select = c(1, 15, 5, 6, 8, 9))
# Remove Ber
genes.from.path <- droplevels(subset(genes.from.path,
                                     sample_2 != "Ber_6_uM"))
genes.from.path

# Keep only certain comparisons
genes.from.path <- droplevels(subset(genes.from.path, 
                                     (sample_2 == "HG" & sample_1 == "LG") |
                                       (sample_2 == "MIC_1_5_uM" & sample_1 == "HG") |
                                       (sample_2 == "FX_1_uM" & sample_1 == "HG")))
genes.from.path

# CHECKPOINT: do the genes I found map to TNFa correctly?
m1 <- enrichPathway(gene = get.genes,
                    pvalueCutoff = 0.01, 
                    readable = TRUE, 
                    organism = "mouse")
t1 <- as.data.table(m1)
t1
barplot(m1, 
        showCategory = 10)

# Create a variable of contrast labels----
dt2 <- genes.from.path
dt2$contr <- paste(dt2$sample_2,
                   "-",
                   dt2$sample_1)
dt2$contr <- factor(dt2$contr,
                    levels = unique(dt2$contr))
write.csv(dt2,
          file = "mes13/tmp/mes13_rna_genes.from.path.csv")

# Separate genes with |log2| > 2
genes.keep <- unique(subset(dt2,
                     log2diff > 2 |
                       log2diff < -2)$gene)
genes.keep
tmp <- subset(dt2,
              gene %in% genes.keep)

p1 <- ggplot(data = tmp) +
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
  # ggtitle("MES13 RNA-Seq: Genes With At Least 2-Fold Change\nAnd Mapped To Known Pathways") +
  theme(axis.text.x = element_text(angle = 20,
                                   hjust = 1),
        plot.title = element_text(hjust = 0.5))
p1

tiff(filename = "mes13/tmp/mes13_rna_heatmap.tiff",
     height = 7,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# sink()