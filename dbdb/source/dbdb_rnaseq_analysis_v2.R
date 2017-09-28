# |----------------------------------------------------------------------------------|
# | Project: DB/DB                                                                   |
# | Script: RNA-seq data analysis and visualization                                  |
# | Author: Davit Sargsyan                                                           |
# | Created: 09/19/2017                                                              |
# |----------------------------------------------------------------------------------|
# Header----
# Save consol output to a log file
# sink(file = "tmp/log_dbdb_rnaseq_analysis_v2.txt")

require(data.table)
require(ggplot2)
require(ReactomePA)

# Part I: Data----
# Data Script: .../kidney.ro1/dbdb/source/dbdb_rnaseq_deseq2_v1.R
load("dbdb/data/rna_seq/dbdb.rnaseq.deseq2.RData")
dt1

# Keep DB vs. Control only----
unique(dt1$contr)
dt2 <- droplevels(subset(dt1,
                         contr == "DB vs. Control"))
dt2

# Part II: get specific pathways----
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

# David's selection of pathways, 09/20/2017----
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

genes.from.path <- unique(dt2[dt2$entraz_id%in% get.genes, ])
write.csv(genes.from.path,
          file = "dbdb/tmpdbdb_rna_genes.from.path.csv")

# CHECKPOINT: do the genes I found map to pathways correctly?
m1 <- enrichPathway(gene = unique(genes.from.path$entraz_id),
                    universe = unique(dt1$entraz_id),
                    pvalueCutoff = 0.01, 
                    readable = TRUE, 
                    organism = "mouse")
t1 <- as.data.table(m1)
t1
barplot(m1, 
        showCategory = 50)

tiff(filename = "dbdb/tmp/dbdb_rna_genes.from.path_pathways.tiff",
     height = 6,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
barplot(m1, 
        showCategory = 30)
graphics.off()

# Heatmap of the selected genes----
dt2 <- dt2[order(dt2$log2diff,
                 decreasing = TRUE), ]
dt2 <- dt2[order(dt2$time), ]
dt2
# Top 20 Overexpressed in DB at 16 weeks
tmp <- dt2[dt2$gene %in% dt2$gene[c(1:20,
                                    (nrow(dt2) - 19):(nrow(dt2)))], ]
tmp

p1 <- ggplot(data = tmp) +
  geom_tile(aes(x =  time,
                y = gene,
                fill = log2diff),
            color = "black") +
  scale_fill_gradient2(high = "green",
                       mid = "black",
                       low = "red",
                       midpoint = 0,
                       name = "Difference of\nLog2(RNA)") +
  scale_x_discrete("Time",
                   expand = c(0, 0)) +
  scale_y_discrete("Gene",
                     expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 0,
                                   hjust = 1),
        # legend.position = "top",
        plot.title = element_text(hjust = 0.5))
p1

tiff(filename = "dbdb/tmp/dbdb_rna_genes.from.pathways_heatmap.tiff",
     height = 6,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()
write.csv(genes.from.path,
          file = "dbdb/tmp/dbdb_rna_genes.from.pathways_heatmap.csv")

# sink()