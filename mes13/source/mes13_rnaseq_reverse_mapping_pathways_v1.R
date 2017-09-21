# |----------------------------------------------------------------------------------|
# | Project: Study of Diabetes in MES13 cells                                        |
# | Script: RNA-seq reverse mathching to diabetes-related pathways                   |
# | Author: Davit Sargsyan                                                           |
# | Created: 09/20/2017                                                              |
# |----------------------------------------------------------------------------------|
# Header----
# References:
# 1. https://www.bioconductor.org/packages/devel/bioc/vignettes/AnnotationDbi/inst/doc/IntroToAnnotationPackages.pdf
# 2. https://www.biostars.org/p/70264/

# Save consol output to a log file
# sink(file = "tmp/log_mes13_rnaseq_figures_v1.txt")

require(data.table)
require(ggplot2)
require(org.Mm.eg.db)
require(ReactomePA)

# Load data----
dt1 <- fread("mes13/data/rna_seq/gene_exp.csv")

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
dt1 <- dt1[, c(1, 15, 5, 6, 8, 9)]
dt1

# # Pathways databases----
# # 1. org.Mm.eg.db----
# ls(org.Mm.eg.db)
# org.Mm.eg.db
# 
# columns(org.Mm.eg.db)
# help("PATH")
# help("SYMBOL")
# 
# keytypes(org.Mm.eg.db)
# keys(org.Mm.eg.db)
# 
# ll <- select(org.Mm.eg.db,
#              keys = keys(org.Mm.eg.db),
#              columns = c("SYMBOL","PATH"))
# ll
#
# # 2. KEGG----
# require(KEGG.db)
# path.names <- do.call("rbind",
#                       as.list(KEGG.db::KEGGPATHID2NAME))

# 3. REACTOME----
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
path.mm[path == "TNF signaling"]

# Map path IDs to Entrez gene IDs----
path2gene <- as.list(reactome.db::reactomePATHID2EXTID)
path2gene <- path2gene[names(path2gene) %in% path.mm$reactome_id]

# Get TNF signaling gene Entrez IDs----
get.genes <- path2gene[names(path2gene) == path.mm$reactome_id[path.mm$path == "TNF signaling"]][[1]]
dt1[dt1$entraz_id %in% get.genes, ]

# CHECKPOINT: do the genes I found map to TNFa correctly?
m1 <- enrichPathway(gene = get.genes,
                    pvalueCutoff = 0.01, 
                    readable = TRUE, 
                    organism = "mouse")
t1 <- as.data.table(m1)
t1
barplot(m1, 
        showCategory = 10)
# Yes

# Find all expressed genes from the list of TNFa genes, HG vs LG----
tmp <- dt1[dt1$entraz_id %in% get.genes &
             dt1$sample_1 == "LG" &
             dt1$sample_2 == "HG", ]
# tmp <- dt1[dt1$entraz_id %in% get.genes &
#              dt1$sample_1 == "HG" &
#              dt1$sample_2 == "FX_1_uM", ]
length(unique(tmp$gene))
tmp$log2diff <- log2(tmp$value_2 + 1) - log2(tmp$value_1 + 1)
tmp

# Prepare data for viewPathway function----
geneList <- tmp$log2diff[!duplicated(tmp$entraz_id)]
names(geneList) <- tmp$entraz_id[!duplicated(tmp$entraz_id)]
length(geneList)
geneList

viewPathway(pathName = "TNF signaling",
            readable = TRUE,
            organism = "mouse",
            foldChange = geneList)

tiff(filename = "mes13/tmp/rnaseq_hg-lg_tnfa_pathway.tiff",
     height = 15,
     width = 20,
     units = 'in',
     res = 300,
     compression = "lzw+p")
# tiff(filename = "mes13/tmp/rnaseq_fx-hg_tnfa_pathway.tiff",
#      height = 15,
#      width = 20,
#      units = 'in',
#      res = 300,
#      compression = "lzw+p")
viewPathway(pathName = "TNF signaling",
            readable = TRUE,
            organism = "mouse",
            foldChange = geneList)
graphics.off()