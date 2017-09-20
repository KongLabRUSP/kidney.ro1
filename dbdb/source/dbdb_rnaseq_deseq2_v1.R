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
require(org.Mm.eg.db)
require(DESeq2)
require(BiocParallel)

# Load data----
# ORIGINAL DATA SOURCE: /home/administrator/Documents/David_RNAseq_bgi/combraw.count
# Created by David and Renyi, 07/27/2017
dt1 <- read.table("dbdb/data/rna_seq/combraw.count",
                  sep = "\t",
                  header = TRUE,
                  stringsAsFactors = FALSE)

# Rename the samples
colnames(dt1)[7:14] <- paste("DR", 
                             1:8, 
                             sep = "")

# Part I: Run DESeq2----
# Specify sample groups
grp <- data.frame(trt = rep(rep(c("Control",
                                  "DB"),
                                each = 2),
                            2),
                  time = rep(c("16w",
                               "21w"),
                             each = 4),
                  repl = rep(paste("Repl",
                                   1:2,
                                   sep = ""),
                             4))
grp

# Bild a DESeq2 model with treatment*time interaction----
# NOTE: ADD 1 TO EVERY OBSERVATION TO CALCULATE LOG2!
dds <- DESeqDataSetFromMatrix(dt1[, 7:14] + 1, 
                              grp,
                              ~ trt + time + trt:time)
dds

# Set cores for parallel processing of DESeq----
# snowparam <- SnowParam(workers = snowWorkers(), 
#                        type = "SOCK")
# register(snowparam, 
#          default = TRUE)

# Run DESeq----
dds <- DESeq(dds,
             fitType = "local",
             parallel = TRUE)
resultsNames(dds)



# Contrasts----
# b. DB vs. Control at Week 16----
res.db.vs.ctrl.16w <- results(dds,
                              name = "trt_DB_vs_Control")
# # SAME AS:
# res.db.vs.ctrl <- results(dds,
#                           contrast = c("trt",
#                                        "DB",
#                                        "Control"))

res.db.vs.ctrl.16w <- data.table(gene = dt1$Geneid,
                                 contr = "DB vs. Control",
                                 time = "16 Weeks",
                                 log2diff = res.db.vs.ctrl.16w$log2FoldChange,
                                 padj = res.db.vs.ctrl.16w$padj)
res.db.vs.ctrl.16w

# c. DB vs. Control at Week 21----
res.db.vs.ctrl.21w <- results(dds,
                              contrast = list(c("trt_DB_vs_Control" ,
                                                "trtDB.time21w")))
res.db.vs.ctrl.21w <- data.table(gene = dt1$Geneid,
                                 contr = "DB vs. Control",
                                 time = "21 Weeks",
                                 log2diff = res.db.vs.ctrl.21w$log2FoldChange,
                                 padj = res.db.vs.ctrl.21w$padj)

res.db.vs.ctrl.21w

# d. Week 16 vs. Week 21, Control----
w21.vs.w16.ctrl <- results(dds,
                           name = "time_21w_vs_16w")
w21.vs.w16.ctrl <- data.table(gene = dt1$Geneid,
                                 contr = "Control",
                                 time = "21 vs. 16 Weeks",
                                 log2diff = w21.vs.w16.ctrl$log2FoldChange,
                                 padj = w21.vs.w16.ctrl$padj)
w21.vs.w16.ctrl

# e. Week 16 vs. Week 21, DB----
w21.vs.w16.db <- results(dds,
                         contrast = list(c("time_21w_vs_16w" ,
                                           "trtDB.time21w")))
w21.vs.w16.db <- data.table(gene = dt1$Geneid,
                            contr = "DB",
                            time = "21 vs. 16 Weeks",
                            log2diff = w21.vs.w16.db$log2FoldChange,
                            padj = w21.vs.w16.db$padj)
w21.vs.w16.db 

# Merge all tables
dbdb.rnaseq.deseq2 <- rbindlist(list(res.db.vs.ctrl.16w,
                                     res.db.vs.ctrl.21w,
                                     w21.vs.w16.ctrl,
                                     w21.vs.w16.db))
dbdb.rnaseq.deseq2

# Part II: Map  gene names to Entrez Gene ID (SLOW!)----
gene_names <- unique(dt1$Geneid)

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
dt1 <- merge(gene_map,
             dbdb.rnaseq.deseq2,
             by = "gene")
dt1

# Save all results----
save(dt1,
     file = "dbdb/data/rna_seq/dbdb.rnaseq.deseq2.RData")
write.csv(dt1,
          file = "dbdb/data/rna_seq/dbdb.rnaseq.deseq2.csv")

# sink()