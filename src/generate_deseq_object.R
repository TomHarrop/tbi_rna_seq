#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log,
     type = "message")
sink(log,
     append = TRUE,
     type = "output")

library(data.table)
library(DESeq2)

star_files <- snakemake@input[["quant_files"]]
dds_file <- snakemake@output[["dds"]]

# dev
# star_files <- list.files("output/025_star",
#            pattern = "ReadsPerGene",
#            recursive = TRUE,
#            full.names = TRUE)
# dds_file <- "tmp/dds.Rds"

names(star_files) <- sub(".ReadsPerGene.out.tab",
                         "",
                         basename(star_files))

star_list <- lapply(star_files, fread)
star_counts_raw <- rbindlist(star_list, idcol = "sample")

# remove unmapped reads
star_counts <- star_counts_raw[!startsWith(V1, "N_")]

# take a peek at the counts
# column 1:  gene ID
# column 2:  counts for unstranded RNA-seq
# column 3:  counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
# column 4:  counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)
x <- star_counts[, lapply(.SD, sum), .SDcols = c("V2", "V3", "V4"),
                 by = sample]
# head(x)

# generate count matrix
# illumina truseq stranded == column 4
sample_order <- star_counts[, gtools::mixedsort(unique(sample))]
count_dt <- dcast(star_counts, V1 ~ sample, value.var = "V4")
setcolorder(count_dt, sample_order)

# generate col_data
col_data <- data.table(
    names = sample_order
)
col_data[, splitname := gsub("^(n[[:digit:]]+)([d|f])([[:digit:]]+)",
                             "\\1_\\2_\\3",
                             names)]
col_data[, c("hemisphere", "sex", "injury", "indiv") := tstrsplit(splitname, "_")]
col_data[, indiv := as.factor(gsub("[^[:digit:]]+", "", indiv))]
col_data[, splitname := NULL]

dds <- DESeqDataSetFromMatrix(
    countData = as.matrix(data.frame(count_dt, row.names = "V1")),
    colData = data.frame(col_data, row.names = "names"),
    design = ~ sex + hemisphere + injury)

# write output
saveRDS(dds, dds_file)

# log
sessionInfo()
