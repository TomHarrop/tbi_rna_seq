#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log,
     type = "message")
sink(log,
     append = TRUE,
     type = "output")

library(DESeq2)
library(data.table)
library(ggplot2)

dds_file <- snakemake@input[["dds"]]
inj_res_file <- snakemake@output[["injury_results"]]
group1_res_file <- snakemake@output[["group1_results"]]
group2_res_file <- snakemake@output[["group2_results"]]
threads <- snakemake@threads

# dev
# dds_file <- "tmp/dds.Rds"
# threads <- 8

BiocParallel::register(
    BiocParallel::MulticoreParam(workers = 8))

#############
# FILTERING #
#############

dds <- readRDS(dds_file)
dds <- estimateSizeFactors(dds)

# remove samples that were re-run
keep_samples <- colnames(dds)[!grepl("[b|c]$", colnames(dds))]

# remove genes with low expression
keep_genes <- rowSums(counts(dds, normalized = TRUE) >= 10) >= 3

dds_filtered <- dds[keep_genes, keep_samples]

# add group factors for later
dds_filtered$group1 <- factor(paste(
    dds_filtered$sex,
    dds_filtered$injury,
    sep = "_"))
dds_filtered$group2 <- factor(paste(
    dds_filtered$hemisphere,
    dds_filtered$injury,
    sep = "_"))

##########
# TEST 1 #
##########

# Do any genes correlate with injury, across sex and hemisphere?
design(dds_filtered) <- ~ sex + hemisphere + injury
dds_filtered <- DESeq(dds_filtered, parallel = TRUE)

# extract the results for contrast of interest
# can play with lfcThreshold and alpha (both arbitrary)
injury_res <- data.table(results(dds_filtered,
                             name = "injury_T_vs_S",
                             tidy = TRUE,
                             lfcThreshold = log(1.5, 2),
                             alpha = 0.1))

# 
injury_res[, al2fc := abs(log2FoldChange)]
setorder(injury_res, -al2fc)
injury_res[, al2fc := NULL]

##########
# TEST 2 #
##########

# Do any genes correlate with injury, across hemispheres but specifically in
# females?
design(dds_filtered) <- ~ hemisphere + group1
dds_filtered <- DESeq(dds_filtered, parallel = TRUE)

# extract the results for contrast of interest
# can play with lfcThreshold and alpha (both arbitrary)
group1_res <- data.table(results(dds_filtered,
                                 name = "group1_F_T_vs_F_S",
                                 tidy = TRUE,
                                 lfcThreshold = log(1.5, 2),
                                 alpha = 0.1))

group1_res[, al2fc := abs(log2FoldChange)]
setorder(group1_res, -al2fc)
group1_res[, al2fc := NULL]

# uncomment next line to take a peek at a gene
# plotCounts(dds_filtered, "Pmch", intgroup = c("injury", "sex"))


##########
# TEST 3 #
##########

# Do any genes correlate with injury, across sexes but specifically in
# the ipsilateral hemisphere?
design(dds_filtered) <- ~ sex + group2
dds_filtered <- DESeq(dds_filtered, parallel = TRUE)

resultsNames(dds_filtered)

# extract the results for contrast of interest
# can play with lfcThreshold and alpha (both arbitrary)
group2_res <- data.table(results(dds_filtered,
                                 contrast = c("group2",
                                              "I_T",
                                              "I_S"),
                                 tidy = TRUE,
                                 lfcThreshold = log(1.5, 2),
                                 alpha = 0.1))

group2_res[, al2fc := abs(log2FoldChange)]
setorder(group2_res, -al2fc)
group2_res[, al2fc := NULL]

# uncomment next line to take a peek at a gene
# plotCounts(dds_filtered, "Hcrt", intgroup = c("injury", "hemisphere"))


##########
# OUTPUT #
##########

fwrite(injury_res, inj_res_file)
fwrite(group1_res, group1_res_file)
fwrite(group2_res, group2_res_file)

sessionInfo()
