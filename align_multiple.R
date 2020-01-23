#!/usr/bin/env Rscript
#compatible with R version 3.5.1 & bowtie2/2.3.3
#align multiple fastq files to taxanomic fasta reference
#Brett Youtsey
# DADA2 PIPELINE: script 2 after `orient.py`

library(optparse)
#Saving inputs as variables with optparse
option_list <- list(
  make_option(c("-U", "--taxonRef"), help = "Taxonomic reference in fasta format with taxonomic strata in the header"),
  make_option(c("-r", "--reads"), help = "Path to directory containing demultiplexed reads (_R1.fastq & _R2.fastq).[Required]"),
  make_option(c("-t", "--threads"), default = 1, type = "integer", help = "Number of threads (default 1). Must be less than or equal to the number of cores assigned on SLURM")

)
parser = OptionParser(usage = "Builds bowtie2 database & aligns directory of paired end fastq files with a specified bowtie2 index. R version 3.5.1 & bowtie2/2.3.3 MUST BE LOADED",
                      option_list = option_list)

inputs <- parse_args(parser)

taxonRef <- inputs$taxonRef
readsPath <- inputs$reads
numThreads <- inputs$threads
index <- inputs$indexName


forwardReads <- sort(list.files(readsPath, pattern = "_R1.fastq", full.names = TRUE))
reverseReads <- sort(list.files(readsPath, pattern = "_R2.fastq", full.names = TRUE))

system2("bowtie2", args = "--version")

system2("bowtie2-build", args = c("-f", taxonRef, "--threads", numThreads, "align_multiple_database"))

for(fileIndex in seq_along(forwardReads)){

  #remove extension and file path from sample name. This is so the output SAM will be {sampleName}.sam
  outputWPath <- unlist(strsplit(forwardReads[fileIndex], "_R"))[1]
  outputName <- gsub(readsPath, "", outputWPath)
  outputName <- paste(gsub("/", "", outputName), ".sam", sep = "")

  cat("Aligning sample: ", outputName, "\n")
  system2("bowtie2", args = c("-x", "align_multiple_database", "-1", forwardReads[fileIndex], "-2", reverseReads[fileIndex], "-p", numThreads, "-S", outputName))
}
