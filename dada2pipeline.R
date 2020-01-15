#!/usr/bin/env Rscript

# McGuire Lab, University of Oregon. Created by Ben Cosgrove, Brett Youtsey, Josh Sakai, & Sam Velazquez
# Adapted from tutorials created by Devin Dinwiddie & DADA2 tutorial 1.12 (https://benjjneb.github.io/dada2/tutorial.html)
# Compatible with R version 3.5.1

exit <- function() {
# Ends the script. Will be used in the event of an error.
  .Internal(.invokeRestart(list(NULL, NULL), NULL))
}




########################COMMAND LINE INPUTS
library(optparse)
#Saving inputs as variables with optparse
option_list <- list(
  make_option(c("-r", "--reads"), help = "Path to directory containing demultiplexed reads (_R1.fastq.gz & _R2.fastq.gz). MUST BE COMPRESSED (gz)"),
  make_option(c("-o", "--outputDirectory"), help = "Will create directory where outputs are stored"),
  make_option(c("-U", "--taxonRef"), help = "Taxon reference in fasta format with taxanomic strata in the header"),
  make_option(c("-m", "--metadata"), help = "Metadata csv containing same sample names as input reads in the first column"),
  make_option(c("-f", "--formula"), help = "Formula for experimental design in deseq2 (used in variance transformation) FORMAT: ~var1+var2+varx (Variables must exactly match colnames in the metadata file)")
)
parser = OptionParser(usage = "This script filters reads and runs DADA2 for assessment of amplicon sequencing data across multiple samples. The resulting ASV tables are then paired with metadata and variance transformed with deseq2.  \n\n\t OUTPUTS \n\t ASV_table.tsv       : ASV+taxonomy table DADA2 outputs \n\t ASV_summary.tsv   : Sum counts of all ASV's per sample \n\t track_reads.csv     : Table counting the number of reads kept at each step of the pipeline \n\t ASV-VT.tsv          : A variance transformed version of ASV_table.tsv \n\t phyloseq_object.RDS : phyloseq object of the variance transformed ASV with metadata \n\n\t *Will use all available threads",
                      option_list = option_list)

inputs <- parse_args(parser)

#Save file path to reads
path <- inputs$reads
#End of the path must have "/". The next line will add one just in case. Ending with "//" is the same as "/"
path <- paste(path,"/",sep = "")

#Save name of output Directory
outDir <- inputs$outputDirectory
outDir <- paste(outDir,"/",sep = "")

#Save path to taxonomy reference:
taxonReference <- inputs$taxonRef

#Save name of metadata
metaPath <- inputs$metadata

#Save formula as a formula object
design <- as.formula(inputs$formula)

#make new output directory in current directory
system(paste("mkdir", outDir, sep = " "))


#check to see if all inputs are in
if (!(all(c("reads", "outputDirectory", "taxonRef") %in% names(inputs)))) {
  print("You do not have all the required inputs. Use --help for more info")
  exit()
}




########################LOAD REQUIRED R PACKAGES
if (!"BiocManager"%in% installed.packages()){
  install.packages("BiocManager") # dada2 requires BioManager to install for R 3.5 +
}

#load dada2
if (!"dada2"%in% installed.packages()){
  BiocManager::install("dada2", version = "3.8")
}
library(dada2)

#load Biostrings
if(!"Biostrings" %in% installed.packages()){
  BiocManager::install("Biostrings")
}
library(Biostrings)

#load ShortRead
if (!"ShortRead"%in% installed.packages()){
  BiocManager::install("ShortRead")
}
library(ShortRead)

if(!"dplyr" %in% installed.packages()){
  install.packages("dplyr")
}
library(dplyr)


if (!"phyloseq"%in% installed.packages()){
  install.packages("phyloseq")
}
library(phyloseq)

if (!"DESeq2"%in% installed.packages()){
  install.packages("DESeq2")
}
library(DESeq2)

if (!"vegan"%in% installed.packages()){
  install.packages("vegan")
}
library(vegan)

if (!"rowr"%in% installed.packages()){
  install.packages("rowr")
}
library(rowr)




######################## Filter Reads

#sort the names of the forward and reverse reads so they are in the same order
#order R1
name.fS <- sort(list.files(path, pattern = "_R1.fastq.gz", full.names = TRUE)) #CHANGE "_R1.fastq" if your files do not have that ending pattern
#order R2
name.rS <- sort(list.files(path, pattern = "_R2.fastq.gz", full.names = TRUE)) #CHANGE "_R1.fastq" if your files do not have that ending pattern

#check if R can read files
if (length(name.fS) == 0){
  print("DADA2PIPELINE ERROR: filenames for input reads do not follow correct pattern. Refer to the description --reads")
  exit()
}

#create path for
names.fS.filtN = file.path(outDir, "filtN", basename(name.fS)) # path to put N filtered reads path+filtN --forward reads
names.rS.filtN = file.path(outDir, "filtN", basename(name.rS)) # path to put N filtered reads path+filtN --Reverse reads

#filter N sequences from file
filterAndTrim(name.fS, names.fS.filtN, name.rS, names.rS.filtN, maxN = 0, multithread = TRUE,compress=T)

#in case samples got tossed out during this step lets reset file names
pathFiltN=paste(outDir, "filtN/", sep="")

name.fS <- sort(list.files(pathFiltN, pattern = "_R1.fastq.gz", full.names = TRUE))
name.rS <- sort(list.files(pathFiltN, pattern = "_R2.fastq.gz", full.names = TRUE))
names.fS.filtN = file.path(outDir, "filtN", basename(name.fS))
names.rS.filtN = file.path(outDir, "filtN", basename(name.rS))


cutFs <- sort(list.files(pathFiltN, pattern = "_R1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(pathFiltN, pattern = "_R2.fastq.gz", full.names = TRUE))

# Extract sample names, assuming filenames have format: McG.L2NLitter_R2.fastq
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))


filtFWD <- file.path(pathFiltN, "filtered", basename(cutFs))
filtREV <- file.path(pathFiltN, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFWD, cutRs, filtREV, maxN = 0, maxEE = c(2, 2),
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = F, multithread = TRUE)  #filter the reads

filtPath=file.path(pathFiltN,"filtered")

filtFs <- sort(list.files(filtPath, pattern = "_R1.fastq.gz", full.names = TRUE))
filtRs <- sort(list.files(filtPath, pattern = "_R2.fastq.gz", full.names = TRUE))

filtFWD <- file.path(filtPath, basename(filtFs))
filtREV <- file.path(filtPath, basename(filtRs))
sample.names <- unname(sapply(filtFs, get.sample.name))




######################## Run DADA2
set.seed(111)
errF <- learnErrors(filtFWD, multithread=TRUE)
errR <- learnErrors(filtREV, multithread=TRUE)

derepFWD <- derepFastq(filtFWD, verbose=TRUE)
derepREV <- derepFastq(filtREV, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFWD) <- sample.names
names(derepREV) <- sample.names

dadaFWD <- dada(derepFWD, err=errF, multithread=TRUE)
dadaREV <- dada(derepREV, err=errR, multithread=TRUE)

merge <- mergePairs(dadaFWD, derepFWD, dadaREV, derepREV, verbose=TRUE)

ASVtab <- makeSequenceTable(merge)

# remove chimeric sequences
ASV.nochim <- removeBimeraDenovo(ASVtab, method="consensus", multithread=TRUE, verbose=TRUE)




######################## Track Reads
getN <- function(x) sum(getUniques(x))
track <- cbind.fill(out, sapply(dadaFWD, getN), sapply(dadaREV, getN), sapply(merge,
                                                                         getN), rowSums(ASV.nochim), fill = 0)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged",
                     "nonchim")
rownames(track) <- sample.names
trackOut <- paste(outDir, "track_reads.csv")
write.csv(track, trackOut)




######################## Assign Taxonomy & Save ASV+taxonomy as tsv
taxa <- assignTaxonomy(ASV.nochim, taxonReference, multithread = TRUE, tryRC = TRUE) #tryRC=try revese compliment

#clear all non-essential variables to recover memory
rm(list= setdiff(ls(), c("ASV.nochim", "taxa", "outDir", "metaPath", "design")))

ASV.nochim <- t(ASV.nochim)
tab=merge(taxa,ASV.nochim,by = "row.names",sort=F)
tab=subset(tab,Kingdom != "k__unidentified") # remove unidentified

combinedOut <- paste(outDir, "ASV_table.tsv", sep = "")
write.table(tab, combinedOut, sep='\t', row.names=F, quote=FALSE) # combined table

sumo=data.frame(X=colSums(tab[9:ncol(tab)])) # get sample sums for rareification
sum1=sumo[order(sumo$X),,drop=F] # order low to high

summaryOut <- paste(outDir, "ASV_summary.tsv", sep = "")
write.table(sum1, summaryOut, sep="\t") #change file name --sample sum table




######################## Variance Stabilization
rm(list= setdiff(ls(), c("tab", "outDir", "metaPath", "design")))
otu <- tab
#change name of first column containing sequences to "X"
colnames(otu)[1] <- "X"

meta=read.csv(metaPath) # read in you metadata file
colnames(meta)[1] <- "SampleID"
meta$SampleID=make.names(meta$SampleID, unique=TRUE)

# `-` in sample names need to be replaced with `_` to match metadata
colnames(otu)=gsub("[-]","_",colnames(otu))

#detect taxa columns:
for(colNumber in seq(1,dim(otu)[2])){
  if(grepl("[a-z]__",otu[1,colNumber])){
    lastTaxaColumn <- colNumber
  }}

otuMap=otu[,-c(1:lastTaxaColumn)]

otuMap1=otuMap %>% select(one_of(as.character(meta$SampleID)))

samp=data.frame(SampleID=colnames(otuMap1))
meta2=merge(samp,meta,by='SampleID')

#make phyloseq object

taxa=as.matrix(otu[,2:lastTaxaColumn]) # assign taxonomy data as a matrix make sure you are selecting the correct columns [2:7 in this example]

OTU = otu_table(otuMap1, taxa_are_rows = TRUE)  # create a phyloseq otu table object
TAX = tax_table(taxa) # create a phyloseq taxa table object

physeq = phyloseq(OTU,TAX) #merge the otu and the taxa into a phyloseq object

sampledata = sample_data(data.frame(meta2,row.names=sample_names(physeq),stringsAsFactors=FALSE)) # make phyloseq sample_data object.
#make sure the row names of this object match the sample _ID

#combine phylo objects --add sample_data
physeq1 = phyloseq(OTU,TAX, sampledata)

#some samples may have zero counts, which will cause errors in estimateDispersions. This line will remove the samples
physeq1 <- prune_samples(sample_sums(physeq1) > 0, physeq1)

#first convert your phyloseq object to a deseq2 object
diagdds = phyloseq_to_deseq2(physeq1, design) # change variable (~tree_species) to variable of interest.
#for some reason it didn't change the values when I did this with different variables. I don't think it matters.
#this estimates size factors based on the mean of the sample
diagdds = estimateSizeFactors(diagdds,type="poscounts") # need to estimate the size factors *iterate acounts for 0's in the data

diagdds = estimateDispersions(diagdds) # estimate the dispersion parameters
diagvst = getVarianceStabilizedData(diagdds) # do variance stabilization

# the new variance table can be written to a new file.
# all otu with count 0 will now have a negative number associated with them so you need to conver them back to 0

taxa1=data.frame(tax_table(physeq1)) # pull out the taxonomy data
otu_table(physeq1) <- otu_table(diagvst, taxa_are_rows = TRUE) # place variance stabilized table as otu table in phyloseq object
ot=data.frame(otu_table(physeq1)) # pull out OTU table from phyloseq object
ot[ot < 0.0] <- 0.0 #convert negative numnbers to 0
deseq=cbind(taxa1,ot) #combine taxa and ASV table
deseq.X=cbind(otu$X,deseq)

# save variance transformed ASV
asvOut <- paste(outDir, "ASV-VT.tsv", sep = "")
write.table(deseq.X, asvOut, sep='\t', row.names=F, quote=FALSE) #write out the new table

# Repeat steps to make phyloseq object this time for the variance transformed ASV
otuMap=deseq.X[,-c(1:lastTaxaColumn)]
otuMap1=otuMap %>% select(one_of(as.character(meta$SampleID)))
taxa=as.matrix(deseq.X[,2:lastTaxaColumn])
OTU = otu_table(otuMap1, taxa_are_rows = TRUE)  # create a phyloseq otu table object
TAX = tax_table(taxa) # create a phyloseq taxa table object

physeq = phyloseq(OTU,TAX) #merge the otu and the taxa into a phyloseq object

sampledata = sample_data(data.frame(meta2,row.names=sample_names(physeq),stringsAsFactors=FALSE)) # make phyloseq sample_data object.
#make sure the row names of this object match the sample _ID

# combine phylo objects --add sample_data
physeq1 = phyloseq(OTU,TAX, sampledata)

#save phyloseq object
phyOut <- paste(outDir, "phyloseq_object.R", sep = "")
saveRDS(physeq1, phyOut)
