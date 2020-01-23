#!/bin/bash
#Requires samtools/1.5 & umi_tools 1.0.0
#Converts multiple SAM files to sorted BAM
#Performs umi_tools dedup on BAM files
#Converts deduped output to paired end fastq files
#Brett Youtsey
#DADA2 PIPELINE: script 3 after `align_multiple.R`


mkdir samFiles
mv *sam samFiles

mkdir bamFiles

#convert SAM to BAM
for sample in `ls samFiles`; do \
	samtools view -S -b samFiles/$sample > bamFiles/$sample.bam; done

#sort BAM files by left-most position
for sample in `ls bamFiles/*bam`; do \
	samtools sort $sample -o $sample.sorted.bam; done

#deduplicate sorted BAM files
for sample in `ls bamFiles/*sorted*`; do \
	samtools index $sample; umi_tools dedup --stdin $sample --stdout $sample.deduped.bam --paired --method unique; done

#sort deduped BAM files by read name
for sample in `ls bamFiles/*deduped*`; do \
	samtools sort $sample -n -o $sample.sorted.bam; done

#convert SAM to FASTQ
mkdir deduped_reads
for sample in `ls bamFiles/*deduped*sorted*bam`; do \
	output="$(basename $sample)"; samtools fastq -1 deduped_reads/$output.R1.fastq -2 deduped_reads/$output.R2.fastq -s junk.fq $sample; done

#remove all singletons
rm junk.fq

#properly name fastq files
rename .sam.bam.sorted.bam.deduped.bam.sorted.bam.R1.fastq _R1.fastq deduped_reads/*R1*
rename .sam.bam.sorted.bam.deduped.bam.sorted.bam.R2.fastq _R2.fastq deduped_reads/*R2*
