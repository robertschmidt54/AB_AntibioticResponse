#!usr/bin/bash


# Align files using Bowtie2, compress output, remove sam, sort, and index bam files.

for fq in $( ls data/* )
do
	ID=$( basename ${fq%%.fastq} )
	bowtie2 -x RefGenome/Index/AB17978_ -U ${fq} -S BowtieOut/${ID}.sam
       #samtools view -b BowtieOut/${ID}.sam > BowtieOut/${ID}.bam 
       samtools sort BowtieOut/${ID}.sam > BowtieOut/${ID}.bam  
       samtools index -b BowtieOut/${ID}.bam
       rm BowtieOut/${ID}.sam
done

#Run HTSeq to count reads in gene.
