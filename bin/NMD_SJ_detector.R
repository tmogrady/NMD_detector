#to assess likelihood of NMD targeting by STAR-reported splice junctions

library("dplyr")

#read in STAR SJ.out.tab file
sj <- read.table("input_data/MC1_truncated_test-SJ.out.tab")

#name columns for my own sanity
colnames(sj) <- c("chr","int_start","int_end","strand","motif","ann","unique","multi","overhang")

#first, get annotated genes (transcripts?) for each SJ
#check if coding
#check if SJ coordinates here are annotated (check donor & acceptor separately)
#if so, look at exons between and divide by 3 to check for frameshift
    #how to deal with multiple isoforms?
#if no frameshift, assume no NMD
#if frameshift, check for PTC
#start with just junctions in annotated genes, with annotated donors & acceptors

