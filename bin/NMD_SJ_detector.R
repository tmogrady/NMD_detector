#to assess likelihood of NMD targeting by STAR-reported splice junctions

library("dplyr")
library("GenomicRanges")

#read in STAR SJ.out.tab file
sj <- read.table("input_data/MC1_truncated_test-SJ.out.tab")

#name columns for GRanges (and my own sanity)
colnames(sj) <- c("chr","start","end","strand","motif","ann","unique","multi","overhang")
#change strand info to +/- for GRanges
sj <- sj %>%
  mutate(strand = ifelse(strand == 1, "+",
                         ifelse(strand == 2, "-",
                                strand)))

#first, get annotated genes (transcripts?) for each SJ
sj_gr <- makeGRangesFromDataFrame(sj, keep.extra.columns = TRUE)


#check if coding
#check if SJ coordinates here are annotated (check donor & acceptor separately)
#if so, check if there are annotated exons between
#if so, look at exon(s) and divide by 3 to check for frameshift
    #how to deal with multiple isoforms?
#if no frameshift, assume no NMD
#if frameshift, check for PTC
#start with just junctions in annotated genes, with annotated donors & acceptors

