#to assess likelihood of NMD targeting by STAR-reported splice junctions

library("dplyr")
library("GenomicRanges")
library("rtracklayer")
#library("ensembldb")
#library("EnsDb.Hsapiens.v86") #most recent version available

#read in STAR SJ.out.tab file
sj <- read.table("input_data/MC1_truncated_test-SJ.out.tab")

#name columns for GRanges (and my own sanity)
colnames(sj) <- c("chr","start","end","strand","motif","ann","unique","multi","overhang")
#change strand info to +/- for GRanges
sj <- sj %>%
  mutate(strand = ifelse(strand == 1, "+",
                         ifelse(strand == 2, "-",
                                strand)))
#make SJ data a GRanges object
sj_gr <- makeGRangesFromDataFrame(sj, keep.extra.columns = TRUE)

#get annotation (gtf)
#try gff
ann_gtf <- import("/Applications/Genomics_applications/Genomes_and_transcriptomes/hg38_plus_Akata_inverted.gtf")
ann_gtf_gene <- ann_gtf[mcols(ann_gtf)$type == "gene"]
  
ovrlp <- findOverlaps(sj_gr, ann_gtf_gene, type = "within", select = "first")
#later consider how to handle SJs within multiple genes (not common)
sj_ann <- sj_gr
sj_ann$gene_id <- ann_gtf_gene$gene_id[ovrlp]
sj_ann$gene_biotype <- ann_gtf$gene_biotype[ovrlp]
#here filter to only protein-coding (maybe count and report others later)
#first remove SJs not in annotated genes
#because they don't have ORF information, and to avoid later error with NAs
sj_ann <- sj_ann[!is.na(sj_ann$gene_id)]
sj_pc <- sj_ann[sj_ann$gene_biotype == "protein_coding"]

#next, for each SJ & gene combo,
#run through the possible transcripts and check for NMD triggering
#may need to subset ann_gtf by gene's location, 
#or otherwise work around metadata NA errors
#this is the first (annotated, protein-coding) example SJ's gene:
gene_gr <- ann_gtf[ann_gtf$gene_id == "ENSG00000122435"]
#move to dataframes because they're easier and avoid GRanges metadata NA problem
gene_df <- data.frame(gene_gr)
transcripts <- unique(gene_df$transcript_id)
transcript_df <- gene_df %>%
  filter(transcript_id == "ENST00000370141") %>%
  filter(type == "CDS")

for (transcript in transcripts) {
  if (is.na(transcript)) {
    next
  }
  else {
    print(transcript)
    transcript_df <- gene_df %>%
      filter(transcript_id == transcript) %>%
      filter(type == "CDS")
    if (nrow(transcript_df) == 0) {
      print(paste(transcript, "has no CDS"))
      next
    }
    else {
      relevant <- data.frame()
      for (i in 1:nrow(transcript_df)) {
        if (nrow(relevant) == 0) {
          if (start(sj_pc[1])-1 == transcript_df[i,3]) {
            relevant <- transcript_df[i,]
          }
          else {
            next
          }
        }
        else { #if there is content in "relevant"
          #print(i)
          if (end(sj_pc[1])+1 == transcript_df[i,2]) {
            relevant <- rbind(relevant, transcript_df[i,])
            break
          }
          else {
            relevant <- rbind(relevant, transcript_df[i,])
          }
        }
      }
    }

    if (nrow(relevant) == 2) {
      print("no SE")
    }
    else {
      se <- relevant[2:(nrow(relevant) - 1), ]
      se$exon_length <- se$end - se$start + 1
      if (sum(se$exon_length) %% 3 == 0) {
        print("in frame SE")
      }
      else {
        print(paste("out-of-frame SE:", sum(se$exon_length), "nt"))
      }
    }
  }
}


ann_gtf[ann_gtf$transcript_id == "ENST00000370143"]
#Error: logical subscript contains NAs
#could possibly get around this error by subsetting based on genomic coordinates of gene

#get (only) the CDS exons
#if it matches the end of one, check the beginning of the next
#if it matches, exit: it's an annotated splice junction that doesn't skip an exon in that isoform
#if it doesn't match, store that exon and check the next one
#continue till it matches
#then take the skipped exon(s) and divide by 3
#if divisible, it's an in-frame skip

#to assess things at the end, will need to somehow compare depth (not just number I think)
#of NMD-targeting splice junctions to other splice junctions
#becuase eveerything is going to have some of both

#at end of transcripts for that gene, put SJ into NMD list or not-NMD list


#previous attempts with bed file:
ann <- import("/Applications/Genomics_applications/Genomes_and_transcriptomes/hg38_plus_Akata_inverted.bed.converted.bed.reconfigured.bed")
ovrlp <- findOverlaps(sj_gr, ann, type = "within", select = "first")
#will need to deal with multiple isoforms
sj_ann <- sj_gr
sj_ann$name <- ann$name[ovrlp]
#Error: subscript contains NAs when I try to bring over thick and blocks
#not a problem with name (or score or itemRgb, though all have at least 1 NA)
sj_ann$thick <- ann$thick[ovrlp]
sj_ann$blocks <- ann$blocks[ovrlp]

#check if SJ coordinates here are annotated (check donor & acceptor separately)
#if so, check if there are annotated exons between
#if so, look at exon(s) and divide by 3 to check for frameshift
    #how to deal with multiple isoforms?
#if no frameshift, assume no NMD
#if frameshift, check for PTC
#start with just junctions in annotated genes, with annotated donors & acceptors

