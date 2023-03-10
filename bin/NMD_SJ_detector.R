#to assess likelihood of NMD targeting by STAR-reported splice junctions

library("dplyr")
library("GenomicRanges")
library("GenomicFeatures")
library("rtracklayer")
library("BSgenome.Hsapiens.UCSC.hg38")
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
#this is the first (annotated, protein-coding) example SJ's gene:
gene_gr <- ann_gtf[ann_gtf$gene_id == "ENSG00000122435"]
#get only CDS. That's what we need
#and this will avoid NA errors later when subsetting on transcript
gene_cds_gr <- gene_gr[gene_gr$type == "CDS"]
#move to dataframes because they're easier and avoid GRanges metadata NA problem
gene_df <- data.frame(gene_cds_gr)
transcripts <- unique(gene_df$transcript_id)
transcript_df <- gene_df %>%
  filter(transcript_id == "ENST00000370141") %>%
  filter(type == "CDS")

for (transcript in transcripts) { #go through each transcript
  if (is.na(transcript)) { #ignore gene lines (NA in transcript field)
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
    else { #if the transcript is coding, look further
      relevant <- data.frame() #set up a df for relevant exons
      for (i in 1:nrow(transcript_df)) {
        if (nrow(relevant) == 0) { #look for matching sj donor
          if (start(sj_pc[1])-1 == transcript_df[i,3]) {
            relevant <- transcript_df[i,] #SJ donor is a relevant exon. Add it.
          }
          else {
            next
          }
        }
        else { #if there is content in "relevant" (i.e. we've found a matching donor)
          if (end(sj_pc[1])+1 == transcript_df[i,2]) { #look for matching acceptor
            relevant <- rbind(relevant, transcript_df[i,])
            break #if matching acceptor is found, we have all the exons we need
          }
          else {
            relevant <- rbind(relevant, transcript_df[i,])
          }
        }
      }
    } #end of transcript
    if (nrow(relevant) == 0 | (end(sj_pc[1])+1 != relevant[nrow(relevant),2])) {
      print("splice sites not annotated in this transcript") #ignore for now
    }  
    else {
      if (nrow(relevant) == 2) {
        print("annotated SJ: no new SE") #assume no NMD
      }
      else {
        se <- relevant[2:(nrow(relevant) - 1), ]
        se$exon_length <- se$end - se$start + 1 #could use width column here
        if (sum(se$exon_length) %% 3 == 0) { #in-frame: assume no NMD
          print("new in frame SE")
        }
        else {
          print(paste("new out-of-frame SE:", sum(se$exon_length), "nt"))
          #get sequence of this transcript
          #remove SE
          #and look for stop codons
        }
      }
    }
  }
}


trans_gr <- gene_cds_gr[gene_cds_gr$transcript_id == "ENST00000370141"]

#make a df to figure out what I'm doing:
trans_df <- data.frame(trans_gr)

test_seq <- getSeq(Hsapiens, trans_gr)
#worked! But is a list of exons rather than a complete sequence
#warning message but ok
#note that stop codon is not included

#combine exons in one sequence:
#(got to be a better way to do this :/)
test_seq_trans <- ""
for (i in 1:length(test_seq)) {
  test_seq_trans <- paste(test_seq_trans, as.character(test_seq[i]), sep = "")
}
test_seq_trans <- DNAString(test_seq_trans)
test_whole_trans_aa <- translate(test_seq_trans)
countPattern("*", test_whole_trans_aa)
#0. Good.

#try to remove an exon:
#for transcript ENST00000370141 it's exon 5
test_seq_se <- test_seq[c(1:4,6:length(test_seq))]
test_seq_se_trans <- ""
for (i in 1:length(test_seq_se)) {
  test_seq_se_trans <- paste(test_seq_se_trans, as.character(test_seq_se[i]), sep = "")
}
test_seq_se_trans <- DNAString(test_seq_se_trans)
test_whole_se_trans_aa <- translate(test_seq_se_trans)
#warning that last two bases were ignored. Out-of-frame! Hooray!
countPattern("*", test_whole_se_trans_aa)
#17! Good!

trans_gr <- ann_gtf[ann_gtf$transcript_id == "ENST00000370143"]
#Error: logical subscript contains NAs

#to assess things at the end, will need to somehow compare depth (not just number I think)
#of NMD-targeting splice junctions to other splice junctions
#because everything is going to have some of both

