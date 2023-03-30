#to assess likelihood of NMD targeting by STAR-reported splice junctions

library("dplyr")
library("GenomicRanges")
library("GenomicFeatures")
library("rtracklayer")
library("BSgenome.Hsapiens.UCSC.hg38")

#functions ####
check_NMD <- function(transcript, exons, gene_gr) {
  #get GRanges object for transcript CDS exons
  trans_gr <- gene_cds_gr[gene_cds_gr$transcript_id == transcript]
  #get sequence of the (annotated) transcript
  trans_full_seq <- getSeq(Hsapiens, trans_gr)
  #remove skipped exons from transcript's set of exons
  used_exons <- c(1:length(trans_full_seq))
  used_exons <- used_exons[!used_exons %in% exons]
  trans_se_seq <- trans_full_seq[used_exons]
  #put together the exons to get the AS isoform sequence:
  trans_se_seq_uni <- ""
  for (i in 1:length(trans_se_seq)) {
    trans_se_seq_uni <- paste(trans_se_seq_uni, as.character(trans_se_seq[i]), sep = "")
  }
  trans_se_seq_uni <- DNAString(trans_se_seq_uni)
  trans_se_seq_aa <- translate(trans_se_seq_uni)
  #for now assume that the last CDS exon is the last exon
  #(otherwise it might be NMD-triggering anyways)
  #get aa position of first *
  #reverse translate that to nucleotide position
  ptc_pos <- (unlist(gregexpr("\\*", trans_se_seq_aa))[1])*3
  #use the list of exons to see if it's >50 nt upstream of a junction
  exon_sizes <- width(trans_se_seq)
  exon_starts <- 1
  for (i in 1:length(exon_sizes)-1) {
    exon_starts <- append(exon_starts, exon_starts[i] + exon_sizes[i])
  }
  if (any(ptc_pos < (exon_starts - 50))) {
    return("yes")
  } else {
    return("no")
  }
}


#read in data ####
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
ann_gtf <- import("/Applications/Genomics_applications/Genomes_and_transcriptomes/hg38_plus_Akata_inverted.gtf")
ann_gtf_gene <- ann_gtf[mcols(ann_gtf)$type == "gene"]

#get genes that contain detected SJs ####
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

#input going into this for loop:
     #transcripts: list of transcripts (in one gene that has been identified as containing an SJ; e.g. one gene from sj_pc)
     #gene_df: data frame of gtf rows for that gene
     #sj_pc[1]: a detected splice junction and the gene it's in
#eventually, put this all in a big for loop based on sj_pc
#and generate gene_df and transcripts from the info in each line of sj_pc

sj_NMD <- data.frame()
sj_no_NMD <- data.frame() #not doing anything with this yet. Worth keeping?

for (i in 1:length(sj_pc)) {
  print(i)
  print(start(sj_pc[i]))
  print(sj_pc[i]$gene_id)
  gene_gr <- ann_gtf[ann_gtf$gene_id == sj_pc[i]$gene_id]
  gene_cds_gr <- gene_gr[gene_gr$type == "CDS"]
  gene_df <- data.frame(gene_cds_gr)
  transcripts <- unique(gene_df$transcript_id)
}


for (i in 1:length(sj_pc)) {
  #need to get gene_df and transcripts list here
  print(i)
  print(sj_pc[i]$gene_id)
  gene_gr <- ann_gtf[ann_gtf$gene_id == sj_pc[i]$gene_id]
  gene_cds_gr <- gene_gr[gene_gr$type == "CDS"]
  gene_df <- data.frame(gene_cds_gr)
  transcripts <- unique(gene_df$transcript_id)
  for (transcript in transcripts) { #go through each transcript
    if (is.na(transcript)) { #ignore gene lines (NA in transcript field). Though actually shouldn't be any
      next
    } else {
      print(transcript)
      transcript_df <- gene_df %>%
        filter(transcript_id == transcript) %>%
        filter(type == "CDS")
      if (nrow(transcript_df) == 0) {
        print(paste(transcript, "has no CDS"))
        next
      } else { #if the transcript is coding, look further
        relevant <- data.frame() #set up a df for relevant exons
        for (j in 1:nrow(transcript_df)) {
          if (nrow(relevant) == 0) { #look for matching sj donor
            if (start(sj_pc[i])-1 == transcript_df[j,3]) {
              relevant <- transcript_df[j,] #SJ donor is a relevant exon. Add it.
            } else {
              next
            }
          } else { #if there is content in "relevant" (i.e. we've found a matching donor)
            if (end(sj_pc[i])+1 == transcript_df[j,2]) { #look for matching acceptor
              relevant <- rbind(relevant, transcript_df[j,])
              break #if matching acceptor is found, we have all the exons we need
            } else {
              relevant <- rbind(relevant, transcript_df[j,])
            }
          }
        }
      } #end of transcript
      if (nrow(relevant) == 0) {
        print("splice junction either fully annotated or splice sites unannotated: assume no NMD") #ignore for now
      }  
      else if (end(sj_pc[i])+1 != relevant[nrow(relevant),2]) { 
        print("splice sites not annotated in this transcript") #ignore for now
      }
      else {
        if (nrow(relevant) == 2) {
          print("annotated SJ: no new SE") #assume no NMD
        } else {
          se <- relevant[2:(nrow(relevant) - 1), ]
          se$exon_length <- se$end - se$start + 1 #could use width column here
          if (sum(se$exon_length) %% 3 == 0) { #in-frame: assume no NMD
            print("new in frame SE")
          } else {
            print(paste("new out-of-frame SE:", sum(se$exon_length), "nt"))
            NMD <- check_NMD(transcript, se$exon_number, gene_cds_gr)
            if (NMD == "yes") {
              print("NMD!") #should break here for efficiency: go to next SJ
              if (nrow(sj_NMD) > 0) {
                print("adding to df of NMD SJs")
                print(paste("which already contains", nrow(sj_NMD), "SJ(s)") )
                sj_NMD <- rbind(sj_NMD, data.frame(sj_pc[i]))
                print(paste("now it contains", nrow(sj_NMD), "SJs") )
                break
                }
              else {
                print("starting df of NMD SJs")
                sj_NMD <- data.frame(sj_pc[i])
                break
                }
            }
            else {print("no NMD!")}
          }
        }
      }
    }
  }
}

for (transcript in transcripts) { #go through each transcript
  if (is.na(transcript)) { #ignore gene lines (NA in transcript field). Though actually shouldn't be any
    next
  } else {
    print(transcript)
    transcript_df <- gene_df %>%
      filter(transcript_id == transcript) %>%
      filter(type == "CDS")
    if (nrow(transcript_df) == 0) {
      print(paste(transcript, "has no CDS"))
      next
    } else { #if the transcript is coding, look further
      relevant <- data.frame() #set up a df for relevant exons
      for (i in 1:nrow(transcript_df)) {
        if (nrow(relevant) == 0) { #look for matching sj donor
          if (start(sj_pc[1])-1 == transcript_df[i,3]) {
            relevant <- transcript_df[i,] #SJ donor is a relevant exon. Add it.
          } else {
            next
          }
        } else { #if there is content in "relevant" (i.e. we've found a matching donor)
          if (end(sj_pc[1])+1 == transcript_df[i,2]) { #look for matching acceptor
            relevant <- rbind(relevant, transcript_df[i,])
            break #if matching acceptor is found, we have all the exons we need
          } else {
            relevant <- rbind(relevant, transcript_df[i,])
          }
        }
      }
    } #end of transcript
    if (nrow(relevant) == 0 | (end(sj_pc[1])+1 != relevant[nrow(relevant),2])) {
      print("splice sites not annotated in this transcript") #ignore for now
    }  else {
      if (nrow(relevant) == 2) {
        print("annotated SJ: no new SE") #assume no NMD
      } else {
        se <- relevant[2:(nrow(relevant) - 1), ]
        se$exon_length <- se$end - se$start + 1 #could use width column here
        if (sum(se$exon_length) %% 3 == 0) { #in-frame: assume no NMD
          print("new in frame SE")
        } else {
          print(paste("new out-of-frame SE:", sum(se$exon_length), "nt"))
          NMD <- check_NMD(transcript, se$exon_number, gene_cds_gr)
          if (NMD == "yes") {
            print("NMD!")
            if (nrow(sj_NMD) > 0) {rbind(sj_NMD, data.frame(sj_pc[1]))}
            else {sj_NMD <- data.frame(sj_pc[1])}
            }
          else {print("no NMD!")}
        }
      }
    }
  }
}


#code contributing to check_NMD function:
#to delete once function is better tested
#make a function given transcript name, exon number(s) and gene_cds_gr
#transcript = "ENST00000482437"
transcript = "ENST00000370143"
ses = se$exon_number

#get GRanges object for transcript CDS exons
trans_gr <- gene_cds_gr[gene_cds_gr$transcript_id == transcript]
#get sequence of the (annotated) transcript
trans_full_seq <- getSeq(Hsapiens, trans_gr)

#remove skipped exons from transcript's set of exons
used_exons <- c(1:length(trans_full_seq))
used_exons <- used_exons[!used_exons %in% ses ]
trans_se_seq <- trans_full_seq[used_exons]

#put together the exons to get the AS isoform sequence:
trans_se_seq_uni <- ""
for (i in 1:length(trans_se_seq)) {
  trans_se_seq_uni <- paste(trans_se_seq_uni, as.character(trans_se_seq[i]), sep = "")
}
trans_se_seq_uni <- DNAString(trans_se_seq_uni)
trans_se_seq_aa <- translate(trans_se_seq_uni)
#warning that last two bases were ignored. Out-of-frame! Hooray!
countPattern("*", trans_se_seq_aa)
#17! Good!

#for now assume that the last CDS exon is the last exon
#(otherwise it might be NMD-triggering anyways)

#get aa position of first *
#reverse translate that to nucleotide position
ptc_pos <- (unlist(gregexpr("\\*", trans_se_seq_aa))[1])*3
# nucleotide 330 in this example

#use the list of exons to see if it's >50 nt upstream of a junction
exon_sizes <- width(trans_se_seq)
exon_starts <- 1
for (i in 1:length(exon_sizes)-1) {
  exon_starts <- append(exon_starts, exon_starts[i] + exon_sizes[i])
}
print(exon_starts)

if (any(ptc_pos < (exon_starts - 50))) {
  print("NMD!")
} else {
  print("No NMD!")
}

#to assess things at the end, will need to somehow compare depth (not just number I think)
#of NMD-targeting splice junctions to other splice junctions
#because everything is going to have some of both

