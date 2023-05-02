#to assess likelihood of NMD targeting by STAR-reported splice junctions

library("dplyr")
library("GenomicRanges")
library("GenomicFeatures")
library("rtracklayer")
library("BSgenome.Hsapiens.UCSC.hg38")
library("ggplot2")
library("foreach")
library("doParallel")

#functions ####
check_NMD <- function(transcript, exons, gene_gr) {
  #get GRanges object for transcript CDS exons
  trans_gr <- gene_cds_gr[gene_cds_gr$transcript_id == transcript]
  #get sequence of the (annotated) transcript
  trans_full_seq <- suppressWarnings(getSeq(Hsapiens, trans_gr))
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
  trans_se_seq_aa <- suppressWarnings(translate(trans_se_seq_uni))
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

#set parameters ####
#number of unique reads required to consider junctions
#later maybe normalize to read depth (sj reads per million?)
sj_thres <- 0

#read in data ####
#read in STAR SJ.out.tab file
sj <- read.table("input_data/MC1_truncated_test-SJ.out.tab")
sj <- read.table("../temp/MC1_S34_L004_R1_001_MC1_S34_L004_R2_001.hg38plusAkata_inverted-SJ.out.tab")

#name columns for GRanges (and my own sanity)
colnames(sj) <- c("chr","start","end","strand","motif","ann","unique","multi","overhang")
#change strand info to +/- for GRanges
sj <- sj %>%
  mutate(strand = ifelse(strand == 1, "+",
                         ifelse(strand == 2, "-",
                                "*")))
sj_filt <- sj %>%
  filter(unique >= sj_thres)

#make SJ data a GRanges object
sj_gr <- makeGRangesFromDataFrame(sj_filt, keep.extra.columns = TRUE)

#get annotation (gtf)
#should get this from a package instead of a gtf file (or have option)
ann_gtf <- import("/Applications/Genomics_applications/Genomes_and_transcriptomes/hg38_plus_Akata_inverted.gtf")
ann_gtf_gene <- ann_gtf[mcols(ann_gtf)$type == "gene"]

#get genes that contain detected SJs ####
ovrlp <- findOverlaps(sj_gr, ann_gtf_gene, type = "within", select = "first")
#gives index in ann_gtf_gene for each SJ
#has NAs for SJs in unannotated regions
#later consider how to handle SJs within multiple genes (not common)

#get the annotated SJs:
sj_ann <- sj_gr
sj_ann$gene_id <- ann_gtf_gene$gene_id[ovrlp]
sj_ann$gene_name <- ann_gtf_gene$gene_name[ovrlp]
sj_ann$gene_biotype <- ann_gtf$gene_biotype[ovrlp]
#here filter to only protein-coding
#first remove SJs not in annotated genes
#because they don't have ORF information, and to avoid later error with NAs
sj_ann <- sj_ann[!is.na(sj_ann$gene_id)]
sj_pc <- sj_ann[sj_ann$gene_biotype == "protein_coding"]

#get intergenic SJs to report:
intrgn <- which(is.na(ovrlp))
sj_notann <- sj_gr[intrgn]
sj_notann_df <- data.frame(sj_notann) #maybe not necessary

#get SJs not in coding genes to report:
sj_notPC <- sj_ann[sj_ann$gene_biotype != "protein_coding"]
sj_notPC_df <- data.frame(sj_notPC) #maybe not necessary to df

#test each SJ for NMD ####
#first set up parallelization
#trivial check:
n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
some_numbers <- sample(1:1000, 100, replace=TRUE)
classified_numbers <- foreach(
  i = 1:length(some_numbers), 
  .combine = 'rbind'
) %dopar% {
  if (some_numbers[i] < 500) {
    squareroot <- sqrt(some_numbers[i])
    data.frame("number" = squareroot, "type" = "lower")
  }
  else {
    squareroot <- sqrt(some_numbers[i])
    data.frame("number" = squareroot, "type" = "higher")
  }
}
parallel::stopCluster(cl = my.cluster)
#working

n.cores <- parallel::detectCores() - 1
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)

#for each SJ & gene combo,
#run through the possible transcripts and check for NMD triggering
sj_NMD_or_no <- foreach(
  i = 1:length(sj_pc),
  .combine = 'rbind'
) %dopar% {
  NMD = "unknown"
  gene_gr <- ann_gtf[ann_gtf$gene_id == sj_pc[i]$gene_id]
  gene_cds_gr <- gene_gr[gene_gr$type == "CDS"] #get coding transcripts of the gene
  gene_df <- data.frame(gene_cds_gr)
  transcripts <- unique(gene_df$transcript_id)
  for (transcript in transcripts) { #check each coding transcript in this gene
    transcript_df <- gene_df %>%
      filter(transcript_id == transcript) %>%
      filter(type == "CDS")
    relevant <- data.frame() #set up a df for relevant exons
    for (j in 1:nrow(transcript_df)) { #check each exon in this transcript
      if (nrow(relevant) == 0) { #no donor found yet; look for one
        if (start(sj_pc[i])-1 == transcript_df[j,3]) {
          relevant <- transcript_df[j,] #SJ donor is a relevant exon. Add it.
        } else { next }
      } else { #if there is content in "relevant" (i.e. we've found a matching donor)
        if (end(sj_pc[i])+1 == transcript_df[j,2]) { #look for matching acceptor
          relevant <- rbind(relevant, transcript_df[j,])
          break #if matching acceptor is found, we have all the exons we need
        } else {
          relevant <- rbind(relevant, transcript_df[j,]) #no acceptor yet: might be a SE
        }
      }
    } #end of this transcript
    #if there are skipped exons, check annotation status:
    if (nrow(relevant) < 2) { next } #splice site not annotated; ignore for now
    else if (end(sj_pc[i])+1 != relevant[nrow(relevant),2]) { next } #splice site not annotated; ignore for now
    else {
      if (nrow(relevant) == 2) { next } #fully annotated; ignore for now
      else {
        se <- relevant[2:(nrow(relevant) - 1), ] #extract skipped exons
        se$exon_length <- se$end - se$start + 1 #get exon length. Could use width column instead
        #splice sites are annotated but junction isn't, so check frame:
        if (sum(se$exon_length) %% 3 == 0) { next } #in-frame: assume no NMD
        else {
          #test: remove the call to check_NMD()
          #NMD <- check_NMD(transcript, se$exon_number, gene_cds_gr)
          #put in a meaningless test
          if (i %% 2) { NMD = "yes" }
          else {NMD = "no"}
          if (NMD == "yes") { #NMD. Add to output and move to next SJ
            new_row <- data.frame(sj_pc[i])
            new_row$NMD <- "yes"
            #return(new_row)
            break
          }
          else { next } #no NMD in this transcript, move to the next one
        }
      }
    }
  }
  #once all the transcripts are checked, if no NMD is found add the SJ to the no-NMD list
  if (NMD != "yes") {
    new_row <- data.frame(sj_pc[i])
    new_row$NMD <- "no"
    #return(new_row)
  }
}
#doesn't finish, even with small example. NEED TO TROUBLESHOOT
parallel::stopCluster(cl = my.cluster)

#old nonparallel:
sj_NMD_or_no <- data.frame()

for (i in 1:length(sj_pc)) {
  #need to get gene_df and transcripts list here
  NMD = "unknown"
  #print(i)
  #print(sj_pc[i]$gene_id)
  gene_gr <- ann_gtf[ann_gtf$gene_id == sj_pc[i]$gene_id]
  gene_cds_gr <- gene_gr[gene_gr$type == "CDS"]
  gene_df <- data.frame(gene_cds_gr)
  transcripts <- unique(gene_df$transcript_id)
  for (transcript in transcripts) { #go through each transcript
    if (is.na(transcript)) { #ignore gene lines (NA in transcript field). Though actually shouldn't be any
      next
    } else {
      #print(transcript)
      transcript_df <- gene_df %>%
        filter(transcript_id == transcript) %>%
        filter(type == "CDS")
      
      #if the transcript is coding, get skipped exons:
      if (nrow(transcript_df) == 0) {
        #print(paste(transcript, "has no CDS"))
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
      
      #if there are skipped exons, check annotation status:
      if (nrow(relevant) < 2) {
        #print("one or both splice sites unannotated: not assessed") #ignore for now
      }  #check if acceptor is annotated
      else if (end(sj_pc[i])+1 != relevant[nrow(relevant),2]) { 
        #print("one splice site unannotated: not assessed") #ignore for now
      }
      else {
        if (nrow(relevant) == 2) {
          #print("No skipped exon in this transcript: assume no NMD") #assume no NMD
        } else {
          se <- relevant[2:(nrow(relevant) - 1), ] #extract skipped exons
          se$exon_length <- se$end - se$start + 1 #could use width column here
          
          #splice sites are annotated but junction isn't, so check frame:
          if (sum(se$exon_length) %% 3 == 0) { #in-frame: assume no NMD
            #print("new in frame SE")
          } else {
            #print(paste("new out-of-frame SE:", sum(se$exon_length), "nt"))
            NMD <- check_NMD(transcript, se$exon_number, gene_cds_gr)
            
            #if SJ triggers NMD, print it. If not, check next transcript
            if (NMD == "yes") {
              new_row <- data.frame(sj_pc[i])
              new_row$NMD <- "yes"
              sj_NMD_or_no <- rbind(sj_NMD_or_no, new_row)
              #print(paste("now it contains", nrow(sj_NMD), "SJs") )
              break
            }
            else {
              #print("no NMD!")
              }
          }
        }
      }
    }
  }
  
  #once all the transcripts are checked, if no NMD is found add the SJ to the no-NMD list
  if (NMD != "yes") {
    new_row <- data.frame(sj_pc[i])
    new_row$NMD <- "no"
    sj_NMD_or_no <- rbind(sj_NMD_or_no, new_row)
  }
}


#output of this loop: 
#     sj_NMD_or_no (df of SJs in protein-coding genes, and whether or not they 
#          are likely NMD-targeting)
#maybe split this into two dfs to print out?

#should produce other lists as well:
#    in-frame SJs
#    out-of-frame SJs that don't meet NMD rules
#    SJs that aren't assessed for NMD (intergenic or with unannotated sites)
# and maybe a chart to compare numbers? Yes, probably
# one with numbers of junctions, one with junction read depth (maybe distribution?)
# also allow filtering by read number or maybe TPM (would require more input)
# plots at both SJ level and gene level

#simple bargraph for single sample
#later will do a version with multiple samples
#NEED TO UPDATE TO REFLECT COMBINED NMD/NO NMD DATAFRAME!!!
type_sum <- data.frame(
  dataset = c(1,1,1,1),
  type = c("intergenic", "noncoding", "NMD", "no_NMD"),
  count = c(nrow(sj_notann_df), nrow(sj_notPC_df), nrow(sj_NMD), nrow(sj_no_NMD)),
  reads_unique = c(sum(sj_notann_df$unique), sum(sj_notPC_df$unique), sum(sj_NMD$unique), sum(sj_no_NMD$unique)),
  reads_multi = c(sum(sj_notann_df$multi), sum(sj_notPC_df$multi), sum(sj_NMD$multi), sum(sj_no_NMD$multi))
)

type_sum$type <- factor(type_sum$type, levels = c("no_NMD", "intergenic", "noncoding","NMD"))

ggplot(type_sum, aes(x=dataset, y=count, fill = type)) +
  geom_bar(stat = "identity") +
  ylab("Number of splice junctions") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank())

ggplot(type_sum, aes(x=dataset, y=reads_unique, fill = type)) +
  geom_bar(stat = "identity") +
  ylab("Number of uniquely-mapped splice junction reads") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank())


#test code ####
#to delete once function/loop is better tested
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

#transcript = "ENST00000482437"
transcript = "ENST00000370143"
ses = se$exon_number

#get GRanges object for transcript CDS exons
trans_gr <- gene_cds_gr[gene_cds_gr$transcript_id == transcript]
#get sequence of the (annotated) transcript
trans_full_seq <- getSeq(Hsapiens, trans_gr)
trans_full_seq <- suppressWarnings(getSeq(Hsapiens, trans_gr))

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

