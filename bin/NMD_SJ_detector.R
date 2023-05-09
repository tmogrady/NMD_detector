#to assess likelihood of NMD targeting by STAR-reported splice junctions
#IN DEVELOPMENT

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
#number of unique reads per million required to consider junctions
sj_thres <- 1

#get annotation (gtf) ####
#should get this from a package instead of a gtf file (or have option)
ann_gtf <- import("/Applications/Genomics_applications/Genomes_and_transcriptomes/hg38_plus_Akata_inverted.gtf")
ann_gtf_gene <- ann_gtf[mcols(ann_gtf)$type == "gene"]

#read in data ####
#read in STAR SJ.out.tab file
#sj <- read.table("input_data/MC1_truncated_test-SJ.out.tab")
#sj <- read.table("../temp/MC1_S34_L004_R1_001_MC1_S34_L004_R2_001.hg38plusAkata_inverted-SJ.out.tab")
sj <- read.table("../temp/MZ1_S38_L004_R1_001_MZ1_S38_L004_R2_001.hg38plusAkata_inverted-SJ.out.tab")

#name columns for GRanges (and my own sanity)
colnames(sj) <- c("chr","start","end","strand","motif","ann","unique","multi","overhang")
#change strand info to +/- for GRanges
sj <- sj %>%
  mutate(strand = ifelse(strand == 1, "+",
                         ifelse(strand == 2, "-",
                                "*")))

# #explore sj read depth a bit ####
# #using MC1_S34_L004_R1_001_MC1_S34_L004_R2_001.hg38plusAkata_inverted-SJ.out.tab
# #as an example
# ggplot(sj, aes(x=unique)) +
#   geom_histogram()
# 
# sum(sj$unique)
# #76,115,116 for MC1
# #29,220,512 for MZ1
# quantile(sj$unique)
# #MC1:
# # 0%   25%   50%   75%  100%
# # 0     1     4    66 86704
# #MZ1:
# # 0%   25%   50%   75%  100% 
# # 0     1     2    11 33376 
# nrow(sj)
# #377,437 for MC1
# #483,998 for MZ1
# sj %>%
#   filter(unique > 1) %>%
#   nrow()
# #255,235 for MC1
# #271,847 for MZ1
# sj %>%
#   filter(unique > 10) %>%
#   nrow()
# #138,704 for MC1
# #121,571 for MZ1
# 
# #calculate unique reads per million unique splice-junction mapped reads
# sj$uniquepm <- (sj$unique/sum(sj$unique))*1e06
# quantile(sj$uniquepm)
# #MC1
# # 0%          25%          50%          75%         100%
# # 0.000000e+00 1.313799e-02 5.255198e-02 8.671077e-01 1.139117e+03
# #MZ1 
# # 0%          25%          50%          75%         100% 
# # 0.000000e+00 3.422254e-02 6.844507e-02 3.764479e-01 1.142211e+03 
# sj %>%
#   filter(uniquepm > 0.5) %>%
#   nrow()
# #105375 for MC1
# #111629 for MZ1
# sj %>%
#   filter(uniquepm > 1) %>%
#   nrow()
# #90,947 for MC1
# #93,094 for MZ1
# sj %>%
#   filter(uniquepm > 10) %>%
#   nrow()
# #22407 for MC1
# #20404 for MZ1

#what is the best way to pick a threshold?
#possible improvement: machine learning component to select best threshold

#filter & sort SJs ####
sj$uniquepm <- (sj$unique/sum(sj$unique))*1e06
sj_filt <- sj %>%
  filter(uniquepm >= sj_thres)

#make SJ data a GRanges object:
sj_gr <- makeGRangesFromDataFrame(sj_filt, keep.extra.columns = TRUE)

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
#first set up parallelization:
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
  .combine = 'rbind',
  .packages = c("dplyr", "GenomicRanges", "GenomicFeatures", "BSgenome.Hsapiens.UCSC.hg38")
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
          NMD <- check_NMD(transcript, se$exon_number, gene_cds_gr)
          if (NMD == "yes") { #NMD. Add to output and move to next SJ
            new_row <- data.frame(sj_pc[i])
            new_row$NMD <- "yes"
            return(new_row)
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
    new_row$NMD <- "no" #try new_row$NMD <- NMD. Could be "no" or "unknown" i.e. unassessed. Right?
    return(new_row)
  }
}

parallel::stopCluster(cl = my.cluster)

#output of this loop: 
#     sj_NMD_or_no (df of SJs in protein-coding genes, and whether or not they 
#          are likely NMD-targeting)

#write some files to play with
# write.table(sj_NMD_or_no, "../temp/MC1_S34_L004_R1_001_MC1_S34_L004_R2_001.hg38plusAkata_inverted-SJ_NMD_or_no.txt",
#             row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

sj_NMD <- sj_NMD_or_no %>%
  filter(NMD == "yes")
sj_no_NMD <- sj_NMD_or_no %>%
  filter(NMD == "no")

#should eventually produce other lists as well:
#    in-frame SJs
#    out-of-frame SJs that don't meet NMD rules
#    SJs that aren't assessed for NMD (unannotated sites)
# also allow filtering by read number or maybe TPM (would require more input)
# plots at both SJ level and gene level

#simple bargraph for single sample
#later will do a version with multiple samples
type_sum <- data.frame(
  dataset = c(1,1,1,1),
  type = c("intergenic", "noncoding", "NMD", "no_NMD"),
  count = c(nrow(sj_notann_df), nrow(sj_notPC_df), nrow(sj_NMD), nrow(sj_no_NMD)),
  reads_unique = c(sum(sj_notann_df$unique), sum(sj_notPC_df$unique), sum(sj_NMD$unique), sum(sj_no_NMD$unique)),
  reads_multi = c(sum(sj_notann_df$multi), sum(sj_notPC_df$multi), sum(sj_NMD$multi), sum(sj_no_NMD$multi))
)

type_sum$type <- factor(type_sum$type, levels = c("no_NMD", "intergenic", "noncoding","NMD"))

#plot number of splice junctions
ggplot(type_sum, aes(x=dataset, y=count, fill = type)) +
  geom_bar(stat = "identity") +
  ylab("Number of splice junctions") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank())

#plot number of splice junction-spanning reads
ggplot(type_sum, aes(x=dataset, y=reads_unique, fill = type)) +
  geom_bar(stat = "identity") +
  ylab("Number of uniquely-mapped splice junction reads") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank())

#will probably need to somehow compare depth (not just number I think)
#of NMD-targeting splice junctions to other splice junctions
#because everything is going to have some of both

