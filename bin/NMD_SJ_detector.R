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
gene_test_gr <- ann_gtf[ann_gtf$gene_id == "ENSG00000122435"]
gene_test_df <- data.frame(gene_test_gr)
ann_gtf[ann_gtf$transcript_id == "ENST00000370143"]
#Error: logical subscript contains NAs
#could possibly get around this error by subsetting based on genomic coordinates of gene


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

