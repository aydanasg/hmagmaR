#' Gene level analysis carried out by magma (https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004219)
#'
#' Generates an annotation file for magma input
#' @param hic Dataframe of chromatin interaction data with 6 columns in the format: "chrom1", "start1", "end1", "chrom2", "start2", "end2". Chromosome columns should have "chr" before the number
#' @param regulatoryRegions Dataframe of regulatory regions (H3K27ac) in the .bed format (chr, start, end). Chromosome columns should have "chr" before the number
#' @param snps Dataframe of reference snps. This could be the ".bim" file from g1000 reference genome in the format: "chr", "position", "rsid"
#' @param annotated_genes Dataframe of reference genes. This should contain columns with names "chr", "start", "end",	"ensg". Chromosome columns should have "chr" before the number
#' @param snpgeneexon Dataframe of genes and snps within it. This should contain columns with names "rsid", "ensg"
#' @param loopNumber Number to sample down the loop number to. This is useful when comparing more than one sample to ensure that the difference is not down to the difference in loop number
#' @param sampleDownNumber Number of times sample down the loop number randomly
#' @param AnnotationFile Path where to save the sampled down annotation files
#' @return Function generates sampled down gene level analysis files (.genes.raw, .genes.out, .log.suppl, .log)
#' @examples
#' data("test_dataset.microglia.hg19")
#' out <- tempfile()
#' SampledDownAnnotation(
#'   hic = hic, regulatoryRegions = regulatoryRegions, snps = snps,
#'   annotated_genes = annotated_genes, snpgeneexon = snpgeneexon,
#'   loopNumber = 5000, sampleDownNumber = 1, AnnotationFile = out
#' )
#' @export

## PLAC-seq filtered to regulatoryRegions promoter interactions
SampledDownAnnotation <- function(hic, regulatoryRegions, snps, annotated_genes, snpgeneexon, loopNumber, sampleDownNumber, AnnotationFile){

  #Save and restore the caller's RNG state so the set.seed() calls below
  #(needed for reproducible sample-down iterations) don't leak into the
  #caller's session
  if (exists(".Random.seed", envir = .GlobalEnv)) {
    oldseed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", oldseed, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  }

  #reading snps file -  g1000 reference genome from European ancestry (.bim) and selecting only relevant columns (chr, rsid, position)
  snps<-GRanges(snps$chr, IRanges(snps$Position, snps$Position), rsid=snps$SNP)

  #Reading regulatoryRegions data (Nott et al., 2019 paper)
  colnames(regulatoryRegions)<-c("chr", "start", "end")
  regulatoryRegions_ranges<-GRanges(regulatoryRegions$chr, IRanges(as.numeric(regulatoryRegions$start), as.numeric(regulatoryRegions$end)))

  #Creating gene coordinate information (chr:start:end)
  annotated_genes$index <- paste(annotated_genes$chr, annotated_genes$start, annotated_genes$end, sep=":")

#Original plac-seq interactions, sampled down fresh on each iteration below
hic_full <- hic

for(s in 1:sampleDownNumber) {
   #Sampling 60000 random rows from the plac-seq data table using set.seed() function in data.table package
   #(https://stackoverflow.com/questions/8273313/sample-random-rows-in-dataframe)
   set.seed(s)
   hic <- hic_full[sample(nrow(hic_full), size=loopNumber),]
   
  #Doubling the hic and flipping the start1-end1 with int1-int2 (this is for GRanges function and IRanges)
  hic.int1 <- hic[,1:6]
  hic.int2 <- hic[,c(4:6,1:3)]
  colnames(hic.int1) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2")
  colnames(hic.int2) <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2")
  hic <- rbind(hic.int1, hic.int2)
  
  #Creating GRanges for findOverlaps() function
  hicranges <- GRanges(hic$chrom1, IRanges(as.numeric(hic$start1), as.numeric(hic$end1)),
                       int1=hic$start2,int2=hic$end2)
  
  #Selecting plac-seq interactions that overlap with promoters
  olap<-findOverlaps(query = hicranges, subject = regulatoryRegions_ranges) #default maxgap=-1L, when one range's start/end strictly inside the other, #the gap is considered to be -1.
  placranges<-hicranges[queryHits(olap)]
  mcols(placranges)<-cbind(mcols(placranges), mcols(regulatoryRegions_ranges[subjectHits(olap)]))
  
  #annotating consensus peaks to promoters/exons
  annotation<-annotatePeak(peak = placranges, tssRegion = c(-2000, 500), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb = "org.Hs.eg.db")
  annotation_df<-as.data.frame(annotation@anno)
  #selecting only promoter and exon annotations peaks 
  pattern<-paste(c("Promoter", "Exon"), collapse = "|")
  genes_df<- annotation_df[grepl(pattern = pattern, x = annotation_df$annotation), ]
  genes_df<-unique(genes_df)
  #creating GRanges for those peaks
  genesnpranges<-GRanges(seqnames = genes_df$seqnames, ranges = IRanges(start = genes_df$int1, end = genes_df$int2), snp.start=genes_df$int1, snp.end=genes_df$int2, gene.start=genes_df$start, gene.end=genes_df$end, ensg=genes_df$ENSEMBL)
  
  #Selecting SNPs that fall within regions that interact with exons and promoters
  olap <- findOverlaps(query = snps, subject = genesnpranges) #default maxgap=-1L, when one range's start/end strictly inside the other, the gap is considered to be -1.
  snpint <- snps[queryHits(olap)]
  mcols(snpint) <- cbind(mcols(snpint), mcols(genesnpranges[subjectHits(olap)]))
  
  #Making a dataframe with plac-ranges interacting regions and ONLY THEN removing row duplicates
  #(unique() on a GRanges directly only compares seqnames/ranges/strand, not mcols like ensg -
  #that would silently collapse a SNP that interacts with multiple different genes down to one)
  snpint_data<-as.data.frame(snpint)
  snpint_data <- unique(snpint_data)

  #Making GRanges for the SNPs located within exon/promoter-interacting regions
  snpint<-GRanges(seqnames = snpint_data$seqnames, ranges = IRanges(snpint_data$start, snpint_data$end), strand = snpint_data$strand,
                  rsid=snpint_data$rsid, ensg=snpint_data$ensg, snp.start=snpint_data$snp.start, snp.end=snpint_data$snp.end, gene.start=snpint_data$gene.start, gene.end=snpint_data$gene.end)

  #Selecting exonic/promoter-interacting SNPS overlapping with regulatoryRegions
  olap<-findOverlaps(query = snpint, subject = regulatoryRegions_ranges) #default maxgap=-1L, when one range's start/end strictly inside the other, the gap is considered to be -1.
  regulatoryRegions_snp_ranges<-snpint[queryHits(olap)]
  mcols(regulatoryRegions_snp_ranges)<-cbind(mcols(regulatoryRegions_snp_ranges), mcols(regulatoryRegions_ranges[subjectHits(olap)]))
  
  #Making a dataframe with exonic/promoter-interacting regulatoryRegions SNPs and ONLY THEN removing row duplicates 
  regulatoryRegions_snp_ranges_data<-as.data.frame(regulatoryRegions_snp_ranges) #•
  regulatoryRegions_snp_ranges_data<-unique(regulatoryRegions_snp_ranges_data) #• 
  
  #Making GRanges for the regulatoryRegions SNPs located within exon/promoter-interacting regions 
  regulatoryRegions_snp_ranges<-GRanges(seqnames = regulatoryRegions_snp_ranges_data$seqnames, ranges = IRanges(regulatoryRegions_snp_ranges_data$start, regulatoryRegions_snp_ranges_data$end), strand = regulatoryRegions_snp_ranges_data$strand,
                                        rsid=regulatoryRegions_snp_ranges_data$rsid, ensg=regulatoryRegions_snp_ranges_data$ensg)
  
  #Saving exonic/promoter interacting regulatoryRegions SNPs and associated genes 
  snpintmat<-unique(data.frame(rsid=regulatoryRegions_snp_ranges$rsid, ensg=regulatoryRegions_snp_ranges$ensg))
  snpintmat<-na.omit(snpintmat)

  #Selecting genes and all the associated SNPs that have interactions with regulatoryRegions:
  #the regulatoryRegions SNPs plus the exonic/promoter SNPs of those same genes
  exonmat <- snpgeneexon[snpgeneexon$ensg %in% snpintmat$ensg, c("rsid", "ensg")]
  snpcomb_loop_melt <- na.omit(unique(rbind(snpintmat, exonmat)))
  
  #Making MAGMA compatible annotation file
  #Aggregate rsids by ensg
  aggregated_df <- aggregate(rsid ~ ensg, data = snpcomb_loop_melt,
                              FUN = function(x) paste(x, collapse = "\t"))

  #Adding the index information to the gene-snp file
  aggregated_df$index<-annotated_genes[match(aggregated_df$ensg, annotated_genes$ensg),"index"]

  #Removing any NA
  aggregated_df <- aggregated_df[!is.na(aggregated_df$index),]

  writeLines(paste(aggregated_df$ensg, aggregated_df$index, aggregated_df$rsid, sep = "\t"),
             paste0(AnnotationFile, ".transcript.annot.", s))
}
}
  