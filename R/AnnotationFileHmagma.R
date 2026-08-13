#' Gene level analysis carried out by magma (https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004219)
#'
#' Generates an annotation file for magma input
#' @param hic Dataframe of chromatin interaction data with 6 columns in the format: "chrom1", "start1", "end1", "chrom2", "start2", "end2". Chromosome columns should have "chr" before the number
#' @param regulatoryRegions Dataframe of regulatory regions (H3K27ac) in the .bed format (chr, start, end). Chromosome columns should have "chr" before the number
#' @param promoterRegions Dataframe of promoters in the .bed format (chr, start, end). Chromosome columns should have "chr" before the number
#' @param enhancerRegions Dataframe of enhancers in the .bed format (chr, start, end). Chromosome columns should have "chr" before the number
#' @param snps Dataframe of reference snps. This could be the ".bim" file from g1000 reference genome in the format: "chr", "position", "rsid"
#' @param annotated_genes Dataframe of reference genes. This should contain columns with names "chr", "start", "end",	"ensg". Chromosome columns should have "chr" before the number
#' @param snpgeneexon Dataframe of genes and snps within it. This should contain columns with names "rsid", "ensg"
#' @param AnnotationFile Path prefix for the output file; ".transcript.annot" is appended
#' @return Function generates gene level analysis files (.genes.raw, .genes.out, .log.suppl, .log)
#' @examples
#' data("test_dataset.microglia.hg19")
#' out <- tempfile()
#' AnnotationFileHmagma(
#'   hic = hic, promoterRegions = promoterRegions, enhancerRegions = enhancerRegions,
#'   snps = snps, annotated_genes = annotated_genes, snpgeneexon = snpgeneexon,
#'   AnnotationFile = out
#' )
#' @export

## PLAC-seq filtered to regulatoryRegions promoter interactions
AnnotationFileHmagma <- function(hic, regulatoryRegions = NULL, promoterRegions = NULL, enhancerRegions = NULL, snps, annotated_genes, snpgeneexon, AnnotationFile) {

  # Two supported modes. They share the same pipeline and differ only in which
  # regions anchor the hic interactions and which regions the SNPs must fall in:
  #   - regulatoryRegions: the same set is used for both steps
  #   - promoterRegions + enhancerRegions: promoters anchor, enhancers filter SNPs
  if (!is.null(regulatoryRegions)) {
    message("Using regulatoryRegions...")
    anchorRegions <- regulatoryRegions
    snpRegions <- regulatoryRegions
  } else if (!is.null(promoterRegions) && !is.null(enhancerRegions)) {
    message("Using promoterRegions and enhancerRegions")
    anchorRegions <- unique(promoterRegions)
    snpRegions <- enhancerRegions
  } else {
    stop("Either regulatoryRegions or all three (promoterRegions and enhancerRegions) must be provided.")
  }

  #reading snps file -  g1000 reference genome from European ancestry (.bim) and selecting only relevant columns (chr, rsid, position)
  snpranges <- GRanges(snps$chr, IRanges(snps$Position, snps$Position), rsid = snps$SNP)

  #Doubling the hic and flipping the start1-end1 with int1-int2 (this is for GRanges function and IRanges)
  bedCols <- c("chrom1", "start1", "end1", "chrom2", "start2", "end2")
  hic.int1 <- stats::setNames(hic[, 1:6], bedCols)
  hic.int2 <- stats::setNames(hic[, c(4:6, 1:3)], bedCols)
  hic <- rbind(hic.int1, hic.int2)

  #Creating GRanges for findOverlaps() function
  hicranges <- GRanges(hic$chrom1, IRanges(as.numeric(hic$start1), as.numeric(hic$end1)),
                       int1 = hic$start2, int2 = hic$end2)

  #Selecting plac-seq interactions that overlap with promoters
  #default maxgap=-1L, when one range's start/end strictly inside the other, the gap is considered to be -1.
  anchorRanges <- .bedToGRanges(anchorRegions)
  placranges <- IRanges::subsetByOverlaps(hicranges, anchorRanges)

  #annotating consensus peaks to promoters/exons
  annotation <- annotatePeak(peak = placranges, tssRegion = c(-2000, 500),
                             TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb = "org.Hs.eg.db")
  annotation_df <- as.data.frame(annotation@anno)
  #selecting only promoter and exon annotations peaks
  genes_df <- unique(annotation_df[grepl("Promoter|Exon", annotation_df$annotation), ])

  #creating GRanges for those peaks: the interacting (distal) window, tagged with the gene it reaches
  genesnpranges <- GRanges(seqnames = genes_df$seqnames,
                           ranges = IRanges(start = genes_df$int1, end = genes_df$int2),
                           ensg = genes_df$ENSEMBL)

  #Selecting SNPs that fall within regions that interact with exons and promoters
  olap <- findOverlaps(query = snpranges, subject = genesnpranges)
  snpint <- snpranges[queryHits(olap)]
  mcols(snpint) <- cbind(mcols(snpint), mcols(genesnpranges[subjectHits(olap)]))

  #Collapsing duplicate SNP-gene pairs here keeps the second overlap step small
  #(unique() on a GRanges only compares seqnames/ranges/strand, not mcols like ensg - that would
  #silently collapse a SNP interacting with multiple different genes down to one, so dedupe on
  #the full row instead)
  snpint_df <- unique(as.data.frame(snpint))
  snpint <- GRanges(seqnames = snpint_df$seqnames, ranges = IRanges(snpint_df$start, snpint_df$end),
                     rsid = snpint_df$rsid, ensg = snpint_df$ensg)

  #Selecting exonic/promoter-interacting SNPs overlapping with the regulatory/enhancer regions
  snpRanges <- .bedToGRanges(snpRegions)
  snpint <- IRanges::subsetByOverlaps(snpint, snpRanges)

  #Saving exonic/promoter interacting regulatoryRegions SNPs and associated genes
  snpintmat <- unique(data.frame(rsid = snpint$rsid, ensg = snpint$ensg, stringsAsFactors = FALSE))
  snpintmat <- na.omit(snpintmat)

  #Selecting genes and all the associated SNPs that have interactions with the regulatory regions:
  #the interacting SNPs plus the exonic/promoter SNPs of those same genes
  exonmat <- snpgeneexon[snpgeneexon$ensg %in% snpintmat$ensg, c("rsid", "ensg")]
  snpcomb <- na.omit(unique(rbind(snpintmat, exonmat)))

  #Making MAGMA compatible annotation file
  # Aggregate rsids by ensg
  aggregated_df <- aggregate(rsid ~ ensg, data = snpcomb,
                              FUN = function(x) paste(x, collapse = "\t"))

  #Creating gene coordinate information (chr:start:end) and attaching it to each gene
  index <- paste(annotated_genes$chr, annotated_genes$start, annotated_genes$end, sep = ":")
  aggregated_df$index <- index[match(aggregated_df$ensg, annotated_genes$ensg)]
  aggregated_df <- aggregated_df[!is.na(aggregated_df$index), ]

  writeLines(paste(aggregated_df$ensg, aggregated_df$index, aggregated_df$rsid, sep = "\t"),
             paste0(AnnotationFile, ".transcript.annot"))
}

# Coerce a 3-column bed-style data frame (chr, start, end) to GRanges
.bedToGRanges <- function(bed) {
  GRanges(bed[[1]], IRanges(as.numeric(bed[[2]]), as.numeric(bed[[3]])))
}
