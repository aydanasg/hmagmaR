#' hmagmaR: Incorporation of chromatin interaction and epigenetic profiles to predict risk genes
#'
#' @import GenomicRanges
#' @importFrom IRanges IRanges subsetByOverlaps
#' @importFrom S4Vectors mcols mcols<- queryHits subjectHits
#' @importFrom ChIPseeker annotatePeak
#' @importFrom TxDb.Hsapiens.UCSC.hg19.knownGene TxDb.Hsapiens.UCSC.hg19.knownGene
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom stats aggregate na.omit
#' @keywords internal
"_PACKAGE"
