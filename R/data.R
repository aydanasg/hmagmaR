#' Example microglia PLAC-seq dataset
#'
#' A small, randomly subsampled set of objects used to demonstrate the
#' \code{hmagmaR} workflow in the package vignette. Because \code{LazyData}
#' is enabled, all of the objects below become available directly by name
#' as soon as \code{library(hmagmaR)} is called. Because the rows are drawn
#' independently at random from each source table (rather than hand-picked
#' around loci known to overlap), results produced from this dataset are
#' illustrative only and not biologically meaningful.
#'
#' @format Seven objects are provided:
#' \describe{
#'   \item{hic}{A data frame/tibble of 20,000 microglia PLAC-seq chromatin
#'     interactions (Nott et al., 2019), with columns \code{chr1}, \code{start1},
#'     \code{end1}, \code{chr2}, \code{start2}, \code{end2}, plus interaction-calling
#'     metadata (\code{count}, \code{expected}, \code{fdr}, \code{ClusterLabel},
#'     \code{ClusterSize}, \code{ClusterType}, \code{ClusterNegLog10P},
#'     \code{ClusterSummit}).}
#'   \item{promoterRegions}{A data frame of 12,709 microglia promoter regions
#'     (Nott et al., 2019) in BED-like format: \code{chr}, \code{start}, \code{end}.}
#'   \item{enhancerRegions}{A data frame of 44,090 microglia enhancer regions
#'     (Nott et al., 2019) in BED-like format: \code{chr}, \code{start}, \code{end}.}
#'   \item{snps}{A data table of 20,000 reference SNPs subsampled from 1000
#'     Genomes (European ancestry), with columns \code{chr}, \code{Position},
#'     \code{SNP}.}
#'   \item{annotated_genes}{A data table of 20,000 reference gene coordinates,
#'     with columns \code{chr}, \code{start}, \code{end}, \code{ensg}.}
#'   \item{snpgeneexon}{A data table of 20,000 SNPs falling within exons/promoters
#'     of a gene, with columns \code{rsid}, \code{ensg}.}
#'   \item{gwas}{A data table of 20,000 example GWAS summary statistics rows,
#'     with columns \code{SNP}, \code{CHR}, \code{BP}, \code{A1}, \code{A2},
#'     \code{FRQ}, \code{BETA}, \code{SE}, \code{P}, \code{N}.}
#' }
#' @source Subsampled from microglia PLAC-seq interactions and regulatory
#'   regions (Nott et al., 2019, \url{https://doi.org/10.1126/science.aay0793}),
#'   1000 Genomes reference SNPs, and example GWAS summary statistics.
#' @docType data
#' @keywords datasets
#' @name test_dataset.microglia.hg19
#' @aliases hic promoterRegions enhancerRegions snps annotated_genes snpgeneexon gwas
NULL

#' @rdname test_dataset.microglia.hg19
"hic"

#' @rdname test_dataset.microglia.hg19
"promoterRegions"

#' @rdname test_dataset.microglia.hg19
"enhancerRegions"

#' @rdname test_dataset.microglia.hg19
"snps"

#' @rdname test_dataset.microglia.hg19
"annotated_genes"

#' @rdname test_dataset.microglia.hg19
"snpgeneexon"

#' @rdname test_dataset.microglia.hg19
"gwas"
