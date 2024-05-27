#' Gene level analysis carried out by magma (https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004219)
#'
#' Generates an annotation file for magma input
#' @param g1000 Path to 1,000 Genomes reference data files (.bin, .bed, .fam, .synonyms). Include the name of the of the files in the path (e.g. /path_to_g1000/g1000_eur). Obtain from https://cncr.nl/research/magma/ under Auxiliary files.
#' @param gwas Path to the GWAS summary statistics to be used where rsid column name is "SNP", p-value column name is "P", and sample size column name is "N"
#' @param AnnotationFile Annotation file generated using AnnotationFileHmagma() or SampledDownAnnotation()  
#' @param output Name of the output file 
#' @return Function generates gene level analysis files (.genes.raw, .genes.out, .log.suppl, .log)
#' @export
GeneLevelAnalysis_hmagma <- function(g1000, gwas, AnnotationFile, output) {
    system(paste("HMAGMA_system/magma", "--bfile", g1000, "--pval", gwas, "use=SNP,P ncol=N", "--gene-annot", AnnotationFile, "--out", output))
}
