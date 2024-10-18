remotes::install_github("aydanasg/hmagmaR", force=TRUE)
library(hmagmaR)
??hmagmaR::AnnotationFileHmagma
hmagmaR::AnnotationFileHmagma()

library(roxygen2)
roxygenise()



library(hmagmaR)
hmagmaR::GeneLevelAnalysis_hmagma(
  magma="/Volumes/aa19618/projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/magma/HMAGMA_system/magma_v1.10_mac/magma",
  g1000="/Volumes/aa19618/home/HMAGMA_Protocol/required_files/g1000/g1000_files/g1000_eur",
  gwas="/Volumes/aa19618/projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/gwas_studies/munged_files/AD_Jansen2019_munged.tsv",
  AnnotationFile="/Volumes/aa19618/projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/magma/annotation_files_SD/ADvas_AAA_20240422_ac/SNP_aggregate_transcriptpu1.1.chipSeeker.transcript.annot",
  output="/Volumes/aa19618/projects/epinott/live/user_analysed_data/Aydan/vasculature_disease_epi/magma/magma_output/pu1.1")
