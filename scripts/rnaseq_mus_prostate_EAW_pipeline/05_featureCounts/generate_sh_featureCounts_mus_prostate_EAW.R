setwd("/Users/chris/Desktop/20240308_mmus_prostate_EAW")

##### module load subread/2.0.3

library(tidyverse)
library(knitr)

fastq_list_filename <- "rnaseq_mmus_prostate_EAW_fastq_list.txt"
df <- read_tsv(file.path("input", "rnaseq_mmus_prostate_EAW", fastq_list_filename))
output_pipeline_dir <- "rna-pipeline_mmus_prostate_EAW-GRCh38_PE"
script_pipeline_dir <- "rnaseq_mus_prostate_EAW_pipeline"
workdir <- "/home/chris11/projects/def-stbil30/chris11/20240308_mmus_prostate_EAW"

strandness <- 2
strandness_chr <- paste0("s", strandness)

alignment_dir <- file.path("output", output_pipeline_dir, "alignment")
featureCount_dir <- file.path("output", output_pipeline_dir, paste(sep = "_", "featureCounts", strandness_chr))
ensembl_path <- "input/ensembl/Homo_sapiens.GRCh38.104.gtf"
# gencode_path <- file.path(workdir, "input", "gencode", "gencode.v29.annotation.gtf")

message("mkdir -p ", file.path(workdir, featureCount_dir))

header_sh <- c("#!/bin/sh",
               "#SBATCH --time=3:00:00",
               "#SBATCH --nodes=1",
               "#SBATCH --ntasks-per-node=1",
               "#SBATCH --cpus-per-task=4",
               "#SBATCH --mem-per-cpu=8G",
               "#SBATCH --account=def-stbil30",
               "#SBATCH --mail-user=christophe.tav@gmail.com",
               "#SBATCH --mail-type=ALL")

for (i in 1:nrow(df)) {
  # for (i in 1:5) {
  # timepoint_i <- df$time_point[i]
  # rep_i <- df$isogenic_replicate[i]
  # basename <- paste(sep = "_", "A549", "RNA", timepoint_i, rep_i)
  # message("# ", i, " | ", basename)
  
  sample_name_i_tmp <- df$sample_name[i]
  sample_name_i <- gsub("\\+", "x", sample_name_i_tmp)
  basename <- sample_name_i
  # message("# ", i, " | ", basename)
  
  
  # bam
  input_filename <- paste0(basename, ".QueryNameSorted.bam")
  input_filepath <- file.path(workdir, alignment_dir, basename, input_filename)
  
  # output
  output_filename <- paste0(basename, ".featureCounts.tsv")
  output_filepath <- file.path(workdir, featureCount_dir, output_filename)
  
  call_featureCounts <- paste("featureCounts",
                              "-T", "1", "-p", "--countReadPairs", "-B", "-t", "exon", "-g", "gene_id",
                              "-s", strandness,
                              "-a", ensembl_path ,
                              "-o", output_filepath,
                              input_filepath)
  # message(call_featureCounts)
  
  file_sh <- file.path("scripts", script_pipeline_dir, "05_featureCounts/batch_sh",
                       paste0("featureCounts_", basename, "_", strandness_chr, ".sh"))
  message("sbatch ", file_sh)
  fileConn <- file(file_sh)
  writeLines(c(header_sh, "\n", call_featureCounts), fileConn)
  close(fileConn)
}
