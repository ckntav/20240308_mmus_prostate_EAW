setwd("/Users/chris/Desktop/20240308_mmus_prostate_EAW")

library(tidyverse)

#
fastq_list_filename <- "rnaseq_mmus_prostate_EAW_fastq_list.txt"
df <- read_tsv(file.path("input", "rnaseq_mmus_prostate_EAW", fastq_list_filename))
output_pipeline_dir <- "rna-pipeline_mmus_prostate_EAW-GRCh38_PE"

sparam <- "s0"

df_list <- list()
for (i in 1:nrow(df)) {
# for (i in 1:5) {
  sample_id_i_tmp <- df$sample_name[i]
  sample_i <- gsub("\\+", "x", sample_id_i_tmp)
    message("# ", i, " ", sample_i)
    
    sample_i_filename <- paste0(sample_i, ".featureCounts.tsv")
    
    rawData <- read_tsv(file.path("output", output_pipeline_dir, paste(sep = "_", "featureCounts", sparam), sample_i_filename),
                        col_names = TRUE, show_col_types = FALSE, skip = 1) %>%
        set_names("gene_id", "chr", "start", "end", "strand", "length", "count")
    
    count_sample <- rawData %>% dplyr::select(gene_id, count) %>% 
        set_names("gene_id", sample_i)
    
    df_list[[sample_i]] <- count_sample
}

raw_counts <- plyr::join_all(df_list, type = "full")

head(raw_counts)

raw_counts_filepath <- file.path("output", output_pipeline_dir, paste(sep = "_", "featureCounts", sparam), "RNA_mus_prostate_EAW_raw_counts.csv")
write_csv(raw_counts, file = raw_counts_filepath)
message(" > Raw counts saved in : ", raw_counts_filepath)
