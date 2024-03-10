setwd("/Users/chris/Desktop/20240308_mmus_prostate_EAW")

library(tidyverse)
library(knitr)

#
fastq_list_filename <- "rnaseq_mmus_prostate_EAW_fastq_list.txt"
df <- read_tsv(file.path("input", "rnaseq_mmus_prostate_EAW", fastq_list_filename))
output_pipeline_dir <- "rna-pipeline_mmus_prostate_EAW-GRCh38_PE"

#
slist <- c("s0", "s1", "s2")

df_sumcount <- data.frame()
df_raw_counts <- list()
df_list <- list()
for (i in 1:nrow(df)) {
  sample_id_i_tmp <- df$sample_name[i]
  sample_id_i <- gsub("\\+", "x", sample_id_i_tmp)
    res <- c(sample_id_i)
    
    ### featureCounts
    for (sparam in slist) {
        message("# ", i, " | ", sample_id_i, " | ", sparam)
        sample_i_filename <- paste0(sample_id_i, ".featureCounts.tsv")
        sample_i_filepath <- file.path("output", output_pipeline_dir, paste(sep = "_", "featureCounts", sparam), sample_i_filename)
        rawData_s <- read_tsv(sample_i_filepath, col_names = TRUE, show_col_types = FALSE, skip = 1) %>% 
            set_names("gene_id", "chr", "start", "end", "strand", "length", "count")
        
        # sumcount
        sumcount_sample_s <- rawData_s %>% pull(count) %>% sum
        message("\t> ", sparam, " : ", sumcount_sample_s)
        res <- c(res, sumcount_sample_s)
        
        # df_count 
        df_count_i <- rawData_s %>% dplyr::select(gene_id, count) %>% 
            set_names("gene_id", sparam)
        df_list[[sparam]] <- df_count_i
        
    }
    
    ### htseq_count
        # sumcount
    sample_i_filename <- paste0(sample_id_i, ".readcounts.tsv")
    sample_i_filepath <- file.path("output", output_pipeline_dir, "raw_counts", sample_i_filename)
    rawData_s_tmp <- read_tsv(sample_i_filepath, col_names = FALSE, show_col_types = FALSE) %>%
        set_names("gene_id", "count")
    
    rawData_s <- rawData_s_tmp %>% dplyr::filter(grepl(pattern = "ENSG", x = gene_id))

    sumcount_sample_s <- rawData_s %>% pull(count) %>% sum
    message("\t> ", "htseq_reverse", " : ", sumcount_sample_s)
    rawData_s_tmp %>% dplyr::filter(!grepl(pattern = "ENSG", x = gene_id)) %>% kable %>% print
    res <- c(res, sumcount_sample_s)

        # sumcount gather everything
    df_sumcount <- rbind(df_sumcount, res)
    
        # df_count
    df_count_i <- rawData_s %>% dplyr::select(gene_id, count) %>% 
        set_names("gene_id", "htseq_reverse")
    df_list[["htseq_reverse"]] <- df_count_i
    
    # all_raw_counts
    all_raw_counts <- plyr::join_all(df_list, type = "full")
    df_raw_counts[[sample_id_i]] <- all_raw_counts
}

# sumcount
colnames(df_sumcount) <- c("sample_id", "s0", "s1", "s2", "htseq_reverse")
df_sumcount %>% 
    mutate(s0 = as.numeric(s0) %>% format(big.mark=" "),
           s1 = as.numeric(s1) %>% format(big.mark=" "),
           s2 = as.numeric(s2) %>% format(big.mark=" "),
           htseq_reverse = as.numeric(htseq_reverse) %>% format(big.mark=" ")) %>% 
    kable

# df_raw_counts
names(df_raw_counts)
df_raw_counts[[1]] %>% head(50)
