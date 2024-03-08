setwd("/Users/chris/Desktop/20240308_mmus_prostate_EAW")

library(tidyverse)
library(kableExtra)

#
fastq_list_filename <- "rnaseq_mmus_prostate_EAW_fastq_list.txt"
df <- read_tsv(file.path("input", "rnaseq_mmus_prostate_EAW", fastq_list_filename))
fastq_folder <- "rnaseq_mmus_prostate_EAW"
workdir <- "/home/chris11/scratch/20240308_mmus_prostate_EAW"
output_dir <- file.path(workdir, "raw", fastq_folder, "raw_fastq")

#
df %>% dplyr::select(-fastq_R1_filepath, -fastq_R2_filepath) %>% 
  kbl %>% kable_classic_2(full_width = F, font_size = 12)

#
link_part1 <- "https://www.dropbox.com/scl/fo/wxqg2kd57yi1jtv615byh/h/GEO_MUS_RNASeq/"
link_part2 <- "?rlkey=0jkzmfq0qcbli40jl7tgh8ybk&dl=0"

#
all_cmd_line <- list()
# for (i in 1:nrow(df)) {
for (i in 1:nrow(df)) {  
  message("########################################")
  sample_name <- df$sample_name[i]
  message("# ", i, " | ", sample_name)
  
  ##### R1
  df[i, ] %>% dplyr::select(-fastq_R1_filepath, -fastq_R2_filepath, -fastq_R2_filename) %>% 
    t %>% kbl %>% kable_classic_2(full_width = F, font_size = 20) %>% print
  fastq_R1_filename <- df$fastq_R1_filename[i]
  message("\t> R1 : ", fastq_R1_filename)
  fastq_R1_filepath <- file.path(output_dir, fastq_R1_filename)
  link_R1 <- paste0(link_part1, fastq_R1_filename, link_part2)
  log_filename_R1 <- paste0("dl_fastq_", sample_name, "_R1.log")
  log_filepath_R1 <- file.path(output_dir, log_filename_R1)
  
  call_wget_R1 <- paste("wget",
                        # "--no-verbose",
                        "-O", fastq_R1_filepath,
                        "-o", log_filepath_R1,
                        link_R1, "&")
  message(call_wget_R1)
  # system(call_wget_R1)
  all_cmd_line[[paste(sep = "_", sample_name, "R1")]] <- call_wget_R1
  
  ##### R2
  df[i, ] %>% dplyr::select(-fastq_R1_filepath, -fastq_R2_filepath, -fastq_R1_filename) %>% 
    t %>% kbl %>% kable_classic_2(full_width = F, font_size = 20) %>% print
  fastq_R2_filename <- df$fastq_R2_filename[i]
  message("\t> R2 : ", fastq_R2_filename)
  fastq_R2_filepath <- file.path(output_dir, fastq_R2_filename)
  link_R2 <- paste0(link_part1, fastq_R2_filename, link_part2)
  log_filename_R2 <- paste0("dl_fastq_", sample_name, "_R2.log")
  log_filepath_R2 <- file.path(output_dir, log_filename_R2)
  
  call_wget_R2 <- paste("wget",
                        # "--no-verbose",
                        "-O", fastq_R2_filepath,
                        "-o", log_filepath_R2,
                        link_R2, "&")
  message(call_wget_R2)
  # system(call_wget_R2, wait = TRUE)
  all_cmd_line[[paste(sep = "_", sample_name, "R2")]] <- call_wget_R2
  
}

all_cmd_line
unlist(all_cmd_line)

bash_file <- "scripts/rnaseq_mus_prostate_EAW_pipeline/01b_download_fastq_onNarval/wget_fastq_cmdline.sh"
writeLines(unlist(all_cmd_line), bash_file)

