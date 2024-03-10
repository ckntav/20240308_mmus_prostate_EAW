setwd("/Users/chris/Desktop/20240308_mmus_prostate_EAW")

library(tidyverse)
library(DESeq2)
library(RColorBrewer)
library(ComplexHeatmap)
library(wesanderson)

#
sparam <- "s2"

#
condition_list <- c("CTL", "E2", "ET", "T")
type_list <- c("Normal", "PCa")

# Prepare inputs
fastq_list_filename <- "rnaseq_mmus_prostate_EAW_fastq_list.txt"
design <- read_tsv(file.path("input", "rnaseq_mmus_prostate_EAW", fastq_list_filename)) %>%
  mutate(sample_id = sample_name) %>% 
  dplyr::select(sample_id, condition)
design$type <- factor(design$type, levels = type_list)
design$condition <- factor(design$condition, levels = condition_list)

bad_replicates_list <- c()

# 
design <- design %>%
  dplyr::filter(!sample_id %in% bad_replicates_list)

#
raw_counts <- read_csv(file.path("output/rna-pipeline_mmus_prostate_EAW-GRCh38_PE",
                                 paste(sep = "_", "featureCounts", sparam),
                                 "RNA_mus_prostate_EAW_raw_counts.csv"))

#
current_design <- design
current_samples <- design$sample_id
myCountData <- raw_counts %>% dplyr::select(all_of(current_samples))

#
dds <- DESeqDataSetFromMatrix(countData = myCountData,
                              colData = current_design,
                              design = ~ condition)

rownames(dds) <- raw_counts$gene_id

#
vsd <- vst(dds, blind = FALSE)
colnames(vsd) <- gsub("RNA_", "", current_design$sample_id)
head(assay(vsd), 10)

##### Heatmap
#
sampleDists <- dist(t(assay(vsd)))
sampleDists

#
sampleDistMatrix <- as.matrix(sampleDists)
sampleDistMatrix

print(max(sampleDistMatrix))

annotMatrix <- data.frame(sample_id = rownames(sampleDistMatrix)) %>% 
  mutate(condition = str_split(sample_id, pattern = "_") %>% map(1) %>% unlist)

tpMat <- annotMatrix %>% dplyr::select(type, condition) %>% as.matrix %>% t

cols_ht <- wes_palette("Zissou1", 4, type = "continuous") %>% as.vector
# names(cols_ht) <- tpMat[1, ] %>% unique # c("0h", time_point_list)

rowAnnot <- rowAnnotation("type" = tpMat[1, ],
                          "condition" = tpMat[1, ],
                          annotation_legend_param = list(
                            type = list(labels = type_list,
                                        at = type_list),
                            condition = list(labels = condition_list,
                                             at = condition_list)),
                          col = list(type = c("Normal" = "#8E6E53", "PCa" = "#5D3FD3"),
                                     condition = c("CTL" = "#3A9AB2", "E2" = "#ADC397",
                                                   "ET" = "#E5A208", "T" = "#F11B00")))

#
colorBlues <- colorRampPalette(rev(brewer.pal(9, "Blues")))(9)
colorBlues

#
ht <- Heatmap(sampleDistMatrix, col = colorBlues, name = "distance",
              row_dend_width = unit(25, "mm"), show_column_dend = FALSE,
              column_dend_side = "bottom", column_names_side = "top", column_names_rot = 45,
              right_annotation = rowAnnot)

ht

# per time_point
for (t in condition_list[2:4]) {
  message("# ", t)
  col_index <- grepl(paste0("^CTL|^", t), colnames(sampleDistMatrix))
  row_index <- grepl(paste0("^CTL|^", t), rownames(sampleDistMatrix))
  
  sampleDistMatrix_t <- sampleDistMatrix[col_index, row_index]
  
  annotMatrix <- data.frame(sample_id = rownames(sampleDistMatrix_t)) %>% 
    mutate(time_point = str_split(sample_id, pattern = "_") %>% map(1) %>% unlist)
  
  tpMat <- annotMatrix %>% dplyr::select(time_point) %>% as.matrix %>% t
  
  rowAnnot <- rowAnnotation("condition" = tpMat[1, ])
  
  #
  colorBlues <- colorRampPalette(rev(brewer.pal(9, "Blues")))(9)
  colorBlues
  
  #
  ht <- Heatmap(sampleDistMatrix_t, col = colorBlues, name = "distance",
                row_dend_width = unit(25, "mm"), show_column_dend = FALSE,
                column_dend_side = "bottom", column_names_side = "top", column_names_rot = 45,
                right_annotation = rowAnnot)
  
  print(ht)
  
  readline(prompt="Press [enter] to continue")
}

