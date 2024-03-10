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
  dplyr::select(sample_id, type, condition)
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


##### PCA
pca_plot <- plotPCA(vsd, intgroup = "condition")
pca_plot + geom_text(label =  gsub("RNA_", "", current_design$sample_id),
                     size = 2,
                     vjust = 2) +
  theme_minimal()

pca_plot + theme_minimal()

#
cols_pca <- wes_palette("Zissou1", 4, type = "continuous") %>% as.vector
pca_plot +
  theme_minimal() +
  scale_colour_manual(values = cols_pca)

#
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)

pca_data %>% 
  ggplot(aes(x = PC1, y = PC2, color = group)) +
  geom_point() +
  geom_text(label = pca_data$name, vjust = 2) +
  theme_minimal() +
  scale_colour_manual(values = cols_pca)

pca_data %>% 
  rownames_to_column("sample_id") %>% 
  mutate(type = str_split(pattern = "_", sample_id) %>% map(1) %>% unlist) %>% 
  ggplot(aes(x = PC1, y = PC2, color = group)) +
  geom_point(aes(shape = type), size = 3) +
  geom_text(label = pca_data$name, size = 2, vjust = 3) +
  theme_minimal() +
  scale_colour_manual(values = cols_pca) +
  coord_fixed()
 
# 
# 
# for (t in condition_list) {
#   message("# ", t)
#   
#   pca_data_t <- pca_data %>% mutate(labelt = ifelse(group == t, name, ""),
#                                     sizet = ifelse(group == t, 4, 2))
#   
#   plot_t <- 
#   pca_data_t %>% 
#     ggplot(aes(x = PC1, y = PC2, color = group)) +
#     geom_point(size = pca_data_t$sizet) +
#     geom_text(label = pca_data_t$labelt, vjust = 2) +
#     theme_minimal() +
#     scale_colour_manual(values = cols_pca)
#   
#   print(plot_t)
#   
#   readline(prompt="Press [enter] to continue")
# }

