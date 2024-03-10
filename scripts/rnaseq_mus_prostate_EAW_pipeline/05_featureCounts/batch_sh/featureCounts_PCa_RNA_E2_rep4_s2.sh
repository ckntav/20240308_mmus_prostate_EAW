#!/bin/sh
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --account=def-stbil30
#SBATCH --mail-user=christophe.tav@gmail.com
#SBATCH --mail-type=ALL


featureCounts -T 1 -p --countReadPairs -B -t exon -g gene_id -s 2 -a input/ensembl/Homo_sapiens.GRCh38.104.gtf -o /home/chris11/projects/def-stbil30/chris11/20240308_mmus_prostate_EAW/output/rna-pipeline_mmus_prostate_EAW-GRCh38_PE/featureCounts_s2/PCa_RNA_E2_rep4.featureCounts.tsv /home/chris11/projects/def-stbil30/chris11/20240308_mmus_prostate_EAW/output/rna-pipeline_mmus_prostate_EAW-GRCh38_PE/alignment/PCa_RNA_E2_rep4/PCa_RNA_E2_rep4.QueryNameSorted.bam
