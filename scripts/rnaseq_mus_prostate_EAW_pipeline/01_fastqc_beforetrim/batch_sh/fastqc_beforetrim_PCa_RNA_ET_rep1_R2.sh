#!/bin/sh
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --account=def-stbil30
#SBATCH --mail-user=christophe.tav@gmail.com
#SBATCH --mail-type=ALL


mkdir -p /home/chris11/projects/def-stbil30/chris11/20240308_mmus_prostate_EAW/output/rna-pipeline_mmus_prostate_EAW-GRCh38_PE/fastqc_beforetrim_output/PCa_RNA_ET_rep1_R2


/cvmfs/soft.computecanada.ca/easybuild/software/2020/Core/fastqc/0.11.9/fastqc --outdir /home/chris11/projects/def-stbil30/chris11/20240308_mmus_prostate_EAW/output/rna-pipeline_mmus_prostate_EAW-GRCh38_PE/fastqc_beforetrim_output/PCa_RNA_ET_rep1_R2 --format fastq /home/chris11/projects/def-stbil30/chris11/20240308_mmus_prostate_EAW/raw/rnaseq_mmus_prostate_EAW/raw_fastq/PCa_ET-1_R2.fastq.gz
