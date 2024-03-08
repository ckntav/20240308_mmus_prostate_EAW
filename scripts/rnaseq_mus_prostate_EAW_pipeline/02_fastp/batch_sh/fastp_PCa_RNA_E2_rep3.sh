#!/bin/sh
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --account=def-stbil30
#SBATCH --mail-user=christophe.tav@gmail.com
#SBATCH --mail-type=ALL


mkdir -p /home/chris11/projects/def-stbil30/chris11/20240308_mmus_prostate_EAW/output/rna-pipeline_mmus_prostate_EAW-GRCh38_PE/fastp_output/PCa_RNA_E2_rep3


mkdir -p /home/chris11/projects/def-stbil30/chris11/20240308_mmus_prostate_EAW/raw/rnaseq_mmus_prostate_EAW/fastp_output


/cvmfs/soft.computecanada.ca/easybuild/software/2020/avx2/Core/fastp/0.23.4/bin/fastp --in1 /home/chris11/projects/def-stbil30/chris11/20240308_mmus_prostate_EAW/raw/rnaseq_mmus_prostate_EAW/raw_fastq/PCa_E2-3_R1.fastq.gz --in2 /home/chris11/projects/def-stbil30/chris11/20240308_mmus_prostate_EAW/raw/rnaseq_mmus_prostate_EAW/raw_fastq/PCa_E2-3_R2.fastq.gz --detect_adapter_for_pe --overrepresentation_analysis --overrepresentation_sampling 10 --thread 8 --length_required 10 --out1 /home/chris11/projects/def-stbil30/chris11/20240308_mmus_prostate_EAW/raw/rnaseq_mmus_prostate_EAW/fastp_output/PCa_RNA_E2_rep3_1.fastq.gz --out2 /home/chris11/projects/def-stbil30/chris11/20240308_mmus_prostate_EAW/raw/rnaseq_mmus_prostate_EAW/fastp_output/PCa_RNA_E2_rep3_2.fastq.gz --html /home/chris11/projects/def-stbil30/chris11/20240308_mmus_prostate_EAW/output/rna-pipeline_mmus_prostate_EAW-GRCh38_PE/fastp_output/PCa_RNA_E2_rep3/PCa_RNA_E2_rep3_fastp_report.html --json /home/chris11/projects/def-stbil30/chris11/20240308_mmus_prostate_EAW/output/rna-pipeline_mmus_prostate_EAW-GRCh38_PE/fastp_output/PCa_RNA_E2_rep3/PCa_RNA_E2_rep3_fastp_report.json
