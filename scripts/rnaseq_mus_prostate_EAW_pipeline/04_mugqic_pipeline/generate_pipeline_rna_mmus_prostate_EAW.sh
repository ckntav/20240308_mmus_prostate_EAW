export RAP_ID=def-stbil30
export JOB_MAIL="christophe.tav@gmail.com"

mkdir -p $HOME/projects/def-stbil30/chris11/20240308_mmus_prostate_EAW/output/rna-pipeline_mmus_prostate_EAW-GRCh38_PE

rnaseq.py --job-scheduler slurm -s 1-13 \
  --log debug \
  --readsets raw/rnaseq_mmus_prostate_EAW/readset_rnaseq_mmus_prostate_EAW_20240308.txt \
  -o output/rna-pipeline_mmus_prostate_EAW-GRCh38_PE \
  --no-json \
  --config $MUGQIC_PIPELINES_HOME/pipelines/rnaseq/rnaseq.base.ini \
      $MUGQIC_PIPELINES_HOME/pipelines/common_ini/narval.ini \
      $MUGQIC_INSTALL_HOME/genomes/species/Mus_musculus.GRCm38/Mus_musculus.GRCm38.ini