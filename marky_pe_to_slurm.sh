#!/bin/bash				
#SBATCH -A snic2021-5-53				
#SBATCH -n 12				
#SBATCH --time=10:00:00	
chmod +x workflow/genesearch.sh

sample="$1"

module load conda
source conda_init.sh
conda activate coco
snakemake --cores 12 ${sample}_tmp/${sample}_megahit
snakemake --cores 12 ${sample}_tmp/${sample}.index.1.bt2
snakemake --cores 12 ${sample}_outputs/${sample}_hgcA_final.txt
