#!/bin/bash
#SBATCH -o job.pyscenic.%j.out
#SBATCH --partition=C032M0128G
#SBATCH --qos=low
#SBATCH -J pyscenic
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=30

cancers=("BRCA" "BRCA_valid1" "BRCA_valid2" "LUAD" "LUAD_valid1" "LUAD_valid2" "COAD" "PAAD" "STAD" "PRAD")

# load singularity modulefile
module load singularity/3.5.2

for i in {0..9};
do
    echo "Begin ${cancers[$i]} at `date`:"

    f_loom_path_scenic="data/${cancers[$i]}_filtered_whole.loom"
    f_TF_reg_path="data/${cancers[$i]}_reg_metab.gmt"
    f_TF_aucell_path="res/${cancers[$i]}_TF_aucell.csv"

    f_metab_reg_path="data/gaudeMetab.gmt"
    f_metab_aucell_path="res/${cancers[$i]}_metab_aucell.csv"
    
    singularity run aertslab-pyscenic-0.11.2.sif \
        pyscenic aucell \
            ${f_loom_path_scenic} \
            ${f_TF_reg_path} \
            --output ${f_TF_aucell_path} \
            --num_workers 20
    
    singularity run aertslab-pyscenic-0.11.2.sif \
        pyscenic aucell \
            ${f_loom_path_scenic} \
            ${f_metab_reg_path} \
            --output ${f_metab_aucell_path} \
            --num_workers 20

    echo "Over ${cancers[$i]} at `date`:"
    echo ""
    echo ""
    echo ""
done