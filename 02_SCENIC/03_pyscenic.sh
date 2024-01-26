#!/bin/bash
#SBATCH -o job.pyscenic.%j.out
#SBATCH --partition=C032M0128G
#SBATCH --qos=low
#SBATCH -J pyscenic
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=30

proj=("EMBO_2021_Pal_BC" "NG_2021_Wu_Spatial" "WangShu_p_EMBO_BRCA" "NC_2020_Kim_LUAD" "NatMed_2018_Lam_lung" "WangGX_Lung" "Cell_2021_Pelka_colon" "CellRes_2019_Peng_PDAC" "CancerDiscov_2022_Kumar_STAD" "NC_2022_Song_PRAD")
cancers=("BRCA" "BRCA_valid1" "BRCA_valid2" "LUAD" "LUAD_valid1" "LUAD_valid2" "COAD" "PAAD" "STAD" "PRAD")

# load singularity modulefile
module load singularity/3.5.2

for i in {0..8};
do
        # cd ../${proj[${i}]}
        echo "Begin ${proj[${i}]} at `date`:"
        pwd

        # resource files
        f_tfs="hs_hgnc_tfs_symbols_updated.txt"
        f_motif="motifs-v9-nr.hgnc-m0.001-o0.0-symbols-updated.tbl"
        f_db=`ls hg19-*-symbols-updated.mc9nr.feather`

        # input files
        f_loom_expr_mat="${cancers[${i}]}_filtered_scenic.loom"
        if [ -f ${f_loom_expr_mat} ];
      then
                echo "${f_loom_expr_mat} exist."
        fi

        # output files
        f_adj="${cancers[${i}]}_adj.csv"
        f_reg="${cancers[${i}]}_reg.csv"

        singularity run aertslab-pyscenic-0.11.2.sif \
                pyscenic grn \
                    --num_workers 20 \
                -o ${f_adj} \
                ${f_loom_expr_mat} \
                ${f_tfs}

        singularity run ../aertslab-pyscenic-0.11.2.sif \
            pyscenic ctx ${f_adj} \
            ${f_db} \
            --annotations_fname ${f_motif} \
            --expression_mtx_fname ${f_loom_expr_mat} \
            --output ${f_reg} \
            --mask_dropouts \
            --num_workers 20

        echo "Over ${proj[${i}]} at `date`:"
        echo ""
        echo ""
        echo ""
done
