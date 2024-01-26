#!/bin/bash
# python virtual enviroment
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate wes

# specify the working dirs
maindir=/store/Pancancer_Metab/Public/PanCan/
bamdir=${maindir}/res/3_epi/RNAseq/STAR/
qualimap_dir=${maindir}/res/3_epi/RNAseq/qualimap
logdir=${maindir}/res/log

# creating dirs
mkdir -p ${qualimap_dir}

cut -f1 config_PANC1 | while read sample
do
    echo -e "______________________________\nQualimap of ${sample} begin at `date`\n------------------------------\n"
	qualimap rnaseq \
		-outdir ${qualimap_dir}/${sample} \
		-a proportional \
		-bam ${bamdir}/${sample}_Aligned.out.bam \
		-p non-strand-specific \
		-gtf /data/genome/GRCh38.GENCODE/gencode.v32.primary_assembly.annotation.gtf \
		--java-mem-size=32G \
		-pe
    echo -e "------------------------------\nQualimap of ${sample} Over at `date`\n______________________________\n"
done

multiqc ${qualimap_dir}/*/rnaseq_qc_results.txt -o ${qualimap_dir}/multiqc_qualimap