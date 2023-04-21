#!/bin/bash
# python virtual enviroment
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate wes

# specify the working dirs
maindir=/store/Pancancer_Metab/Public/PanCan/
bamdir=${maindir}/res/3_epi/RNAseq/STAR/
counts_dir=${maindir}/res/3_epi/RNAseq/counts/
logdir=${maindir}/res/3_epi/RNAseq/log/featureCounts

# creating dirs
mkdir -p ${counts_dir} ${logdir}

cut -f1 config_PANC1 | while read sample
do
    echo -e "______________________________\nfeatureCounts of ${sample} begin at `date`\n------------------------------\n"
	featureCounts -T 12 \
		-a /data/genome/GRCh38.GENCODE/gencode.v32.primary_assembly.annotation.gtf \
		-o ${counts_dir}/${sample}_featurecounts.txt \
		-p -g gene_id --extraAttributes gene_name \
		${bamdir}/${sample}_Aligned.out.bam \
		2> ${logdir}/${sample}_featurecounts.log
    echo -e "------------------------------\nfeatureCounts of ${sample} Over at `date`\n______________________________\n"
done

multiqc ${counts_dir}/*.summary -o ${counts_dir}/multiqc_featureCounts
