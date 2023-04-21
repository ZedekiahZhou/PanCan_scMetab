#!/bin/bash
# python virtual enviroment
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate wes

# specify the working dirs
maindir=/store/Pancancer_Metab/Public/PanCan/
indir=${maindir}/data/TF_RNAseq/clean/
outdir=${maindir}/res/3_epi/RNAseq/salmon
# logdir=${maindir}/res/log

# creating dirs
mkdir -p ${outdir}

cut -f1 config_PANC1 | while read sample
do
	echo -e "______________________________\nSalmon of ${sample} begin at `date`\n------------------------------\n"
    salmon quant -i /data/genome/GRCh38.GENCODE/transcriptome/transcripts_index \
	-l A -1 ${indir}/${sample}/${sample}_R1.fq.gz -2 ${indir}/${sample}/${sample}_R2.fq.gz \
	-p 10 -o ${outdir}/${sample}_quant \
	--seqBias --gcBias
	echo -e "------------------------------\nSalmon of ${sample} Over at `date`\n______________________________\n"
done