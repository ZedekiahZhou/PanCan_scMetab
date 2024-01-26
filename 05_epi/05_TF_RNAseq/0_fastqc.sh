#!/bin/bash
# python virtual enviroment
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate wes

# specify the working dirs
maindir=/store/Pancancer_Metab/Public/PanCan/
indir=${maindir}/data/TF_RNAseq/
fastqc_raw=${maindir}/res/3_epi/RNAseq/fastqc/raw
fastqc_clean=${maindir}/res/3_epi/RNAseq/fastqc/clean

# creating dirs
mkdir -p ${fastqc_raw} ${fastqc_clean}

cut -f1 config_PANC1 | while read sample
do
	# fastqc on raw seq file
	echo -e "______________________________\nFastqc on Raw Data for ${sample} begin at `date`\n------------------------------\n"
	fastqc ${indir}/raw/${sample}/${sample}_R1.fq.gz ${indir}/raw/${sample}/${sample}_R2.fq.gz \
		--threads 12 -o ${fastqc_raw} 2>${fastqc_raw}/${sample}_fastqc.log 2>&1
	echo -e "------------------------------\nFastqc Over at `date`\n______________________________\n"

	# fastqc on clean seq file
    echo -e "______________________________\nFastqc on Clean Data for ${sample} begin at `date`\n------------------------------\n"
    fastqc ${indir}/clean/${sample}/${sample}_R1.fq.gz ${indir}/clean/${sample}/${sample}_R2.fq.gz \
		--threads 12 -o ${fastqc_clean} > ${fastqc_clean}/${sample}_fastqc.log 2>&1
    echo -e "------------------------------\nFastqc Over at `date`\n______________________________\n"
done

multiqc ${fastqc_raw}/*.zip -o ${fastqc_raw}/multiqc
multiqc ${fastqc_clean}/*.zip -o ${fastqc_clean}/multiqc

