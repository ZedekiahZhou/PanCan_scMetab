#!/bin/bash
# python virtual enviroment
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate wes

# specify the working dirs
maindir=/store/Pancancer_Metab/Public/PanCan/
indir=${maindir}/data/TF_RNAseq/clean/
outdir=${maindir}/res/3_epi/RNAseq/STAR/
# logdir=${maindir}/res/3_epi/RNAseq/log

# creating dirs
mkdir -p ${outdir}

cut -f1 config_PANC1 | while read sample
do
	echo -e "______________________________\nSTAR aligning of ${sample} begin at `date`\n------------------------------\n"
    STAR --runThreadN 10 --genomeDir /data/genome/GRCh38.GENCODE/STAR.index \
	--readFilesIn ${indir}/${sample}/${sample}_R1.fq.gz ${indir}/${sample}/${sample}_R2.fq.gz \
    --readFilesCommand zcat --outFileNamePrefix ${outdir}/${sample}_ \
    --outSAMtype BAM Unsorted \
    --outBAMsortingThreadN 10 \
    --outSAMunmapped Within \
	--quantMode GeneCounts \
	--sjdbGTFfile /data/genome/GRCh38.GENCODE/gencode.v32.primary_assembly.annotation.gtf	  
    # SortedByCoordinate bam可以取消
	echo -e "------------------------------\nSTAR aligning of ${sample} Over at `date`\n______________________________\n"
done

multiqc ${outdir}/*Log.final.out -o ${outdir}/multiqc_STAR