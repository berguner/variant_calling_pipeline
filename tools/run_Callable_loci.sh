#!/bin/bash
#SBATCH --partition=mediumq
#SBATCH --qos=mediumq
#SBATCH -c 2
#SBATCH --mem=16000
#SBATCH --time=2-00:00:00

module unload java
module load Java/1.8

JAVATMP="/nobackup/lab_bsf/users/berguener/tmpdir/"
SAMPLE_NAME=$1
SAMPLE_BAM="$SAMPLE_NAME.bam"
INTERVAL_LIST="/research/lab_bsf/resources/intervals/Ensembl_e99_Twist_RefSeq_exome_targets_hg38.interval_list"
BED="/research/lab_bsf/resources/intervals/Ensembl_e99_Twist_RefSeq_exome_targets_hg38.bed"
REF_FASTA="/nobackup/lab_bsf/resources/GATK/hg38/v0/Homo_sapiens_assembly38.fasta"
GTF="/research/lab_bsf/resources/variant_calling_resources/Homo_sapiens.GRCh38.100.gtf.gz"
GATK_3_JAR="/research/lab_bsf/resources/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.jar"

date

java -Djava.io.tmpdir=$JAVATMP -Xmx3G \
	-jar $GATK_3_JAR --analysis_type CallableLoci \
	--input_file $SAMPLE_BAM \
	--intervals $INTERVAL_LIST \
	--minDepth 10 \
	--out ${SAMPLE_NAME}_callable_loci.bed \
	--reference_sequence $REF_FASTA \
	--summary ${SAMPLE_NAME}_callable_loci.txt \
	> ${SAMPLE_NAME}_callable_loci.log 2>&1

bedSort ${SAMPLE_NAME}_callable_loci.bed  ${SAMPLE_NAME}_callable_sorted.bed

bedToBigBed -type=bed4 ${SAMPLE_NAME}_callable_sorted.bed ${REF_FASTA}.fai ${SAMPLE_NAME}_callable_loci.bb

~/src/variant_calling_pipeline/tools/bsf_variant_calling_coverage.R --callable-loci ${SAMPLE_NAME}_callable_loci.bed \
	--exon-flanks 15 \
	--exons $GTF \
	--targets $BED

date
