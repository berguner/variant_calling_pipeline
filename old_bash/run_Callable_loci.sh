#!/bin/bash
#SBATCH --partition=longq
#SBATCH --ntasks=1
#SBATCH --mem=8000
#SBATCH --time=7-00:00:00


JAVATMP="/scratch/users/berguener/tmpdir/"
SAMPLE_NAME=$1
SAMPLE_BAM="$SAMPLE_NAME.bam"
INTERVAL_LIST="/data/groups/lab_bsf/resources/interval_lists/Ensembl_e99_Twist_RefSeq_exome_targets_hg38.interval_list"
BED="/data/groups/lab_bsf/resources/interval_lists/Ensembl_e99_Twist_RefSeq_exome_targets_hg38.bed"
REF_FASTA="/scratch/lab_bsf/resources/GATK/hg38/v0/Homo_sapiens_assembly38.fasta"
GTF="/data/groups/lab_bsf/resources/variant_calling_resources/Homo_sapiens.GRCh38.100.gtf.gz"
GATK_3_JAR="/scratch/lab_bsf/modules/applications/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar"

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

bedToBigBed -type=bed4 ${SAMPLE_NAME}_callable_sorted.bed $REF_FASTA.fai ${SAMPLE_NAME}_callable_loci.bb

~/src/variant_calling_pipeline/tools/bsf_variant_calling_coverage.R --callable-loci ${SAMPLE_NAME}_callable_loci.bed \
	--exon-flanks 15 \
	--exons $GTF \
	--targets $BED

date
