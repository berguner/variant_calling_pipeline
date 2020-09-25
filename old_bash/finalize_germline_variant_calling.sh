#!/bin/bash
#SBATCH --partition=mediumq
#SBATCH --cpus=4
#SBATCH --mem=16000
#SBATCH --time=2-00:00:00

#source ~mschuster/src/bsfbash/bsf_software.bash
#module load java/jdk/1.8.0/102 ensembl/95

JAVATMP="/scratch/users/berguener/tmpdir/"
VCF_FOLDER="./vcf/"
REF_FASTA="/scratch/lab_bsf/resources/GATK/hg38/v0/Homo_sapiens_assembly38.fasta"
DBSNP_VCF="/scratch/lab_bsf/resources/GATK/hg38/v0/dbsnp_146.hg38.vcf.gz"

PROJECT_NAME=$1
COHORT_VCF="${PROJECT_NAME}_cohort.vcf.gz"
COHORT_VEP_VCF="${PROJECT_NAME}_cohort.vep.vcf.gz"

cd $VCF_FOLDER

INPUT_VCFS=""
INPUT_VEP_VCFS=""
for i in {1..30};
do
	INPUT_VCFS="$INPUT_VCFS $i.vcf.gz";
	bcftools index -t $i.vep.vcf.gz
	INPUT_VEP_VCFS="$INPUT_VEP_VCFS $i.vep.vcf.gz";
done

bcftools concat -a -Oz -o ${COHORT_VCF} ${INPUT_VCFS}
bcftools index -t ${COHORT_VCF}

bcftools concat -a -Oz -o ${COHORT_VEP_VCF} ${INPUT_VEP_VCFS}
bcftools index -t ${COHORT_VEP_VCF}

SAMPLE_LIST=`bcftools query -l ${COHORT_VCF}`
for sample in ${SAMPLE_LIST};
do
	bcftools view -Oz -o ${sample}.vcf.gz -s ${sample} ${COHORT_VCF}
	bcftools index -t ${sample}.vcf.gz
	bcftools view -Oz -o ${sample}.vep.vcf.gz -s ${sample} ${COHORT_VEP_VCF}
	bcftools index -t ${sample}.vep.vcf.gz
	sbatch -J ${sample}.VCM \
		--error=${sample}.vcf.gz.variant_calling_metrics.log \
		--output=${sample}.vcf.gz.variant_calling_metrics.log \
		--partition=shortq --mem=8000 --cpus=2 \
		--wrap="gatk --java-options \"-Xmx4g -Xms128m -Djava.io.tmpdir:$JAVATMP -d64\" \
		CollectVariantCallingMetrics \
		--INPUT ${sample}.vcf.gz \
		--OUTPUT ${sample} \
		--DBSNP ${DBSNP_VCF}"
done

date
