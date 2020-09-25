#!/bin/bash
#SBATCH --partition=longq
#SBATCH --ntasks=1
#SBATCH --mem=8000
#SBATCH --time=7-00:00:00
#SBATCH --exclude=i0[01-22]

#source ~mschuster/src/bsfbash/bsf_software.bash
#module load java/jdk/1.8.0/102

JAVATMP="/scratch/users/berguener/tmpdir/"
SAMPLE_NAME=$1
NORMAL_NAME=$2
CHUNK=$3
SAMPLE_BAM="./bam/$SAMPLE_NAME.bam"
NORMAL_BAM="./bam/$NORMAL_NAME.bam"
SAMPLE_MUTECT_FOLDER="./mutect2/$SAMPLE_NAME"
SAMPLE_MUTECT_VCF="$SAMPLE_MUTECT_FOLDER/somatic_${CHUNK}.vcf.gz"
INTERVAL_LIST="/data/groups/lab_bsf/resources/interval_lists/hg38_ACGTmer/hg38_parts/hg38_${CHUNK}.interval_list"
REF_FASTA="/scratch/lab_bsf/resources/GATK/hg38/v0/Homo_sapiens_assembly38.fasta"
GATK_LOCAL_JAR="/scratch/users/berguener/bin/gatk-4.1.8.1/gatk-package-4.1.8.1-local.jar"
GERMLINE_RESOURCE="/scratch/lab_bsf/resources/GATK/hg38/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz"
VARIANTS_FOR_CONTAMINATION="/scratch/lab_bsf/resources/GATK/hg38/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz"


# create necessary output folders
if [ ! -d "$SAMPLE_MUTECT_FOLDER" ]; then
  mkdir $SAMPLE_MUTECT_FOLDER
fi

date
gatk --java-options "-Xmx6g -Xms128m -Djava.io.tmpdir:$JAVATMP -d64" \
    Mutect2 \
    -L $INTERVAL_LIST \
    -O $SAMPLE_MUTECT_VCF \
    -R $REF_FASTA \
	--genotype-germline-sites true \
	--f1r2-tar-gz $SAMPLE_MUTECT_FOLDER/somatic_${CHUNK}.f1r2.tar.gz \
	--germline-resource $GERMLINE_RESOURCE \
    -I $SAMPLE_BAM -tumor $SAMPLE_NAME \
	-I $NORMAL_BAM -normal $NORMAL_NAME > $SAMPLE_MUTECT_VCF.log 2>&1

gatk --java-options "-Xmx6g -Xms128m -Djava.io.tmpdir:$JAVATMP -d64" \
	GetPileupSummaries \
	-R $REF_FASTA \
	-I $SAMPLE_BAM \
	--interval-set-rule INTERSECTION -L $INTERVAL_LIST \
	-V $VARIANTS_FOR_CONTAMINATION \
	-L $VARIANTS_FOR_CONTAMINATION \
	-O $SAMPLE_MUTECT_FOLDER/tumor-pileups_${CHUNK}.table > $SAMPLE_MUTECT_FOLDER/tumor-pileups_${CHUNK}.table.log 2>&1

gatk --java-options "-Xmx6g -Xms128m -Djava.io.tmpdir:$JAVATMP -d64" \
	GetPileupSummaries \
	-R $REF_FASTA \
	-I $NORMAL_BAM \
	--interval-set-rule INTERSECTION -L $INTERVAL_LIST \
	-V $VARIANTS_FOR_CONTAMINATION \
	-L $VARIANTS_FOR_CONTAMINATION \
	-O $SAMPLE_MUTECT_FOLDER/normal-pileups_${CHUNK}.table > $SAMPLE_MUTECT_FOLDER/normal-pileups_${CHUNK}.table.log 2>&1

date
