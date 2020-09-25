#!/bin/bash
#SBATCH --partition=longq
#SBATCH --cpus=8
#SBATCH --mem=64000
#SBATCH --time=7-00:00:00
#SBATCH --exclude=i0[01-22]

#source ~mschuster/src/bsfbash/bsf_software.bash
#module load java/jdk/1.8.0/102 parallel/20140122

SAMPLE_NAME=$1
HC_threads=8
JAVATMP="/scratch/users/berguener/tmpdir/"
DBSNP_VCF="/scratch/lab_bsf/resources/GATK/hg38/v0/dbsnp_146.hg38.vcf.gz"
SAMPLE_BAM="./bam/$SAMPLE_NAME.bam"
SAMPLE_GVCF_FOLDER="./gvcf/$SAMPLE_NAME"
REF_FASTA="/scratch/lab_bsf/resources/GATK/hg38/v0/Homo_sapiens_assembly38.fasta"

# create necessary output folders
if [ ! -d "./gvcf" ]; then
  mkdir ./gvcf
fi

# create necessary output folders
if [ ! -d "$SAMPLE_GVCF_FOLDER" ]; then
  mkdir $SAMPLE_GVCF_FOLDER
fi

for CHUNK in {1..30}; do
	SAMPLE_INTERVAL_GVCF="./gvcf/$SAMPLE_NAME/${SAMPLE_NAME}_${CHUNK}.g.vcf.gz"
	INTERVAL_LIST="/data/groups/lab_bsf/resources/interval_lists/hg38_ACGTmer/hg38_parts/hg38_${CHUNK}.interval_list"
    echo -e "gatk \
    --java-options \"-Xmx6g -Xms128m -Djava.io.tmpdir:$JAVATMP -d64\" \
    HaplotypeCaller \
    --dbsnp $DBSNP_VCF \
    -L $INTERVAL_LIST \
    -O $SAMPLE_INTERVAL_GVCF \
    -R $REF_FASTA \
    -ERC GVCF \
    -I $SAMPLE_BAM > $SAMPLE_INTERVAL_GVCF.log 2>&1";
done | parallel --no-notice -j ${HC_threads} ::: > $SAMPLE_GVCF_FOLDER/HC.parallel.log 2>&1;
date
