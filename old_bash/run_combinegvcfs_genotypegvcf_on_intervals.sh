#!/bin/bash
#SBATCH --partition=longq
#SBATCH --cpus=2
#SBATCH --mem=16000
#SBATCH --time=7-00:00:00
#SBATCH --exclude=i0[01-22]

#source ~mschuster/src/bsfbash/bsf_software.bash
#module load java/jdk/1.8.0/102 ensembl/95

JAVATMP="/scratch/users/berguener/tmpdir/"
DBSNP_VCF="/scratch/lab_bsf/resources/GATK/hg38/v0/dbsnp_146.hg38.vcf.gz"
CHUNK=$1
GVCF_FOLDER="./gvcf/"
VCF_FOLDER="./vcf/"
OUT_GVCF_PATH="$GVCF_FOLDER/combined_gvcfs/$CHUNK.g.vcf.gz"
INTERVAL_LIST="/data/groups/lab_bsf/resources/interval_lists/hg38_ACGTmer/hg38_parts/hg38_${CHUNK}.interval_list"
REF_FASTA="/scratch/lab_bsf/resources/GATK/hg38/v0/Homo_sapiens_assembly38.fasta"


INPUT_GVCF="./gvcf/combined_gvcfs/$CHUNK.g.vcf.gz"
VCF_FOLDER="./vcf/"
OUT_VCF_PATH="$VCF_FOLDER/$CHUNK.vcf.gz"

# create necessary output folders
if [ ! -d "$GVCF_FOLDER/combined_gvcfs/" ]; then
  mkdir $GVCF_FOLDER/combined_gvcfs/
fi

if [ ! -d "$VCF_FOLDER" ]; then
  mkdir $VCF_FOLDER
fi

INPUT_GVCFS=""
for gvcf in $GVCF_FOLDER/*/*_${CHUNK}.g.vcf.gz;
do
	INPUT_GVCFS="$INPUT_GVCFS --variant $gvcf";
done

date
gatk --java-options "-Xmx12g -Xms128m -Djava.io.tmpdir:$JAVATMP -d64" \
    CombineGVCFs \
    -L $INTERVAL_LIST \
    -O $OUT_GVCF_PATH \
    -R $REF_FASTA \
    $INPUT_GVCFS > $OUT_GVCF_PATH.log 2>&1

gatk --java-options "-Xmx12g -Xms128m -Djava.io.tmpdir:$JAVATMP -d64" \
    GenotypeGVCFs \
    --dbsnp $DBSNP_VCF \
    -L $INTERVAL_LIST \
    -O $OUT_VCF_PATH \
    -R $REF_FASTA \
    -V $INPUT_GVCF > $OUT_VCF_PATH.log 2>&1

cd $VCF_FOLDER

sbatch --cpus=16 --mem=16000 --partition=shortq --error=${CHUNK}.vep.log --output=${CHUNK}.vep.log --job-name=${CHUNK}.vep \
	$ENSEMBL_HOME/ensembl-vep/vep \
	--allele_number --allow_non_variant --assembly GRCh38 \
	--cache --dir_cache $ENSEMBL_VEP_CACHE --offline \
	--dir_plugins $ENSEMBL_VEP_PLUGINS \
	--dont_skip --everything --failed 1 \
	--fasta $ENSEMBL_VEP_CACHE/homo_sapiens_merged/100_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
	--flag_pick_allele_gene --force_overwrite --format vcf \
	--gencode_basic --hgvsg --exclude_predicted \
	--plugin CADD,/scratch/lab_bsf/resources/CADD/hg38/1.6/whole_genome_SNVs.tsv.gz,/scratch/lab_bsf/resources/CADD/hg38/1.6/gnomad.genomes.r3.0.indel.tsv.gz \
	--species homo_sapiens --merged \
	--tmpdir ./ \
	--vcf --compress_output bgzip --fork 16 \
	--input_file ${CHUNK}.vcf.gz \
	--output_file ${CHUNK}.vep.vcf.gz \
	--stats_file ${CHUNK}.vep.html \
	--warning_file ${CHUNK}.vep.warning.txt

date
