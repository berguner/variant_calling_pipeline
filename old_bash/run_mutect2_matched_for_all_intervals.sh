#!/bin/bash
#SBATCH --partition=longq
#SBATCH --cpus=8
#SBATCH --mem=64000
#SBATCH --time=7-00:00:00
#SBATCH --exclude=i0[01-22]

#source ~mschuster/src/bsfbash/bsf_software.bash
#module load java/jdk/1.8.0/102 parallel/20140122 


SAMPLE_NAME=$1
NORMAL_NAME=$2
MERGE_ONLY=$3
MUTECT_FOLDER="./mutect2"
SAMPLE_MUTECT_FOLDER="./mutect2/$SAMPLE_NAME"
HC_INTERVAL_SCRIPT="~/workspace/run_Mutect2_matched_on_intervals.sh"
GATK_LOCAL_JAR="/scratch/users/berguener/bin/gatk-4.1.8.1/gatk-package-4.1.8.1-local.jar"
JAVATMP="/scratch/users/berguener/tmpdir/"
REF_DICT="/scratch/lab_bsf/resources/GATK/hg38/v0/Homo_sapiens_assembly38.dict"
REF_FASTA="/scratch/lab_bsf/resources/GATK/hg38/v0/Homo_sapiens_assembly38.fasta"

# create necessary output folders
if [ ! -d "$MUTECT_FOLDER" ]; then
  mkdir $MUTECT_FOLDER
fi

# create necessary output folders
if [ ! -d "$SAMPLE_MUTECT_FOLDER" ]; then
  mkdir $SAMPLE_MUTECT_FOLDER
fi

for i in {1..30}; do
	STATS="$STATS -stats $SAMPLE_MUTECT_FOLDER/somatic_${i}.vcf.gz.stats"
	TUMOR_PILEUPS="$TUMOR_PILEUPS -I $SAMPLE_MUTECT_FOLDER/tumor-pileups_${i}.table"
	NORMAL_PILEUPS="$NORMAL_PILEUPS -I $SAMPLE_MUTECT_FOLDER/normal-pileups_${i}.table"
	F1R2="$F1R2 -I $SAMPLE_MUTECT_FOLDER/somatic_${i}.f1r2.tar.gz"
done

if [ $MERGE_ONLY != "merge" ]; then
date
echo {1..30} | tr ' ' '\n' | parallel -j 8 $HC_INTERVAL_SCRIPT $SAMPLE_NAME $NORMAL_NAME
date
fi

COUNT=$(grep "SUCCESS" $SAMPLE_MUTECT_FOLDER/*vcf.gz.log | wc -l)
if [ "$COUNT" -eq "30" ]; then
  echo "All 30 intervals were completed successfully, now concatenating the VCF files"
  VARIANTS=""
  for i in {1..30}; do VARIANTS="$VARIANTS $SAMPLE_MUTECT_FOLDER/somatic_${i}.vcf.gz"; done
  bcftools concat --threads 8 -o $MUTECT_FOLDER/$SAMPLE_NAME.raw.vcf.gz -O z $VARIANTS
  bcftools index -t $MUTECT_FOLDER/$SAMPLE_NAME.raw.vcf.gz
  echo "Finished concatenating vcf files"
  
  echo "Merging stats"
  gatk --java-options "-Xmx6g -Xms128m -Djava.io.tmpdir:$JAVATMP -d64" \
	MergeMutectStats $STATS -O $SAMPLE_MUTECT_FOLDER/mutect2.stats > $SAMPLE_MUTECT_FOLDER/mutect2.stats.log 2>&1
  
  echo "Merging pile-ups"
  gatk --java-options "-Xmx6g -Xms128m -Djava.io.tmpdir:$JAVATMP -d64" \
	GatherPileupSummaries --sequence-dictionary $REF_DICT \
	$TUMOR_PILEUPS \
	-O $SAMPLE_MUTECT_FOLDER/tumor_pileup.tsv > $SAMPLE_MUTECT_FOLDER/tumor_pileup.tsv.log 2>&1 
  gatk --java-options "-Xmx6g -Xms128m -Djava.io.tmpdir:$JAVATMP -d64" \
	GatherPileupSummaries --sequence-dictionary $REF_DICT \
	$NORMAL_PILEUPS \
	-O $SAMPLE_MUTECT_FOLDER/normal_pileup.tsv > $SAMPLE_MUTECT_FOLDER/normal_pileup.tsv.log 2>&1 
  
  echo "Running CalculateContamination"
  gatk --java-options "-Xmx6g -Xms128m -Djava.io.tmpdir:$JAVATMP -d64" \
	CalculateContamination \
	-I $SAMPLE_MUTECT_FOLDER/tumor_pileup.tsv \
	-matched $SAMPLE_MUTECT_FOLDER/normal_pileup.tsv \
	--tumor-segmentation $SAMPLE_MUTECT_FOLDER/segments.table \
	-O $SAMPLE_MUTECT_FOLDER/contamination.table > $SAMPLE_MUTECT_FOLDER/contamination.table.log 2>&1
  
  echo "Running LearnReadOrientationModel"
  gatk --java-options "-Xmx6g -Xms128m -Djava.io.tmpdir:$JAVATMP -d64" \
	LearnReadOrientationModel $F1R2 \
	-O $SAMPLE_MUTECT_FOLDER/artifact-priors.tar.gz > $SAMPLE_MUTECT_FOLDER/artifact-priors.tar.gz.log 2>&1 

  echo "Filtering"
  gatk --java-options "-Xmx6g -Xms128m -Djava.io.tmpdir:$JAVATMP -d64" \
	FilterMutectCalls \
	-V $MUTECT_FOLDER/$SAMPLE_NAME.raw.vcf.gz \
	-R $REF_FASTA \
	-O $MUTECT_FOLDER/$SAMPLE_NAME.vcf.gz \
	--contamination-table $SAMPLE_MUTECT_FOLDER/contamination.table \
	--tumor-segmentation $SAMPLE_MUTECT_FOLDER/segments.table \
	--ob-priors $SAMPLE_MUTECT_FOLDER/artifact-priors.tar.gz \
	-stats $SAMPLE_MUTECT_FOLDER/mutect2.stats \
	--filtering-stats $MUTECT_FOLDER/$SAMPLE_NAME.filtering.stats > $MUTECT_FOLDER/$SAMPLE_NAME.filtering.stats.log 2>&1
else
  echo "Only $COUNT of the 30 intervals were completed, exiting";
  exit 1;
fi

perl $ENSEMBL_HOME/ensembl-vep/vep \
	--allele_number --allow_non_variant --assembly GRCh38 \
	--cache --dir_cache $ENSEMBL_VEP_CACHE --offline \
	--dir_plugins $ENSEMBL_VEP_PLUGINS \
	--dont_skip --everything --failed 1 \
	--fasta $ENSEMBL_VEP_CACHE/homo_sapiens_merged/100_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
	--flag_pick_allele_gene --force_overwrite --format vcf \
	--gencode_basic --hgvsg --exclude_predicted \
	--plugin CADD,/scratch/lab_bsf/resources/CADD/hg38/1.6/whole_genome_SNVs.tsv.gz,/scratch/lab_bsf/resources/CADD/hg38/1.6/gnomad.genomes.r3.0.indel.tsv.gz \
	--species homo_sapiens --merged \
	--tmpdir $SAMPLE_MUTECT_FOLDER \
	--vcf --compress_output bgzip --fork 8 \
	--input_file $MUTECT_FOLDER/$SAMPLE_NAME.vcf.gz \
	--output_file $MUTECT_FOLDER/$SAMPLE_NAME.vep.vcf.gz \
	--stats_file $MUTECT_FOLDER/$SAMPLE_NAME.vep.html \
	--warning_file $MUTECT_FOLDER/$SAMPLE_NAME.vep.warning_file.txt > $MUTECT_FOLDER/$SAMPLE_NAME.vep.log 2>&1

bcftools index -t $MUTECT_FOLDER/$SAMPLE_NAME.vep.vcf.gz

date
