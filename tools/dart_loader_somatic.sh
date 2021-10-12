#!/bin/bash
#SBATCH --job-name=VCF_Loader
#SBATCH --output %j.log
#SBATCH --error %j.err
#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2G
#SBATCH --time=12:00:00
#SBATCH --nodes=1

set -e
set -x

#module load java/jdk/1.8.0/102

for i in "DART_PROJECT_NAME" "DART_USER" "DART_PWD" "DART_GRP"; do
  if [ -z "${!i}" ]; then
    echo "${i} name not set" 
    exit 1
  fi
done

#Type can be GERMLINE (cohort load), SOMATIC (multiple file load), SINGLE_FILE.
sample_name=$1
uid=$2
project_name="$DART_PROJECT_NAME"

#User setup
username="$DART_USER"
password="$DART_PWD"
usergroup="$DART_GRP"

# Fix the header of the VCF file, this is necessary because AS_FilterStatus and AS_SB_TABLE have wrong number values.
VCF_FILE="${BSF_PROJECTS}/${DART_PROJECT_NAME}/hg38/mutect2/${sample_name}.vep.vcf.gz"
bcftools view ${VCF_FILE} | sed 's/##INFO=<ID=AS_FilterStatus,Number=A/##INFO=<ID=AS_FilterStatus,Number=./' | sed  's/##INFO=<ID=AS_SB_TABLE,Number=1/##INFO=<ID=AS_SB_TABLE,Number=./' | bcftools view -O z -o ${VCF_FILE/.vep.vcf.gz/.tmp.vep.vcf.gz} -
bcftools index -t ${VCF_FILE/.vep.vcf.gz/.tmp.vep.vcf.gz}
[ ! -d "${BSF_PROJECTS}/${DART_PROJECT_NAME}/hg38/mutect2/broken_header" ] && mkdir "${BSF_PROJECTS}/${DART_PROJECT_NAME}/hg38/mutect2/broken_header"
mv "${BSF_PROJECTS}/${DART_PROJECT_NAME}/hg38/mutect2/${sample_name}.vep.vcf.gz" "${BSF_PROJECTS}/${DART_PROJECT_NAME}/hg38/mutect2/${sample_name}.vep.vcf.gz.tbi" "${BSF_PROJECTS}/${DART_PROJECT_NAME}/hg38/mutect2/broken_header"
mv ${VCF_FILE/.vep.vcf.gz/.tmp.vep.vcf.gz} ${VCF_FILE}
mv ${VCF_FILE/.vep.vcf.gz/.tmp.vep.vcf.gz.tbi} "${VCF_FILE}.tbi"

#Server setup
dart_server="pub-dart.int.cemm.at:8080"
config_file="${BSF_PROJECTS}/${DART_PROJECT_NAME}/hg38/mutect2/${sample_name}.DART_config.tsv"
echo -e "SAMPLE\tSAMPLE_ALIAS\tINCLUDE\tVCF_URL\tBAM_URL\tCOVERAGE_TRACK_URL\tCOVERAGE_FILE" > ${config_file}
VCF_URL="https://biomedical-sequencing.at/projects/${DART_PROJECT_NAME}_${uid}/hg38/mutect2/${sample_name}.vep.vcf.gz"
BAM_URL="https://biomedical-sequencing.at/projects/${DART_PROJECT_NAME}_${uid}/hg38/bam/${sample_name}.bam"
COVERAGE_TRACK_URL="https://biomedical-sequencing.at/projects/${DART_PROJECT_NAME}_${uid}/hg38/bam/${sample_name}_callable_loci.bb"
COVERAGE_FILE="${BSF_PROJECTS}/${DART_PROJECT_NAME}/hg38/bam/${sample_name}_non_callable_regions.tsv"
echo -e "${sample_name}\t${sample_name}\tTRUE\t${VCF_URL}\t${BAM_URL}\t${COVERAGE_TRACK_URL}\t${COVERAGE_FILE}" >> ${config_file}

chromosomes="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"

#Environment setup
jar_path="/nobackup/lab_bsf/applications/software/"
jar_name="DART-loader-1.1.1-SNAPSHOT-jar-with-dependencies.jar"
main_class="org.open.medgen.dart.loader.job.VCFLoader"

resource_path="/nobackup/lab_bsf/projects/${DART_PROJECT_NAME}/hg38/"
#vep_enums_kbl="${resource_path}VEPEnums.kbl_additions.tsv"
#vep_fields_kbl="${resource_path}VEPOutputFields.kbl_addition.tsv"

#-------------------------------------------------------------------------------

vcf_type="SOMATIC"
file=${VCF_FILE}

#-------------------------------------------------------------------------------


options=" -type "${vcf_type}
options=${options}" -config "${config_file}
options=${options}" -VCF "${file}
options=${options}" -group ${usergroup}"
options=${options}" -server "${username}":"${password}"@"${dart_server}
## add additions
options=${options}" -chrom ${chromosomes}"
#options=${options}" -vepOutputFields "${vep_fields_kbl}

classpath=${jar_path}"/"${jar_name}

java_opt="-Xmx20000m"
java_opt=${java_opt}" -cp "${classpath}

echo ""
echo "-------------LOADING ${file}-----------------------"
echo "JAVA PARAMS: "${java_opt}
echo "class: "${main_class}
echo "OPTS: "${options}
echo ""

java ${java_opt} ${main_class} ${options}

