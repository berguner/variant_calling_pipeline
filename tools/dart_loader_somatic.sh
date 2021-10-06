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
project_type="GERMLINE"
project_name="$DART_PROJECT_NAME"

#User setup
username="$DART_USER"
password="$DART_PWD"
usergroup="$DART_GRP"

#Server setup
dart_server="pub-dart.int.cemm.at:8080"
config_file="${HOME}/src/bsfconfiguration/variant_calling/${DART_PROJECT_NAME}_variant_calling_DART_config.tsv"
chromosomes="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM"

#Environment setup
jar_path="/nobackup/lab_bsf/applications/software/"
jar_name="DART-loader-1.1.1-SNAPSHOT-jar-with-dependencies.jar"
main_class="org.open.medgen.dart.loader.job.VCFLoader"

resource_path="/nobackup/lab_bsf/projects/${DART_PROJECT_NAME}/hg38/"
#vep_enums_kbl="${resource_path}VEPEnums.kbl_additions.tsv"
#vep_fields_kbl="${resource_path}VEPOutputFields.kbl_addition.tsv"

#-------------------------------------------------------------------------------

vcf_type="GERMLINE"
file="/nobackup/lab_bsf/projects/${DART_PROJECT_NAME}/hg38/vcf/${project_name}_cohort.vep.vcf.gz"

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

