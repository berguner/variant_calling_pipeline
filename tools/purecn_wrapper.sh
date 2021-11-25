#!/usr/bin/env bash
#SBATCH --partition=shortq
#SBATCH --qos=shortq
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16000
#SBATCH --time=12:00:00
#SBATCH --error ./purecn_%j.log
#SBATCH --output ./purecn_%j.log
#SBATCH --nodes=1
#SBATCH --ntasks=1


PROJECT_PATH=$1
SAMPLE_NAME=$2
CNVKIT_SIF="/nobackup/lab_bsf/users/berguener/singularity/cnvkit_0.9.8.sif"
BIOCONDUCTOR_SIF="/nobackup/lab_bsf/users/berguener/singularity/bsf_bioconductor_3.12.sif"
PURECN=`singularity exec -e --bind /nobackup:/nobackup $BIOCONDUCTOR_SIF Rscript -e 'cat(system.file("extdata", package="PureCN"))'`

date
echo "Sample: $SAMPLE_NAME"

singularity exec -e \
	--bind /nobackup:/nobackup \
	${CNVKIT_SIF} \
	cnvkit.py export seg ${PROJECT_PATH}/cnvkit/results/${SAMPLE_NAME}.cns \
	-o ${PROJECT_PATH}/cnvkit/results/${SAMPLE_NAME}.seg \
	--enumerate-chroms

singularity exec -e \
	--bind /nobackup:/nobackup \
	${BIOCONDUCTOR_SIF} \
	Rscript ${PURECN}/PureCN.R \
	--out ${PROJECT_PATH}/purecn/${SAMPLE_NAME} \
	--sampleid ${SAMPLE_NAME} \
	--tumor ${PROJECT_PATH}/cnvkit/results/${SAMPLE_NAME}.cnr \
	--segfile ${PROJECT_PATH}/cnvkit/results/${SAMPLE_NAME}.seg \
	--vcf ${PROJECT_PATH}/mutect2/${SAMPLE_NAME}.vep.vcf.gz \
	--genome hg38 --force --postoptimize --seed 123

date
