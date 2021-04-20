#!/bin/bash
#SBATCH --partition=mediumq
#SBATCH --qos=mediumq
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=128000
#SBATCH --time=2-00:00:00
#SBATCH --nodes=1

#source ~mschuster/src/bsfbash/bsf_software.bash
#module unload java
#module load java/jdk/1.8.0/102 parallel

date
THREADS=24
ARGS=("$@")
PICARD="/nobackup/lab_bsf/applications/software/picard/2.19.2_cemm/picard.jar"
REF="/nobackup/lab_bsf/resources/GATK/hg38/v0/Homo_sapiens_assembly38.fasta"
JAVATMP="/nobackup/lab_bsf/users/berguener/tmpdir/"
BAM="./bam"
RAW="./raw"
SAMPLE_NAME=${1}

if [ ! -f "$BAM/$SAMPLE_NAME.bam.bai" ] ;
then
RAW_BAMS=$(ls $RAW/*$SAMPLE_NAME*.bam)
echo "Raw bam files: $RAW_BAMS"
echo "Trimming adapters"
RG="@RG\tID:$SAMPLE_NAME\tSM:$SAMPLE_NAME\tPL:illumina\tLB:$SAMPLE_NAME"
echo "Running BWA mapping"
for i in $RAW_BAMS;
do
	samtools fastq $i;
done | fastp --stdin --interleaved_in --thread $THREADS --stdout --html $BAM/$SAMPLE_NAME.fastp.html --json $BAM/$SAMPLE_NAME.fastp.json --verbose 2> $BAM/$SAMPLE_NAME.fastp.err | \
	bwa mem -t $THREADS -R $RG -p $REF - 2> $BAM/$SAMPLE_NAME.bwa.log | \
	samtools fixmate -O SAM - - 2> $BAM/$SAMPLE_NAME.samtools_fixmate.log | \
	samblaster 2> $BAM/$SAMPLE_NAME.samblaster.log | \
	samtools sort -m 512m -o $BAM/$SAMPLE_NAME.bam -@ 4 - 2> $BAM/$SAMPLE_NAME.samtools_sort.log
samtools index -@ $THREADS $BAM/$SAMPLE_NAME.bam

echo "Plotting fragment size distribution"
samtools view -f66 -F3840 $BAM/$SAMPLE_NAME.bam chr21 | awk '{ if($9<0) dist=-$9; else dist=$9; if(dist<1000){ counts[dist] = counts[dist] + 1; sum = sum + 1; } } END {for (word in counts) {print word, (counts[word]/sum)*100;} }' | sort -n > $BAM/$SAMPLE_NAME.fragsize_hist
gnuplot <<EOF
set terminal png size 800,600
set xlabel "Fragment Size"
set ylabel "Percentage"
set output "$BAM/$SAMPLE_NAME.fragsize.png"
plot "$BAM/$SAMPLE_NAME.fragsize_hist" with lines
EOF
fi

PICARD_CMD=""
echo "Collecting Picard insert size metrics"
PICARD_CMD="${PICARD_CMD}java -Xmx2g -Djava.io.tmpdir=$JAVATMP -d64 -server -jar $PICARD CollectInsertSizeMetrics I=$BAM/$SAMPLE_NAME.bam O=$BAM/$SAMPLE_NAME.insert_size_metrics.tsv H=$BAM/$SAMPLE_NAME.insert_size_metrics.pdf > $BAM/$SAMPLE_NAME.insert_size_metrics.log 2>&1\n"
echo "Collecting Picard alignment summary metrics"
PICARD_CMD="${PICARD_CMD}java -Xmx2g -Djava.io.tmpdir=$JAVATMP -d64 -server -jar $PICARD CollectAlignmentSummaryMetrics R=$REF I=$BAM/$SAMPLE_NAME.bam O=$BAM/$SAMPLE_NAME.alignment_summary_metrics.tsv > $BAM/$SAMPLE_NAME.alignment_summary_metrics.log 2>&1\n"
echo "Collecting Picard WGS metrics"
PICARD_CMD="${PICARD_CMD}java -Xmx2g -Djava.io.tmpdir=$JAVATMP -d64 -server -jar $PICARD CollectWgsMetrics R=$REF I=$BAM/$SAMPLE_NAME.bam O=$BAM/$SAMPLE_NAME.wgs_metrics.tsv > $BAM/$SAMPLE_NAME.wgs_metrics.log 2>&1\n"
echo "Collecting Picard GC bias metrics"
PICARD_CMD="${PICARD_CMD}java -Xmx2g -Djava.io.tmpdir=$JAVATMP -d64 -server -jar $PICARD CollectGcBiasMetrics R=$REF I=$BAM/$SAMPLE_NAME.bam O=$BAM/$SAMPLE_NAME.gc_bias_metrics.tsv CHART=$BAM/$SAMPLE_NAME.gc_bias_metrics.pdf S=$BAM/$SAMPLE_NAME.summary_metrics.tsv > $BAM/$SAMPLE_NAME.gc_metrics.log 2>&1\n"
echo "Running picard markdup"
PICARD_CMD="${PICARD_CMD}java -Xmx2g -Djava.io.tmpdir=$JAVATMP -d64 -server -jar $PICARD MarkDuplicates TMP_DIR=$JAVATMP MAX_RECORDS_IN_RAM=4000000 OPTICAL_DUPLICATE_PIXEL_DISTANCE=10000 I=$BAM/$SAMPLE_NAME.bam O=/dev/null M=$BAM/$SAMPLE_NAME.picard_markdup.tsv > $BAM/$SAMPLE_NAME.picard_markdup.log 2>&1\n"

echo -e $PICARD_CMD | parallel --no-notice -j 5

date
