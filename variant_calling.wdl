version 1.0

struct Sample {
    String sample_name
    String raw_bams
    String library
    String sample_type
    String target_intervals
}


task bwa_align_ubam {
    input {
        String output_dir
        String ref_fasta
        String? gatk_jar
        Sample sample

        Int cpus = 8
        Int memory = 16000
        String partition = "mediumq"
        String time = "2-00:00:00"
    }

    String bam_dir = "~{output_dir}/bam"
    String raw_bams = sample.raw_bams
    String RG = "@RG\\tID:~{sample.sample_name}\\tSM:~{sample.sample_name}\\tPL:illumina\\tLB:~{sample.library}"

    command <<<
        if [ ! -f "~{bam_dir}/~{sample.sample_name}.bam" ]; then

            [ ! -d "~{output_dir}" ] && mkdir -p ~{output_dir};
            [ ! -d "~{bam_dir}" ] && mkdir -p ~{bam_dir};

            for i in ~{raw_bams}; do samtools fastq $i 2>> "~{bam_dir}/~{sample.sample_name}.samtools.log" ; done | \
                bwa mem -t ~{cpus} -R "~{RG}" -p ~{ref_fasta} - 2> "~{bam_dir}/~{sample.sample_name}.bwa.log" | \
                samtools fixmate -O SAM - - 2>> "~{bam_dir}/~{sample.sample_name}.samtools.log" | \
                samtools sort -m 512m -@ ~{cpus / 2} -o "~{bam_dir}/~{sample.sample_name}.aligned.bam" - 2>> "~{bam_dir}/~{sample.sample_name}.samtools.log"

            set -e
            export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_jar}

            gatk --java-options "-Djava.io.tmpdir=~{output_dir} -Xmx8g" MarkDuplicates \
                -I "~{bam_dir}/~{sample.sample_name}.aligned.bam" \
                -M "~{bam_dir}/~{sample.sample_name}.duplicate_metrics.tsv" \
                -O "~{bam_dir}/~{sample.sample_name}.bam" \
                --OPTICAL_DUPLICATE_PIXEL_DISTANCE 10000 \
                --ASSUME_SORTED true \
                --CREATE_INDEX true \
                > "~{bam_dir}/~{sample.sample_name}.markduplicates.log" 2>&1;

            if [ $? -eq "0" ]; then
                rm "~{bam_dir}/~{sample.sample_name}.aligned.bam";
                mv "~{bam_dir}/~{sample.sample_name}.bai" "~{bam_dir}/~{sample.sample_name}.bam.bai";
            fi

        fi
    >>>

    runtime {
        rt_cpus: cpus
        rt_mem: memory
        rt_queue: partition
        rt_time: time
    }

    output {
        File output_bam = "~{bam_dir}/~{sample.sample_name}.bam"
        File output_bai = "~{bam_dir}/~{sample.sample_name}.bam.bai"
    }
}

task collect_wes_metrics {
    input{
        String output_dir
        Sample sample
        String ref_fasta
        File bam
        File bai
        String? gatk_jar

        Int cpus = 2
        Int memory = 8000
        String partition = "mediumq"
        String time = "2-00:00:00"
    }

    String sample_name = "~{sample.sample_name}"
    String target_intervals = "~{sample.target_intervals}"

    command{
        set -e

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_jar}

        gatk --java-options "-Xmx4g" CollectInsertSizeMetrics \
            -R ~{ref_fasta} \
            -I ~{bam} \
            -H ~{output_dir}/~{sample_name}.insert_size_metrics.pdf \
            -O ~{output_dir}/~{sample_name}.insert_size_metrics.tsv \
            > ~{output_dir}/~{sample_name}.insert_size_metrics.log 2>&1;

        gatk --java-options "-Xmx4g" CollectAlignmentSummaryMetrics \
            -R ~{ref_fasta} \
            -I ~{bam} \
            -O ~{output_dir}/~{sample_name}.alignment_summary_metrics.tsv \
            > ~{output_dir}/~{sample_name}.alignment_summary_metrics.log 2>&1;

        gatk --java-options "-Xmx4g" CollectHsMetrics \
            -R ~{ref_fasta} \
            -I ~{bam} \
            --BAIT_INTERVALS ~{target_intervals} \
            --TARGET_INTERVALS ~{target_intervals} \
            -O ~{output_dir}/~{sample_name}.HS_metrics.tsv \
            > ~{output_dir}/~{sample_name}.HS_metrics.log 2>&1;

    }
    runtime {
        rt_cpus: cpus
        rt_mem: memory
        rt_queue: partition
        rt_time: time
    }
    output{
        File insert_size_metrics = "~{output_dir}/~{sample_name}.insert_size_metrics.tsv"
        File alignment_summary_metrics = "~{output_dir}/~{sample_name}.alignment_summary_metrics.tsv"
        File HS_metrics = "~{output_dir}/~{sample_name}.HS_metrics.tsv"
        File duplicate_metrics = "~{output_dir}/~{sample_name}.duplicate_metrics.tsv"
    }
}

workflow variant_calling {
    input {
        String project_name
        String project_path
        String variant_calling_intervals_folder
        String ref_fasta
        String ref_fai
        String genome
        String? ref_dict
        String dbsnp_vcf
        String dbsnp_idx
        String? gatk_jar
        Array[String] sample_list
    }

    String output_dir = "~{project_path}/~{genome}"
    String config_dir = "~{project_path}/config_files"
    String bam_dir = "~{project_path}/~{genome}/bam"
    scatter(sample_name in sample_list) {
        File sample_tsv  = "~{config_dir}/~{sample_name}.tsv"
        Map[String, String] sample_map = read_map(sample_tsv)
        Sample sample = { "sample_name": sample_name,
                            "raw_bams": sample_map["raw_bams"],
                            "library": sample_map["library"],
                            "sample_type": sample_map["sample_type"],
                            "target_intervals": sample_map["target_intervals"]
                        }

        call bwa_align_ubam {
            input:
                sample = sample,
                ref_fasta = ref_fasta,
                output_dir = output_dir,
                gatk_jar = gatk_jar
        }
        call collect_wes_metrics {
            input:
                sample = sample,
                ref_fasta = ref_fasta,
                output_dir = bam_dir,
                bam = bwa_align_ubam.output_bam,
                bai = bwa_align_ubam.output_bai,
                gatk_jar = gatk_jar
        }
    }

}