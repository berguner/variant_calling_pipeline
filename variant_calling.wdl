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

        # runtime parameters
        Int cpus = 8
        Int memory = 16000
        String partition = "mediumq"
        String time = "2-00:00:00"
        String rt_additional_parameters = ""
    }

    # bam_dir would be /project_path/genome/bam
    String bam_dir = "~{output_dir}/bam"
    String raw_bams = sample.raw_bams
    # Read group ID and various info which is needed for GATK tools, BWA adds it to the aligned bam
    String RG = "@RG\\tID:~{sample.sample_name}\\tSM:~{sample.sample_name}\\tPL:illumina\\tLB:~{sample.library}"

    command <<<
        [ ! -d "~{output_dir}" ] && mkdir -p ~{output_dir};
        [ ! -d "~{bam_dir}" ] && mkdir -p ~{bam_dir};

        # Run alignment if the output_bam file does not exist
        if [ ! -f "~{bam_dir}/~{sample.sample_name}.bam" ]; then
            # Set this for enabling summation of return codes from the piped commands
            set -o pipefail

            for i in ~{raw_bams}; do samtools fastq $i 2>> "~{bam_dir}/~{sample.sample_name}.samtools.log" ; done | \
                bwa mem -t ~{cpus} -R "~{RG}" -p ~{ref_fasta} - 2> "~{bam_dir}/~{sample.sample_name}.bwa.log" | \
                samtools fixmate -O SAM - - 2>> "~{bam_dir}/~{sample.sample_name}.samtools.log" | \
                samtools sort -m 1024m -@ ~{cpus / 2} -o "~{bam_dir}/~{sample.sample_name}.aligned.bam" - 2>> "~{bam_dir}/~{sample.sample_name}.samtools.log"

            # Continue if nothing had failed during the alignment step
            if [ $? -eq "0" ]; then
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

                # Remove temporary aligned.bam file
                if [ $? -eq "0" ]; then
                    rm "~{bam_dir}/~{sample.sample_name}.aligned.bam";
                    mv "~{bam_dir}/~{sample.sample_name}.bai" "~{bam_dir}/~{sample.sample_name}.bam.bai";
                fi
            fi
        fi
    >>>

    runtime {
        rt_cpus: cpus
        rt_mem: memory
        rt_queue: partition
        rt_time: time
        rt_additional_parameters: rt_additional_parameters
    }

    output {
        File output_bam = "~{bam_dir}/~{sample.sample_name}.bam"
        File output_bai = "~{bam_dir}/~{sample.sample_name}.bam.bai"
        Sample processed_sample = sample
    }
}

task collect_wes_metrics {
    input{
        String output_dir
        Sample sample
        String ref_fasta
        String? gatk_jar

        # runtime parameters
        Int cpus = 2
        Int memory = 8000
        String partition = "mediumq"
        String time = "2-00:00:00"
        String rt_additional_parameters = ""
    }

    # output_dir should be /project_path/genom/bam/
    File bam = "~{output_dir}/~{sample.sample_name}.bam"
    File bai = "~{output_dir}/~{sample.sample_name}.bam.bai"
    String sample_name = "~{sample.sample_name}"
    String target_intervals = "~{sample.target_intervals}"

    command {
        # Run each step unless the output files are there
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_jar}
        if [ ! -f "~{output_dir}/~{sample_name}.insert_size_metrics.tsv" ]; then
            gatk --java-options "-Xmx4g" CollectInsertSizeMetrics \
                -R ~{ref_fasta} \
                -I ~{bam} \
                -H ~{output_dir}/~{sample_name}.insert_size_metrics.pdf \
                -O ~{output_dir}/~{sample_name}.insert_size_metrics.tsv \
                > ~{output_dir}/~{sample_name}.insert_size_metrics.log 2>&1;
        fi

        if [ ! -f "~{output_dir}/~{sample_name}.alignment_summary_metrics.tsv" ]; then
            gatk --java-options "-Xmx4g" CollectAlignmentSummaryMetrics \
                -R ~{ref_fasta} \
                -I ~{bam} \
                -O ~{output_dir}/~{sample_name}.alignment_summary_metrics.tsv \
                > ~{output_dir}/~{sample_name}.alignment_summary_metrics.log 2>&1;
        fi

        if [ ! -f "~{output_dir}/~{sample_name}.HS_metrics.tsv" ]; then
            gatk --java-options "-Xmx4g" CollectHsMetrics \
                -R ~{ref_fasta} \
                -I ~{bam} \
                --BAIT_INTERVALS ~{target_intervals} \
                --TARGET_INTERVALS ~{target_intervals} \
                -O ~{output_dir}/~{sample_name}.HS_metrics.tsv \
                > ~{output_dir}/~{sample_name}.HS_metrics.log 2>&1;
        fi

    }
    runtime {
        rt_cpus: cpus
        rt_mem: memory
        rt_queue: partition
        rt_time: time
        rt_additional_parameters: rt_additional_parameters
    }
    output{
        File insert_size_metrics = "~{output_dir}/~{sample_name}.insert_size_metrics.tsv"
        File alignment_summary_metrics = "~{output_dir}/~{sample_name}.alignment_summary_metrics.tsv"
        File HS_metrics = "~{output_dir}/~{sample_name}.HS_metrics.tsv"
        File duplicate_metrics = "~{output_dir}/~{sample_name}.duplicate_metrics.tsv"
    }
}

task generate_sample_gvcf {
    input {
        String output_dir
        Sample sample
        String ref_fasta
        String ref_dict
        String calling_intervals
        Int? interval_padding
        String? gatk_jar

        # runtime parameters
        Int cpus = 8
        Int memory = 16000
        String partition = "mediumq"
        String time = "2-00:00:00"
        String rt_additional_parameters = ""
    }

    String parent_dir = "~{output_dir}/gvcf"
    String gvcf_dir = "~{parent_dir}/~{sample.sample_name}"
    File bam = "~{output_dir}/bam/~{sample.sample_name}.bam"
    File bai = "~{output_dir}/bam/~{sample.sample_name}.bam.bai"

    command <<<
        [ ! -d "~{parent_dir}" ] && mkdir -p ~{parent_dir};
        [ ! -d "~{gvcf_dir}" ] && mkdir -p ~{gvcf_dir};

        if [ ! -f "~{output_dir}/~{sample.sample_name}.g.vcf.gz" ]; then
            gatk --java-options "-Djava.io.tmpdir=~{gvcf_dir} -Xmx1g" \
                SplitIntervals -R ~{ref_fasta} \
                -L ~{calling_intervals} \
                --scatter-count ~{cpus} \
                -O ~{gvcf_dir} \
                > ~{gvcf_dir}/SplitIntervals.log 2>&1;

            for interval_file in ~{gvcf_dir}/*.interval_list; do
                BN=$(basename $interval_file);
                interval_name=${BN/.interval_list/};
                echo -e "gatk --java-options \"-Xmx3g -Xms128m -Djava.io.tmpdir:~{gvcf_dir} -d64\" \
                    HaplotypeCaller \
                    -L ${interval_file} \
                    ~{"--interval-padding " + interval_padding } \
                    -O ~{gvcf_dir}/${interval_name}.g.vcf.gz \
                    -R ~{ref_fasta} \
                    -ERC GVCF \
                    -I ~{bam} > ~{gvcf_dir}/${interval_name}.g.vcf.gz.log 2>&1";
            done | parallel --no-notice -j ~{cpus} ::: > ~{gvcf_dir}/parallel.log 2>&1;

            GVCF_LIST=""
            for interval_file in ~{gvcf_dir}/*.interval_list; do
                BN=$(basename $interval_file);
                interval_name=${BN/.interval_list/};
                GVCF_LIST="${GVCF_LIST} -I ~{gvcf_dir}/${interval_name}.g.vcf.gz";
            done

            gatk --java-options "-Djava.io.tmpdir=~{gvcf_dir} -Xmx8g" \
                MergeVcfs \
                -D ~{ref_dict} \
                ${GVCF_LIST} \
                -O ~{parent_dir}/~{sample.sample_name}.g.vcf.gz \
                > ~{parent_dir}/~{sample.sample_name}.MergeVcfs.log 2>&1;

            if [ $? -eq "0" ]; then
                rm ~{gvcf_dir}/*.g.vcf.gz ~{gvcf_dir}/*.interval_list;
            fi
        fi
    >>>

    runtime {
        rt_cpus: cpus
        rt_mem: memory
        rt_queue: partition
        rt_time: time
        rt_additional_parameters: rt_additional_parameters
    }

    output {
        File output_gvcf = "~{parent_dir}/~{sample.sample_name}.g.vcf.gz"
        File output_gvcf_tbi = "~{parent_dir}/~{sample.sample_name}.g.vcf.gz.tbi"
        Sample gvcf_sample = sample
    }
}
workflow variant_calling {
    input {
        String project_name
        String project_path
        String variant_calling_intervals
        Int? interval_padding
        String ref_fasta
        String ref_fai
        String genome
        String ref_dict
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
    }

    scatter(sample in bwa_align_ubam.processed_sample) {
        call collect_wes_metrics {
            input:
                sample = sample,
                ref_fasta = ref_fasta,
                output_dir = bam_dir,
                gatk_jar = gatk_jar
        }
    }

    scatter(sample in bwa_align_ubam.processed_sample) {
        if(sample.sample_type == "germline") {
            call generate_sample_gvcf {
                input:
                    sample = sample,
                    ref_fasta = ref_fasta,
                    ref_dict = ref_dict,
                    calling_intervals = variant_calling_intervals,
                    interval_padding = interval_padding,
                    output_dir = output_dir,
                    gatk_jar = gatk_jar
            }
        }
    }

}