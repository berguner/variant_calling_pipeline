version 1.0

import "mutect2_wdl/mutect2.wdl" as m2

struct Sample {
    String sample_name
    String raw_bams
    String library
    String sample_type
    String matched_normal
    String? target_intervals
    String? UMI
}

workflow variant_calling {
    input {
        String project_name
        String project_path
        File? variant_calling_intervals
        Int? interval_padding
        File ref_fasta
        File ref_fai
        File? ref_alt
        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_pac
        File ref_sa
        String genome_version
        File ref_dict
        File dbsnp_vcf
        File dbsnp_idx
        String? gatk_jar
        Array[String] sample_list
        Array[String]? germline_samples
        Array[String]? tumor_samples
        Int genotype_scatter_count = 10
        String? split_intervals_extra_args
        String perform_germline_variant_calling = "yes"
        String perform_somatic_variant_calling = "yes"
        String generate_germline_sample_vcfs = "yes"
        String? rt_image
        String tools_folder
        File? empty_file
    }

    String output_dir = "~{project_path}/~{genome_version}"
    String config_dir = "~{project_path}/config_files"
    String bam_dir = "~{project_path}/~{genome_version}/bam"
    scatter(sample_name in sample_list) {
        File sample_tsv  = "~{config_dir}/~{sample_name}.tsv"
        Map[String, String] sample_map = read_map(sample_tsv)
        Sample sample = { "sample_name": sample_name,
                            "raw_bams": sample_map["raw_bams"],
                            "library": sample_map["library"],
                            "sample_type": sample_map["sample_type"],
                            "target_intervals": sample_map["target_intervals"],
                            "UMI": sample_map["UMI"],
                            "matched_normal": sample_map["matched_normal"]
                        }

        call bwa_align_ubam {
            input:
                sample = sample,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_alt = ref_alt,
                ref_amb = ref_amb,
                ref_ann = ref_ann,
                ref_bwt = ref_bwt,
                ref_pac = ref_pac,
                ref_sa = ref_sa,
                output_dir = output_dir,
                rt_image = rt_image,
                tools_folder = tools_folder
        }

        call markduplicates {
            input:
                sample = bwa_align_ubam.processed_sample,
                aligned_bam = bwa_align_ubam.output_bam,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                genome_version = genome_version,
                output_dir = output_dir,
                rt_image = rt_image
        }
    }

    scatter(sample in markduplicates.processed_sample) {
        call collect_wes_metrics {
            input:
                sample = sample,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                output_dir = bam_dir,
                gatk_jar = gatk_jar,
                rt_image = rt_image
        }
    }

    if(perform_germline_variant_calling == "yes") {
        scatter(sample in markduplicates.processed_sample) {
            if(sample.sample_type == "germline") {
                call generate_sample_gvcf {
                    input:
                        sample = sample,
                        ref_fasta = ref_fasta,
                        ref_fai = ref_fai,
                        ref_dict = ref_dict,
                        calling_intervals = variant_calling_intervals,
                        interval_padding = interval_padding,
                        output_dir = output_dir,
                        gatk_jar = gatk_jar,
                        rt_image = rt_image
                }
            }
        }

        call split_intervals {
            input:
                scatter_count = genotype_scatter_count,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                intervals = variant_calling_intervals,
                gatk_jar = gatk_jar,
                split_intervals_extra_args = split_intervals_extra_args,
                rt_image = rt_image
        }

            Array[String?] my_gvcfs = generate_sample_gvcf.output_gvcf
            scatter(interval_file in split_intervals.interval_files) {
                call combine_genotype_gvcfs {
                    input:
                        interval_file = interval_file,
                        ref_fasta = ref_fasta,
                        ref_fai = ref_fai,
                        ref_dict = ref_dict,
                        interval_padding = interval_padding,
                        dbsnp_vcf = dbsnp_vcf,
                        dbsnp_idx = dbsnp_idx,
                        gatk_jar = gatk_jar,
                        output_dir = output_dir,
                        input_gvcfs =  generate_sample_gvcf.output_gvcf,
                        rt_image = rt_image
                }
            }


            call merge_combined_gvcfs {
                input:
                    output_dir = output_dir,
                    project_name = project_name,
                    ref_dict = ref_dict,
                    gatk_jar = gatk_jar,
                    input_gvcfs = combine_genotype_gvcfs.combined_gvcf,
                    input_gvcf_tbis = combine_genotype_gvcfs.combined_gvcf_tbi,
                    input_vcfs = combine_genotype_gvcfs.genotyped_vcf,
                    input_vcf_tbis = combine_genotype_gvcfs.genotyped_vcf_tbi,
                    rt_image = rt_image
            }

            String output_vcf_dir = "~{output_dir}/vcf"
            call annotate_vcf_vep {
                input:
                    output_vcf_dir = output_vcf_dir,
                    input_vcf = merge_combined_gvcfs.cohort_vcf,
                    input_vcf_tbi = merge_combined_gvcfs.cohort_vcf_tbi
            }

            if(generate_germline_sample_vcfs == "yes") {
                call generate_sample_vcfs {
                    input:
                        output_vcf_dir = output_vcf_dir,
                        input_vcf = annotate_vcf_vep.annotated_vcf,
                        input_vcf_tbi = annotate_vcf_vep.annotated_vcf_tbi,
                        dbsnp_vcf = dbsnp_vcf,
                        dbsnp_vcf_tbi = dbsnp_idx,
                        gatk_jar = gatk_jar,
                        rt_image = rt_image
                }
            }

    }
    if(perform_somatic_variant_calling == "yes") {
        scatter(sample in markduplicates.processed_sample) {
            if(sample.sample_type == "tumor") {
                File tumor_bam = "~{output_dir}/bam/~{sample.sample_name}.bam"
                File tumor_bai = "~{output_dir}/bam/~{sample.sample_name}.bam.bai"
                File? normal_bam = if(sample.matched_normal != "NULL") then "~{output_dir}/bam/~{sample.matched_normal}.bam" else empty_file
                File? normal_bai = if(sample.matched_normal != "NULL") then "~{output_dir}/bam/~{sample.matched_normal}.bam.bai" else empty_file
                call m2.Mutect2 {
                    input:
                        tumor_reads = tumor_bam,
                        tumor_reads_index = tumor_bai,
                        normal_reads = normal_bam,
                        normal_reads_index = normal_bai,
                        ref_fasta = ref_fasta,
                        ref_fai = ref_fai,
                        ref_dict = ref_dict,
                        intervals = variant_calling_intervals,
                        compress_vcfs = true
                }

                String output_vcf_dir = "~{output_dir}/mutect2"
                call annotate_vcf_vep as somatic_vcf_vep{
                    input:
                        output_vcf_dir = output_vcf_dir,
                        input_vcf = Mutect2.filtered_vcf,
                        input_vcf_tbi = Mutect2.filtered_vcf_idx
                }

                call collect_variant_calling_metrics {
                    input:
                        sample = sample,
                        output_vcf_dir = output_vcf_dir,
                        input_vcf = Mutect2.filtered_vcf,
                        input_vcf_tbi = Mutect2.filtered_vcf_idx,
                        dbsnp_vcf = dbsnp_vcf,
                        dbsnp_vcf_tbi = dbsnp_idx,
                        gatk_jar = gatk_jar,
                        rt_image = rt_image
                }
                call copy_mutect2_results {
                    input:
                        sample = sample,
                        output_vcf_dir = output_vcf_dir,
                        input_vcf = Mutect2.filtered_vcf,
                        input_vcf_tbi = Mutect2.filtered_vcf_idx,
                        filtering_stats = Mutect2.filtering_stats,
                        mutect_stats = Mutect2.mutect_stats,
                        contamination_table = Mutect2.contamination_table,
                        rt_image = rt_image
                }
            }
        }
    }
}

task bwa_align_ubam {
    input {
        String output_dir
        File ref_fasta
        File ref_fai
        File? ref_alt
        File ref_amb
        File ref_ann
        File ref_bwt
        File ref_pac
        File ref_sa
        Sample sample
        String tools_folder

        # runtime parameters
        Int cpus = 8
        Int memory = 16000
        String partition = "mediumq"
        String time = "2-00:00:00"
        String? rt_additional_parameters
        String? rt_image
    }

    # bam_dir would be /project_path/genome_version/bam
    String bam_dir = "~{output_dir}/bam"
    String raw_bams = sample.raw_bams
    # Read group ID and various info which is needed for GATK tools, BWA adds it to the aligned bam
    String RG = "@RG\\tID:~{sample.sample_name}\\tSM:~{sample.sample_name}\\tPL:illumina\\tLB:~{sample.library}"
    String samtools_fastq = if (sample.UMI == "yes") then "samtools fastq -n -T 'RX,QX' " else "samtools fastq -n "
    String move_umi = if (sample.UMI == "yes") then "python3 ~{tools_folder}/move_umi_to_tag.py |" else ""
    command <<<
        [ ! -d "~{output_dir}" ] && mkdir -p ~{output_dir};
        [ ! -d "~{bam_dir}" ] && mkdir -p ~{bam_dir};

        # Set this for enabling summation of return codes from the piped commands
        set -o pipefail
        pipe_returned=0
        if [[ ! -f "~{bam_dir}/~{sample.sample_name}.aligned.bam" && ! -f "~{bam_dir}/~{sample.sample_name}.bam.bai" ]]; then
            for i in ~{raw_bams}; do ~{samtools_fastq} $i 2>> "~{bam_dir}/~{sample.sample_name}.samtools.log" ; done | \
                awk '{ if(NR % 4 == 1) gsub("\t", "", $0); print $0}' | \
                bwa mem -t ~{cpus} -R "~{RG}" -p ~{ref_fasta} - 2> "~{bam_dir}/~{sample.sample_name}.bwa.log" | ~{move_umi} \
                samtools fixmate -O SAM - - 2>> "~{bam_dir}/~{sample.sample_name}.samtools.log" | \
                samtools sort -m 1024m -@ ~{cpus / 2} -o "~{bam_dir}/~{sample.sample_name}.aligned.bam" - 2>> "~{bam_dir}/~{sample.sample_name}.samtools.log"
            pipe_returned=$?
        fi

        if [ ! -f "~{bam_dir}/~{sample.sample_name}.aligned.bam" ]; then
            touch "~{bam_dir}/~{sample.sample_name}.aligned.bam"
        fi

        exit $pipe_returned
    >>>

    runtime {
        rt_cpus: cpus
        rt_mem: memory
        rt_queue: partition
        rt_time: time
        rt_additional_parameters: rt_additional_parameters
        rt_image: rt_image
    }

    output {
        File output_bam = "~{bam_dir}/~{sample.sample_name}.aligned.bam"
        Sample processed_sample = sample
    }
}

task markduplicates {
    input {
        String output_dir
        File aligned_bam
        File ref_fasta
        File ref_fai
        String? gatk_jar
        Sample sample
        String genome_version

        # runtime parameters
        Int cpus = 2
        Int memory = 8000
        String partition = "mediumq"
        String time = "2-00:00:00"
        String? rt_additional_parameters
        String? rt_image
    }

    # bam_dir would be /project_path/genome_version/bam
    String bam_dir = "~{output_dir}/bam"
    String mark_duplicates_method = if (sample.UMI == "yes") then "UmiAwareMarkDuplicatesWithMateCigar --UMI_METRICS ~{bam_dir}/~{sample.sample_name}.umi_metrics.tsv" else "MarkDuplicates"
    Int mark_duplicates_memory = memory - 2000

    command <<<
        [ ! -d "~{output_dir}" ] && mkdir -p ~{output_dir};
        [ ! -d "~{bam_dir}" ] && mkdir -p ~{bam_dir};

        set -e
        export GATK_LOCAL_JAR=~{default="/gatk/gatk.jar" gatk_jar}

        if [ ! -f "~{bam_dir}/~{sample.sample_name}.bam.bai" ]; then
            gatk --java-options "-Djava.io.tmpdir=~{output_dir} -Xmx~{mark_duplicates_memory}m" ~{mark_duplicates_method} \
                -I ~{aligned_bam} \
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

        if [[ -f "~{bam_dir}/~{sample.sample_name}.bam.bai" && -f "~{bam_dir}/~{sample.sample_name}.aligned.bam" ]]; then
            rm "~{bam_dir}/~{sample.sample_name}.aligned.bam";
        fi

        # Infer sex
        BAM="~{bam_dir}/~{sample.sample_name}.bam"
        REF=~{genome_version}
        SEX_FILE="~{bam_dir}/~{sample.sample_name}.sex"

        if [ $REF == "hg38" ]; then
            REG="chrY:2786847-2787698";
        elif [ $REF == "b37" ]; then
            REG="Y:2654887-2655798";
        elif [ $REF == "hg19" ]; then
            REG="chrY:2654887-2655798";
        else
            echo "Invalid genome";
        fi

        SAMERR=$((samtools view -c ${BAM} ${REG}) 2>&1 )
        SRY=`samtools view -c ${BAM} ${REG}`;

        if [[ "$SAMERR" == *"main"* ]];then
            echo "Invalid region"
        elif [ $SRY -gt 10 ]; then
            echo "~{sample.sample_name}: MALE" > $SEX_FILE;
        else
            echo "~{sample.sample_name}: FEMALE" > $SEX_FILE;
        fi

        >>>

    runtime {
        rt_cpus: cpus
        rt_mem: memory
        rt_queue: partition
        rt_time: time
        rt_additional_parameters: rt_additional_parameters
        rt_image: rt_image
    }

    output {
        File output_bam = "~{bam_dir}/~{sample.sample_name}.bam"
        File output_bai = "~{bam_dir}/~{sample.sample_name}.bam.bai"
        File duplicate_metrics = "~{bam_dir}/~{sample.sample_name}.duplicate_metrics.tsv"
        Sample processed_sample = sample
    }
}

task collect_wes_metrics {
    input{
        String output_dir
        Sample sample
        File ref_fasta
        File ref_fai
        String? gatk_jar

        # runtime parameters
        Int cpus = 4
        Int memory = 16000
        String partition = "mediumq"
        String time = "2-00:00:00"
        String? rt_additional_parameters
        String? rt_image
    }

    # output_dir should be /project_path/genom/bam/
    File bam = "~{output_dir}/~{sample.sample_name}.bam"
    File bai = "~{output_dir}/~{sample.sample_name}.bam.bai"
    String sample_name = "~{sample.sample_name}"

    command <<<
        # Run each step unless the output files are there
        set -e
        # Assume WGS if target intervals file was not provided
        export TARGET_INTERVALS=~{default="NULL" sample.target_intervals}
        export GATK_LOCAL_JAR=~{default="/gatk/gatk.jar" gatk_jar}
        COMMANDS=""
        if [ ! -f "~{output_dir}/~{sample_name}.insert_size_metrics.tsv" ]; then
            COMMANDS="${COMMANDS}gatk --java-options \"-Xmx3g\" CollectInsertSizeMetrics \
                -R ~{ref_fasta} \
                -I ~{bam} \
                -H ~{output_dir}/~{sample_name}.insert_size_metrics.pdf \
                -O ~{output_dir}/~{sample_name}.insert_size_metrics.tsv \
                > ~{output_dir}/~{sample_name}.insert_size_metrics.log 2>&1;\n"
        fi

        if [ ! -f "~{output_dir}/~{sample_name}.alignment_summary_metrics.tsv" ]; then
            COMMANDS="${COMMANDS}gatk --java-options \"-Xmx3g\" CollectAlignmentSummaryMetrics \
                -R ~{ref_fasta} \
                -I ~{bam} \
                -O ~{output_dir}/~{sample_name}.alignment_summary_metrics.tsv \
                > ~{output_dir}/~{sample_name}.alignment_summary_metrics.log 2>&1;\n"
        fi

        if [[ ! -f "~{output_dir}/~{sample_name}.HS_metrics.tsv" && "$TARGET_INTERVALS" != "NULL" ]]; then
            COMMANDS="${COMMANDS}gatk --java-options \"-Xmx3g\" CollectHsMetrics \
                -R ~{ref_fasta} \
                -I ~{bam} \
                --BAIT_INTERVALS ~{sample.target_intervals} \
                --TARGET_INTERVALS ~{sample.target_intervals} \
                -O ~{output_dir}/~{sample_name}.HS_metrics.tsv \
                > ~{output_dir}/~{sample_name}.HS_metrics.log 2>&1;\n"
        fi

        if [[ ! -f "~{output_dir}/~{sample_name}.gc_bias_metrics.tsv" && "$TARGET_INTERVALS" == "NULL" ]]; then
            COMMANDS="${COMMANDS}gatk --java-options \"-Xmx3g\" CollectGcBiasMetrics \
                -R ~{ref_fasta} \
                -I ~{bam} \
                --CHART_OUTPUT ~{output_dir}/~{sample_name}.gc_bias_metrics.pdf \
                --SUMMARY_OUTPUT ~{output_dir}/~{sample_name}.gc_summary_metrics.tsv \
                -O ~{output_dir}/~{sample_name}.gc_bias_metrics.tsv \
                > ~{output_dir}/~{sample_name}.gc_bias_metrics.log 2>&1;\n"
        fi

        if [[ ! -f "~{output_dir}/~{sample_name}.wgs_metrics.tsv" && "$TARGET_INTERVALS" == "NULL" ]]; then
            COMMANDS="${COMMANDS}gatk --java-options \"-Xmx3g\" CollectWgsMetrics \
                -R ~{ref_fasta} \
                -I ~{bam} \
                -O ~{output_dir}/~{sample_name}.wgs_metrics.tsv \
                > ~{output_dir}/~{sample_name}.wgs_metrics.log 2>&1;\n"
        fi

        echo -e ${COMMANDS} | parallel --no-notice -j ~{cpus}

    >>>
    runtime {
        rt_cpus: cpus
        rt_mem: memory
        rt_queue: partition
        rt_time: time
        rt_additional_parameters: rt_additional_parameters
        rt_image: rt_image
    }
    output{
        File insert_size_metrics = "~{output_dir}/~{sample_name}.insert_size_metrics.tsv"
        File alignment_summary_metrics = "~{output_dir}/~{sample_name}.alignment_summary_metrics.tsv"
        File? HS_metrics = "~{output_dir}/~{sample_name}.HS_metrics.tsv"
        File? GC_metrics = "~{output_dir}/~{sample_name}.gc_bias_metrics.tsv"
        File? WGS_metrics = "~{output_dir}/~{sample_name}.wgs_metrics.tsv"
    }
}

task generate_sample_gvcf {
    input {
        String output_dir
        Sample sample
        File ref_fasta
        File ref_fai
        File ref_dict
        File? calling_intervals
        Int? interval_padding
        String? gatk_jar

        # runtime parameters
        Int cpus = 8
        Int memory = 16000
        String partition = "mediumq"
        String time = "2-00:00:00"
        String? rt_additional_parameters
        String? rt_image
    }

    String parent_dir = "~{output_dir}/gvcf"
    String gvcf_dir = "~{parent_dir}/~{sample.sample_name}"
    File bam = "~{output_dir}/bam/~{sample.sample_name}.bam"
    File bai = "~{output_dir}/bam/~{sample.sample_name}.bam.bai"

    command <<<
        [ ! -d "~{parent_dir}" ] && mkdir -p ~{parent_dir};

        if [ ! -f "~{parent_dir}/~{sample.sample_name}.g.vcf.gz" ]; then
            [ ! -d "~{gvcf_dir}" ] && mkdir -p ~{gvcf_dir};

            gatk --java-options "-Djava.io.tmpdir=~{gvcf_dir} -Xmx1g" \
                SplitIntervals -R ~{ref_fasta} \
                ~{"-L " + calling_intervals} \
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
                rm -r ~{gvcf_dir};
            fi
        fi
    >>>

    runtime {
        rt_cpus: cpus
        rt_mem: memory
        rt_queue: partition
        rt_time: time
        rt_additional_parameters: rt_additional_parameters
        rt_image: rt_image
    }

    output {
        String output_gvcf = "~{parent_dir}/~{sample.sample_name}.g.vcf.gz"
        String output_gvcf_tbi = "~{parent_dir}/~{sample.sample_name}.g.vcf.gz.tbi"
    }
}

task split_intervals {
    input {
        File? intervals
        File ref_fasta
        File ref_fai
        File ref_dict
        Int scatter_count
        String? split_intervals_extra_args
        String? gatk_jar

        # runtime parameters
        Int cpus = 1
        Int memory = 4000
        String partition = "shortq"
        String time = "1:00:00"
        String? rt_additional_parameters
        String? rt_image
    }

    command {
        set -e
        export GATK_LOCAL_JAR=~{default="/gatk/gatk.jar" gatk_jar}

        mkdir interval-files
        gatk --java-options "-Xmx2g" SplitIntervals \
        -R ~{ref_fasta} \
        ~{"-L " + intervals} \
        -scatter ~{scatter_count} \
        -O interval-files
        cp interval-files/*.interval_list .
    }

    runtime {
        rt_cpus: cpus
        rt_mem: memory
        rt_queue: partition
        rt_time: time
        rt_additional_parameters: rt_additional_parameters
        rt_image: rt_image
    }

    output {
        Array[File] interval_files = glob("*.interval_list")
    }
}

task combine_genotype_gvcfs {
    input {
        String output_dir
        File interval_file
        Int? interval_padding
        File ref_fasta
        File ref_fai
        File ref_dict
        Array[String?] input_gvcfs
        File dbsnp_vcf
        File dbsnp_idx
        String? gatk_jar

        # runtime parameters
        Int cpus = 4
        Int memory = 16000
        String partition = "mediumq"
        String time = "2-00:00:00"
        String? rt_additional_parameters
        String? rt_image
    }

    # output_dir should be /project_path/genome_version/
    String combined_gvcf_dir = "~{output_dir}/gvcf/combined_gvcf"
    String BN = basename("~{interval_file}")
    String interval_name = sub(BN, "\\.interval_list$", "")

    command <<<
        [ ! -d "~{combined_gvcf_dir}" ] && mkdir -p ~{combined_gvcf_dir};

        set -e
        export GATK_LOCAL_JAR=~{default="/gatk/gatk.jar" gatk_jar}

        ### This is a temporary solution to avaoid the cromwell error:
        ### cannot interpolate Array[File?] into a command string with attribute set
        VAR_FILES=""
        for i in ~{output_dir}/gvcf/*.g.vcf.gz; do
            VAR_FILES=" -V $i $VAR_FILES"
        done

        gatk --java-options "-Xmx12g -Xms128m -Djava.io.tmpdir:~{combined_gvcf_dir} -d64" \
            CombineGVCFs \
            -L ~{interval_file} \
            ~{"--interval-padding " + interval_padding } \
            -O "~{combined_gvcf_dir}/~{interval_name}.g.vcf.gz"  \
            -R ~{ref_fasta} \
            --sequence-dictionary ~{ref_dict} \
            $VAR_FILES > "~{combined_gvcf_dir}/~{interval_name}.g.vcf.gz.log" 2>&1

        gatk --java-options "-Xmx12g -Xms128m -Djava.io.tmpdir:~{combined_gvcf_dir} -d64" \
            GenotypeGVCFs \
            --dbsnp ~{dbsnp_vcf} \
            -L ~{interval_file} \
            ~{"--interval-padding " + interval_padding } \
            -O "~{combined_gvcf_dir}/~{interval_name}.vcf.gz" \
            -R ~{ref_fasta} \
            --sequence-dictionary ~{ref_dict} \
            -V "~{combined_gvcf_dir}/~{interval_name}.g.vcf.gz" > "~{combined_gvcf_dir}/~{interval_name}.vcf.gz.log" 2>&1

    >>>

    runtime {
        rt_cpus: cpus
        rt_mem: memory
        rt_queue: partition
        rt_time: time
        rt_additional_parameters: rt_additional_parameters
        rt_image: rt_image
    }

    output {
        File genotyped_vcf = "~{combined_gvcf_dir}/~{interval_name}.vcf.gz"
        File genotyped_vcf_tbi = "~{combined_gvcf_dir}/~{interval_name}.vcf.gz.tbi"
        File combined_gvcf = "~{combined_gvcf_dir}/~{interval_name}.g.vcf.gz"
        File combined_gvcf_tbi = "~{combined_gvcf_dir}/~{interval_name}.g.vcf.gz.tbi"
    }
}

task merge_combined_gvcfs {
    # This task will merge the scattered combined (multisample) gvcf and vcf files to create cohort level gvcf and vcf files
    input {
        String output_dir
        String project_name
        File ref_dict
        Array[File]? input_gvcfs
        Array[File]? input_gvcf_tbis
        Array[File]? input_vcfs
        Array[File]? input_vcf_tbis
        String? gatk_jar

        # runtime parameters
        Int cpus = 2
        Int memory = 8000
        String partition = "mediumq"
        String time = "2-00:00:00"
        String? rt_additional_parameters
        String? rt_image
    }

    # output_dir should be /project_path/genome_version/
    String cohort_gvcf_dir = "~{output_dir}/gvcf"
    String cohort_vcf_dir = "~{output_dir}/vcf"
    String cohort_gvcf = "~{cohort_gvcf_dir}/~{project_name}_cohort.g.vcf.gz"
    String cohort_vcf = "~{cohort_vcf_dir}/~{project_name}_cohort.vcf.gz"

    command <<<
        # Creat necessary folders
        [ ! -d "~{cohort_gvcf_dir}" ] && mkdir -p ~{cohort_gvcf_dir};
        [ ! -d "~{cohort_vcf_dir}" ] && mkdir -p ~{cohort_vcf_dir};

        set -e
        export GATK_LOCAL_JAR=~{default="/gatk/gatk.jar" gatk_jar}
        # Rename the existing cohort files
        if [ -f "~{cohort_gvcf}" ]; then
            mv "~{cohort_gvcf}" "~{cohort_gvcf}.old";
            mv "~{cohort_gvcf}.tbi" "~{cohort_gvcf}.tbi.old";
        fi

        if [ -f "~{cohort_vcf}" ]; then
            mv "~{cohort_vcf}" "~{cohort_vcf}.old";
            mv "~{cohort_vcf}.tbi" "~{cohort_vcf}.tbi.old";
        fi

        # Merge the scattered gvcf files
        gatk --java-options "-Djava.io.tmpdir=~{cohort_gvcf_dir} -Xmx6g" \
                MergeVcfs \
                -D ~{ref_dict} \
                -I ~{sep=' -I ' input_gvcfs} \
                -O ~{cohort_gvcf} \
                > ~{cohort_gvcf}.MergeVcfs.log 2>&1;

        # Remove the old cohort gvcf file and scattered gvcf files
        if [ $? -eq "0" ]; then
            if [ -f "~{cohort_gvcf}.old" ]; then
                rm "~{cohort_gvcf}.old" "~{cohort_gvcf}.tbi.old";
            fi
            for i in ~{sep=' ' input_gvcfs}; do rm $(readlink -f $i); done;
            for i in ~{sep=' ' input_gvcf_tbis}; do rm $(readlink -f $i); done;
        fi

        # Merge the scattered vcf files
        gatk --java-options "-Djava.io.tmpdir=~{cohort_vcf_dir} -Xmx6g" \
                MergeVcfs \
                -D ~{ref_dict} \
                -I ~{sep=' -I ' input_vcfs} \
                -O ~{cohort_vcf} \
                > ~{cohort_vcf}.MergeVcfs.log 2>&1;

        # Remove the old cohort vcf file and scattered vcf files
        if [ $? -eq "0" ]; then
            if [ -f "~{cohort_vcf}.old" ]; then
                rm "~{cohort_vcf}.old" "~{cohort_vcf}.tbi.old";
            fi
            for i in ~{sep=' ' input_vcfs}; do rm $(readlink -f $i); done;
            for i in ~{sep=' ' input_vcf_tbis}; do rm $(readlink -f $i); done;
        fi

    >>>

    runtime {
        rt_cpus: cpus
        rt_mem: memory
        rt_queue: partition
        rt_time: time
        rt_additional_parameters: rt_additional_parameters
        rt_image: rt_image
    }

    output {
        File cohort_gvcf = "~{output_dir}/gvcf/~{project_name}_cohort.g.vcf.gz"
        File cohort_gvcf_tbi = "~{output_dir}/gvcf/~{project_name}_cohort.g.vcf.gz.tbi"
        File cohort_vcf = "~{output_dir}/vcf/~{project_name}_cohort.vcf.gz"
        File cohort_vcf_tbi = "~{output_dir}/vcf/~{project_name}_cohort.vcf.gz.tbi"
    }
}

task annotate_vcf_vep {
    input {
        String vep_executable
        String dir_cache
        String dir_plugins
        String CADD_SNVs
        String CADD_indels
        String fasta
        String assembly
        String output_vcf_dir
        String vcf_clinvar
        File input_vcf
        File input_vcf_tbi

        # runtime parameters
        Int cpus = 16
        Int memory = 32000
        String partition = "shortq"
        String time = "12:00:00"
        String? rt_additional_parameters
        String? rt_image
    }

    String input_basename = basename("~{input_vcf}")
    # remove the -filtered suffix added by Mutect2
    String pre_prefix = sub(input_basename, "-filtered", "")
    String output_prefix = sub(pre_prefix, "\\.vcf.gz$", "")

    command {
        [ ! -d "~{output_vcf_dir}" ] && mkdir -p ~{output_vcf_dir};
        ~{vep_executable} \
            --allele_number --allow_non_variant \
            --assembly ~{assembly} \
            --cache --dir_cache ~{dir_cache} --offline \
            --dir_plugins ~{dir_plugins} \
            --dont_skip --everything --failed 1 \
            --fasta ~{fasta} \
            --flag_pick_allele_gene --force_overwrite --format vcf \
            --gencode_basic --hgvsg --exclude_predicted \
            --plugin CADD,~{CADD_SNVs},~{CADD_indels} \
            --custom ~{vcf_clinvar},ClnV,vcf,exact,0,ALLELEID,ORIGIN,CLNSIG,CLNREVSTAT,CLNDN \
            --species homo_sapiens --merged \
            --tmpdir ~{output_vcf_dir} \
            --vcf --compress_output bgzip --fork ~{cpus} \
            --input_file ~{input_vcf} \
            --output_file ~{output_vcf_dir}/~{output_prefix}.vep.vcf.gz \
            --stats_file ~{output_vcf_dir}/~{output_prefix}.vep.html \
            --warning_file ~{output_vcf_dir}/~{output_prefix}.vep.warning.txt;

        tabix ~{output_vcf_dir}/~{output_prefix}.vep.vcf.gz;
    }

    runtime {
        rt_cpus: cpus
        rt_mem: memory
        rt_queue: partition
        rt_time: time
        rt_additional_parameters: rt_additional_parameters
        rt_image: rt_image
    }

    output {
        File annotated_vcf = "~{output_vcf_dir}/~{output_prefix}.vep.vcf.gz"
        File annotated_vcf_tbi = "~{output_vcf_dir}/~{output_prefix}.vep.vcf.gz.tbi"
    }
}

task generate_sample_vcfs {
    input {
        File input_vcf
        File input_vcf_tbi
        File dbsnp_vcf
        File dbsnp_vcf_tbi
        String output_vcf_dir
        String? gatk_jar

        # runtime parameters
        Int cpus = 16
        Int memory = 32000
        String partition = "mediumq"
        String time = "2-00:00:00"
        String? rt_additional_parameters
        String? rt_image
    }

    command <<<
        set -e
        export GATK_LOCAL_JAR=~{default="/gatk/gatk.jar" gatk_jar}

        SAMPLE_LIST=`bcftools query -l ~{input_vcf}`
        for sample in ${SAMPLE_LIST};
        do
            echo -e "bcftools view -s ${sample} ~{input_vcf} | bcftools view -e 'GT=\"RR\" || GT=\"mis\"' -Oz -o ~{output_vcf_dir}/${sample}.vcf.gz"
        done | parallel --no-notice -j ~{cpus} :::

        for sample in ${SAMPLE_LIST};
        do
            echo -e "bcftools index -t ~{output_vcf_dir}/${sample}.vcf.gz"
        done | parallel --no-notice -j ~{cpus} :::

        for sample in ${SAMPLE_LIST};
        do
            echo -e "gatk --java-options \"-Xmx1g -Xms128m -Djava.io.tmpdir:~{output_vcf_dir} -d64\" \
                CollectVariantCallingMetrics \
                --INPUT ~{output_vcf_dir}/${sample}.vcf.gz \
                --OUTPUT ~{output_vcf_dir}/${sample} \
                --DBSNP ~{dbsnp_vcf} > ~{output_vcf_dir}/${sample}.CollectVariantCallingMetrics.log 2>&1"
        done | parallel --no-notice -j ~{cpus} :::
    >>>

    runtime {
        rt_cpus: cpus
        rt_mem: memory
        rt_queue: partition
        rt_time: time
        rt_additional_parameters: rt_additional_parameters
        rt_image: rt_image
    }

    output {
        Array[File] all_vcf_files = glob("~{output_vcf_dir}/*.vcf.gz")
    }
}

task collect_variant_calling_metrics {
    input {
        Sample sample
        File input_vcf
        File input_vcf_tbi
        File dbsnp_vcf
        File dbsnp_vcf_tbi
        String output_vcf_dir
        String? gatk_jar

        # runtime parameters
        Int cpus = 2
        Int memory = 8000
        String partition = "shortq"
        String time = "12:00:00"
        String? rt_additional_parameters
        String? rt_image
    }

    command <<<
        [ ! -d "~{output_vcf_dir}" ] && mkdir -p ~{output_vcf_dir};
        set -e
        export GATK_LOCAL_JAR=~{default="/gatk/gatk.jar" gatk_jar}

        gatk --java-options "-Xmx4g -Xms128m -Djava.io.tmpdir:~{output_vcf_dir} -d64" \
                CollectVariantCallingMetrics \
                --INPUT ~{input_vcf} \
                --OUTPUT ~{output_vcf_dir}/~{sample.sample_name} \
                --DBSNP ~{dbsnp_vcf} > ~{output_vcf_dir}/~{sample.sample_name}.CollectVariantCallingMetrics.log 2>&1
    >>>

    runtime {
        rt_cpus: cpus
        rt_mem: memory
        rt_queue: partition
        rt_time: time
        rt_additional_parameters: rt_additional_parameters
        rt_image: rt_image
    }

    output {
        File summary_metrics_file = "~{output_vcf_dir}/~{sample.sample_name}.variant_calling_summary_metrics"
        File detail_metrics_file = "~{output_vcf_dir}/~{sample.sample_name}.variant_calling_detail_metrics"
    }
}

task copy_mutect2_results {
    input {
        Sample sample
        File input_vcf
        File input_vcf_tbi
        File filtering_stats
        File mutect_stats
        File? contamination_table
        String output_vcf_dir

        # runtime parameters
        Int cpus = 1
        Int memory = 1000
        String partition = "tinyq"
        String time = "2:00:00"
        String? rt_additional_parameters
        String? rt_image
    }

    String copy_contamination_table = if defined(contamination_table) then "cp ~{contamination_table} ~{output_vcf_dir}/~{sample.sample_name}.mutect2_contamination_table ;" else ""

    command <<<
        [ ! -d "~{output_vcf_dir}" ] && mkdir -p ~{output_vcf_dir};
        cp ~{input_vcf} ~{output_vcf_dir}/~{sample.sample_name}.vcf.gz ;
        cp ~{input_vcf_tbi} ~{output_vcf_dir}/~{sample.sample_name}.vcf.gz.tbi ;
        cp ~{filtering_stats} ~{output_vcf_dir}/~{sample.sample_name}.mutect2_filtering_stats ;
        cp ~{mutect_stats} ~{output_vcf_dir}/~{sample.sample_name}.mutect2_calling_stats ;
        ~{copy_contamination_table}
    >>>

    runtime {
        rt_cpus: cpus
        rt_mem: memory
        rt_queue: partition
        rt_time: time
        rt_additional_parameters: rt_additional_parameters
        rt_image: rt_image
    }

    output {
        File copied_vcf = "~{output_vcf_dir}/~{sample.sample_name}.vcf.gz"
        File copied_vcf_tbi = "~{output_vcf_dir}/~{sample.sample_name}.vcf.gz.tbi"
        File copied_filtering_stats = "~{output_vcf_dir}/~{sample.sample_name}.mutect2_filtering_stats"
        File copied_mutect2_stats = "~{output_vcf_dir}/~{sample.sample_name}.mutect2_calling_stats"
        File? copied_contamination_table = "~{output_vcf_dir}/~{sample.sample_name}.mutect2_contamination_table"
    }
}
