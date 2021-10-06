version 1.0
import "variant_calling.wdl" as variant_calling

workflow mutect2_multiple_tumor {
    input {
        String patient_name
        String root_dir
        File calling_intervals
        Int? interval_padding
        File ref_fasta
        File ref_fai
        File ref_dict
        String tumor_samples
        String normal_sample
        File? pon
        File? pon_idx
        File? gnomad
        File? gnomad_idx
        String? m2_extra_args
        String? getpileupsummaries_extra_args
        File? gga_vcf
        File? gga_vcf_idx
        File? variants_for_contamination
        File? variants_for_contamination_idx
        String? gatk_jar

        # runtime parameters
        Int cpus = 8
        Int memory = 32000
        String partition = "mediumq"
        String time = "2-00:00:00"
        String? rt_additional_parameters
        String? rt_image
    }

    call mutect2 {
        input:
            patient_name = patient_name,
            root_dir = root_dir,
            calling_intervals = calling_intervals,
            interval_padding = interval_padding,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            tumor_samples = tumor_samples,
            normal_sample = normal_sample,
            pon = pon,
            pon_idx = pon_idx,
            gnomad = gnomad,
            gnomad_idx = gnomad_idx,
            m2_extra_args = m2_extra_args,
            getpileupsummaries_extra_args = getpileupsummaries_extra_args,
            gga_vcf = gga_vcf,
            gga_vcf_idx = gga_vcf_idx,
            variants_for_contamination = variants_for_contamination,
            variants_for_contamination_idx = variants_for_contamination_idx,
            gatk_jar = gatk_jar,

            cpus = cpus,
            memory = memory,
            partition = partition,
            time = time,
            rt_image = rt_image,
            rt_additional_parameters = rt_additional_parameters
    }

    call variant_calling.annotate_vcf_vep {
        input:
            input_vcf = mutect2.output_vcf,
            input_vcf_tbi = mutect2.output_vcf_idx,
            output_vcf_dir = "~{root_dir}/mutect2"
    }
}

task mutect2 {
    input {
        String patient_name
        String root_dir
        File calling_intervals
        Int? interval_padding
        File ref_fasta
        File ref_fai
        File ref_dict
        String tumor_samples
        String normal_sample
        File? pon
        File? pon_idx
        File? gnomad
        File? gnomad_idx
        String? m2_extra_args
        String? getpileupsummaries_extra_args
        File? gga_vcf
        File? gga_vcf_idx
        File? variants_for_contamination
        File? variants_for_contamination_idx
        String? gatk_jar

        # runtime parameters
        Int cpus = 8
        Int memory = 32000
        String partition = "mediumq"
        String time = "2-00:00:00"
        String? rt_additional_parameters
        String? rt_image
    }

    String output_dir = "~{root_dir}/mutect2/~{patient_name}"

    String normal_input = if defined(normal_sample) then "-I ~{root_dir}/bam/~{normal_sample}.bam -normal ~{normal_sample}" else ""

    command <<<

        [ ! -d "~{root_dir}" ] && mkdir -p ~{root_dir}
        [ ! -d "~{output_dir}" ] && mkdir -p ~{output_dir}
        [ ! -d "~{root_dir}/mutect2" ] && mkdir -p "~{root_dir}/mutect2"

        export GATK_LOCAL_JAR=~{default="/gatk/gatk.jar" gatk_jar}

        gatk --java-options "-Djava.io.tmpdir=~{output_dir} -Xmx1g" \
            SplitIntervals -R ~{ref_fasta} \
            ~{"-L " + calling_intervals} \
            --scatter-count ~{cpus} \
            -O ~{output_dir} \
            > ~{output_dir}/SplitIntervals.log 2>&1;

        IFS=',' read -ra TUMOR_LIST <<< "~{tumor_samples}"
        TUMOR_INPUT=""
        for i in "${TUMOR_LIST[@]}"; do TUMOR_INPUT="${TUMOR_INPUT} -I ~{root_dir}/bam/${i}.bam" ; done

        NORMAL_INPUT=""
        if [[ ! -z "~{normal_sample}" ]]; then
            NORMAL_INPUT="-I ~{root_dir}/bam/~{normal_sample}.bam -normal ~{normal_sample}" ;
        fi

        for interval_file in ~{output_dir}/*.interval_list; do
            BN=$(basename $interval_file);
            interval_name=${BN/.interval_list/};
            echo -e "gatk --java-options \"-Xmx3g -Xms128m -Djava.io.tmpdir:~{output_dir} -d64\" Mutect2 \
                -R ~{ref_fasta} \
                $TUMOR_INPUT \
                $NORMAL_INPUT \
                ~{"--germline-resource " + gnomad} \
                ~{"-pon " + pon} \
                ~{"--alleles " + gga_vcf} \
                --f1r2-tar-gz ~{output_dir}/${interval_name}.f1r2.tar.gz \
                ~{m2_extra_args} \
                -L ${interval_file} \
                ~{"--interval-padding " + interval_padding } \
                -O ~{output_dir}/${interval_name}.vcf.gz \
                > ~{output_dir}/${interval_name}.vcf.gz.log 2>&1";
        done | parallel --no-notice -j ~{cpus} ::: > ~{output_dir}/parallel.log 2>&1;

        ### Merge VCF files
        VCF_LIST=""
        for interval_file in ~{output_dir}/*.interval_list; do
            BN=$(basename $interval_file);
            interval_name=${BN/.interval_list/};
            VCF_LIST="${VCF_LIST} -I ~{output_dir}/${interval_name}.vcf.gz";
        done

        gatk --java-options "-Djava.io.tmpdir=~{output_dir} -Xmx3g" \
            MergeVcfs \
            -D ~{ref_dict} \
            ${VCF_LIST} \
            -O ~{output_dir}/~{patient_name}.raw.vcf.gz \
            > ~{output_dir}/~{patient_name}.MergeVcfs.log 2>&1;

        ### Merge Mutect2 stats
        STATS=""
        for s in ~{output_dir}/*.stats; do
            STATS="${STATS} -stats ${s}";
        done

        gatk --java-options "-Djava.io.tmpdir=~{output_dir} -Xmx3g" \
            MergeMutectStats \
            ${STATS} \
            -O ~{output_dir}/~{patient_name}.merged.stats \
            > ~{output_dir}/~{patient_name}.MergeVcfs.log 2>&1;

        ### GetPileupSummaries
        for T in "${TUMOR_LIST[@]}"; do
            for interval_file in ~{output_dir}/*.interval_list; do
                BN=$(basename $interval_file);
                interval_name=${BN/.interval_list/};
                echo -e "gatk --java-options \"-Xmx3g -Xms128m -Djava.io.tmpdir:~{output_dir} -d64\" GetPileupSummaries \
                    -R ~{ref_fasta} \
                    -I ~{root_dir}/bam/${T}.bam \
                    --interval-set-rule INTERSECTION -L ${interval_file} \
                    -L ~{variants_for_contamination} \
                    -V ~{variants_for_contamination} \
                    -O ~{output_dir}/${T}.${interval_name}.pileups.table \
                    ~{getpileupsummaries_extra_args} \
                    > ~{output_dir}/${T}.${interval_name}.pileups.table.log 2>&1";
            done | parallel --no-notice -j ~{cpus} ::: >> ~{output_dir}/parallel.log 2>&1;

            PILEUPS=""
            for interval_file in ~{output_dir}/*.interval_list; do
                BN=$(basename $interval_file);
                interval_name=${BN/.interval_list/};
                PILEUPS="${PILEUPS} -I ~{output_dir}/${T}.${interval_name}.pileups.table"
            done

            gatk --java-options "-Xmx3g -Xms128m -Djava.io.tmpdir:~{output_dir} -d64" GatherPileupSummaries \
                --sequence-dictionary ~{ref_dict} \
                ${PILEUPS}\
                -O ~{output_dir}/${T}.pileups.table \
                > ~{output_dir}/${T}.pileups.table.log 2>&1;
        done

        if [[ ! -z "~{normal_sample}" ]]; then
            T="~{normal_sample}"
            for interval_file in ~{output_dir}/*.interval_list; do
                BN=$(basename $interval_file);
                interval_name=${BN/.interval_list/};
                echo -e "gatk --java-options \"-Xmx3g -Xms128m -Djava.io.tmpdir:~{output_dir} -d64\" GetPileupSummaries \
                    -R ~{ref_fasta} \
                    -I ~{root_dir}/bam/${T}.bam \
                    --interval-set-rule INTERSECTION -L ${interval_file} \
                    -L ~{variants_for_contamination} \
                    -V ~{variants_for_contamination} \
                    -O ~{output_dir}/${T}.${interval_name}.pileups.table \
                    ~{getpileupsummaries_extra_args} \
                    > ~{output_dir}/${T}.${interval_name}.pileups.table.log 2>&1";
            done | parallel --no-notice -j ~{cpus} ::: >> ~{output_dir}/parallel.log 2>&1;

            PILEUPS=""
            for interval_file in ~{output_dir}/*.interval_list; do
                BN=$(basename $interval_file);
                interval_name=${BN/.interval_list/};
                PILEUPS="${PILEUPS} -I ~{output_dir}/${T}.${interval_name}.pileups.table";
            done

            gatk --java-options "-Xmx3g -Xms128m -Djava.io.tmpdir:~{output_dir} -d64" GatherPileupSummaries \
                --sequence-dictionary ~{ref_dict} \
                ${PILEUPS}\
                -O ~{output_dir}/${T}.pileups.table \
                > ~{output_dir}/${T}.pileups.table.log 2>&1;
        fi

        ### CalculateContamination
        NORMAL_PILEUPS=""
        if [[ ! -z "~{normal_sample}" ]]; then
            NORMAL_PILEUPS="-matched ~{output_dir}/~{normal_sample}.pileups.table"
        fi
        CONTAMINTAION_TABLES=""
        for T in "${TUMOR_LIST[@]}"; do
            gatk --java-options "-Xmx3g -Xms128m -Djava.io.tmpdir:~{output_dir} -d64" CalculateContamination  \
                -I ~{output_dir}/${T}.pileups.table \
                ${NORMAL_PILEUPS} \
                -O ~{output_dir}/${T}.contamination.table \
                > ~{output_dir}/${T}.contamination.table.log 2>&1;
            CONTAMINTAION_TABLES="${CONTAMINTAION_TABLES} --contamination-table ~{output_dir}/${T}.contamination.table"
        done

        ### LearnReadOrientationModel
        F1R2_COUNTS=""
        for i in ~{output_dir}/*.f1r2.tar.gz; do F1R2_COUNTS="${F1R2_COUNTS} -I ${i}"; done
        gatk --java-options "-Xmx3g -Xms128m -Djava.io.tmpdir:~{output_dir} -d64" LearnReadOrientationModel \
            ${F1R2_COUNTS} \
            -O "~{output_dir}/artifact-priors.tar.gz" \
            > "~{output_dir}/artifact-priors.tar.gz.log" 2>&1;

        ### FilterMutectCalls
        gatk --java-options "-Xmx3g -Xms128m -Djava.io.tmpdir:~{output_dir} -d64" FilterMutectCalls \
            -V ~{output_dir}/~{patient_name}.raw.vcf.gz \
            -R ~{ref_fasta} \
            -O ~{output_dir}/~{patient_name}.vcf.gz \
            ${CONTAMINTAION_TABLES} \
            --ob-priors ~{output_dir}/artifact-priors.tar.gz \
            -stats ~{output_dir}/~{patient_name}.merged.stats \
            --filtering-stats ~{output_dir}/~{patient_name}.filtering.stats
            > ~{output_dir}/~{patient_name}.filtering.log 2>&1;

        ## Clean up
        rm ~{output_dir}/*-scattered*
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
        File output_vcf = "~{output_dir}/~{patient_name}.vcf.gz"
        File output_vcf_idx = "~{output_dir}/~{patient_name}.vcf.gz.tbi"
    }
}
