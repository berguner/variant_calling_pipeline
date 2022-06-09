# BSF Variant Calling Pipeline

## 1. Introduction

This repository contains the code and configurations required to run the genomic SNV/indel detection pipeline.
The pipeline is developed by the [Biomedical Sequencing Facility](https://www.biomedical-sequencing.org/) at [CeMM](https://cemm.at/) 
in order to process whole genome, whole exome and hybrid captured panel data.
The pipeline takes unaligned BAM files as raw sequnece data and generates germline and/or somatic variant calls.
In summary, [BWA-MEM](https://github.com/lh3/bwa) is used for alignment,
[GATK4](https://gatk.broadinstitute.org/hc/en-us) HaplotypeCaller is used for germline variant calling and
[GATK4](https://gatk.broadinstitute.org/hc/en-us) Mutect2 is used for somatic variant calling.
Detected variants are annotated using [Ensembl VEP](https://www.ensembl.org/info/docs/tools/vep/index.html) and stored in VCF format.
Finally, a [MultiQC](https://multiqc.info/) report can be generated that contains a rich set of quality metrics.

## 2. Setting up the pipeline

### 2.1. Requirements

- Your favorite flavor of Linux
- python3.7 or greater
- Java 8 (JDK)
- One of the options below:
    - Docker or Singularity
    - BWA, Samtools, Bcftools, Samblaster, GATK4, Ensembl VEP
- [Cromwell-59](https://github.com/broadinstitute/cromwell/releases/download/59/cromwell-59.jar)

### 2.2. Installation
- Clone this repository
- Install the required software and/or pull the docker images:
```
docker pull berguner/bsf_variant_calling:0.3
docker pull ensemblorg/ensembl-vep:release_103.0
```
or
```
singularity pull docker://berguner/bsf_variant_calling:0.3
singularity pull docker://ensemblorg/ensembl-vep:release_103.0
```
- Install the custom MultiQC plugin locally since it isn't included in the Docker image:
```
python3 /path/to/variant_calling_pipeline/multiqc_variant_calling/setup.py install
```
- Collect the necessary reference files listed in the pipeline configuration file `variant_calling_config_hg38.yaml` and update the paths in this file.
  Most of the files are from [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle).
- Prepare the VEP cache and plugins:
```
singularity exec -e --bind /nobackup:/nobackup \
  /path/to/ensembl-vep_release_103.0.sif \
  INSTALL.pl \
  --AUTO cfp \
  --CACHEDIR /nobackup/lab_bsf/resources/VEP/ \
  --PLUGINS all \
  --PLUGINSDIR /nobackup/lab_bsf/resources/VEP/Plugins_103 \
  --SPECIES homo_sapiens \
  --ASSEMBLY GRCh38
```
- Download and use the [latest](https://github.com/broadinstitute/gatk/releases) GATK4 local jar file by changing the `variant_calling.Mutect2.gatk_override` option in the `variant_calling_config_hg38.yaml` file.
  The `bsf_variant_calling` Docker image contains a copy of GATK4, however this is a quick way to switch to a more recent version of GATK4.
  Major updates of GATK4 might require updating the Dockerfile and building the image again.

## 3. Using the Pipeline

### 3.1. Creating a New Project
For every project you need to prepare a sample annotation sheet in CSV format and a project configuration file in YAML format.

**Project Configuration File:**

- The project configuration file should include the project specific parameters such as the output path, project name etc.
  You can simply copy the `test/BSA_0000_variant_calling.yaml` file and update the parameters accordingly.
- You can override any configuration defined in the pipeline configuration file by setting the same parameter in the project configuration file.
- `data_sources` is an important configuration which should be included in either the pipeline or the project configuration file.
  Under this parameter, there should be a template string that describes the full path of each `data source`.
  This will be used to locate the raw data files of the samples.
- `exploratory_columns`: The list of fields from the sample annotation sheet that you want to see in the MultiQC report.
- `include_sex`: `yes` or `no`. The pipeline checks if there are some reads mapped to the human SRY gene to infer the sex.
  You can set this parameter to `yes` if you think that this method could infer the sex of your samples (i.e. WGS and WES) and you want to report the result in the MultiQC report. 
- `bypass_parser`: You can specify any configuration that you want to pass to the `inputs.json` file which will be used to execute the pipeline.

**Sample Annotation Sheet:**

- Sample annotation sheet should include all attributes of a sample necessary to configure the analysis tasks for that sample.
  You can simply copy the `test/BSA_0000_test_variant_calling_sample_annotations.csv` file and update the fields accordingly.
- If there are multiple data files (i.e. from multiple lanes) for a sample, just add multiple lines for that sample changing only the data source field(s).
- `sample_name`: Required.
- `data_source`: Required together with the fields to define the data source such as `flowcell`, `lane`, `BSF_name`.
- `library`: Required but only used in the BAM metadata and MultiQC report.
- `sample_type`: Required, should be either `germline` or `tumor`. Determines if the sample will go through germline or somatic variant calling.
- `matched_normal`: Sample name of the matching normal sample. It can be empty.
- `target_intervals`: Required. The path of the `.interval_list` file of targeted regions, this is used **only** for generating the hybrid selection metrics.
  Put `NULL` for whole genome samples, WGS metrics will be generated in that case.
- `sex`: The reported sex of the sample which can be included in the MultiQC report. It can be empty.
- `UMI`: `yes` or `no`. The UMI tags (`RX:Z` and `QX:Z`) in the unaligned BAM file will be transferred to the aligned BAM file using the `move_umi_to_tag.py` script.
  Samples with UMIs will be processed using `UmiAwareMarkDuplicatesWithMateCigar`. 
- You can add additional exploratory columns to populate the MultiQC report.
  Don't forget to include them in the `exploratory_columns` list in your project configuration file.

### 3.2 Running the Pipeline
- Cromwell is the workflow execution engine choise.
  You can define and use different backend configurations for various compute environments. 
- The default backend configuration used on the CeMM HPC cluster is `slurm_singularity.conf`.
  If you want to run the pipeline with one of the local backend configurations,
  adjust the `concurrent-job-limit` according to the size of the machine. 
- First, run the configurator scrtipt to prepare the project folder and `my_project.inputs.json` file.
  This script also checks if the data source files are available and their total sizes for each sample.
  Samples without any available data source will be skipped.
```
python3 /path/to/variant_calling_pipeline/configurator.py \
  -p /path/to/variant_calling_pipeline/variant_calling_config_hg38.yaml \
  -c /path/to/my_project_config.yaml
```
- Then change into the project output folder and submit the Cromwell command to Slurm:
```
cd /path/to/my_project_folder
sbatch -J my_project \
  --partition longq --qos longq \
  --time 20-00:00:00 \
  --mem 8000 --cpus-per-task 2 \
  --wrap "java -Xmx6g \
  -Dconfig.file=/path/to/variant_calling_pipeline/backends/slurm_singularity.conf \
  -Xmx24g -jar /path/to/cromwell-59.jar \
  run /path/to/variant_calling_pipeline/variant_calling.wdl \
  --inputs config_files/my_project.inputs.json"
```

- You can also run the pipeline locally:
```
java -Xmx6g \
  -Dconfig.file=/path/to/variant_calling_pipeline/backends/local_singularity.conf \
  -Xmx24g -jar /path/to/cromwell-59.jar \
  run /path/to/variant_calling_pipeline/variant_calling.wdl \
  --inputs config_files/my_project.inputs.json
```
- Cromwell `call-caching` is activated in the `slurm_singularity.conf` and `local_singularity.conf` backends.
  The files required for call caching will be stored under `/path/to/my_project_folder/cromwell-executions` folder.
  If an instance of pipeline fails to complete, the pipeline can be resubmitted and cached calls will be skipped.
  You can delete this folder after confirming that the pipeline has completed and all the files were generated.

### 3.3 Preparing the MultiQC Report
- Make sure that the `variant-calling-report`package is installed by running `pip list`.
- Once the pipeline has finished running, you can run MultiQC to generate the project report.
```
multiqc -f -c /path/to/my_project_config.yaml /path/to/my_project_folder
```
- You should see the info message `Running variant_calling_report MultiQC Plugin` if everything was configured correctly.

### 3.4 Adding New Samples to Your Project
You can simply add your new samples to the sample annotation sheet and run the pipeline.
Here are some of the diffences how various outputs are updated.

**Aligned BAM files:**
The pipeline will submit alignment jobs for all samples, but the jobs will simply exit if the sample was aligned before.

**Germline variants:**
Sample gVCF files won't be generated if they already exist.
Cohort merging and genotyping is performed from scratch for all the samples, therefore this pipeline is not ideal for more than ~1000 samples.

**Somatic variants:**
Currently, somatic variant calling is performed for all the samples.
Ideally, it should be performed only for the new samples and this remains as a future improvement.

**MultiQC report:**
Just run the `multiqc -f` command as usual. It will regenerate the report with all the samples and overwrite the old report.
