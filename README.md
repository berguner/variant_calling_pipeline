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
- Install the required software and/or pull the docker image:
`docker pull berguner/bsf_variant_calling:0.3`
- Install the custom MultiQC plugin:
`python3 /path/to/multiqc_variant_calling/setup.py install`
- Collect the necessary reference files listed in the pipeline configuration file `variant_calling_config_hg38.yaml` and update the paths in this file.
Most of the files are from [GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle).

## 3. Using the Pipeline

### 3.1. Creating a New Project
For every project you need to prepare a sample annotation sheet (.csv) and project configuration file (.yaml).

**Project Configuration File:**

- The project configuration file should include the project specific parameters such as the output path.
  You can simply copy the `test/BSA_0000_variant_calling.yaml` file and update the parameters accordingly.
- You can override any configuration defined in the pipeline configuration file by setting the same parameter in the project configuration file.
- `data_sources` is an important configuration which should be included in either the pipeline configuration file.
  Under this parameter, there should be a template string that describes the full path for each `data source`.
  This will be used to locate the raw data files of the samples. 

**Sample Annotation Sheet:**

- Sample annotation sheet should include all the attributes of a sample to be able to configure the analysis tasks for that sample.
  You can simply copy the `test/BSA_0000_test_variant_calling_sample_annotations.csv` file and update the fields accordingly.
- If there are multiple data files (i.e. from multiple lanes) for a sample, just add multiple lines for that sample changing only the data source field(s). 

### 3.2 Running the Pipeline
- We recommend using Cromwell as the workflow execution engine.
  You can define and use different backend configurations for various compute environments.
  The default backend configuration used on the CeMM HPC cluster is `slurm_singularity.conf`.
- First run the configurator scrtipt to prepare the project folder and inputs json file.
```
python3 /path/to/variant_calling_pipeline/configurator.py \
  -p /path/to/variant_calling_pipeline/variant_calling_config_hg38.yaml \
  -c /path/to/my_project_config.yaml
```
- Then change into the project output folder and start/submit cromwell.
```
cd /path/to/my_project_folder
sbatch -J my_project \
  --partition longq --qos longq \
  --time 20-00:00:00 \
  --mem 32000 --cpus-per-task 2 \
  --wrap "java -Xmx24g -Dconfig.file=/path/to/variant_calling_pipeline/backends/slurm_singularity.conf -Xmx24g -jar /path/to/cromwell-59.jar run /path/to/variant_calling_pipeline/variant_calling.wdl  --inputs config_files/my_project.inputs.json"
```

### 3.3 Preparing the MultiQC Report
- Make sure that the `variant-calling-report`package is installed by running `pip list`.
- Once the pipeline has finished running, you can run MultiQC to generate the project report.
```
multiqc -f -c /path/to/my_project_config.yaml /path/to/my_project_folder
```
- You should see the info message `Running variant_calling_report MultiQC Plugin` if everything was configured correctly.