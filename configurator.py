#!/usr/bin/env python3

import argparse
import yaml, json
import csv
import os
import sys
from collections.abc import Mapping
from string import Template

def update(d, u):
    """
    Recursively updates the entries in a given dictionary
    :param d: The dictionary to be updated
    :param u: The values which will be added to the input dictionary
    :return: Updated dictionary
    """
    for k, v in u.items():
        if isinstance(v, Mapping):
            d[k] = update(d.get(k, {}), v)
        else:
            d[k] = v
    return d


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="JSON config generator for the Variant Calling Pipeline")
    parser.add_argument('--project-config', '-c',
                        dest='project_config',
                        help='Project specific .yaml config file',
                        type=str)
    parser.add_argument('--pipeline-config', '-p',
                        dest='pipeline_config',
                        help='Pipeline specific .yaml config file',
                        type=str)
    args = parser.parse_args()

    # Read config files
    # default will contain pipeline configurations
    default = {}
    # specific will contain the project configurations
    specific = {}
    with open(args.pipeline_config, 'r') as stream:
        try:
            default = yaml.safe_load(stream)
        except yaml.YAMLError as exception:
            sys.stderr.write(str(exception))
    with open(args.project_config, 'r') as stream:
        try:
            specific = yaml.safe_load(stream)
        except yaml.YAMLError as exception:
            sys.stderr.write(str(exception))
    update(default, specific)
    # Rename the updated dictionary
    config = default

    # Create the output folders
    project_path = config['project_path']
    output_path = os.path.join(project_path, config['genome'])
    json_path = os.path.join(project_path, 'config_files')
    if not os.path.exists(project_path):
        os.mkdir(project_path)
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    if not os.path.exists(json_path):
        os.mkdir(json_path)

    # Parse necessary configurations
    project_genome = config['genome']
    inputs_dict = {'variant_calling.project_name': config['project_name'],
                   'variant_calling.project_path': config['project_path'],
                   'variant_calling.genome': project_genome
                   }

    # Collect direct configurations for WDL tasks
    for key in config['bypass_parser']:
        inputs_dict[key] = config['bypass_parser'][key]

    # Parse sample specific configurations
    sas_file = config['sample_annotation']
    sas_dict = {}
    germline_samples = []
    tumor_samples = []
    with open(sas_file, 'r') as sas:
        reader = csv.DictReader(sas, dialect='excel')
        for row in reader:
            if 'sample_name' in row:
                if row['sample_name'] in sas_dict:
                    sas_dict[row['sample_name']].append(row)
                else:
                    sas_dict[row['sample_name']] = [row]
                    if row['sample_type'] == 'germline':
                        germline_samples.append(row['sample_name'])
                    if row['sample_type'] == 'tumor':
                        tumor_samples.append(row['sample_name'])

    # WDL workflow will iterate over these samples
    inputs_dict['variant_calling.sample_list'] = list(sas_dict.keys())
    if len(germline_samples) > 0:
        inputs_dict['variant_calling.germline_samples'] = germline_samples
    if len(tumor_samples) > 0:
        inputs_dict['variant_calling.tumor_samples'] = tumor_samples

    # Dump the inputs.json file
    project_json = os.path.join(json_path, '{}.inputs.json'.format(config['project_name']))
    with open(project_json, 'w') as output:
        json.dump(inputs_dict, output, indent=2)

    # Read sample details
    sample_dicts = []
    for sample in sas_dict:
        sample_dict = {'sample_name': sample,
                       'target_intervals': sas_dict[sample][0]['target_intervals'],
                       'sample_type': sas_dict[sample][0]['sample_type'],
                       'library': sas_dict[sample][0]['library'],
                       'raw_bams': ''}

        row_list = sas_dict[sample]
        number_of_rows = len(row_list)
        bam_sources = []
        raw_size_mb = 0
        # Collect the list of raw data input files for each sample
        for i in range(number_of_rows):
            if 'data_source' in row_list[i] and row_list[i]['data_source'] != '':
                source_template = config['data_sources'][row_list[i]['data_source']]
                source = source_template.format(**row_list[i])
                if os.path.exists(source):
                    bam_sources.append(source)
                    if os.path.exists(source):
                        source_stats = os.stat(source)
                        raw_size_mb += int(source_stats.st_size / (1024 * 1024))
                else:
                    print('WARNING: Could not locate {} for sample {}'.format(source))
        # Skip the sample in case the input files weren't located
        if len(bam_sources) == 0:
            print('WARNING: Could not locate any raw data files for sample {}, skipping.'.format(sample))
        else:
            print('Sample {} has {} data sources, total size is {}MB'.format(sample, len(bam_sources), raw_size_mb))
            sample_dict['raw_bams'] = ' '.join(bam_sources)
            sample_dict['raw_size_mb'] = raw_size_mb
            sample_tsv = os.path.join(json_path, '{}.tsv'.format(sample))
            with open(sample_tsv, 'w') as output:
                for key in sample_dict:
                    output.write('{}\t{}\n'.format(key, sample_dict[key]))
