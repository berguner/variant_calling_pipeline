#!/usr/bin/env python
""" MultiQC BSF plugin functions
We can add any custom Python functions here and call them
using the setuptools plugin hooks.
"""


from __future__ import print_function
from pkg_resources import get_distribution
import logging
import os, csv

from multiqc.utils import report, util_functions, config

# Initialise the main MultiQC logger
log = logging.getLogger('multiqc')

# Save this plugin's version number (defined in setup.py) to the MultiQC config
config.variant_calling_report_version = get_distribution("variant_calling_report").version


# Add default config options for the things that are used in bsf_reports
def variant_calling_report_execution_start():
    """
    Code to execute after the config files and
    command line flags have been parsed self.
    this setuptools hook is the earliest that will be able
    to use custom command line flags.
    """
    # Halt execution if we've disabled the plugin
    if config.kwargs.get('disable_variant_calling_report', True):
        return None

    log.info("Running variant_calling_report MultiQC Plugin v{}".format(config.variant_calling_report_version))

    # Create symlink for the web server
    if hasattr(config, 'base_url') and hasattr(config, 'project_uuid') and hasattr(config, 'public_html_folder'):
        project_url = os.path.join(config.base_url, config.project_uuid)
        os.chdir(config.public_html_folder)
        if not os.path.exists(os.path.join(config.public_html_folder, config.project_uuid)):
            # The symlink has to be relative so that the web server can locate the project folder
            relative_path = os.path.relpath(config.project_path)
            os.symlink(relative_path, config.project_uuid)
        log.info('## You can access the project report from: ##\n{}\n'.format(os.path.join(project_url,
                                                                                           config.genome_version,
                                                                                           'report',
                                                                                           'multiqc_report.html')))
    else:
        log.error('Please provide base_url, project_uuid and public_html_folder in the configuration file')
        exit(1)

    report_dir = os.path.join(config.project_path, config.genome_version, 'report')
    if not os.path.exists(report_dir):
        os.mkdir(report_dir)
    config.output_dir = report_dir
    config.analysis_dir = [report_dir]
    os.chdir(report_dir)
    for folder in ['bam', 'vcf', 'relatedness', 'mutect2']:
        if os.path.exists(os.path.join(config.project_path, config.genome_version, folder)) and not os.path.islink(os.path.join(report_dir, folder)):
            relative_path = os.path.relpath(os.path.join(config.project_path, config.genome_version, folder), os.curdir)
            os.symlink(relative_path, os.path.join(report_dir, folder))
