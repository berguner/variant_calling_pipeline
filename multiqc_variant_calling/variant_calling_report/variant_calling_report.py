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
