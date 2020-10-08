#!/usr/bin/env python
"""
MultiQC plugin for reporting results of pipelines used at BSF
"""

from setuptools import setup, find_packages

version = '0.1'

setup(
    name = 'variant_calling_report',
    version = 0.1,
    author = 'Bekir Erguener',
    author_email = 'berguener@cemm.at',
    description = "MultiQC plugin for reporting results of the variant calling pipeline used at BSF",
    long_description = __doc__,
    keywords = 'bioinformatics',
    url = 'https://github.com/berguner/multiqc_variant_calling',
    download_url = 'https://github.com/berguner/multiqc_variant_calling/releases',
    license = 'MIT',
    packages = find_packages(),
    include_package_data = True,
    install_requires =[
        'multiqc',
        'click'
    ],
    entry_points = {
        'multiqc.modules.v1': [
            'variant_calling = variant_calling_report.modules.variant_calling:MultiqcModule'
        ],
        'multiqc.cli_options.v1': [
            'disable_variant_calling_report = variant_calling_report.cli:disable_variant_calling_report'
        ],
        'multiqc.hooks.v1': [
            'execution_start = variant_calling_report.variant_calling_report:variant_calling_report_execution_start'
        ]
    },
    classifiers = [
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Environment :: Web Environment',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX',
        'Operating System :: Unix',
        'Programming Language :: Python',
        'Programming Language :: JavaScript',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Visualization',
    ]
)
