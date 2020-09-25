#!/usr/bin/env python

"""
This script transfers the UMI tags in the read names to BAM tags.

Coder: Bekir Erguner, 2018
"""

import pysam
from argparse import ArgumentParser
import os


parser = ArgumentParser(description='Add UMI to BAM')
parser.add_argument('-b', '--bam-file', dest="b", help="Input BAM file which has UMI tags in read name",
                    required=True, type=str)
parser.add_argument('-o', '--output', dest="o", help="Output BAM file with UMI tags",
                    required=True, type=str)
args = parser.parse_args()

bam_file = pysam.AlignmentFile(args.b, "rb")
out_bam = pysam.AlignmentFile(args.o, "wb", template=bam_file)
count = 0
for r in bam_file:
    count += 1
    if count % 100000 == 0:
        print("{} reads processed".format(count))
    readname = r.query_name.split(":")
    for i in range(len(readname)):
        if readname[i] == "RX" or readname[i] == "QX":
            break
        if i == len(readname) - 3:
            print("No UMI tags were detected in the read name")
            exit(1)
    outname = ":".join(readname[:i])
    r.query_name = outname.split("#")[0]
    tag1 = (readname[i],readname[i+2])
    tag2 = (readname[i+3],readname[i+5])
    tmp_tags = r.get_tags()
    tmp_tags.extend([tag1, tag2])
    r.tags = tmp_tags
    out_bam.write(r)
