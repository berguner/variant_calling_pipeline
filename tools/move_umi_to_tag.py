#!/usr/bin/env python3

"""
This script transfers the UMI tags in the read names to BAM tags.

Coder: Bekir Erguner, 2018
"""

from argparse import ArgumentParser
import sys


parser = ArgumentParser(description='Move UMI tag in the read name to BAM tags, reads SAM from STDIN and writes to STDOUT')
parser.add_argument('--barcode-tag', dest="bt", default="RX:Z:", help="UMI barcode tag",
                    required=False, type=str)
parser.add_argument('--quality-tag', dest="qt", default="QX:Z:", help="UMI quality tag",
                    required=False, type=str)
args = parser.parse_args()

line = sys.stdin.readline().rstrip("\n")
while line:
    if line[0] == '#':
        sys.stdout.write(line + "\n")
        line = sys.stdin.readline().rstrip("\n")
        continue

    line_arr = line.split("\t")
    readname = line_arr[0]
    new_readname = ''
    bt_str = ''
    qt_str = ''
    if readname.find(args.bt) > -1 and readname.find(args.qt) > -1:
        new_readname = readname[:readname.find(args.bt)]
        bt_str = readname[readname.find(args.bt):readname.find(args.qt)]
        qt_str = readname[readname.find(args.qt):]
        line_arr.append(bt_str)
        line_arr.append(qt_str)
        line_arr[0] = new_readname
        sys.stdout.write("\t".join(line_arr) + "\n")
    else:
        sys.stdout.write(line + "\n")

    line = sys.stdin.readline().rstrip("\n")
