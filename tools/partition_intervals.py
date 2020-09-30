#!/usr/bin/env python

import sys
import os
from argparse import ArgumentParser


parser = ArgumentParser(description='Partition an interval_list file into smaller intervals')
parser.add_argument('-i', '--interval-file', dest="i", help="Input interval_list file including all of the intervals",
                    required=True, type=str)
parser.add_argument('-o', '--output-folder', dest="o", help="Output folder path to write interval_list files",
                    required=True, type=str)
parser.add_argument('-p', '--prefix', dest="p", help="Prefix for generated interval_list files", default="",
                    required=False, type=str)
parser.add_argument('-s', '--size', dest="s", help="Maximum target size of the generated interval_lists",
                    required=True, type=int)
args = parser.parse_args()


intervalFile = open(args.i, 'r')
outFolder = args.o
partitionSize = args.s

if not os.path.exists(outFolder):
    os.mkdir(outFolder)

header = ''
line = intervalFile.readline()
while line[0] == '@':
    header += line
    line = intervalFile.readline()

prefix = args.p
suffix = '.interval_list'
counter = 1
currentSize = 0
currentIntervals = []
while line:
    currentIntervals.append(line)
    currentSize += (int(line.split('\t')[2]) - int(line.split('\t')[1]))
    if currentSize >= partitionSize and len(currentIntervals) > 1:
        tmpInterval = currentIntervals.pop(len(currentIntervals)-1)
        outName = prefix + str(counter) + suffix
        outPath = os.path.join(outFolder, outName)
        outFile = open(outPath, 'w')
        outFile.write(header)
        print('Size of the interval {} is {}'.format(counter, currentSize))
        for interval in currentIntervals:
            outFile.write(interval)
        outFile.close()
        counter += 1
        currentIntervals = [tmpInterval]
        currentSize = (int(line.split('\t')[2]) - int(line.split('\t')[1]))
    line = intervalFile.readline()

if currentSize >= 1:
    outName = prefix + str(counter) + suffix
    outPath = os.path.join(outFolder, outName)
    outFile = open(outPath, 'w')
    outFile.write(header)
    print('Size of the interval {} is {}'.format(counter, currentSize))
    for interval in currentIntervals:
        outFile.write(interval)
    outFile.close()

