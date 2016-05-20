#!/usr/bin/env python2.7

import os
import sys
import HTSeq
import random
import itertools


def usage():
    print 'reduceFastq.py <1.fq> <2.fq> <prob>'
    sys.exit()


def main():
    if len(sys.argv) != 4:
        usage()

    fq1 = HTSeq.FastqReader(sys.argv[1])
    fq2 = HTSeq.FastqReader(sys.argv[2])

    fq1o = open('output.1.fq', 'w')
    fq2o = open('output.2.fq', 'w')

    count = 0
    written = 0
    for r1, r2 in itertools.izip(fq1, fq2):
        count += 1
        if random.random() <= float(sys.argv[3]):
            written += 1
            r1.write_to_fastq_file(fq1o)
            r2.write_to_fastq_file(fq2o)

    print '{} Reads total; Wrote {}'.format(count, written)


if __name__ == '__main__':
    main()
