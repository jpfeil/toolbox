#!/usr/bin/env python2.7

import os
import sys
import HTSeq
import re
import string
import random
import itertools

def usage():
    print 'mergeAndSample.py <merge.fq> <sample.fq> <num sample>'
    sys.exit()

def main():
    if len(sys.argv) != 4:
        usage()

    merge_fqname, sample_fqname, num_sample = sys.argv[1:4]
    merge_fqname1 = merge_fqname
    _, ext = os.path.splitext(merge_fqname)
    merge_fqname2 = re.sub('1\.(?:fq|fastq)$', '2{}'.format(ext), merge_fqname)

    merge_fq1 = HTSeq.FastqReader(merge_fqname1)
    merge_fq2 = HTSeq.FastqReader(merge_fqname2)

    merge_read_names = set()

    outfile1 = open('merge_and_sample_1.fq', 'w')
    outfile2 = open('merge_and_sample_2.fq', 'w')
    for read1, read2 in itertools.izip(merge_fq1, merge_fq2):
        assert read1.name.split('/')[0] == read2.name.split('/')[0]
        if read1.name not in merge_read_names:
            read1.write_to_fastq_file(outfile1)
            read2.write_to_fastq_file(outfile2)
            merge_read_names.add(read1.name)
        else:
            print "DUPLICATE READS"

    sample_fqname1 = sample_fqname
    _, ext = os.path.splitext(sample_fqname)
    sample_fqname2 = re.sub('1\.(?:fq|fastq)$', '2{}'.format(ext), sample_fqname)

    sample_fq1 = [read1 for read1 in HTSeq.FastqReader(sample_fqname1)]
    sample_fq2 = {read2.name: read2 for read2 in HTSeq.FastqReader(sample_fqname2)}

    prev = -1
    count = 0
    num_sample = int(num_sample)
    while count < num_sample:
        read1 = random.choice(sample_fq1)
        if read1.name not in merge_read_names:
            read2_name = string.replace(read1.name, '/1', '/2')
            read2 = sample_fq2[read2_name]
            assert read1.name.split('/')[0] == read2.name.split('/')[0]
            read1.write_to_fastq_file(outfile1)
            read2.write_to_fastq_file(outfile2)
            count += 1
        else:
            print "DUPLICATE READS"









if __name__ == '__main__':
    main()