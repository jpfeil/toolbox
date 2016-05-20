#!/usr/bin/env python2.7

"""
Simple script for making pc_translation file from gencode reference. Used
for making a slimmed down reference data set for CI.

Removes fasta entries that are not in the gtf file.
"""

import sys
import HTSeq
import dill

def usage():
    print 'makeGencodeTranlation <gtf> <fasta file>'
    sys.exit()

def main():
    if len(sys.argv) != 3:
        usage()

    gtf, fa = sys.argv[1:3]

    try:
        with open('gencode_transcript_ids.pkl', 'r') as pklf:
            transcript_ids = dill.load(pklf)
    except IOError:
        transcript_ids = set()
        gtf_file = HTSeq.GFF_Reader(gtf)
        for feature in gtf_file:
            if feature.type == 'transcript':
                transcript_ids.add(feature.attr['transcript_id'])
        with open('gencode_transcript_ids.pkl', 'w') as pklf:
            dill.dump(transcript_ids, pklf)

    print '\n\n### Finished reading in GTF ###\n\n'

    fa_file = HTSeq.FastaReader(fa)
    fa_fout = open('output.fa', 'w')
    for fa in fa_file:
        ids = fa.name.split('|')
        if ids[0] in transcript_ids:
            fa.write_to_fasta_file(fa_fout)
    fa_fout.close()

if __name__ == '__main__':
    main()