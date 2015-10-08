#!/usr/bin/env python

import sys
import os
import re
import screed
from screed import fasta
from screed import fastq

CUTOFF = 0.2
K = 4
def count_uniq_kmer(seq):
    length = len(seq)
    kmer_set = set()
    for i in xrange(length-K+1):
        kmer = seq[i:i+K-1]
        kmer_set.add(kmer)

    return len(kmer_set)
        
def main():
    '''
    Usage: python <thisfile> <infile> <outfile>
    '''
    if len(sys.argv) != 3:
        mes = ('Usage: python {} <infile> <outfile>')
        print >> sys.stderr, mes.format(os.path.basename(sys.argv[0]))
        sys.exit(1)

    infile = sys.argv[1]
    outfile = sys.argv[2]

    try:
        if outfile == '-':
            fw = sys.stdout
        else:
            fw = open(outfile, 'wb')

        lowcomp_fw = open('low_complexity.fa', 'wb')

        for n, record in enumerate(screed.open(infile)):
            name = record['name']
            seq = record['sequence']
            uniq_kmer_count = count_uniq_kmer(seq)
            if uniq_kmer_count * 1.0/(len(seq) - K + 1) < CUTOFF:
                lowcomp_fw.write('>{}\n{}\n'.format(name, seq)) #fasta output
                continue

            fw.write('>{}\n{}\n'.format(name, seq)) #fasta output

        try:
            n
        except NameError:
            print >> sys.stderr, '*** No seqs are in seqfile'

    except IOError as err:
        if outfile == '-':
            pass
        else:
            print >> sys.stderr, '*** {}'.format(err)
            sys.exit(1)

if __name__ == '__main__':
    main()
