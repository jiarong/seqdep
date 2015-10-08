#!/usr/bin/env python

import sys
import os
import re
import screed
from screed import fasta
from screed import fastq


def main():
    '''
    Usage: python <thisfile> <infile> length numseq2keep tag <outfile>
    '''
    if len(sys.argv) != 6:
        mes = ('Usage: python {} <infile> length numseq2keep tag <outfile>\n'
                '*** tag can be used screen OUT seq names\n')
        print >> sys.stderr, mes.format(os.path.basename(sys.argv[0]))
        sys.exit(1)

    infile = sys.argv[1]
    length = int(sys.argv[2])
    num = int(sys.argv[3])
    tag = sys.argv[4]
    outfile = sys.argv[5]

    try:
        if infile == '-':
            fp = sys.stdin
        else:
            fp = open(infile)

        if outfile == '-':
            fw = sys.stdout
        else:
            fw = open(outfile, 'wb')

        for n, record in enumerate(fasta.fasta_iter(fp)):
            if n == num:
                break
            name = record['name']
            seq = record['sequence']

            if len(seq) < length:
                continue

            if tag in name:
                continue

            new_seq = seq[:length]

            fw.write('>{}\n{}\n'.format(name, new_seq)) #fasta output

        try:
            n
        except NameError:
            print >> sys.stderr, '*** No seqs are in seqfile'

        if n < num:
            mes = '*** Not enough seqs in {} ({} < {}), only {} subsampled'
            print >> sys.stderr, mes.format(os.path.basename(infile),
                                              n, num, n)

    except IOError as err:
        print >> sys.stderr, err
        fw.close()
        sys.exit(1)

if __name__ == '__main__':
    main()
