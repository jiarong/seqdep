#! /usr/bin/env python
# subsample lines
# by gjr; 012315

# python <thisfile> <file.medcount> num <outfile>

import sys
import os
import random

def get_file_line_num(f):
    with open(f) as fp:
        i = -1
        for i, l in enumerate(fp):
            pass

    assert i > -1, "{} is empty".format(os.path.basename(f))
    return i + 1

def main():
    if len(sys.argv) != 4:
        mes = ('Usage: python {} <file.medcount> num <outfile>\n'
                'use "-" for <outfile> if to stdout')
        print >> sys.stderr, mes.format(sys.argv[0])
        sys.exit(1)

    infile = sys.argv[1]
    sub_num = int(sys.argv[2])
    outfile = sys.argv[3]

    total = get_file_line_num(infile)
    if sub_num > total:
        print >> sys.stderr, '*** Subsample size ({}) is larger than total ({}), subsample size changed to {}'.format(sub_num, total, total)
        sub_num = total

    to_pick_list = random.sample(xrange(total), sub_num)
    to_pick_set = set(to_pick_list)

    with open(infile) as fp, open(outfile, 'wb') as fw:
        for n, line in enumerate(fp):
            if n in to_pick_set:
                fw.write(line)

if __name__ == '__main__':
   main()
    


