#! /usr/bin/env python
# take an abundance table and get the redundancy curve
# by gjr; 06/26/2014

import sys,os
import fileinput
import numpy

STEP_NUM=100
REP_NUM=100
START=0
END=0.05
SINGLETON=1
REDUN_RATIO_CUTOFF=0.2

def main():
    if len(sys.argv) != 6:
        mes = 'Usage: python {} <abundance.table> start end steps reps'
        print >> sys.stderr, mes.format(os.path.basename(sys.argv[0]))
        sys.exit(1)

    table_file = sys.argv[1]
    START=float(sys.argv[2])
    END=float(sys.argv[3])
    STEP_NUM=int(sys.argv[4])
    REP_NUM=int(sys.argv[5])

    # split at linear space
    #subsample_rate_list = numpy.linspace(START, END, STEP_NUM)
 
    # split at log space
    if START == 0:  # avoid log(0)
        START = 0.0001
        subsample_rate_list = numpy.logspace(
                           numpy.log10(START),
                           numpy.log10(END),
                           STEP_NUM,
                           base=10,
                                                )
    total_read_num = 0
    redun_read_num = 0
    counting_array_list = [numpy.zeros(REP_NUM)] * STEP_NUM

    line_num = 0
    for line in fileinput.input(table_file):
        line_num += 1
        if line.startswith('#'):
            continue
        line = line.rstrip()
        #name, cnt = line.split('\t')[:2]
        name, cnt = line.split()[:2]
        cnt = int(cnt)
        total_read_num += 1
        if cnt > SINGLETON:
            redun_read_num += 1
        else:
            continue

        if total_read_num % 100000 ==0:
            print >> sys.stderr, '{} reads processed..'.format(total_read_num)
            ratio = redun_read_num*1.0/total_read_num
            if  ratio < REDUN_RATIO_CUTOFF:
                _m = ('** Redundancy ratio {} falls below 20%, '
                      'more sequencing is needed for sequencing depth '
                      'estimation'
                     )
                print >> sys.stderr, _m.format(ratio)
                sys.exit(1)

        # interate throught sample rates
        for ind, p in enumerate(subsample_rate_list):
            if cnt*p <=1: 
                continue                 
            else:
                # if redundant, update the array of a sample rate
                ys = numpy.random.binomial(1, p, REP_NUM)
                counting_array_list[ind] = counting_array_list[ind] + ys

    subsample_num_list = subsample_rate_list * total_read_num
    _temp_list = zip(subsample_num_list, counting_array_list)
    for subsample_num, counting_array in _temp_list:
        if subsample_num == 0:
            assert sum(counting_array) == 0
            redun_ratio_array = counting_array
        else:
            redun_ratio_array = counting_array/subsample_num
        redun_ratio_array[redun_ratio_array > 1] = 1
        mean = numpy.mean(redun_ratio_array)
        std = numpy.std(redun_ratio_array)
        q1, q2, q3 = numpy.percentile(redun_ratio_array, [25, 50, 75])
        print '{}\t{}\t{}\t{}\t{}\t{}'.format(
                   subsample_num, mean, std, q1, q2, q3,
                                             )

if __name__ == '__main__':
    main()
