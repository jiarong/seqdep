#! /usr/bin/env python
# plot gamma cdf curve for one or more samples
# by gjr; 070314

import sys, os
import matplotlib
matplotlib.use('Pdf')
import matplotlib.pylab as plt
import numpy

def main():
    if len(sys.argv) != 4:
        mes = ('Usage: python {} "<file1.curve.model>,<fil2.curve.model>.."'
                                '"label1,label2.." out.pdf'
              )
        print >> sys.stderr, mes.format(os.path.basename(sys.argv[0]))
        sys.exit(1)

    infiles=[f.strip() for f in sys.argv[1].split(',')]
    labels =[l.strip() for l in sys.argv[2].split(',')]
    outfile=sys.argv[3]
    
    fig, ax = plt.subplots()
    for f, label in zip(infiles, labels):
        arr = numpy.loadtxt(f)
        total_xs = arr[:,0]
        x_c = total_xs[-1]
        xs = total_xs[:-1]
        total_ys = arr[:,1]
        y_c = total_ys[-1]
        ys = total_ys[:-1]
        ax.plot(xs, ys, label=label)
        ax.scatter(x_c, y_c, marker='s', facecolor='w')

    ax.legend(loc=2)
    ax.axhline(y=0.95, linestyle='--')
    ax.grid(False)
    ax.set_ylabel('Coverage')
    ax.set_xlabel('Number of reads')
    ax.set_xscale('log')
    #outfile='{}.py.pdf'.format(infiles[0])
    plt.savefig(outfile)

if __name__ == '__main__':
    main()
