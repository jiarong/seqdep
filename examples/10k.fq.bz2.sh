#! /bin/bash
#PBS -q main
#PBS -l nodes=1:ppn=1,walltime=4:00:00
#PBS -l mem=21gb
#PBS -j oe

Readfile=../tests/test-data/10k.fq.bz2

# parameters to samples data for curve
Start=0
End=1
Steps=100
Reps=10

Trimmed_length=90
# Subsample reads from original read data
Subsample_seq_num=1000000
# Subsample representatives from abundance table
Subsample_rep_num=500000

Scriptdir=/mnt/home/guojiaro/Documents/software/gits/seqdep/scripts
Ksize=21

Hashsize=1000000000

#
# change parameters above
#

# big count use 2Byte for each kmer;
Hashmem=$(echo "scale=2; $Hashsize*2*4/1000000000" | bc)
echo "$Hashmem Gb is used for hashtable"


set -e
module load screed
module load NumPy
module load khmer/1.3

Bname=$(basename $Readfile)
Readfile=$(cd "$(dirname $Readfile)" && pwd)/$(basename $Readfile)

Outdir=$Readfile.seqdep.out
mkdir -p $Outdir
cd $Outdir

echo "loading counting table.."
python $Scriptdir/trim-seq-by-len.py $Readfile $Trimmed_length $Subsample_seq_num /2 - |\
tee $Bname.subseq.$Subsample_seq_num | \
load-into-counting.py -k $Ksize -x $Hashsize -N 4 $Bname.subseq.$Subsample_seq_num.ht -
echo "loading counting table finished.."

echo "start counting median kmer.."
count-median.py $Bname.subseq.$Subsample_seq_num.ht $Bname.subseq.$Subsample_seq_num $Bname.subseq.$Subsample_seq_num.K${Ksize}.medcount
echo "finish counting median kmer.."

echo "$Subsample_rep_num reps from $Subsample_seq_num are used for getting redundancy data before curve fitting"

Sweep_ratio=$(echo "scale=1; $Subsample_seq_num/$Subsample_rep_num" | bc)

python $Scriptdir/subsample-line.py $Bname.subseq.$Subsample_seq_num.K${Ksize}.medcount $Subsample_rep_num $Bname.subseq.$Subsample_seq_num.K${Ksize}.medcount.subrep.$Subsample_rep_num

echo "start making .npo file.."
python $Scriptdir/redunancy-curve.py $Bname.subseq.$Subsample_seq_num.K${Ksize}.medcount.subrep.$Subsample_rep_num  $Start $End $Steps $Reps > $Bname.subseq.$Subsample_seq_num.K${Ksize}.medcount.subrep.$Subsample_rep_num.curve
echo "finish making .npo file.."

Rscript $Scriptdir/seqdep.R $Bname.subseq.$Subsample_seq_num.K${Ksize}.medcount.subrep.$Subsample_rep_num.curve $Sweep_ratio > $Bname.subseq.$Subsample_seq_num.K${Ksize}.medcount.subrep.$Subsample_rep_num.curve.R.params

qstat -f ${PBS_JOBID}
