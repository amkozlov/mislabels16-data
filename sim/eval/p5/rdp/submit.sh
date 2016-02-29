#!/bin/bash
#
#$ -cwd -V                              # Shift directories and export variables
#$ -q bridge.q                        # Select the queue
#$ -pe mvapich8 8                          # Set the parallel environment
#$ -l h_rt=24:00:00                     # Request the time for the job
#$ -N sim5_rdp

PCT=5
TAXREP=$1
METHOD=rdp

#INDIR=/hits/sco/zhangje/data/epa_mislable_benchmark/$PCT
SIMDIR=/hits/sco/kozlov/mislabels/paper_v2/sim/full
TESTDIR=/hits/sco/kozlov/mislabels/tax_benchmark
INDIR=$SIMDIR
#OUTDIR=$SIMDIR/mis/p${PCT}/t$TAXREP
OUTDIR=$SIMDIR/mis/p${PCT}/$METHOD/t$TAXREP
TMPDIR=$SIMDIR/tmp

FASTA=$INDIR/ali/sim_full.fa
TAXFILE=$INDIR/tax/p${PCT}/simfull-p${PCT}-t${TAXREP}.tax

NAME=$METHOD-p$PCT-t$TAXREP

module load python

. /hits/sco/kozlov/soft/qiime/activate_cluster.sh

# cd $TESTDIR

python $TESTDIR/script/mis_tests2.py $FASTA $TAXFILE rdp $OUTDIR/$NAME.mis.$2 $2

mv ${JOB_NAME}.o${JOB_ID} $OUTDIR/${NAME}.log.$2
