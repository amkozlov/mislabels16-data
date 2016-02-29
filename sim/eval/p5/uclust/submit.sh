#!/bin/bash
#
#$ -cwd -V                              # Shift directories and export variables
#$ -q bridge.q                        # Select the queue
#$ -pe mvapich8 8                          # Set the parallel environment
#$ -l h_rt=24:00:00                     # Request the time for the job

#source /etc/profile.d/modules.sh

PCT=5
TAXREP=3
METHOD=uclust

#INDIR=/hits/sco/zhangje/data/epa_mislable_benchmark/$PCT
SIMDIR=/hits/sco/kozlov/mislabels/paper_v2/sim/full
TESTDIR=/hits/sco/kozlov/mislabels/tax_benchmark
INDIR=$SIMDIR
#OUTDIR=$SIMDIR/mis/p${PCT}/t$TAXREP
OUTDIR=$SIMDIR/mis/p${PCT}/$METHOD
TMPDIR=$SIMDIR/tmp/uclust

FASTA=$INDIR/ali/sim_full.fa
TAXFILE=$INDIR/tax/p${PCT}/simfull-p${PCT}-t${TAXREP}.tax

NAME=$METHOD-p$PCT-t$TAXREP

#module load python

. /hits/sco/kozlov/soft/qiime/activate_cluster.sh

# cd $TESTDIR

python $TESTDIR/script/mis_tests.py $FASTA $TAXFILE $METHOD $OUTDIR/$NAME.mis

mv ${JOB_NAME}.o${JOB_ID} $OUTDIR/${NAME}.log
