#!/bin/bash
#
#$ -cwd -V                              # Shift directories and export variables
##$ -q bridge.q,bridge_fat.q                        # Select the queue
#$ -q ivyf.q                        # Select the queue
#$ -pe mvapich16 16                          # Set the parallel environment
#$ -l h_rt=01:00:00                     # Request the time for the job
#$ -N sim5_full_rep

DATASET=ltp123_full
#MODE=fast
MODE=thorough
TAXREP=$1
REP=1
RUNS=1
PCT=5
METHOD=sativa

SATDIR=/hits/sco/kozlov/sativa
SIMDIR=/hits/sco/kozlov/mislabels/paper_v2/sim/full
#INDIR=$SIMDIR
INDIR=$SIMDIR/mis/p${PCT}/sativa/t$TAXREP
OUTDIR=$SIMDIR/mis/p${PCT}/sativa5/t$TAXREP
TMPDIR=$SIMDIR/tmp

TAXFILE=$INDIR/tax/p${PCT}/simfull-p$PCT-t${TAXREP}.tax
ALIFILE=$INDIR/ali/sim_full.phy

NAME=$METHOD-p${PCT}-t${TAXREP}-${REP}
#JSON=$OUTDIR/$NAME.refjson
#CFG=$SIMDIR/epac.cfg.cat
CFG=$SATDIR/sativa.cfg

$SATDIR/sativa.py -t $TAXFILE -s $ALIFILE -x BAC -c $CFG -o $OUTDIR -n $NAME -tmpdir $TMPDIR -m $MODE -N $RUNS -v -T 16

