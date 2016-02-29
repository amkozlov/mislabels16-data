#!/bin/bash

#/hits/sco/kozlov/soft/Indelible/INDELibleV1.03/src/indelible > log

INDELIBLE=~/INDELibleV1.03/src/indelible

realphy=../ltp123_tree/ltp123_full.phy

SCRIPTDIR=../script

for part in cons flanks v1 v2 v3 v4 v5 v6 v7 v8 v9;
do
#  ./gen_control.sh $part

  ln -fs control.txt.$part control.txt

  $INDELIBLE > INDEL.log.$part

  echo $part
  
# manual control for gap% in real and simulated alignments
#  $SCRIPTDIR/gaps.py sim_${part}_TRUE.phy
#  $SCRIPTDIR/gaps.py ${realphy}.${part}.phy

echo ""
done

$SCRIPTDIR/merge.py sim_*_TRUE.phy > sim_full.phy

$SCRIPTDIR/merge_fasta.py sim_*.fas > sim_full.fa
