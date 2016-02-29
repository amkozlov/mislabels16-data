#!/bin/bash

PCT=$1
METHOD=$2
CUTOFF=$3

EPATAXDIR=`cd "../script/epatax"; pwd`
TAXDIR=`cd "../tax"; pwd`

TRUETAX=$TAXDIR/p0/true.tax

for TAXREP in 1 2 3;
do

 FULLTAX=$TAXDIR/p$PCT/simfull-p$PCT-t$TAXREP.tax
 TRUEMIS=$TAXDIR/p$PCT/simfull-p$PCT-t$TAXREP.true.tax

 MISFILE=p$PCT/$METHOD/t$TAXREP/$METHOD-p$PCT-t$TAXREP.mis
 FULLASS=p$PCT/$METHOD/t$TAXREP/$METHOD-p$PCT-t$TAXREP.eval$CUTOFF.full
 TPFILE=p$PCT/$METHOD/t$TAXREP/$METHOD-p$PCT-t$TAXREP.eval$CUTOFF.tp

 $EPATAXDIR/calc_stats_mis.py -m $MISFILE -t $FULLTAX -g $TRUEMIS -c $EPATAXDIR/epatax.cfg -x $CUTOFF

 $EPATAXDIR/calc_stats.py -a $FULLASS -t $TRUETAX  -c $EPATAXDIR/epatax.cfg
   
done

cd p$PCT/$METHOD

$EPATAXDIR/summarize_stats.py $PCT $METHOD ident $CUTOFF

$EPATAXDIR/summarize_stats.py $PCT $METHOD corr $CUTOFF

cd ../..
