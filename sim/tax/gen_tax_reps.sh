SCRIPTDIR=../script
SIMDIR=.

for pct in 1 5;
do

  for rep in 1 2 3;
  do

    python $SCRIPTDIR/tax_simulator.py $SIMDIR $pct $rep

  done

done
