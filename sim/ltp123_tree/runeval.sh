
# path to RAxML exectable
RAXML=~/raxml/raxmlHPC-PTHREADS-AVX

# Optimize and print model parameters for each partition
$RAXML -s ltp123_full.phy -t reftree.nw -q ltp123_part.txt -m GTRGAMMA -n evalR -T 16 -f e -M 

# Split alignment to get per-partition subalignments
$RAXML -s ltp123_full.phy -q ltp123_part.txt -m GTRGAMMA -n splitR -T 16 -f s
