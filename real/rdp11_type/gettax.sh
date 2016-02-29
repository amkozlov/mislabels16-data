# STEP 1. Build NCBI accno <-> RDP seqid mapping

grep -e "^ACCESSION " -e "^LOCUS " current_Archaea_unaligned.gb | sed 'N;s/\n/ /' | sed 's/LOCUS[ ]*\([^ ]*\).*ACCESSION[ ]*\([^ ]*\).*/\1\t\2/g' > id2acc_arc.txt

grep -e "^ACCESSION " -e "^LOCUS " current_Bacteria_unaligned.gb | sed 'N;s/\n/ /' | sed 's/LOCUS[ ]*\([^ ]*\).*ACCESSION[ ]*\([^ ]*\).*/\1\t\2/g' > id2acc_bac.txt

cat id2acc_arc.txt id2acc_bac.txt > rdp11_name2acc.txt

rm id2acc_arc.txt id2acc_bac.txt


# STEP 2: Extract taxonomic paths from FASTA

grep "^>" current_Bacteria_unaligned.fa | sed 's/>\([^ ]*\).*Lineage=\(.*\)$/\1\t\2/g' > current_bac.tax.1

sed 's/Root;rootrank;//g; s/;domain;/;/g; s/;phylum;/;/g; s/;class;/;/g; s/;subclass;/;/g; s/;order;/;/g; s/;suborder;/;/g; s/;family;/;/g; s/;genus//g; s/;$//g' current_bac.tax.1 > current_bac.tax

grep "^>" current_Archaea_unaligned.fa | sed 's/>\([^ ]*\).*Lineage=\(.*\)$/\1\t\2/g' > current_arc.tax.1

sed 's/Root;rootrank;//g; s/;domain;/;/g; s/;phylum;/;/g; s/;class;/;/g; s/;subclass;/;/g; s/;order;/;/g; s/;suborder;/;/g; s/;family;/;/g; s/;genus//g; s/;$//g' current_arc.tax.1 > current_arc.tax

cat current_bac.tax current_arc.tax > current_full.tax

rm current_bac.tax.1 current_bac.tax current_bac.tax current_arc.tax


# STEP 3: convert RDP seqids to ARB names

../filter_tax.py ../ltp123_type/ltp123_name2acc.txt rdp11_name2acc.txt current_full.tax | sed 's/;unclassified_.*$//g' > sativa_in.tax

gzip rdp11_name2acc.txt
