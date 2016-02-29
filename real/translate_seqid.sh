ACC2NAME=./acc2name.py

cd gg13_type

$ACC2NAME gg13_type.mis ../ltp123_type/ltp123_name2acc.txt -r > gg13_type.mis.ACCNO

$ACC2NAME gg13_typemis.ACCNO ../gg13_nr99/gg13_name2acc.txt > gg13_type.mis.GGID

cd ../rdp11_type

$ACC2NAME rdp11_type.mis ../ltp123_type/ltp123_name2acc.txt -r > rdp11_type.mis.ACCNO

$ACC2NAME rdp11_type.mis.ACCNO <(zcat ../rdp11_type/rdp11_name2acc.txt.gz) > rdp11_type.mis.RDPID

cd ../ltp123_type

$ACC2NAME ltp123_type.mis ../ltp123_type/ltp123_name2acc.txt -r > ltp123_type.mis.ACCNO

cd ../slv123_type

$ACC2NAME slv123_type.mis ../ltp123_type/ltp123_name2acc.txt -r > slv123_type.mis.ACCNO

cd ../gg13_nr99

$ACC2NAME gg13_nr99.mis gg13_name2acc.txt -r > gg13_nr99.mis.ACCNO

cd ../slv123_nr99 

$ACC2NAME slv123_nr99.mis slv123_name2acc.txt -r > slv123_nr99.mis.ACCNO
