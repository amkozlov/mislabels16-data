grep -v -f <(cut -f1,2,3,4 slv123_nr99.mis | grep -e "Incert" -e "uncultured" -e "unclassified" -e "Unknown" | cut -f1) slv123_nr99.mis > slv123_nr99.mis.corrected
