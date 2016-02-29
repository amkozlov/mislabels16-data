sed 's/unclassified_["A-Za-z0-9 _]*;//g' sativa_in.tax | sed 's/;[^;]*$//g' | sed 's/;uncultured$//g' > sativa_in.tax.genus
