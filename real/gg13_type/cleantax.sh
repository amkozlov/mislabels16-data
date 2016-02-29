grep -v "Unclassified;" sativa_in.tax | sed 's/unclassified_["A-Za-z0-9 _]*;//g' | sed 's/;[^;]*$//g' > sativa_in.tax.genus
