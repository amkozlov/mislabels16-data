#!/bin/bash

INFOFILE=../ltp123_tree/RAxML_info.eval
MODFILE=control

for part in cons flanks v1 v2 v3 v4 v5 v6 v7 v8 v9;
do
  python ../script/gen_control.py $INFOFILE $PART > $MODFILE

  awk 'FNR==NR{s=(!s)?$0:s RS $0;next} /<<MODELS>>/{sub(/<<MODELS>>/, s)} 1' $MODFILE control.templ > control.txt.$PART
done;
