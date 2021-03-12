#!/bin/bash

FILES="/disk01/sklb2/t2kmc/t2kmc_14c/t2k_14a_root/nu-mode/numu_x_numu/root/*.root"
TOTAL_NB=30
STR=" "
_i=0
OUT_FILE="/disk02/usr6/fshaker/numu_kin.root"
for _file in $FILES
do
    if [[ $_i -eq $TOTAL_NB ]]; then
	break
    fi
    _i=$((_i+1))
    STR="$STR $_file"

done

echo "hadd $OUT_FILE $STR"
hadd -f $OUT_FILE $STR
