#!/bin/bash
#NEUT 5.4.0 files
FILES="/disk02/sklb3/t2k-sk.mc.release/t2kmc_19b/root/fhc/numu_x_numu/t2ksk19b.fqv4r0.fhc.numu_x_numu.00*.root"
#NEUT 5.3.2 files
#FILES="/disk01/sklb2/t2kmc/t2kmc_14c/t2k_14a_root/nu-mode/numu_x_numu/root/*.root"
TOTAL_NB=500
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
