#!/bin/bash
#Configuration
#Input files:
#NEUT 5.4.0 NuMu files
#IN_FILES="/disk02/sklb3/t2k-sk.mc.release/t2kmc_19b/root/fhc/numu_x_numu/t2ksk19b.fqv4r0.fhc.numu_x_numu.00*.root"
#NEUT 5.4.0 Nue files
IN_FILES="/disk02/sklb3/t2k-sk.mc.release/t2kmc_19b/root/fhc/nue_x_nue/t2ksk19b.fqv4r0.fhc.nue_x_nue.*.root"
#Output files:
#OUT_FILE="/disk02/usr6/fshaker/radiative_correction_output/root_files/numu_kin.root"
OUT_FILE="/disk02/usr6/fshaker/radiative_correction_output/root_files/nue_kin.root"
#Max number of input files to read (for hadd command)
TOTAL_NB=500
#-------------------------------------------------------------------------------
STR=" "
_i=0

for _file in $IN_FILES
do
    if [[ $_i -eq $TOTAL_NB ]]; then
	break
    fi
    _i=$((_i+1))
    STR="$STR $_file"

done

echo "hadd $OUT_FILE $STR"
hadd -f $OUT_FILE $STR
