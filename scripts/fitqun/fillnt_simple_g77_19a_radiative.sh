#!/bin/sh
#
#   simple shell script to make official ntuple from ZBS file
#   
#   usage : 
#      fillnt_simple.sh [-o (hbook)] [-pc] (zbsfile1) (zbsfile2) ...
#
#   output : 
#      fillnt.hbk  or  specified file by "-o" option
#
#   option :
#      -o       specify output hbook file
#      -pc      set this option in case of reading PC data
#


#
#   uncomment this when you want to change default flag for filling ntuple
#   set "yes"/"no" to active/inactive these environmental variables
#
# FILL_DST="no"         ;   export FILL_DST
# FILL_RING="no"        ;   export FILL_RING
# FILL_RNGSTUDY="no"    ;   export FILL_RNGSTUDY
# FILL_COMPARE="yes"    ;   export FILL_COMPARE
# FILL_MUE="no"         ;   export FILL_MUE
# FILL_EVISCOMP="yes"   ;   export FILL_EVISCOMP
# FILL_LOWF="yes"       ;   export FILL_LOWF
# FILL_PI0FIT="no"      ;   export FILL_PI0FIT
# FILL_PI0FIT2009="yes" ;   export FILL_PI0FIT2009
# FILL_MC_DST="no"      ;   export FILL_MC_DST
# FILL_SCAN="no"        ;   export FILL_SCAN
# FILL_MSFIT="yes"      ;   export FILL_MSFIT
# FILL_STMU="yes"       ;   export FILL_STMU
# FILL_SKHEAD="yes"     ;   export FILL_SKHEAD
# FILL_NUCEFF="no"      ;   export FILL_NUCEFF
# FILL_CCQE="yes"       ;   export FILL_CCQE
FILL_SECONDARY="yes"  ;   export FILL_SECONDARY
#FILL_NEUT="yes"       ;   export FILL_NEUT
# FILL_FSIHIST="yes"    ;   export FILL_FSIHIST
# FILL_CATPCVARS="yes"  ;   export FILL_CATPCVARS
FILL_NEUTSKDET="yes"  ;   export FILL_NEUTSKDET
# FILL_PROTONID="no"    ;   export FILL_PROTONID
# FILL_NTAG="no"    	;   export FILL_NTAG
FILL_FITQUN="1"       ;   export FILL_FITQUN # Set it to 2 to output more detailed fit information 

#
#  for ATMPD extention
#
# FILL_UPMU="no"        ;   export FILL_UPMU
# FILL_OSC_WEIGHT="no"  ;   export FILL_OSC_WEIGHT
# FILL_PC5TH="yes"      ;   export FILL_PC5TH
# FILL_OSC3D="no"       ;   export FILL_OSC3D
# FILL_OSC3D_TEST="yes" ;   export FILL_OSC3D_TEST
# FILL_LIVETIME="no"    ;   export FILL_LIVETIME
# FILL_ASTRO="no"       ;   export FILL_ASTRO
# FILL_SOLARACT="no"    ;   export FILL_SOLARACT
# FILL_LE="yes"         ;   export FILL_LE
# FILL_PROTON="yes"     ;   export FILL_PROTON
# FILL_FC1="yes"        ;   export FILL_FC1
# FILL_TAU="yes"        ;   export FILL_TAU
# FILL_MVFIT="yes"      ;   export FILL_MVFIT
# FILL_GOODRUN="yes"    ;   export FILL_GOODRUN
#
# FLAG_PCUSED="yes"     ;   export FLAG_PCUSED
# FLAG_INPMT="no"       ;   export FLAG_INPMT
# OSWEIGHT_PARAM1="2.5E-3"    ;  export OSWEIGHT_PARAM1 
# OSWEIGHT_PARAM2="1.0"       ;  export OSWEIGHT_PARAM2

#
#  for T2K extention
#
# FILL_K2K="yes"        ;   export FILL_K2K
# FILL_T2K="yes"        ;   export FILL_T2K
# FILL_T2KREDUC="yes"   ;   export FILL_T2KREDUC
# FILL_T2KNUESYS="yes"   ;   export FILL_T2KNUESYS

#
#  Ad-hoc pre-procedure 
#
# DO_PREFILLNT="yes"       ;   export DO_PREFILLNT


#
#  flasher database and PDF file for T2K flasher cut
#  if you want to specify different file, uncomment and modify this
#
# T2KFC4PDFFILE="/home/atmpd/const/fc4/fc4_pdf_sk4.0903.root"  ; export T2KFC4PDFFILE
# T2KFC4FLASHERDB="/home/sklb/reduction/shell/00_mkflasherdb/formc/fc4_makedb.root" ; export T2KFC4FLASHERDB


#
#  used file (good subrun file) for ntuple_module_goodrun
#  
# ATMPDUSEDRUNFILE="/home/atmpd/summary/LIVETIME_SK4/usedrun_mar13.out" ; export ATMPDUSEDRUNFILE

#
# default 
#
output=fillnt.hbk
pc=0
option=
host=

#
# check options
#
nargs=$#
i=1
while [ $i -le $nargs ] ; do 
  case $1 in 
    -o ) 
     shift 
     output=$1
     shift 
     ;;
    -pc ) 
     shift
     pc=1
     ;;
  esac
  i=`expr $i + 1`
done



#
#  environmental valiables
#
FILLNT=/home/skofl/sklib_g77/atmpd_19a/bin//fillnt
#
KAMFLUX=/home/skofl/sklib_g77/neut_5.4.0/src/kamflux
NEUTFLUX=/home/skofl/sklib_g77/atmpd_19a/const/
SACTFILE=/home/skofl/sklib_g77/atmpd_19a/const//neutron.table
export KAMFLUX NEUTFLUX SACTFILE
#
#
INPUTBNK=0
PCFLAG=$pc
AUTOBNK=3
COMPBNK2=3
COMPBNK3=4
export INPUTBNK PCFLAG AUTOBNK COMPBNK2 COMBNK3 KAMFLUX NEUTFLUX

SKPATH=/home/skofl/sklib_g77/atmpd_19a/const/:/home/skofl/sklib_g77/skofl_19a/const/:/home/skofl/sklib_g77/skam/const:${NEUTFLUX}
export SKPATH


RFLIST=rflist.$$.$USER.`hostname`
export RFLIST



#
# make RFLIST
#
#
echo '10{' > $RFLIST 
while [ $# -ge 1 ] ; do
  echo \{\"$1\",LOCAL,,RED,,,\"recl=5670 status=old\"\} >> $RFLIST
  shift
done
echo '}' >> $RFLIST 

cat <<EOF >> $RFLIST
51{
{"$KAMFLUX/kamflx.data/KAMM09",LOCAL,,RED,,,"form=formatted"}
{"$KAMFLUX/kamflx.data/KAMM08",LOCAL,,RED,,,"form=formatted"}
{"$KAMFLUX/kamflx.data/KAMM07",LOCAL,,RED,,,"form=formatted"}
{"$KAMFLUX/kamflx.data/KAMM06",LOCAL,,RED,,,"form=formatted"}
{"$KAMFLUX/kamflx.data/KAMM05",LOCAL,,RED,,,"form=formatted"}
{"$KAMFLUX/kamflx.data/KAMM04",LOCAL,,RED,,,"form=formatted"}
{"$KAMFLUX/kamflx.data/KAMM03",LOCAL,,RED,,,"form=formatted"}
{"$KAMFLUX/kamflx.data/KAMM02",LOCAL,,RED,,,"form=formatted"}
{"$KAMFLUX/kamflx.data/KAMM01",LOCAL,,RED,,,"form=formatted"}
{"$KAMFLUX/kamflx.data/KAM00",LOCAL,,RED,,,"form=formatted"}
{"$KAMFLUX/kamflx.data/KAM01",LOCAL,,RED,,,"form=formatted"}
{"$KAMFLUX/kamflx.data/KAM02",LOCAL,,RED,,,"form=formatted"}
{"$KAMFLUX/kamflx.data/KAM03",LOCAL,,RED,,,"form=formatted"}
{"$KAMFLUX/kamflx.data/KAM04",LOCAL,,RED,,,"form=formatted"}
{"$KAMFLUX/kamflx.data/KAM05",LOCAL,,RED,,,"form=formatted"}
{"$KAMFLUX/kamflx.data/KAM06",LOCAL,,RED,,,"form=formatted"}
{"$KAMFLUX/kamflx.data/KAM07",LOCAL,,RED,,,"form=formatted"}
{"$KAMFLUX/kamflx.data/KAM08",LOCAL,,RED,,,"form=formatted"}
{"$KAMFLUX/kamflx.data/KAM09",LOCAL,,RED,,,"form=formatted"}
{"$KAMFLUX/kamflx.data/KAM10",LOCAL,,RED,,,"form=formatted"}
}
52{
{"$KAMFLUX/kamflx.data/NFLX01",LOCAL,,RED,,,"form=formatted"}
{"$KAMFLUX/kamflx.data/NFLX02",LOCAL,,RED,,,"form=formatted"}
{"$KAMFLUX/kamflx.data/NFLX03",LOCAL,,RED,,,"form=formatted"}
{"$KAMFLUX/kamflx.data/NFLX04",LOCAL,,RED,,,"form=formatted"}
{"$KAMFLUX/kamflx.data/NFLX05",LOCAL,,RED,,,"form=formatted"}
{"$KAMFLUX/kamflx.data/NFLX06",LOCAL,,RED,,,"form=formatted"}
{"$KAMFLUX/kamflx.data/NFLX07",LOCAL,,RED,,,"form=formatted"}
{"$KAMFLUX/kamflx.data/NFLX08",LOCAL,,RED,,,"form=formatted"}
{"$KAMFLUX/kamflx.data/NFLX09",LOCAL,,RED,,,"form=formatted"}
{"$KAMFLUX/kamflx.data/NFLX10",LOCAL,,RED,,,"form=formatted"}
}
61{{"$NEUTFLUX/hkkm06mt.dat",LOCAL,,RED,,,"form=formatted"}}
62{{"$NEUTFLUX/bartol03.dat",LOCAL,,RED,,,"form=formatted"}}
63{{"$NEUTFLUX/fluka03.dat",LOCAL,,RED,,,"form=formatted"}}
64{{"$NEUTFLUX/hkkm03mt.dat",LOCAL,,RED,,,"form=formatted"}}
65{{"$NEUTFLUX/hkkm2.dat",LOCAL,,RED,,,"form=formatted"}}
68{{"$NEUTFLUX/hkkm11mt.dat",LOCAL,,RED,,,"form=formatted"}}
69{{"$NEUTFLUX/hkkm11low.dat",LOCAL,,RED,,,"form=formatted"}}
#
# these file path (ID=66,67,71,72,73,76,77,78) can be found by 
# "findconst routines"
#66{{"$NEUTFLUX/fit_fluka02.dat",LOCAL,,RED,,,"form=formatted"}}
#67{{"$NEUTFLUX/fit_fluka02_mx.dat",LOCAL,,RED,,,"form=formatted"}}
#71{{"$NEUTFLUX/honda96low.dat",LOCAL,,RED,,,"form=formatted"}}
#72{{"$NEUTFLUX/honda97mid.dat",LOCAL,,RED,,,"form=formatted"}}
#73{{"$NEUTFLUX/honda96high.dat",LOCAL,,RED,,,"form=formatted"}}
#76{{"$NEUTFLUX/gaisser96.dat",LOCAL,,RED,,,"form=formatted"}}
#77{{"$NEUTFLUX/honda96low.dat",LOCAL,,RED,,,"form=formatted"}}
# 'hkkm03mt.dat' is for new neut version 
#78{{"$NEUTFLUX/hkkm03mt.dat",LOCAL,,RED,,,"form=formatted"}}
# 'hkkm2.dat' is for SK-I neut version 
#78{{"$NEUTFLUX/hkkm2.dat",LOCAL,,RED,,,"form=formatted"}}
#
87{{"$NEUTFLUX/honda96low.dat",LOCAL,,RED,,,"form=formatted"}}
88{{"$NEUTFLUX/bess71.dat",LOCAL,,RED,,,"form=formatted"}}
89{{"$NEUTFLUX/bess74.dat",LOCAL,,RED,,,"form=formatted"}}
90{{"$output",LOCAL,,WRT,,,"access=direct recl=4096 form=unformatted status=unknown"}}
EOF

echo ' '
echo ' '
cat $RFLIST
echo ' '
echo ' '


$FILLNT

/bin/rm $RFLIST
