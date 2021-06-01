#!/bin/bash

# Usage :  
#  prompt> skdetsim.sh [card name] [output file] [input file] [# of sukipped event]
#

## assume skdetsim exists in the same directory

DIR=$PWD 
if [[ "$DIR" = "$0" ]]
then
   DIR="."
fi
## SUBRUN number is input from enviromment variable.
## This should be editted by user.

export MC_SUBNUM=0

## Random seed file for NEUT.
## This should be controlled by user.

#setenv RANFILE /home/atmpd/skdetsim/skdetsim-v13p90_mar16/random.tbl.000
export RANFILE="/home/atmpd/skdetsim/skdetsim-v13p90_dec17/random.tbl.000"

echo "$#"

if [[ $# -gt 0 ]]
then
    CARD=$1
else
    CARD=supersim.card
fi

if [[ $# -gt 1 ]]
then
    FNAME_OUT=$2
else
    FNAME_OUT=dummy.output
fi

if [[ $# -gt 2 ]]
then
    FNAME_IN=$3
else
    FNAME_IN=dummy.input
fi

if [[ $# -gt 3 ]]
then
    /home/atmpd/skdetsim/skdetsim-v13p90_dec17/bins/skdetsim_high $CARD $FNAME_OUT $FNAME_IN $4
else
    /home/atmpd/skdetsim/skdetsim-v13p90_dec17/bins/skdetsim_high $CARD $FNAME_OUT $FNAME_IN
fi

