#!/bin/bash

export LANG C
export LC_ALL C


export CXX=g++34
export CPP="gcc34 -E"
export CC=gcc34
export F77=g77

export CERN=/home/skofl/sklib_g77/cern
export CERN_LEVEL=2005
export CERN_ROOT=${CERN}/${CERN_LEVEL}

export SKOFL_ROOT=/home/skofl/sklib_g77/skofl_19a
export ATMPD_ROOT=/home/skofl/sklib_g77/atmpd_19a
export NEUT_ROOT=/home/skofl/sklib_g77/neut_5.4.0
export ROOTSYS=/home/skofl/sklib_g77/root_v5.28.00h

export FITQUN_ROOT=/home/skofl/sklib_g77/atmpd_14c/src/recon/fitqun/


#if [ ${LD_LIBRARY_PATH:=${ROOTSYS}/lib} = "${ROOTSYS}/lib" ];
#then
#echo "I am in the LD_LIBRRAY_PATH setup"    
export LD_LIBRARY_PATH=${SKOFL_ROOT}/lib:${ROOTSYS}/lib:$LD_LIBRARY_PATH
#fi
#fsamir replaced that line with another from /home/skofl/sklib_gcc4.8.5/skofl_16c/env.sh
#export PATH=${SKOFL_ROOT}/bin:${ROOTSYS}/bin:$PATH
export PATH=${SKOFL_ROOT}/bin:${ATMPD_ROOT}/bin:${ROOTSYS}/bin:${CERN_ROOT}/bin:$PATH

#fsamir changed the library path as well
#if [ "X${LD_LIBRARY_PATH}" = "X" ];
#then
#    export LD_LIBRARY_PATH=${SKOFL_ROOT}/lib:${ROOTSYS}/lib
#else
#    export LD_LIBRARY_PATH=${SKOFL_ROOT}/lib:${ROOTSYS}/lib:$LD_LIBRARY_PATH
#fi
