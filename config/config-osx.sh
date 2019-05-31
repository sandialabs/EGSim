#!/bin/bash

usage ()
{
  echo "Usage : $0 -d <INSTALL_DIR> -c <gnucompiler> -h"
  exit
}

egsimdir="${PWD}/.."
compilerpath="/opt/local/bin"
installdir="${PWD}/../install"
ctype="gnu8"

while getopts ":c:d:h" opt; do
  case $opt in
    c) ctype="$OPTARG"
    ;;
    d) installdir="$OPTARG"
    ;;
    h) usage
    ;;
    \?) echo "Invalid option -$OPTARG" >&2; usage
    ;;
  esac
done

echo "============================================"
echo "Compiling EGSim with:"
echo " - compiler path: ${compilerpath}"
echo " - compilers:     ${ctype}"
echo " - install:       ${installdir}"
echo "============================================"

GNUROOT=$HOME/local

if [ "${ctype}" == "gnu55" ]; then
  compilerpath="${GNUROOT}/gcc-5.5.0/bin"
  gnuID=5.5.0
elif [ "${ctype}" == "gnu64" ]; then
  compilerpath="${GNUROOT}/gcc-6.4.0/bin"
  gnuID=6.4.0
elif [ "${ctype}" == "gnu72" ]; then
  compilerpath="${GNUROOT}/gcc-7.2.0/bin"
  gnuID=7.2.0
elif [ "${ctype}" == "gnu8" ]; then
  compilerpath="/opt/local/bin"
  gnuID=mp-8
else
  echo "Unknown compiler: ${ctype}"
  exit
fi

cmake -DCMAKE_INSTALL_PREFIX:PATH=${installdir} \
      -DCMAKE_Fortran_COMPILER=${compilerpath}/gfortran-${gnuID} \
      -DCMAKE_CXX_COMPILER=${compilerpath}/g++-${gnuID} \
      ${egsimdir}

