#!/bin/bash

if [ ! -d "$EGSIMBIN" ]; then
  echo "EGSIMBIN variable was not set or directory does not exist"
  exit
fi
flist=("admmat" "admmatI" "admmatR" "linepow" "statsol")

# run 118-bus model
cd case118
echo
echo "--------------------------------------------------------"
echo "Solving the 118 bus model"
echo "--------------------------------------------------------"
${EGSIMBIN}/egsim -d ieee118cdf.txt -a full -p

for fn in ${flist[@]}; do
    if diff "${fn}.dat" "${fn}_118.dat" &> /dev/null ; then
        echo "${fn}.dat matches saved solution"
    else
        echo "${fn}.dat does not match ${fn}_118.dat"
    fi
done

if command -v python &>/dev/null; then
    python mkPlots.py
else
    echo "python is not available! will not run mkPlots.py"
fi
cd ..

# run 300-bus model
cd case300
echo 
echo "--------------------------------------------------------"
echo "Solving the 300 bus model"
echo "--------------------------------------------------------"
${EGSIMBIN}/egsim -d ieee300cdf.txt -a full -p

for fn in ${flist[@]}; do
    if diff "${fn}.dat" "${fn}_300.dat" &> /dev/null ; then
        echo "${fn}.dat matches saved solution"
    else
        echo "${fn}.dat does not match ${fn}_300.dat"
    fi
done

if command -v python &>/dev/null; then
    python mkPlots.py
else
    echo "python is not available! will not run mkPlots.py"
fi
cd ..
