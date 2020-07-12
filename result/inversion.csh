#!/bin/csh -f
#
# inversion.csh
#
# sample script for running geodetic inversion program and plotting the results.
#
# set *.ctp file (setup file for running inversion)
set ctpname = "c0_inversion.ctp"

echo "STARTING INVERSION"
echo " "

# Run the geodetic inversion program
../exe/gdinv.exe ${ctpname}
echo "===== inversion end ====="
echo " "

# Run the script that plots the inversion analysis results
sh ./plot_results.sh
echo "===== plot result end ====="
echo " "

echo "END OF INVERSION"
