#!/bin/sh

#-----------------------------------------------
#       contour_map_PAC.gmt
#             2010/07/22 by fuyuki_H
#-----------------------------------------------

psfile=./contour_map_PAC.ps

trench=./mapdata/trench.dat
vol=./mapdata/voldata.dat
pref=./mapdata/japan.pref.line
#fault=./mapdata/activef.dat

tokai=./mapdata/tokai.region
tonankai=./mapdata/tonankai.region
nankai=./mapdata/nankai.region
tasp=./mapdata/tokai_asperity.data
kasp=./mapdata/kanto_eq.dat

# PAC_Hokkaido (Kita et al., 2010, EPSL)
# PAC_Tohoku (Nakajima and Hasegawa, 2006, GRL)
# PAC_Kanto (Nakajima et al., 2009, JGR)
PACdata=./PAC/plate_combine.dat
grd=pl_plate.grd

#-------------------------------------------
# map view
#-------------------------------------------

proje=-JM17
range=-R139/145/34/43
#range=-R130/147/30/47
tick=-BWeNsa2f1

psbasemap $proje $range $tick -G255/255/255 -P -X2 -Y5 -K > $psfile 
pscoast $proje $range -G230/230/230 -Dh -K -O >> $psfile

#
pscoast $proje $range -W3 -Dh -K -O >> $psfile
psxy $proje $range $trench -: -W4t20_10:0 -M -K -O >> $psfile

# CAMP model
obspos="../jcg.xy"
pltcmp="../camp_pac/gpac.dat"

tmpcmp="./cmpbnd.grd"

xyz2grd ${pltcmp} -G${tmpcmp} -I0.1 -R${rng}

psxy ${obspos} $proje $range -SC0.3 -W1/0/0/0 -G255/255/0 -K -O >> ${psfile} \
  2>/dev/null

grdcontour ${tmpcmp} $proje $range -W3ta -A10ta90  -G16c -L-60/0  -K -O -V >> ${psfile} \
  2>/dev/null
grdcontour ${tmpcmp} $proje $range -W5ta -A0.2ta90 -G16c -L-0.3/0 -K -O -V >> ${psfile} \
  2>/dev/null

rm ${tmpcmp}
# /CAMP model

# PAC
#awk '$3>=100 && $3<=410 {print $1,$2,$3}' $PACdata |\
#blockmedian $range -I0.2 |\
#triangulate -G$grd $range -I0.1

#grdcontour $grd $proje $range -C50 -A50 -W5t20_10:0/255/127/0 -M -K -O >> $psfile

awk '$3>=-10 && $3<=100 {print $1,$2,$3}' $PACdata |\
blockmedian $range -I0.2 |\
triangulate -G$grd $range -I0.1

grdcontour $grd $proje $range -C10 -A10 -W5t20_10:0/255/127/0 -M -K -O >> $psfile

#
# Arrow (Wei & Seno, 1998, Geodynam. Series ed. by M. Flower et al., 27, 337-346)
#psxy $proje $range -Sv0.1/0.3/0.2 -G0 -O -K <<EOF >> $psfile
#140.8 33.5 131 0.6375
#EOF
#echo 140.8 33.3 9 0 0 6 "34 mm/yr" | \
#pstext $proje $range -G0/0/0 -O -K >> $psfile
#
#psxy $proje $range -Sv0.1/0.3/0.2 -G0 -O -K <<EOF >> $psfile
#137 31.1 145 1.05
#EOF
#echo 137.05 31.4 9 0 0 5 "56 mm/yr" | \
#pstext $proje $range -G0/0/0 -O -K >> $psfile
#
#psxy $proje $range -Sv0.1/0.3/0.2 -G0 -N -O -K <<EOF >> $psfile
#145 34.8 156 1.5
#EOF
#echo 144.7 35.2 9 0 0 5 "80 mm/yr" | \
#pstext $proje $range -N -G0/0/0 -O -K >> $psfile
#
#psxy $proje $range -Sv0.1/0.3/0.2 -G0 -N -O -K <<EOF >> $psfile
#146.5 39 155 1.5375
#EOF
#echo 145.9 38.8 9 0 0 6 "82 mm/yr" | \
#pstext $proje $range -N -G0/0/0 -O -K >> $psfile
#
gmtset BASEMAP_TYPE plain
gmtset FRAME_WIDTH 0.2c
gmtset TICK_LENGTH 0.1c

