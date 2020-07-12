#! /usr/bin/sh

# PLOT RESULT => 2 FILES
#  cnt_*.eps : source displacement contour file
#  scd_*.eps : source displacement vector  file
#

pltdt="../plate_boundary/Hirose_gsi/PAC/gphs2.dat"

rng="130/140/30/37"
scl="1:10000000"

gmt makecpt -Cpolar -Ic -T-10/10/0.1 -D > col.cpt
cpt="col.cpt"

cntmax="0.1"
dcnt="1"

# data scale
fwdscl="0.5"
scdscl="1."

length="10cm"

#-------------------#

# FILENAMES
fname="res"

icnt="./data/"${fname}"_scd.cnt"
iscd="./data/"${fname}"_scd.xyz"
ielp="./data/"${fname}"_scd.elp"

ical="./data/"${fname}"_cal.dat"
iobs="./data/"${fname}"_obs.dat"
io_c="./data/"${fname}"_o_c.dat"

ocnt="./cnt_"${fname}".ps"
oscd="./scd_"${fname}".ps"
ofwd="./fwd_"${fname}".ps"
oo_c="./o_c_"${fname}".ps"

pcnt="./cnt_"${fname}".png"
pscd="./scd_"${fname}".png"
pfwd="./fwd_"${fname}".png"
po_c="./o_c_"${fname}".png"

# FILENAMES

mask="../tmp/mask.grd"
nmask="../tmp/nmask.grd"
mmask="../tmp/mmask.grd"
lmask="../tmp/lmask.grd"

trench="../plate_boundary/trench/trench.gmt"
trans="../plate_boundary/trench/transform.gmt"
ridge="../plate_boundary/trench/ridge.gmt"

temp1="../tmp/tmp.grd"
tmpscd="../tmp/tmpscd.xyz"
tmpobs="../tmp/tmpobs.xyz"
tmpcal="../tmp/tmpcal.xyz"
tmpo_c="../tmp/tmpo_c.xyz"
pltmp="../tmp/pltbnd.grd"

qcnt="./"${fname}"_cnt.png"
qscd="./"${fname}"_scd.png"
qfwd="./"${fname}"_fwd.png"
qo_c="./"${fname}"_o_c.png"

gmt gmtset PAGE_ORIENTATION = portrait
gmt gmtset PAPER_MEDIA = a4

#===== CNT =====

gmt triangulate ${pltdt} -G${pltmp} -I0.1 -R${rng} > ../tmp/tri.tmp

gmt blockmean ${icnt} -R${rng} -I0.02 -V > ../tmp/blk.xyz
gmt surface ../tmp/blk.xyz -G${temp1} -R${rng} -I0.02 -T0.5 -V

gmt grdmask ${icnt} -G${mask} -R -I0.02 -V -S20k #-NNaN/1/1
gmt grdmath ${mask} 0 NAN = ${nmask}
gmt grdmath ${nmask} ${temp1} MUL = ${mmask}
gmt grdmath ${mmask} 10 MUL = tempor.grd

gmt psbasemap -Jm${scl} -R${rng} -Bg1f2a2WeSn -K -V > ${ocnt}

gmt grdimage ${mmask} -R${rng} -Jm${scl} -C${cpt} -Q -P -K -O -V >> ${ocnt}
gmt grdimage tempor.grd -R${rng} -Jm${scl} -C${cpt} -Q -P -K -O -V >> ${ocnt}


gmt psscale -C${cpt} -D8.5/2.5/4.0/0.2h -B4:slip -E0.3 -K -O >> ${ocnt}

gmt pstext -R -Jm${scl} -P -K -O -N <<EOF >> ${ocnt}
 143.6 35.0 16 0 0 TC SLIP(m)
EOF

# draw contours of the plate boundary depth (10 km interval)
gmt grdcontour ${pltmp} -R${rng} -Jm${scl} -C10 \
-W0.2,black,. -L-60/-0.1 -K -O -V >> ${ocnt} 2>/dev/null

# draw boundaries
gmt psxy ${trench} -R${rng} -Jm${scl} -W0.5,black  -K -O -V >> ${ocnt} \
  2>/dev/null
gmt psxy ${trans}  -R${rng} -Jm${scl} -W0.5,black -K -O -V >> ${ocnt} \
  2>/dev/null
gmt psxy ${ridge}  -R${rng} -Jm${scl} -W0.5,black -K -O -V >> ${ocnt} \
  2>/dev/null


gmt grdcontour tempor.grd -R${rng} -Jm${scl} -C1 \
-W0.5 -A${dcnt}t+f6 -L-10/10 -K -O -V >> ${ocnt}

gmt pscoast -R -Jm${scl} -Dh -W0.5,black -K -O -V >> ${ocnt}

gmt psbasemap -Jm${scl} -R${rng} -Bg1f2a2WeSn -O -V >> ${ocnt}

#===== SCD =====

gmt psbasemap -Jm${scl} -R${rng} -Bg1f2a2WeSn -K -V > ${oscd}

gmt grdcontour ${pltmp} -R${rng} -Jm${scl} -C10 \
-W0.2,black,. -L-60/-0.1 -K -O -V >> ${oscd} 2>/dev/null
gmt psxy ${trench} -R${rng} -Jm${scl} -W0.5,black -K -O -V >> ${oscd} \
  2>/dev/null
gmt psxy ${trans}  -R${rng} -Jm${scl} -W0.5,black -K -O -V >> ${oscd} \
  2>/dev/null
gmt psxy ${ridge}  -R${rng} -Jm${scl} -W0.5,black -K -O -V >> ${oscd} \
  2>/dev/null

gfortran ../fort/scalexyz.f90 -o ../exe/cgscd.exe
../exe/cgscd.exe ${iscd} ${tmpscd} none none ${scdscl}

gmt psxy ${tmpscd} -Jm${scl} -R${rng} \
-Sv0.007/0.08/0.08n0.2 -W1,red -Gred -P -K -O -V >> ${oscd} \
  2>/dev/null

gmt pscoast -R -Jm${scl} -Dh -W0.5,black -K -O -V >> ${oscd}

gmt psbasemap -Jm${scl} -R${rng} -Bg1f2a2WeSn -O -V >> ${oscd}

#===== FWD =====

../exe/cgscd.exe ${iobs} ${tmpobs} none none ${fwdscl}
../exe/cgscd.exe ${ical} ${tmpcal} none none ${fwdscl}

gmt psbasemap -Jm${scl} -R${rng} -Bg1f2a2WeSn -K -V > ${ofwd}
gmt pscoast -R -Jm${scl} -Dh -W0.5,black -O -V -K >> ${ofwd}

gmt grdcontour ${pltmp} -R${rng} -Jm${scl} -C10 \
-W0.2,black,. -L-60/-0.1 -K -O -V >> ${ofwd} 2>/dev/null
gmt psxy ${trench} -R${rng} -Jm${scl} -M -W1,black  -K -O -V >> ${ofwd} \
  2>/dev/null
gmt psxy ${trans}  -R${rng} -Jm${scl} -M -W1,black   -K -O -V >> ${ofwd} \
  2>/dev/null
gmt psxy ${ridge}  -R${rng} -Jm${scl} -M -W1,black  -K -O -V >> ${ofwd} \
  2>/dev/null

gmt psxy ${tmpobs} -Jm${scl} -R${rng} \
-Sv0.007/0.08/0.08n0.2 -W0.5,black -Gblack -P -O -V -K >> ${ofwd} \
  2>/dev/null
gmt psxy ${tmpcal} -Jm${scl} -R${rng} \
-Sv0.007/0.08/0.08n0.2 -W0.5,red -Gred -P -O -V -K >> ${ofwd} \
  2>/dev/null


gmt psxy -Jm${scl} -R${rng} \
-Sv0.1/0.24/0.2n0.2 -W0/0/0/0 -G0/0/0 -P -O -V -K <<EOF >> ${ofwd}
 143.2 34.8 0.0 1.0
EOF
gmt pstext -R -Jm${scl} -P -K -O -N <<EOF >> ${ofwd}
 143.9 34.8 12 0 0 ML Obs.
EOF

gmt psxy -Jm${scl} -R${rng} \
-Sv0.05/0.24/0.1n0.2 -W0/0/0/0 -G255/255/0 -P -O -V -K <<EOF >> ${ofwd}
 143.2 34.5 0.0 1.0
EOF
gmt pstext -R -Jm${scl} -P -K -O -N <<EOF >> ${ofwd}
 143.9 34.5 12 0 0 ML Cal.
EOF

gmt pstext -R -Jm${scl} -P -K -O -N <<EOF >> ${ofwd}
 143.5 34.2 12 0 0 MC ${length}
EOF

gmt psbasemap -Jm${scl} -R${rng} -Bg1f2a2WeSn -O -V >> ${ofwd}

#===== O-C =====

../exe/cgscd.exe ${io_c} ${tmpo_c} none none ${fwdscl}

gmt psbasemap -Jm${scl} -R${rng} -Bg1f2a2WeSn -K -V > ${oo_c}
gmt pscoast -R -Jm${scl} -Dh -W0.5,black -O -V -K >> ${oo_c}

gmt grdcontour ${pltmp} -R${rng} -Jm${scl} -C10\
-W0.2,black,. -L-60/-0.1 -K -O -V >> ${oo_c} 2>/dev/null
gmt psxy ${trench} -R${rng} -Jm${scl} -W0.5,black -K -O -V >> ${oo_c} \
  2>/dev/null
gmt psxy ${trans}  -R${rng} -Jm${scl} -W0.5,black -K -O -V >> ${oo_c} \
  2>/dev/null
gmt psxy ${ridge}  -R${rng} -Jm${scl} -W0.5,black -K -O -V >> ${oo_c} \
  2>/dev/null

gmt psxy ${tmpo_c} -Jm${scl} -R${rng} \
-Sv0.01/0.07/0.05n0.2 -W0.5,black -Gblack -P -O -V >> ${oo_c} \
  2>/dev/null

#===== DEL TMP =====

gmt ps2raster ./*${fname}.ps -Tg -A -E120 -Qt -Qg
gmt convert -transparent "#ffffff" ${pcnt} ${qcnt}
gmt convert -transparent "#ffffff" ${pscd} ${qscd}
gmt convert -transparent "#ffffff" ${pfwd} ${qfwd}
gmt convert -transparent "#ffffff" ${po_c} ${qo_c}

rm ${pcnt}
rm ${pscd}
rm ${pfwd}
rm ${po_c}

rm -f gmt.conf
rm -f gmt.history

rm ${mask}
rm ${nmask}
rm ${mmask}
rm ${temp1}
rm ${tmpscd}
rm ${tmpobs}
rm ${tmpcal}
rm ${tmpo_c}
rm ${pltmp}
rm tempor.grd
rm -f ../tmp/tri.tmp
