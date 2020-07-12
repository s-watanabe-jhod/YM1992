#!/bin/csh

inp="../gpac.dat"
eps="pac_plate.eps"
rng="139/145/34/43"
scl="1:6000000"

cpt="plate.cpt"
grd1="tmp1.grd"
grd2="tmp2.grd"
mask="mask.grd"
nmask="nmask.grd"
mmask="mmask.grd"

gmtset PAGE_ORIENTATION = portrait
gmtset PAPER_MEDIA = a4+

makecpt -Crainbow -T-105/5/10 -Z > ${cpt}
sed -i -e 's/^B.*/B\   255\     0\     0/' ${cpt} \
       -e 's/^F.*/F\   255\     0\   255/' ${cpt} \
       -e 's/^N.*/N\   255\   255\   255/' ${cpt}

blockmean ${inp} -R${rng} -I0.02 -V | surface -G${grd1} -R${rng} -I0.02 -T0.5 -V

grdmask ${inp} -G${mask} -R${rng} -I0.02 -V -S20k #-NNaN/1/1
grdmath ${mask} 0 NAN = ${nmask}
grdmath ${nmask} ${grd1} MUL = ${grd2}

grdimage ${grd2} -R${rng} -Jm${scl} -C${cpt} -P -K -V > ${eps}

psbasemap -Jm${scl} -R${rng} -B2g1 -K -O -V >> ${eps}

psscale -C${cpt} -D4.0/0.5/1.5/0.1h -Bg100a100 -E0.1 -K -O >> ${eps}
pscoast -R${rng} -Jm${scl} -Dh -W2/0/0/0 -O -V >> ${eps}

rm .gmtdefaults
rm .gmtcommands
rm ${grd1}
rm ${grd2}
rm ${mask}
rm ${nmask}
