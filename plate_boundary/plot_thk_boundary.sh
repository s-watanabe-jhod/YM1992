#! /usr/bin/sh

obspos="./jcg.xy"
pltcmp="./camp_pac/gpac.dat"
pltgsi="./Hirose_gsi/PAC/gsi_pac.dat"

rng="139/145/34/43"
scl="1:6000000"
#-------------------#

# FILENAMES
opac="./PAC.eps"
# FILENAMES

tmpcmp="./cmpbnd.grd"
tmpgsi="./gsibnd.grd"

gmtset PAGE_ORIENTATION = portrait
gmtset PAPER_MEDIA = a4+

#===== CNT =====

xyz2grd ${pltcmp} -G${tmpcmp} -I0.1 -R${rng}

psbasemap -Jm${scl} -R${rng} -B2g1 -K -V > ${opac}

psxy ${obspos} -R -Jm${scl} -SC0.1 -W1/0/0/0 -G255/255/0 -K -O >> ${opac} \
  2>/dev/null

grdcontour ${tmpcmp} -R${rng} -Jm${scl} -W1ta -A10ta90  -G16c -L-60/0  -K -O -V >> ${opac} \
  2>/dev/null
grdcontour ${tmpcmp} -R${rng} -Jm${scl} -W2ta -A0.2ta90 -G16c -L-0.3/0 -K -O -V >> ${opac} \
  2>/dev/null

triangulate ${pltgsi} -G${tmpgsi} -I0.1 -R${rng} \
  1>/dev/null  2>/dev/null
grdcontour ${tmpgsi} -R${rng} -Jm${scl} -W255/0/0 -A10ta90  -G18c -L-60/-0.1  -K -O -V >> ${opac} \
  2>/dev/null
grdcontour ${tmpgsi} -R${rng} -Jm${scl} -W255/0/0 -A0.2ta90 -G18c -L-0.1/0.1 -K -O -V >> ${opac} \
  2>/dev/null

pscoast -R -Jm${scl} -Dh -W2/0/0/0 -K -O -V >> ${opac}

psbasemap -Jm${scl} -R${rng} -B2g1 -O -V >> ${opac}


#===== DEL TMP =====

rm .gmtdefaults*
rm .gmtcommands*
rm ${tmpcmp}
rm ${tmpgsi}

