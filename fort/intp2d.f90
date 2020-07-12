program INTP2D

use mod_ctpin
use mod_sttrans
use mod_starea
use mod_trans
use mod_inp2dp
use mod_intdep

implicit none

integer :: nsource, icheck, ls, lsource, lsd, ipl
real(8) :: ala, alg, x, y, z, x1, y1
character(LEN=24) :: flf, fcof, fctp


call getarg(1,fctp)

open(9,file=fctp)

CALL SETA2D(nsource)

print *, "NSOURCE =", nsource
print *, "ENTER SOURCE NUMBER"
read *, ls

icheck = 0

do lsource = 1, nsource
	lsd = lsource
	CALL SETFC2D(lsd, flf, fcof, ipl)
	if(lsd == ls) then
		icheck = 1
		exit
	end if
end do

if(icheck == 0) stop
if(ipl == 1) stop

CALL STTRANS
CALL STAREA

open(31,file = fcof,status='old')
	CALL INP2DP(31)
close(31)

do
	print *, ' ENTER POSITION(LAT. LON.)'
	print *, ' TO EXIT, ENTER (0.0, 0.0)'
	print *, ' '
	read *, ala, alg

	if(ala == 0.d0 .and. alg == 0.d0) exit

	CALL PLTXY(ala, alg, x, y, 0)
	CALL TRANS(x, y, x1, y1, 1, 1)
	CALL INTDEP(x1, y1, z)

	print '(1X,2(F12.4,1X),1X,3(F10.4,1X))', ala, alg, x, y, z

end do

close(9)

stop
end

