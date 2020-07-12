program GDINV

!-------------------------------------------------------
!    GEODETIC DATA INVERSION WITH abic MINIMUM CRITERION
!
!    @@@ ATTENTION @@@
!          YOU MUST ALREADY PROCESSED PROGRAM GDJAC
!-------------------------------------------------------

use prm_matrix
use prm_inv
use prm_var

use mod_ctpinv
use mod_setinv
use mod_sabicm
use mod_inv
use mod_dissou

implicit none

integer :: nsource, icm, isel, i, j, lname
character(len=24) :: flobspd, cal, fmp, hd, fctp = "                        "
	character(len=24) :: fcalp = "                        "
	character(len=24) :: fobsp = "                        "
	character(len=24) :: fo_cp = "                        "

open(77,file='../tmp/gdinv.out')

call getarg(1,fctp)

open(9,file=fctp)
	CALL SETCP(hd,flobspd,cal,fmp,nsource)
close(9)

open(33,file=flobspd,status='old')
	CALL INPUTD(33)
close(33)

CALL SETUP(flobspd,nsource,icm,0)
CALL SETWT

if(icm /= 0) CALL CHGPM

CALL SABICM(icm,isel)

if(isel /= 2) then

! === CALL DRESID ===
	ysl(1:itotal) = 0.d0
	do j = 1, jsp
		do i = 1, itotal
			ysl(i) = ysl(i) + zz(i,j) * x(j)
		end do
	end do

	do i = 1, itotal
		ysl(i) = ysl(i) * dabs(ee(i))
		yy(i)  = yy(i)  * dabs(ee(i))
	end do
! === /DRESID ===

	lname = len_trim(cal)
	fcalp(1:lname-4) = cal(1:lname-4)
	fobsp(1:lname-4) = cal(1:lname-4)
	fo_cp(1:lname-4) = cal(1:lname-4)
	fcalp(lname-3:lname+4) = "_cal.dat"
	fobsp(lname-3:lname+4) = "_obs.dat"
	fo_cp(lname-3:lname+4) = "_o_c.dat"

	open(60,file=fcalp)
	open(61,file=fobsp)
	open(62,file=fo_cp)

	open(21,file=cal,status='unknown')
		CALL OUTDAT(hd)
	close(21)

	close(60)
	close(61)
	close(62)

	CALL COVAR

	if(icm /= 0) CALL CHGIMP

	open(22,file=fmp,status='unknown')
		CALL OUTSOL(hd,flobspd,nsource)
	close(22)

	CALL DISSOU(hd,nsource)

end if

close(77)

end


