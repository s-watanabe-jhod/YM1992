program GDJACFW
!*************************************************************!
!    COMPUTATION OF JACOBIAN MATRIX FOR THE FORWARD MODELING  !
!                      OF CRUSTAL DEFORMATION                 !
!*************************************************************!

use prm_var

use mod_ctpin
use mod_gdjac

use mod_starea
use mod_subblms

implicit none

integer :: nsource, lsource, ls, ipl, iu
character(len=24) :: fcof,invjac,fwdjac,flobspd,flobspos,hd_jac,fctp="                        "

call getarg(1,fctp)

open(9,file=fctp)

CALL SETA(flobspd,flobspos,nsource) ! ctpin

open(33,file=flobspos,status='old')
	CALL OBSPO(33) ! jaclib
close(33)

do lsource = 1, nsource

	ls = lsource
	print *, ' ---------------------------------------------'
	print *, ' START OF JACOBI MATRIX CALCULATION. SOURCE=', ls

	iu = 40 + lsource

	CALL SETFCP(ls,fcof,invjac,fwdjac,hd_jac,ipl)

	!==== BLMS 断層面の基底関数の設定？ ====

	if(ipl == 0) then ! 断層面が曲面
		open(31,file=fcof,status='old')
			CALL STAREA
			CALL STMSAR
			CALL INP2DP(31)
			CALL MESH
		close(31)
	else if(ipl /= 0) then ! 断層面が平面
		CALL STMSAR2
		CALL MESHSS
	end if

	!==== /BLMS ====

	CALL BIBSP ! BI BSPLINE (DEGREE=NBDEG)

	open(iu,file=fwdjac,form='unformatted')

		CALL OUTPAR(flobspos,hd_jac,ls,iu)
		CALL JACOBIFW(iu)

	close(iu)

end do

close(9)

end


