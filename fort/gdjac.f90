program gdjac

use prm_var

use mod_ctpin
use mod_gdjac

use mod_starea
use mod_inp2dp
use mod_subblms
use mod_sttrans

implicit none

integer :: nsource, lsource, ls, ipl, iu
character(len=24) :: fcof,invjac,fwdjac,flobspd,flobspos,hd_jac,fctp="                        "

call getarg(1,fctp)

open(9,file=fctp)

CALL SETA(flobspd,flobspos,nsource) ! ctpin

open(33,file=flobspd,status='old')
	CALL OBSPO(33) ! jaclib
close(33)


do lsource = 1, nsource

	ls = lsource
	print *, ' ---------------------------------------------'
	print *, ' START OF JACOBI MATRIX CALCULATION. SOURCE=', ls

	iu = 40 + lsource

	CALL SETFCP(ls,fcof,invjac,fwdjac,hd_jac,ipl)

	CALL STTRANS


	!==== BLMS ífëwñ ÇÃäÓíÍä÷êîÇÃê›íËÅH ====

	if(ipl == 0) then ! ífëwñ Ç™ã»ñ 
		open(31,file=fcof,status='old')
			CALL STAREA
			CALL STMSAR
			CALL INP2DP(31)
			CALL MESH
		close(31)
	else if(ipl /= 0) then ! ífëwñ Ç™ïΩñ 
		CALL STMSAR2
		CALL MESHSS
	end if

	!==== /BLMS ====

	CALL BIBSP ! BI BSPLINE (DEGREE=NBDEG)

	open(iu,file=invjac,form='unformatted')

		CALL OUTPAR(flobspd,hd_jac,ls,iu)
		CALL JACOBI(iu)

	close(iu)

end do

close(9)

end

