!=================================================
! SUBROUTINES
!	SETA2D(nsource)
!	SETFC2D(lsource,flf,fcof,ipl)
!	SETA(flobspd,flobspos,nsource)
!	SETFCP(lsource,fcof,invjac,fwdjac,hd_jac,ipl)
!
! ctp file Çì«Ç›çûÇﬁ
!
!=================================================

module mod_ctpin

	implicit none

contains
!==========
subroutine SETA2D(nsource)

	use prm_var

	implicit none

	integer, intent(out) :: nsource

	read(9,*) alat0, alng0, icord
	read(9,*)
	read(9,*)
	read(9,*)
	read(9,*) nsource

	print *, " ORIGIN OF COORDINATE(LAT,LON)=(", alat0, alng0, ")"
	print *, ' SOURCE IS DIVIDED INTO ', nsource,' REGION'

return
end subroutine SETA2D

!==========
subroutine SETFC2D(lsource,flf,fcof,ipl)

	use prm_var

	implicit none

	integer, intent(in) :: lsource
	character(len=24), intent(out) :: flf,fcof
	integer, intent(out) :: ipl

	integer :: lsin

	read(9,*)
	read(9,*) lsin
		if(lsin /= lsource) then
			print *, ' CTR.PR. CHECK ERROR. PROCESS STOP.'
			stop
		end if
	read(9,'(A)') flf
	read(9,'(A)') fcof
	read(9,*) aphi, adlt
	read(9,*) alat, alon
	read(9,*) alen, awid
	read(9,*) ku0, kv0, ndeg
	read(9,*)
	read(9,*)
	read(9,*)
	read(9,*)
	read(9,*)
	read(9,*)
	read(9,*)

	ipl = 0

	if(fcof == "PLANE FAULT             ") ipl = 1

	! == parameter set ==
	nod2 = (ndeg +1) * (ndeg +1)
	ku1 = ku0 + ndeg
	kv1 = kv0 + ndeg

return
end subroutine SETFC2D

!==========
subroutine SETA(flobspd,flobspos,nsource)

	use prm_var

	implicit none

	integer, intent(out) :: nsource
	character(len=24), intent(out) :: flobspd,flobspos

	read(9,*) alat0, alng0, icord
	read(9,'(A)') flobspd
	read(9,'(A)') flobspos
	read(9,*) delt
	read(9,*) nsource

	fde = 1.d0/(1.d0+delt)

	print *, ' ORIGIN OF COODINATE(LAT,LON) =( ', alat0, alng0,')'
	print *, '      GEODETIC DATA FILE NAME =', flobspd
	print *, '     ELASTIC CONST(LAMBDA/MU) =', delt

	print *, '  SOURCE TYPE:  /STR SLIP,DIP SLIP,OPEN C,EXP/)'
	print *, '  LIST ORDER:(STATION / DATA TYPE / SOURCE TYPE / BASE FUNCTION)'

return
end subroutine SETA

!==========
subroutine SETFCP(lsource,fcof,invjac,fwdjac,hd_jac,ipl)

	use prm_var

	implicit none

	integer, intent(in) :: lsource
	integer, intent(out) :: ipl

	character(len=24), intent(out) :: fcof, hd_jac, invjac, fwdjac

	integer :: lsin

	read(9,*) 
	read(9,*) lsin
		if(lsin /= lsource) then
			print *, ' CTR.PR. CHECK ERROR. PROCESS STOP.'
			stop
		end if
	read(9,*) 
	read(9,'(A)') fcof
	read(9,*) aphi, adlt
	read(9,*) alat, alon
	read(9,*) alen, awid
	read(9,*) ku0, kv0, ndeg
	read(9,*) 
	read(9,'(A)') hd_jac
	read(9,'(A)') invjac
	read(9,'(A)') fwdjac
	read(9,*) alaf, alof
	read(9,*) alef, awif, adep
	read(9,*) ks0, kt0, nd2, nbdeg

	ipl = 0
	if(fcof == "PLANE FAULT             ") ipl = 1

!      +++ PARAMETER SET +++

	nod2 = (ndeg+1) * (ndeg+1)
	ku1 = ku0 + ndeg
	kv1 = kv0 + ndeg

	nbod2 = (nbdeg+1) * (nbdeg+1)
	ndd2 = nd2 * nd2
	ms2 = nd2 * ks0
	mt2 = nd2 * kt0

return
end subroutine SETFCP

!===========



end module

