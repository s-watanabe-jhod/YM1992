!=================================================
! SUBROUTINES
!	OBSPO(iu)
!	BIBSP
!	OUTPAR(flobspd,hd_jac,lsource,iu)
!	JACOBI(iu)
!	JACOBIFW(iu)
!
!
!=================================================

module mod_gdjac

	implicit none

contains
!==========
subroutine OBSPO(iu)

	use prm_matrix
	use prm_var
	use mod_trans

	implicit none

	integer, intent(in) :: iu

	integer :: i, ii
	real(8) :: ala, alg
	character(len=4) :: station

	read(iu,*) isn, ih, iv

	if(ih+iv /= isn) stop "ERROR"
	if(isn > kis) stop 'STATION DIMENSION ERROR'

	do i = 1, isn
		read(iu,*) ii, ala, alg, station
		if(ii /= i) stop ' INPUT DATA ERROR'

		CALL PLTXY(ala,alg,st(1,i),st(2,i),0)
	end do

return
end subroutine OBSPO

!==========
subroutine BIBSP

	use prm_matrix
	use prm_var
	use mod_bspline

	integer :: norder, m1, m0, lx, ly, mx, my, j
	real(8) :: bx(4), bs(nd,4)
	real(8) :: px

	norder = nbdeg + 1

	do m1 = 1, nd2

		px = (dble(m1) - 0.5d0) / dble(nd2)
		CALL BSPLINE(px,bx,nbdeg)

		do j = 1, norder
			bs(m1,j) = bx(j)
		end do

	end do

	do m0 = 1, nbod2

		lx = m0 - int((m0-1) / norder) * norder
		ly = 1  + int((m0-1) / norder)

		do m1 = 1, ndd2
			mx = m1 - int((m1-1) / nd2) * nd2
			my = 1  + int((m1-1) / nd2)

			bi(m1,m0) = bs(mx,lx) * bs(my,ly)
		end do

	end do

return
end subroutine BIBSP

!==========
subroutine OUTPAR(flobspd,hd_jac,lsource,iu)

	use prm_var

	implicit none

	integer, intent(in) :: lsource, iu
	character(len=24), intent(in) :: flobspd, hd_jac

	character(len=24) :: uhd = "**** JACOBI MATRIX **** "


	print '(5X,A24,/" ",A24,/," ","SOURCE NUMBER= ",I5)', uhd, hd_jac, lsource
	print '(1X,3(I5,1X),A24," (ISN,IH,IV,OBS DATA FILE)")', isn, ih, iv, flobspd
	print '(1X,3(I5,1X),2X,"SOURCE (BSPL DEG ,KS ,KT)")', nbdeg, ks0, kt0

	print *, '  MEAN SLIP DIRECTION OUT'

	write(iu) uhd, hd_jac, flobspd

return
end subroutine OUTPAR

!==========
subroutine JACOBI(iu)

	use prm_matrix
	use prm_var
	use mod_jaclib
	use mod_trans

	implicit none

	integer, intent(in) :: iu
	integer :: ios, j
	real(8) :: xs, ys

	print *, ' HORIZONTAL DISPLACEMENTS OBS.POINTS= ', ih

	do ios = 1, ih
		CALL TRANS(st(1,ios),st(2,ios),xs,ys,1,1)
		CALL HORDS(xs,ys)
		CALL DSITRNS
		CALL OUTJC(2,iu)
	end do

	print *, ' VERTICAL DISPLACEMENTS OBS.POINTS= ', iv

	do ios = ih+1, isn
		CALL TRANS(st(1,ios),st(2,ios),xs,ys,1,1)
		CALL VERDS(xs,ys)
		CALL OUTJC(1,iu)
	end do

	print *, ' AVERAGED DIRECTION OF SLIP ON THE FAULT'

	CALL SLIPM
	CALL SMITRNS

	write(iu) (ds1(j,1),ds1(j,2),ds1(j,3), ds2(j,1),ds2(j,2),ds2(j,3), j=1,nbi)

return
end subroutine JACOBI

!==========
subroutine JACOBIFW(iu)

	use prm_matrix
	use prm_var
	use mod_jaclib
	use mod_trans

	implicit none

	integer, intent(in) :: iu
	integer :: ios
	real(8) :: xs, ys

	print *, ' HORIZONTAL DISPLACEMENTS OBS.POINTS= ', ih

	do ios = 1, ih
		CALL TRANS(st(1,ios),st(2,ios),xs,ys,1,1)
		CALL HORDS(xs,ys)
		CALL DSITRNS
		CALL OUTJC(2,iu)
	end do

	print *, ' VERTICAL DISPLACEMENTS OBS.POINTS= ', iv

	do ios = ih+1, isn
		CALL TRANS(st(1,ios),st(2,ios),xs,ys,1,1)
		CALL VERDS(xs,ys)
		CALL OUTJC(1,iu)
	end do

return
end subroutine JACOBIFW

!==========


end module mod_gdjac
