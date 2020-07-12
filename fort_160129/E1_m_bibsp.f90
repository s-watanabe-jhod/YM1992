!=================================================
! SUBROUTINES
!	BIBSP
!
!=================================================

module mod_bibsp

	implicit none

contains
!==========
subroutine BIBSP(nd2, nbdeg, ndd2, nbd2, lc_bi)
	
	use prm_matrix
	use mod_bspline
	
	integer, intent(in)  :: nd2, nbdeg, ndd2, nbd2
	real(8), intent(out) :: lc_bi(Ndd,16) ! bi bspline
	
	integer :: norder, m1, m0, lx, ly, mx, my, j
	real(8) :: bx(4), bs(Nd,4)
	real(8) :: px
	
	norder = nbdeg + 1
	
	do m1 = 1, nd2
	
		px = (dble(m1) - 0.5d0) / dble(nd2)
		CALL BSPLINE(px,bx,nbdeg)
		
		do j = 1, norder
			bs(m1,j) = bx(j)
		end do
	
	end do
	
	do m0 = 1, nbd2
		
		lx = m0 - int((m0-1) / norder) * norder
		ly = 1  + int((m0-1) / norder)
		
		do m1 = 1, ndd2
			mx = m1 - int((m1-1) / nd2) * nd2
			my = 1  + int((m1-1) / nd2)
			
			lc_bi(m1,m0) = bs(mx,lx) * bs(my,ly)
		end do
	
	end do

return
end subroutine BIBSP

!==========

end module mod_bibsp

