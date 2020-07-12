!=================================================
! SUBROUTINES
!	INTDEP(xc, yc, zc)
!		intent(in) :: xc, yc
!		intent(out) :: zc
!	AREA2(cx, cy, cz, dx, dy, dz, c2, ar)
!		intent(in) :: cx, cy, cz, dx, dy, dz, c2
!		intent(out) :: ar
!
!=================================================

module mod_intdep

	implicit none

contains
!==========
subroutine INTDEP(xc, yc, zc)

	use prm_var
	use mod_wbb

	implicit none

	real(8), intent(in) :: xc, yc
	real(8), intent(out) :: zc

	integer :: idb(16)
	real(8) :: bsv(16)
	integer :: j, ix, iy, indx, indy
	real(8) :: px, py

	CALL WHERE1(xc,ix,px,xa,xb,indx,ku0)
	CALL WHERE1(yc,iy,py,ya,yb,indy,kv0)

!	write(*,*)'AAAAAAAAAA',indx,indy,xc,xa,xb,yc,ya,yb
	if(indx == 1 .or. indy == 1) then
		print *, ' WARNING SPECIFIED POINT IS OUT OF REGION'
		stop
	end if

	CALL BSVALI(px, py, bsv, ndeg, nod2)
	CALL BSPARI(ix, iy, idb, ndeg, ku1, kv1)

	zc = 0.d0

	do j = 1, nod2
		zc = zc + bsv(j) * xi(idb(j))
	end do

return

end subroutine INTDEP

!==========
subroutine AREA2(cx, cy, cz, dx, dy, dz, c2, ar)

	implicit none

	real(8), intent(in) :: cx, cy, cz, dx, dy, dz, c2
	real(8), intent(out) :: ar
	real(8) :: d2, dd

	d2 = dx * dx + dy * dy + dz * dz
	dd = dx * cx + dy * cy + dz * cz
	ar = dsqrt(d2 * c2 - dd * dd)

return

end subroutine AREA2

end module mod_intdep

