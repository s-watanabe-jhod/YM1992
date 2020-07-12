!=================================================
! SUBROUTINES
!	WHERE1(xx,ix,px,xa,xb,ind,ku0)
!	BSVALI(px,py,bsv,ndeg,nod2)
!	BSPARI(k1,k2,idb,ndeg,kss,ktt)
!
!	no common
!=================================================

module mod_wbb

	implicit none

contains
!=========
subroutine WHERE1(xx,ix,px,xa,xb,ind,ku0)

	implicit none

	integer, intent(in) :: ku0
	real(8), intent(in) :: xx, xa, xb
	integer, intent(out) :: ix, ind
	real(8), intent(out) :: px

	real(8) :: cx, dx, rx

	dx = (xb - xa) / ku0
	ind = 0

	if(xx > xb .or. xx < xa) ind = 1 !OUT OF MODEL REGION

	rx = xx - xa
	ix = int(rx/dx)
	cx = rx - dble(ix)*dx
	px = cx / dx
	ix = ix + 1

return
end subroutine WHERE1

!==========
subroutine BSVALI(px,py,bsv,ndeg,nod2)

	use mod_bspline

	implicit none

	integer, intent(in) :: ndeg, nod2
	real(8), intent(in) :: px, py
	real(8), intent(out) :: bsv(16)

	integer :: ll, m1, m2
	real(8) :: bx(4), by(4)

	CALL BSPLINE(px, bx, ndeg)
	CALL BSPLINE(py, by, ndeg)

	ll = nod2 + 1

	do m2 = 1, ndeg+1
		do m1 = 1, ndeg+1
			ll = ll - 1
			bsv(ll) = bx(m1) * by(m2)
		end do
	end do

return
end subroutine BSVALI

!==========
subroutine BSPARI(k1,k2,idb,ndeg,kss,ktt)

	implicit none

	integer, intent(in) :: k1, k2, ndeg, kss, ktt
	integer, intent(out) :: idb(16)

	integer :: ll, m1, m2

	ll = 0

	do m2 = k2, k2 + ndeg
		do m1 = k1, k1 + ndeg

			ll = ll + 1
			idb(ll) = 0

			if(m2 > ktt .or. m2 < 1) cycle
			if(m1 > kss .or. m1 < 1) cycle

			idb(ll) = m1 + (m2-1) * kss

		end do
	end do

return
end subroutine BSPARI


end module mod_wbb
