!=================================================
! SUBROUTINE
!	MNEWT(x1, y1, z1, x2, y2, z2, x3, y3, z3, di, da2)
!		intent(in) :: x1, y1, z1, x2, y2, z2, x3, di, da2
!		intent(out) :: y3, z3
!
!=================================================
module mod_mnewt

	implicit none

contains
!==========
subroutine MNEWT(x1,y1,z1, x2,y2,z2, x3,y3,z3, di,da2)

	use mod_intdep

	implicit none

	real(8), intent(in) :: x1, y1, z1, x2, y2, z2, x3, di, da2
	real(8), intent(out) :: y3, z3

	real(8) :: ct, cx, cy, cz, c2
	real(8) :: dx, dl, dr, dz, du
	real(8) :: arr, arl, ari

	ct = 1.d-6
	cx = x2 - x1
	cy = y2 - y1
	cz = z2 - z1
	c2 = cx*cx + cy*cy + cz*cz

	dx = x3 - x1
	dl = di
	dr = 0.d0
	arr = 0.d0

	do
		y3 = y1 + dl
		CALL INTDEP(x3,y3,z3)

		dz = z3 - z1
		CALL AREA2(cx, cy, cz, dx, dl, dz, c2, arl)

		if(dabs(arl-da2) < ct) return
		if(arl < da2) then
			dr = dl
			arr = arl
			dl = dl * da2 / arl
			cycle
		end if

		do
			du = dr + (dl - dr) * (da2 - arr) / (arl - arr)
			y3 = y1 + du
			CALL INTDEP(x3,y3,z3)

			dz = z3 - z1
			CALL AREA2(cx, cy, cz, dx, du, dz, c2, ari)

			if(dabs(ari-da2) < ct) return

			if(ari > da2) then
				arl = ari
				dl = du
			else if(ari < da2) then
				arr = ari
				dr = du
			end if

		end do

	end do

end subroutine MNEWT

end module mod_mnewt

