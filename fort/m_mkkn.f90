module mod_mkkn

	implicit none

contains

subroutine MKKN(ds,gf,cnorm,kfe,nbi,nbdeg,ks0,kt0)

	implicit none

	integer, intent(in) :: kfe, nbdeg, ks0, kt0
	real(8), intent(in) :: gf(kfe,16), cnorm

	integer, intent(out) :: nbi
	real(8), intent(inout) :: ds(kfe)

	integer :: is, it, lc, ix, iy, lls, llt, ld
	real(8) :: da

	nbi = 0

	do it = 1, kt0+nbdeg
		do is = 1, ks0+nbdeg

			da = 0.d0
			nbi = nbi + 1
			lc = 0

			do iy = 0, nbdeg

				llt = it - nbdeg +iy

				do ix = 0, nbdeg

					lls = is - nbdeg + ix
					lc = lc + 1

					if(llt < 1 .or. llt > kt0) cycle
					if(lls < 1 .or. lls > ks0) cycle

					ld = lls + (llt-1) * ks0
					da = da + gf(ld,lc)

				end do
			end do

			ds(nbi) = da * cnorm

		end do
	end do

return
end subroutine MKKN

end module mod_mkkn
