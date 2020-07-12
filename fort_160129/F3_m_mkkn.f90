module mod_mkkn

	implicit none

contains

subroutine MKKN(ds,gf,cnorm,kfel,nbi,nbdeg,ks00,kt00)

	implicit none

	integer, intent(in) :: kfel, nbdeg, ks00, kt00
	real(8), intent(in) :: gf(kfel,16), cnorm

	integer, intent(out) :: nbi
	real(8), intent(inout) :: ds(kfel)

	integer :: is, it, lc, ix, iy, lls, llt, ld
	real(8) :: da

	nbi = 0

	do it = 1, kt00+nbdeg
		do is = 1, ks00+nbdeg

			da = 0.d0
			nbi = nbi + 1
			lc = 0

			do iy = 0, nbdeg

				llt = it - nbdeg +iy

				do ix = 0, nbdeg

					lls = is - nbdeg + ix
					lc = lc + 1

					if(llt < 1 .or. llt > kt00) cycle
					if(lls < 1 .or. lls > ks00) cycle

					ld = lls + (llt-1) * ks00
					da = da + gf(ld,lc)

				end do
			end do

			ds(nbi) = da * cnorm

		end do
	end do

return
end subroutine MKKN

end module mod_mkkn
