!*************************************************
! module mod_sabicm
!	SABICM
!
!*************************************************

module mod_sabicm

	implicit none

contains
!=========
subroutine SABICM

	use prm_inv
	use mod_subabic

	implicit none

	real(8), parameter :: pi = 3.14159265358979324d0
	real(8), parameter :: base = 10.0d0

	real(8) :: detc, detf, tn, abic(20), alphalog(20), step, abicmin, offset
	integer :: i, j, imin

	open(89,file="alpha.dat")
	write(89,*) "log(alpha) alpha abic rsq2 sigma" ! detc jct detf jmt jdata"

	step    = 0.5d0
	offset  = -4.d0
	abicmin = 100000.d0
	imin = 0

	do j = 1, 2

		if(Setalpha < 90.d0) then
			imin = 1
			alphalog(imin) = Setalpha
			exit
		end if

		do i = 1, 20
			alphalog(i) = dble(i-10)*step + offset
			Alpha = base ** alphalog(i) ! 2019.10.10 ynakamura log is natural, not base 10

			CALL RESET (Alpha, detc, Jdata) ! Jdata: E_prior/=0 の数 detc = Detca
			CALL FREEMIN(detf, Rsq2, Jdata) !

			tn = dble(Jdata + Jcp - Jsp)

			abic(i) = tn + tn * dlog(2.d0*pi*Rsq2/tn) - 2.d0 * (dble(Jcp) * dlog(Alpha) + detc - detf)
			Sigma = dsqrt(Rsq2/tn)

			write(89,'(F6.2,4X,F10.5,4X,F10.3,4X,F10.5,4X,F10.5)') alphalog(i), Alpha, abic(i), Rsq2, Sigma

			if(abic(i) < abicmin) then
				abicmin = abic(i)
				imin = i
			end if

		end do

		write(89,*)

		if(imin == 1 .or. imin == 20) stop 'OUT OF ALPHA REGION: LOG ALPHA = (-9,+1)'
		if(j == 2) exit

		step = 0.1d0
		offset = alphalog(imin)

	end do

	alphalog(1) = alphalog(imin)
	Alpha = base ** alphalog(1) ! 2019.10.10 ynakamura dlog() is a natural log, not base 10

	CALL RESET (Alpha, detc, Jdata)
	CALL FREEMIN(detf, Rsq2, Jdata)

	tn = dble(Jdata + Jcp - Jsp)

	abic(1) = tn + tn * dlog(2.d0*pi*Rsq2/tn) - 2.d0 * (dble(Jcp) * dlog(Alpha) + detc - detf)
	Sigma = dsqrt(Rsq2/tn)

	print 600, alphalog(1), abic(1), Rsq2, Sigma, detc, Jcp, detf, Jsp, Jdata
		  600 FORMAT(/' ','LOG_ALPHA=',F6.3,' ABIC=',F10.3,' RSQ2=',F10.3, &
             ' SIGMA=',F10.3//' ','DETC= ',E10.4,' JCT= ',I4/' DETF= ',E10.4, &
             ' JMT= ',I4,12X,' JDATA= ',I4)
	write(89,'(F6.2,4X,F10.5,4X,F10.3,4X,F10.5,4X,F10.5)') alphalog(1), Alpha, abic(1), Rsq2, Sigma

	close(89)

return
end subroutine SABICM
!=========

end module mod_sabicm
