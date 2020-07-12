!*************************************************
! module mod_sabicm
!	SABICM(icm, isel)
!
!*************************************************

module mod_sabicm

	implicit none

contains
!=========
subroutine SABICM(icm, isel)

	use prm_inv
	use prm_var
	use mod_subabic

	implicit none

	integer, intent(in) :: icm
	integer, intent(inout) :: isel

	real(8) :: detc, detf, tn, abic(20), alphalog(20), step, abicmin, offset
	integer :: i, j, imin

!	open(88,file="alpha.opt")
!	read (88,*) ch
	open(89,file="alpha.dat")
	write(89,*) "log(alpha) alpha abic rsq2 sigma" ! detc jct detf jmt jdata"

	step = 0.5d0
	offset = - 4.d0
	abicmin = 1000.d0
	imin = 0

	do j = 1, 2

		if(setalpha < 90.d0) then
			imin = 1
			alphalog(imin) = setalpha
			exit
		end if

		do i = 1, 20
!			read(88,*) isel, alphalog
			alphalog(i) = dble(i-10)*step + offset
			alpha = 10.d0 ** alphalog(i)

!			print *, " ENTER ALPHA "
!			read *, alpha

			CALL RESET(alpha, detc, jdata, jct, jmt)
			CALL FREEMIN(detf, rsq2, jdata, jct, jmt)

			if(icm /= 0) CALL INEQMIN(detf)

			tn = dble(jdata + jct - jmt)

			abic(i) = tn + tn * dlog(2.d0*pi*rsq2/tn) - 2.d0 * (dble(jct) * dlog(alpha) + detc - detf)
			sigma = dsqrt(rsq2/tn)

!			print 600, alphalog(i), abic(i), rsq2, sigma, detc, jct, detf, jmt, jdata
			write(89,'(F6.2,4X,F10.5,4X,F10.3,4X,F10.5,4X,F10.5)') alphalog(i), alpha, abic(i), rsq2, sigma
!			print *, ' ENTER ISEL '
!			print *, ' 0/CHANGE ALPHA :  1/COMPUTE RESULT :  2/STOP '
!			read *, isel

!			if(isel == 0) cycle
!			if(isel == 1) exit
!			if(isel /= 0 .and. isel /= 1) stop 'ERROR (in alpha.opt)'

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
	alpha = 10.d0 ** alphalog(1)

	CALL RESET(alpha, detc, jdata, jct, jmt)
	CALL FREEMIN(detf, rsq2, jdata, jct, jmt)

	if(icm /= 0) CALL INEQMIN(detf)

	tn = dble(jdata + jct - jmt)

	abic(1) = tn + tn * dlog(2.d0*pi*rsq2/tn) - 2.d0 * (dble(jct) * dlog(alpha) + detc - detf)
	sigma = dsqrt(rsq2/tn)

	isel = 1

	print 600, alphalog(1), abic(1), rsq2, sigma, detc, jct, detf, jmt, jdata
		  600 FORMAT(/' ','LOG_ALPHA=',F6.3,' ABIC=',F10.3,' RSQ2=',F10.3, &
		     ' SIGMA=',F10.3//' ','DETC= ',E10.4,' JCT= ',I4/' DETF= ',E10.4, &
		     ' JMT= ',I4,12X,' JDATA= ',I4)
	write(89,'(F6.2,4X,F10.5,4X,F10.3,4X,F10.5,4X,F10.5)') alphalog(1), alpha, abic(1), rsq2, sigma

!	close(88)
	close(89)

return
end subroutine SABICM
!=========

end module mod_sabicm

