!*************************************************
! module
!	INSURD(iu)
!	SETKN
!	SETBM2D
!	SABIC2D
!	OUT2DP(flf,iu)
!
!*************************************************

module mod_inv2d

	implicit none

contains
!=========
subroutine INSURD(iu)

	use prm_matrix
	use prm_var
	use mod_trans

	implicit none

	integer, intent(in) :: iu
	integer :: i
	real(8) :: ala, alg, zz, xx, yy, x1, y1

	read(iu,*)
	read(iu,*) isurf

	if(isurf > ki) then
		print *, ' ++WARNING++ INSUFFICIENT ARRAY DIMENSION ', isurf
		print *, ' ONLY ', ki, ' DATA IS USED.'
	end if

	do i = 1, isurf

		read(iu,*,end=999) ala, alg, zz

		CALL PLTXY(ala,alg,xx,yy,0)
		CALL TRANS(xx,yy,x1,y1,1,1)

		gx(i) = x1
		gy(i) = y1
		gz(i) = zz

	end do

	print *, ' DATA INPUT COMPLETE. ISURF= ', isurf

return

999	print *, ' INPUT ERROR SURFACE GEOMETRY OBSERVED DATA'
	stop

end subroutine INSURD

!==========
subroutine SETKN

	use prm_matrix
	use prm_var
	use mod_wbb

	implicit none

	integer :: i, j, ii, ix, iy, indx, indy
	real(8) :: xx, yy, px, py

	real(8) :: bsv(16)
	integer :: idb(16)

	rt(1:ki,1:kmr) = 0.d0

	ii = 0

	do i = 1, isurf
		xx = gx(i)
		yy = gy(i)

		CALL WHERE1(xx,ix,px,xa,xb,indx,ku0)
		CALL WHERE1(yy,iy,py,ya,yb,indy,kv0)
			! ix グリッドを切った時（widthをk0+1個のnode）のnodeの番号
			! px そのnode内からプラス方向にどれだけの位置にあるか（0<px<1）
			! 要は、グリッド空間での位置を算出

		if (indx == 1 .or. indy == 1) then
			print *, ' DATA I= ', i, ' OUT OF MODEL REGION'
			cycle
		end if

		ii = ii + 1

		gz(ii) = gz(i)
		gx(ii) = xx
		gy(ii) = yy

		CALL BSVALI(px,py,bsv,ndeg,nod2)
		CALL BSPARI(ix,iy,idb,ndeg,ku1,kv1)

		do j = 1, nod2
			rt(ii,idb(j)) = bsv(j)
		end do

	end do

	itotal = ii
	jmp = ku1 * kv1

	print *, ' KARNEL MATRIX IS CALCULATED. ITOTAL= ', itotal

return
end subroutine SETKN

!==========
subroutine SETBM2D

	use prm_matrix
	use prm_var

	implicit none

	integer :: i, j, ii, ip, l

	cs2d(1:kmr,1:kmr) = 0.d0

	ii = 0

	do i = 1, ku1
		do j = 1, kv1

			ip = 0

			if(i /= 1 .and. i /= ku1) ip = ip + 1
			if(j /= 1 .and. j /= kv1) ip = ip + 2
			if(ip == 0) cycle

			ii = ii + 1
			l = (j-1) * ku1 + i

			if(ip /= 1) then
				cs2d(ii,l) = cs2d(ii,l) - 2.d0
				cs2d(ii,l-ku1) = 1.d0
				cs2d(ii,l+ku1) = 1.d0
			end if

			if(ip == 2) cycle

			cs2d(ii,l) = cs2d(ii,l) - 2.d0
			cs2d(ii,l-1) = 1.d0
			cs2d(ii,l+1) = 1.d0

		end do
	end do

	jcp = ii
	detca = 0.d0

	print *, ' CONSTRAIN MATRIX IS SET. JCP= ', jcp

return
end subroutine SETBM2D

!==========
subroutine SABIC2D

use prm_var
use mod_rf2d

implicit none

	real(8) :: dlx, crt, abica(7), alph(7), detc, detf, tn, abic
	real(8) :: a1, a2, ew, ek, cp
	integer :: ip, ipcent, istart, ipq, ipr, ipp, i_counter

	ipq = 1
	ipr = 7

	dlx = dlog(10.d0)
	crt = 0.4d0
	ipcent = 4

	do ip = 1, 7
		cp = dble(3-ip)
		alph(ip) = dexp(cp*dlx)
	end do

	istart = 0

loop1: &
	do i_counter = 1, 10000

		if(istart == 0) then
			ipq = 1
			ipr = 7
		else if(istart > 0) then
			ipq = 2
			ipr = 6
		end if

		do ip = ipq, ipr

			if(istart > 0 .and. ip == ipcent) cycle

			alpha = alph(ip)

			CALL RESET2D(alpha, detc, jdata, jct, jmt)
			CALL FREE2D(detf, rsq2, jdata, jct, jmt)

			tn = dble(jdata + jct - jmt)
			abic = tn + tn * dlog(2.d0*pi * rsq2/tn) - 2.d0 * ( dble(jct)*dlog(alpha) + detc - detf)
			sigma = dsqrt(rsq2/tn)

			abica(ip) = abic

			print 610, alpha, abic, rsq2, sigma
			write(77,610) alpha, abic, rsq2, sigma
		  610 FORMAT(//' ','ALPHA=',F10.5,' ABIC=',F10.3,' RSQ2=',F10.3,' SIGMA=',F10.3)

		end do

		do ip = 2, 6

			if(abica(ip-1) > abica(ip) .and. abica(ip) < abica(ip+1)) then

				a1 = alph(ip-1)
				a2 = alph(ip+1)
				ew = (dlog(a2) - dlog(a1)) / dlx

				if(dabs(ew) < crt) then
					alpha = alph(ip)

					CALL RESET2D(alpha, detc, jdata, jct, jmt)
					CALL FREE2D(detf, rsq2, jdata, jct, jmt)

					tn = dble(jdata + jct - jmt)
					abic = tn + tn * dlog(2.d0*pi * rsq2/tn) - 2.d0 * ( dble(jct)*dlog(alpha) + detc - detf)
					sigma = dsqrt(rsq2/tn)

					print *, '------------------------------------------------------'
					print *, ' REPORT : END OF ABIC MINIMUM SEARCH-> COMPUTE RESULT'
					print *, '------------------------------------------------------'

					print 600, alpha, abic, rsq2, sigma, detc, jct, detf, jmt, jdata
					write(77,600) alpha, abic, rsq2, sigma, detc, jct, detf, jmt, jdata
				  600 FORMAT(//' ','ALPHA=',F10.5,' ABIC=',F10.3,' RSQ2=',F10.3,' SIGMA=',F10.3// &
								' ','DETC= ',E10.4,' JCT= ',I4/' DETF= ',E10.4,' JMT= ',I4,12X,' JDATA= ',I4)

					return

				end if

				do ipp = 1, 7
					ek = ew * dlx * dble(ipp-1) / 6.d0
					alph(ipp) = a1 * dexp(ek)
				end do

				abica(1) = abica(ip-1)
				abica(ipcent) = abica(ip)
				abica(7) = abica(ip+1)

				istart = 1

				cycle loop1

			end if

		end do

		print *, 'FATAL ERROR: ABIC MINIMUM SEARCH : NO MINIMUM'
		stop

	end do loop1

end subroutine SABIC2D

!==========
subroutine OUT2DP(flf,iu)

	use prm_matrix
	use prm_var

	integer, intent(in) :: iu
	character(len=24), intent(in) :: flf

	integer :: i

	write(iu,'(A)') '*** 2-DIMENSIONAL FITTING ***'
	write(iu,'(2A)') flf,'/INPUT FILE NAME'

	write(iu,'(1X,I5,1X,I5,1X,I5,1X,I5,1X,I5)') jmt, ku1, kv1, ndeg, jdata
	write(iu,'(A)') 'MODEL PARAMETER, (LENGTH*WIDTH), B-SPL DEGREE, INPUT DATA,'
	write(iu,'(A)') '--------------------------------'

	write(iu,'(5(E12.6,1X))') (x2d(i), i=1,jmt)

return

end subroutine OUT2DP


end module mod_inv2d

