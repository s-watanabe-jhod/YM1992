!=================================================
! SUBROUTINES
!	OUTDAT(hd)
!	COVAR
!	OUTSOL(hd,flobspd,nsource)
!
!=================================================

module mod_inv

	implicit none

contains
!==========
subroutine OUTDAT(hd)

	use prm_matrix
	use prm_inv
	use prm_var

	use mod_trans

	implicit none

	character(len=24), intent(in) :: hd

	integer :: i, k
	character(len= 8) :: cp
	character(len=16) :: cq

	real(8) :: vec_th, vec_r, yoc(3), plong, plati

	write(21,*) ' ***** OBSERVED AND CALCULATED DISPLACEMENTS *****'
	write(21,'(A24,3X,I4)') hd, isn
	write(21,'(1X,3(I5,1X),E12.6," (ISN IH IV IL SIGMA)")') ih, iv, il, sigma
	write(21,'(1X,F10.3," IS ALREADY ADDED TO ALL VERTICAL DISP. DATA")') vrp
	write(21,*) 
	write(21,*) " DATA COVARIANCE"
	write(21,*) "                         HORIZ. DATA     VERTI. DATA    LENGTH CHG"
	write(21,'(A24,7X,F6.3,10X,F6.3,10X,F6.3)') "     FOR ZERO VALUE DATA", (fact(k), k=1,3)
	write(21,'(A24,7X,F6.3,10X,F6.3,10X,F6.3)') "  INCREMENT FOR 1.0 DATA", (fabs(k), k=1,3)

	do i = 1, itotal

		if(ia(i) == 1) cp = "HOR. DX "
		if(ia(i) == 2) cp = "HOR. DY "
		if(ia(i) == 3) cp = "HOR. DZ "
		if(ia(i) == 4) then
			cp = "LGT.CHG."
			write(21,'(I2,1X,A8," ST1",I4," ST2",I4,2X,"OBS ",E12.6,1X,F8.3,1X,"CAL ",E12.6)') &
				ia(i), cp, ib(i), ic(i), yy(i), ee(i), ysl(i)
			cycle
		end if

		if(ic(i) == 0) cq = "NO REF POINT CR."
		if(ic(i) == 1) cq = "   REF POINT CR."

		write(21,'(I2,1X,A8," ST=",I4,1X,"OBS ",E12.6,1X,F8.3,1X,"CAL ",E12.6,1X,A16)') &
			ia(i), cp, ib(i), yy(i), ee(i), ysl(i), cq

		if(ia(i) == 2) then

			CALL PLTXY(plati,plong,st(1,ib(i)),st(2,ib(i)),1)

			!--- OBS DATA ---
			vec_th = datan2(yy(i),yy(i-1)) * 180.d0 / pi
			vec_r  = dsqrt(yy(i)*yy(i) + yy(i-1)*yy(i-1))
			write(61,'(2X,2F10.3,1X,2F10.3)') plong, plati, vec_th, vec_r

			!--- CAL DATA ---
			vec_th = datan2(ysl(i),ysl(i-1)) * 180.d0 / pi
			vec_r  = dsqrt(ysl(i)*ysl(i) + ysl(i-1)*ysl(i-1))
			write(60,'(2X,2F10.3,1X,2F10.3)') plong, plati, vec_th, vec_r

			!--- O-C DATA ---
			yoc(1) = yy(i-1) - ysl(i-1)
			yoc(2) = yy(i) - ysl(i)
			vec_th = datan2(yoc(2),yoc(1)) * 180.d0 / pi
			vec_r  = dsqrt(yoc(2)*yoc(2) + yoc(2)*yoc(2))
			write(62,'(2X,2F10.3,1X,2F10.3)') plong, plati, vec_th, vec_r

		else if(ia(i) == 3) then

			CALL PLTXY(plati,plong,st(1,ib(i)),st(2,ib(i)),1)

			!--- OBS DATA ---
			vec_r = abs(yy(i))
			if(vec_r == yy(i)) vec_th =  90.d0
			if(vec_r /= yy(i)) vec_th = -90.d0
			write(61,'(2X,2F10.3,1X,2F10.3)') plong, plati, vec_th, vec_r

			!--- CAL DATA ---
			vec_r = abs(ysl(i))
			if(vec_r == ysl(i)) vec_th =  90.d0
			if(vec_r /= ysl(i)) vec_th = -90.d0
			write(60,'(2X,2F10.3,1X,2F10.3)') plong, plati, vec_th, vec_r

			!--- O-C DATA ---
			yoc(3) = yy(i) - ysl(i)
			vec_r = abs(yoc(3))
			if(vec_r == yoc(3)) vec_th =  90.d0
			if(vec_r /= yoc(3)) vec_th = -90.d0
			write(62,'(2X,2F10.3,1X,2F10.3)') plong, plati, vec_th, vec_r

		end if

	end do

return
end subroutine OUTDAT

!==========
subroutine COVAR

	use prm_matrix
	use prm_inv
	use prm_var

	use mod_qrdec

	implicit none

	integer :: i, j, k, ii
	real(8) :: sigma2, cc, u(kmp), v(kmp)

!print *, "jmt", jmt, kmp

!	do i = 1, kmp
!		do j = 1, kmp
!			cs(i,j) = 0.d0
!		end do
!	end do
	cs(1:kmp,1:kmp) = 0.d0
	sigma2 = sigma * sigma

	do i = jmt, 1, -1
		ii = i
		u(1:ii) = 0.d0
		u(ii)   = 1.d0

		CALL RINV(v,u,z,kdata,ii,0)

		do j = 1, ii
			z(i,j) = v(j)
		end do
	end do

	do i = 1, jmt
		do j = 1, i

			cc = 0.d0

			do k = 1, jmt
				cc = cc + z(k,i) * z(k,j)
			end do

			cs(i,j) = cc * sigma2
			cs(j,i) = cc * sigma2

		end do
	end do

	print *, ' COVARIANCE MATRIX IS COMPUTED.'

return
end subroutine COVAR

!==========
subroutine OUTSOL(hd,flobspd,nsource)

	use prm_matrix
	use prm_inv
	use prm_var

	implicit none

	character(len=24), intent(in) :: hd, flobspd
	integer, intent(in) :: nsource

	integer :: ls, kss, ktt, k, is, it, j, jj
	real(8) :: er

	write(22,*) "      ***** SOLUTION FOR GEODETIC DATA INVERSION *****"
	write(22,*) hd, "  INPUT DATA=", flobspd
	write(22,'(I3,A25)') nsource, " /NUMBER OF SOURCE REGION"

	write(22,*) "     ALPHA       RSQ2      SIGMA    JMT   JCT JDATA"
	write(22,'(3(E10.4,1X),1X,3(I5,1X))') alpha, rsq2, sigma, jmt, jct, jdata

	do ls = 1, nsource
		kss = js0(ls) + ndgs(ls)
		ktt = jt0(ls) + ndgs(ls)

		write(22,*) "SOURCE NUMBER"
		write(22,'(I5,5X,A24,"/JACOBI MATRIX VERSION")') ls, hd_jac(ls)
		write(22,'(3I5," SOURCE PRESENTATION(KS,KT,DEGREE OF B-SPLINE)")') js0(ls), jt0(ls), ndgs(ls)

		write(22,'(4I3," USED SOURCE TYPE(S-SLIP,D-SLIP,OPC,EXP)")') (just(k,ls), k=1,4)
		write(22,'(4I3," NON-NEGATIVE LEAST SQUARES ",3F8.1," (RSS,ROP,REX)")') (nnls(k,ls), k=1,4), rss(ls), rop(ls), rex(ls)
		write(22,'(4I3," MODEL PARAMETER REGION(S0,S1,T0,T1)")') (ius(k,ls), k=1,4)
		write(22,'(4I3," BOUNDARY CONSTRAIN(0=NOCONSTRAIN,1=CONSTRAIN)")') (iboc(k,ls), k=1,4)

		write(22,*) " SOURCE    IS    IT         M.P.         E.E"

		do k = 1, 4

			if(just(k,ls) == 0) cycle

			do it = 1, ktt
				do is = 1, kss

					j  = is + (it-1) * kss
					jj = idmp(j,k,ls)

					if(jj == 0) cycle

					er = dsqrt(cs(jj,jj))

					write(22,'(2X,I5,1X,I5,1X,I5,1X,E12.6,1X,E12.6)') k, is, it, x(jj), er

				end do
			end do

		end do

	end do


return
end subroutine OUTSOL

!==========

end module mod_inv

