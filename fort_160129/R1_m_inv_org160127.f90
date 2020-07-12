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
	use mod_trans
	
	implicit none
	
	character(len=24), intent(in) :: hd
	
	integer :: i, k
	character(len= 8) :: cp
	character(len=16) :: cq
	
	real(8), parameter :: pi = 3.14159265358979324d0
	real(8) :: vec_th, vec_r, yoc(3), plati, plong
	
	write(21,*) ' ***** OBSERVED AND CALCULATED DISPLACEMENTS *****'
	write(21,'(A24,3X,I4)') hd, isn
	write(21,'(1X,3(I5,1X),E12.6," (ISN IH IV IL SIGMA)")') ih, iv, il, sigma
	write(21,'(1X,F10.3," IS ALREADY ADDED TO ALL VERTICAL DISP. DATA")') vrp
	write(21,*) 
	write(21,*) " DATA COVARIANCE"
	write(21,*) "                         HORIZ. DATA     VERTI. DATA    LENGTH CHG"
	write(21,'(A24,7X,F6.3,10X,F6.3,10X,F6.3)') "     FOR ZERO VALUE DATA", (fact(k), k=1,3)
	write(21,'(A24,7X,F6.3,10X,F6.3,10X,F6.3)') "  INCREMENT FOR 1.0 DATA", (fabs(k), k=1,3)
	
	do i = 1, Itotal
		
		if(ia(i) == 1) cp = "HOR. DX "
		if(ia(i) == 2) cp = "HOR. DY "
		if(ia(i) == 3) cp = "HOR. DZ "
		if(ia(i) == 4) then
			cp = "LGT.CHG."
			write(21,'(I2,1X,A8," ST1",I4," ST2",I4,2X,"OBS ",E12.6,1X,F8.3,1X,"CAL ",E12.6)') &
				ia(i), cp, Ib(i), Ic(i), Y_data(i), E_prior(i), Y_calc(i)
			cycle
		end if
		
		if(Ic(i) == 0) cq = "NO REF POINT CR."
		if(Ic(i) == 1) cq = "   REF POINT CR."
		
		write(21,'(I2,1X,A8," ST=",I4,1X,"OBS ",E12.6,1X,F8.3,1X,"CAL ",E12.6,1X,A16)') &
			Ia(i), cp, Ib(i), Y_data(i), E_prior(i), Y_calc(i), cq
		
		if(ia(i) == 2) then
			
			CALL LL2XY(plati,plong, St(1,Ib(i)),St(2,Ib(i)),1, Alat0,Alng0,Icord)
			
			!--- OBS DATA ---
			vec_th = datan2(Y_data(i),Y_data(i-1)) * 180.d0 / pi
			vec_r  = dsqrt(Y_data(i)*Y_data(i) + Y_data(i-1)*Y_data(i-1))
			write(61,'(2X,2F10.3,1X,2F10.3)') plong, plati, vec_th, vec_r
			
			!--- CAL DATA ---
			vec_th = datan2(Y_calc(i),Y_calc(i-1)) * 180.d0 / pi
			vec_r  = dsqrt(Y_calc(i)*Y_calc(i) + Y_calc(i-1)*Y_calc(i-1))
			write(60,'(2X,2F10.3,1X,2F10.3)') plong, plati, vec_th, vec_r
			
			!--- O-C DATA ---
			yoc(1) = Y_data(i-1) - Y_calc(i-1)
			yoc(2) = Y_data(i) - Y_calc(i)
			vec_th = datan2(yoc(2),yoc(1)) * 180.d0 / pi
			vec_r  = dsqrt(yoc(2)*yoc(2) + yoc(2)*yoc(2))
			write(62,'(2X,2F10.3,1X,2F10.3)') plong, plati, vec_th, vec_r
			
		else if(Ia(i) == 3) then
			
			CALL LL2XY(plati,plong,St(1,Ib(i)),St(2,Ib(i)),1, Alat0,Alng0,Icord)
			
			!--- OBS DATA ---
			vec_r = abs(Y_data(i))
			if(vec_r == Y_data(i)) vec_th =  90.d0
			if(vec_r /= Y_data(i)) vec_th = -90.d0
			write(61,'(2X,2F10.3,1X,2F10.3)') plong, plati, vec_th, vec_r
			
			!--- CAL DATA ---
			vec_r = abs(Y_calc(i))
			if(vec_r == Y_calc(i)) vec_th =  90.d0
			if(vec_r /= Y_calc(i)) vec_th = -90.d0
			write(60,'(2X,2F10.3,1X,2F10.3)') plong, plati, vec_th, vec_r
			
			!--- O-C DATA ---
			yoc(3) = Y_data(i) - Y_calc(i)
			vec_r = abs(yoc(3))
			if(vec_r == yoc(3)) vec_th =  90.d0
			if(vec_r /= yoc(3)) vec_th = -90.d0
			write(62,'(2X,2F10.3,1X,2F10.3)') plong, plati, vec_th, vec_r
			
		end if
		
	end do
	
return
end subroutine OUTDAT

!==========
subroutine COVAR(sgm)

	use prm_matrix
	use prm_inv
	use mod_qrdec
	
	implicit none
	
	real(8), intent(in)  :: sgm
	
	integer :: i, j, k, ii
	real(8) :: cc, u(Kmp), v(Kmp)
	
	Cov(1:Kmp,1:Kmp) = 0.d0
	
	! == 逆行列を計算？
	do i = Jmt, 1, -1
		ii = i
		u(1:ii) = 0.d0
		u(ii)   = 1.d0
		
		CALL RINV(v,u,Z,Kdata,ii,0)
		
		do j = 1, ii
			Z(i,j) = v(j)
		end do
	end do
	! ==
	
	do i = 1, Jmt
		do j = 1, i
			
			cc = 0.d0
			
			do k = 1, Jmt
				cc = cc + Z(k,i) * Z(k,j)
			end do
			
			Cov(i,j) = cc * sgm*sgm
			Cov(j,i) = cc * sgm*sgm
			
		end do
	end do
	
	print *, ' COVARIANCE MATRIX IS COMPUTED.'
	
return
end subroutine COVAR

!==========
subroutine OUTSOL(hd,flobspd,nsource)

	use prm_matrix
	use prm_inv
	
	implicit none
	
	character(len=24), intent(in) :: hd, flobspd
	integer, intent(in) :: nsource
	
	integer :: ls, kss, ktt, k, is, it, j, jj
	real(8) :: er
	
	write(22,*) "      ***** SOLUTION FOR GEODETIC DATA INVERSION *****"
	write(22,*) hd, "  INPUT DATA=", flobspd
	write(22,'(I3,A25)') nsource, " /NUMBER OF SOURCE REGION"
	
	write(22,*) "     ALPHA       RSQ2      SIGMA    JMT   JCT JDATA"
	write(22,'(3(E10.4,1X),1X,3(I5,1X))') Alpha, Rsq2, Sigma, Jmt, Jct, Jdata
	
	do ls = 1, nsource
		kss = Js0(ls) + Ndgs(ls)
		ktt = Jt0(ls) + Ndgs(ls)
		
		write(22,*) "SOURCE NUMBER"
		write(22,'(I5,5X,A24,"/JACOBI MATRIX VERSION")') ls, Hd_jac(ls)
		write(22,'(3I5," SOURCE PRESENTATION(KS,KT,DEGREE OF B-SPLINE)")') Js0(ls), Jt0(ls), Ndgs(ls)
		
		write(22,'(4I3," USED SOURCE TYPE(S-SLIP,D-SLIP,OPC,EXP)")') (Just(k,ls), k=1,4)
		!write(22,'(4I3," NON-NEGATIVE LEAST SQUARES ",3F8.1," (RSS,ROP,REX)")') (0.0, k=1,4), rss(ls), rop(ls), rex(ls)
		write(22,'(4I3," MODEL PARAMETER REGION(S0,S1,T0,T1)")') (Ius(k,ls), k=1,4)
		write(22,'(4I3," BOUNDARY CONSTRAIN(0=NOCONSTRAIN,1=CONSTRAIN)")') (Iboc(k,ls), k=1,4)
		
		write(22,*) " SOURCE    IS    IT         M.P.         E.E"
		
		do k = 1, 4
			
			if(Just(k,ls) == 0) cycle
			
			do it = 1, ktt
				do is = 1, kss
					
					j  = is + (it-1) * kss
					jj = Idmp(j,k,ls)
					
					if(jj == 0) cycle
					
					er = dsqrt(Cov(jj,jj))
					
					write(22,'(2X,I5,1X,I5,1X,I5,1X,E12.6,1X,E12.6)') k, is, it, X_model(jj), er
					
				end do
			end do
			
		end do
	
	end do

return
end subroutine OUTSOL

!==========

end module mod_inv

