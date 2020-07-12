!*************************************************
! module mod_subabic
!	RESET(alpha,detc,jdata_lc,jct_lc,jmt_lc)
!	FREEMIN(detf,rsq2_lc,jdata_lc,jct_lc,jmt_lc)
!	INEQMIN(detf)
!
!
!*************************************************

module mod_subabic

	implicit none

contains
!=========
subroutine RESET(alpha_lc,detc,jdata_lc,jct_lc,jmt_lc)

	use prm_matrix
	use prm_var
	use prm_inv

	implicit none

	real(8), intent(in) :: alpha_lc

	integer, intent(out) :: jdata_lc, jct_lc, jmt_lc
	real(8), intent(out) :: detc

	integer ::i, k

	jdata_lc = 0

	do i = 1, itotal
		if(ee(i) == 0.d0) cycle

		jdata_lc = jdata_lc + 1
		y(jdata_lc) = yy(i)

		do k = 1, jsp
			z(jdata_lc,k) = zz(i,k)
		end do

	end do

	y(jdata_lc+1:jdata_lc+jcp) = 0.d0

	do k = 1, jsp
		do i = 1, jcp
			z(i+jdata_lc,k) = cs(i,k) * alpha_lc
		end do
	end do

	jmt_lc = jsp
	jct_lc = jcp
	detc = detca

return
end subroutine RESET

!=========
subroutine FREEMIN(detf,rsq2_lc,jdata_lc,jct_lc,jmt_lc)

	use prm_matrix
	use prm_inv
	use mod_qrdec

	implicit none

	integer, intent(in) :: jdata_lc, jct_lc, jmt_lc
	real(8), intent(out) :: rsq2_lc, detf

	integer :: jtotal, jmtt, i
	real(8) :: wk1(kdata)

	jtotal = jdata_lc + jct_lc

!	write(*,*)jmt_lc,jmtt
	CALL QRDEC(z,y,detf,kdata,jtotal,jmt_lc,wk1,lu,jmtt,1)
!	write(*,*)'AAAAAAAAAAAAAA'

	if(jmtt /= jmt_lc) then
		print '(" "," ERROR NO SOLUTION FOR BASIC LSQS/ JMT,JMTT",2I5)', jmt_lc, jmtt
		write(77,'(" "," ERROR NO SOLUTION FOR BASIC LSQS/ JMT,JMTT",2I5)') jmt_lc, jmtt
		stop
	end if

	rsq2_lc = 0.d0

	do i = jmt_lc+1, jtotal
		rsq2_lc = rsq2_lc + y(i) * y(i)
	end do

!	write(*,*)'AAAAAAAAAAAAAA'
	CALL RINV(x,y,z,kdata,jmt_lc,0)
!	write(*,*)'AAAAAAAAAAAAAA'

return
end subroutine FREEMIN

!=========
subroutine INEQMIN(detf)

	use prm_matrix
	use prm_var
	use prm_inv
	use mod_qrdec
	use mod_nnls

	implicit none

	real(8), intent(inout) :: detf
	real(8) :: esp

	real(8) :: w(kmpa,kmp), ra(kmpa), eg(kmpa,kmp), fs(kmpa), ux(kmp)
	real(8) :: wk1(kmp), wk2(kmp), xrd
	integer :: la(kmp), i, j, k, ll, jmtt, jmt2

	esp = 1.d-12

!  +++W(TRANSPOSE)= INCC RI, HA-INCC X+++

	do i = 1, inncon

		j = incc(i)

		do k = 1, jmt
			w(k,i) = 0.d0
			wk1(k) = 0.d0
		end do

		wk1(j) = 1.d0

		CALL RINV(fs,wk1,z,kdata,jmt,j)

		do k = j, jmt
			w(k,i) = fs(k)
		end do

		w(jmt+1,i) = ha(j) - x(j)

	end do

! ----- LPDSOL(w,ra,kmpa,jmt+1,inncon,eg,fs,ux,wk1,wk2,la, lu) -----
! ----- LPDSOL(W,RA,KA  ,JMTA, INQ,   EG,FS,UX,WK1,WK2,LA, LB)
	fs(1:jmt) = 0.d0
	fs(jmt+1) = 1.d0

	CALL rtNNLS(w,fs,ux,kmpa,jmt+1,inncon,eg,ra,wk1,wk2,la,lu)

	if(dabs(ra(jmt+1)) < esp) stop ' INEQUALITIES ARE INCONSISTENT'

	do j = 1, jmt
		ra(j) = - ra(j) /ra(jmt+1)
	end do
! ----- /LPDSOL/ -----

!  ++SQUARED RESIDUAL DUE TO INEQ CONSTRAINT++

	do i = 1, jmt
		rsq2 = rsq2 + ra(i) * ra(i)
	end do

!  ++OPTIMUM SOLUTION UNDER INEQ CONSTRAINT++

	CALL RINV(wk1,ra,z,kdata,jmt,0)

	la(1:jmt) = 0

	do i = 1, jmt
		x(i) = x(i) + wk1(i)
	end do

!  +++ SEARCH FREE PARAMETER UNDER INEQ CONSTRAINT +++

	do i = 1, inncon

		j = incc(i)
		xrd = dabs( x(j) - ha(j) )

		if(xrd < esp) then
			la(j) = 1
			x(j) = ha(j)
		end if

	end do

	jmtt = 0

	do i = 1, jmt
		if(la(i) == 1) cycle
		jmtt = jmtt + 1
		lu(jmtt) = i
	end do

!  +++DET FOR MATRIX (A E A- ALPHA2 B)  UNDER INEQ++

	do i = 1, jmtt

		ll = lu(i)

		do j = 1, ll
			z(j,i) = z(j,ll)
		end do

		z(ll+1:jmt,i) = 0.d0

	end do

	CALL QRDEC(z,wk1,detf,kdata,jmt,jmtt,wk2,la,jmt2,0)

	if(jmt2 /= jmtt) then
		print *, "  NO SOL. FOR LSQS UNDER INEQ/ JMT2,JMTT", jmt2, jmtt
		stop
	end if

	jmt = jmtt

	print *, ' Non-NEGATIVE CONSTRAIN, FREE PARAMETER =', jmt
	print *, ' '

return
end subroutine INEQMIN

!=========

end module mod_subabic
