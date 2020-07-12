!*************************************************
! module mod_subabic
!	RESET(alpha,detc,jdata_lc,jct_lc,jmt_lc)
!	FREEMIN(detf,rsq2_lc,jdata_lc,jct_lc,jmt_lc)
!
!*************************************************

module mod_subabic

	implicit none

contains
!=========
subroutine RESET(alpha_lc,detc,jdata_lc)

	use prm_matrix
	use prm_inv

	implicit none

	real(8), intent(in) :: alpha_lc

	integer, intent(out) :: jdata_lc !, jct_lc, jmt_lc
	real(8), intent(out) :: detc

	integer :: i, k

	jdata_lc = 0

	! E_prior がゼロ（＝使わないデータ）を取り除く作業

	do i = 1, Itotal

		if(E_prior(i) == 0.d0) cycle

		jdata_lc = jdata_lc + 1
		Y(jdata_lc) = Y_data(i)

		do k = 1, Jsp
			Z_HaG(jdata_lc,k) = Z_Green(i,k)
		end do

	end do

	Y(jdata_lc+1:jdata_lc+Jcp) = 0.d0

	! Z_HaG の下に Constrain Matrix * alpha を加える

	do k = 1, Jsp
		do i = 1, Jcp
			Z_HaG(i+jdata_lc,k) = Cs(i,k) * alpha_lc
		end do
	end do

	!jmt_lc = Jsp
	!jct_lc = Jcp
	detc = Detca

return
end subroutine RESET

!=========
subroutine FREEMIN(detf,rsq2_lc,jdata_lc)

	use prm_matrix
	use prm_inv
	use mod_qrdec

	implicit none

	integer, intent(in)  :: jdata_lc !, jct_lc, jmt_lc
	real(8), intent(out) :: rsq2_lc, detf

	integer :: jtotal, jmtt, i
	integer :: lu(Kmp)    ! lu  は QRDEC のアウトプット。特に使用されてない。
	real(8) :: wk1(Kdata) ! wk1 は QRDEC のアウトプット。特に使用されてない。

	jtotal = jdata_lc + Jcp

	CALL QRDEC(Z_HaG,Y, detf, Kdata, jtotal, Jsp, wk1,lu, jmtt,1)

	if(jmtt /= Jsp) then
		print    '(" "," ERROR NO SOLUTION FOR BASIC LSQS/ JMT,JMTT",2I5)', Jsp, jmtt
		write(77,'(" "," ERROR NO SOLUTION FOR BASIC LSQS/ JMT,JMTT",2I5)') Jsp, jmtt
		stop
	end if

	rsq2_lc = 0.d0

	do i = Jsp+1, jtotal
		rsq2_lc = rsq2_lc + Y(i) * Y(i)
	end do

	CALL RINV(X_model,Y,Z_HaG,Kdata,Jsp,0)

return
end subroutine FREEMIN
!=========

end module mod_subabic
