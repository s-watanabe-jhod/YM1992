!*************************************************
! module
!	RESET2D(alpha_lc,detc,jdata_lc,jct_lc,jmt_lc)
!	FREE2D(detf,rsq2_lc,jdata_lc,jct_lc,jmt_lc)
!	_lc :: local variable
!*************************************************

module mod_rf2d

	implicit none

contains
!==========
subroutine RESET2D(alpha_lc,detc,jdata_lc,jct_lc,jmt_lc)

	use prm_matrix
	use prm_var

	implicit none

	real(8), intent(in) :: alpha_lc
	real(8), intent(out) :: detc
	integer, intent(out) :: jdata_lc, jct_lc, jmt_lc

	integer :: i, k, ii

	do i = 1, itotal
		y2d(i) = gz(i)
		do k = 1, jmp
			z2d(i,k) = rt(i,k)
		end do
	end do

	do ii = itotal +1, itotal + jcp
		y2d(ii) = 0.d0
	end do

	do k = 1, jmp
		do i = 1, jcp
			z2d(i+itotal,k) = cs2d(i,k) * alpha_lc
		end do
	end do

	jmt_lc = jmp
	jct_lc = jcp

	jdata_lc = itotal
	detc = detca

return

end subroutine RESET2D

!==========
subroutine FREE2D(detf,rsq2_lc,jdata_lc,jct_lc,jmt_lc)

	use prm_matrix
	use prm_var
	use mod_qrdec

	implicit none

	integer, intent(in) :: jdata_lc, jct_lc, jmt_lc
	real(8), intent(out) :: rsq2_lc
	real(8), intent(inout) :: detf

	real(8) :: wk1(kdata2)
	integer :: jtotal, jmtt, i

	jtotal = jdata_lc + jct_lc

	CALL QRDEC(z2d, y2d, detf, kdata2, jtotal, jmt_lc, wk1, lu2d, jmtt, 1)

	if(jmtt /= jmt_lc) then
		print *, ' ERROR NO SOLUTION FOR BASIC LSQS/ JMT,JMTT', jmt_lc, jmtt
		stop
	end if

	rsq2_lc = 0.d0

	do i = jmt_lc+1, jtotal
		rsq2_lc = rsq2_lc + y2d(i) * y2d(i)
	end do

	CALL RINV(x2d, y2d, z2d, kdata2, jmt_lc, 0)

return
end subroutine FREE2D

end module mod_rf2d

