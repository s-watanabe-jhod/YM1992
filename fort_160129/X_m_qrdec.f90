!*************************************************
!	no common
!*************************************************

module mod_qrdec

	implicit none

contains
!=========
subroutine QRDEC(z_lc,y_lc, detl, nkdata, njtotal, njmt, u, lb, jmtt, ind)

	implicit none

	integer, intent(in) :: nkdata, njtotal, njmt, ind
	real(8), intent(inout) :: z_lc(nkdata, njmt), y_lc(njtotal)
	integer, intent(out) :: lb(njmt), jmtt
	real(8), intent(out) :: detl, u(njtotal)
	
	integer :: i, k, kk
	real(8) :: eps, sm, b, sl
	
	eps  = 1.d-12
	detl = 0.d0
	jmtt = 1
	lb   = 0
	
	do k = 1, njmt
		
		sm = 0.d0
		
		do i = jmtt, njtotal
			u(i) = z_lc(i,k)
			sm = sm + u(i) * u(i)
		end do
		
		b = dsqrt(sm)
		
		if(b < eps) cycle
		if(u(jmtt) > 0.d0) b = -b
		
		u(jmtt) = u(jmtt) - b
		
		do kk = k+1, njmt
			
			sl = 0.d0
			
			do i = jmtt, njtotal
				sl = sl + z_lc(i,kk) * u(i)
			end do
			
			sl = sl / u(jmtt) / b
			
			do i = jmtt, njtotal
				z_lc(i,kk) = z_lc(i,kk) + u(i) * sl
			end do
			
		end do
		
		z_lc(jmtt,k) = b
		detl = detl + dlog(dabs(b))
		
		if(ind /= 0) then
			
			sl = 0.d0
			
			do i = jmtt, njtotal
				sl = sl + y_lc(i) * u(i)
			end do
			
			sl = sl / u(jmtt) / b
			
			do i = jmtt, njtotal
				y_lc(i) = y_lc(i) + u(i) * sl
			end do
			
		end if
		
		lb(k) = jmtt
		jmtt = jmtt + 1
		
	end do
	
	jmtt = jmtt - 1
	
return
end subroutine QRDEC

!==========
subroutine RINV(x_lc, y_lc, z_lc, nkdata, jmtt, ind)

	implicit none

	integer, intent(in) :: nkdata, jmtt, ind
	real(8), intent(in) :: y_lc(jmtt), z_lc(nkdata, jmtt)
	real(8), intent(inout) :: x_lc(jmtt)

	integer :: k, j
	real(8) :: uf

	if(ind == 0) then
		do k = jmtt, 1, -1
			
			uf = y_lc(k)
			
			do j = k+1, jmtt
				uf = uf - x_lc(j) * z_lc(k,j)
			end do
			
			x_lc(k) = uf / z_lc(k,k)
			
		end do
		return
	end if

	do k = ind, jmtt
		
		uf = y_lc(k)
		
		do j = ind, k-1
			uf = uf - x_lc(j) * z_lc(j,k)
		end do
		
		x_lc(k) = uf / z_lc(k,k)
		
	end do

return
end subroutine RINV

end module mod_qrdec

