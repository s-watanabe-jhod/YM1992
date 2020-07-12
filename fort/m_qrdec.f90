!*************************************************
! module
!	QRDEC(z,y,detl,n1,n,m,u,lb,j1,ind)
!	RINV(x,y,z,n1,m,ind)
!
!	no common
!*************************************************

module mod_qrdec

	implicit none

contains
!=========
subroutine QRDEC(z_lc,y_lc,detl,n1,n,m,u,lb,j1,ind)

	implicit none

	integer, intent(in) :: n1, n, m, ind
	real(8), intent(inout) :: z_lc(n1,m), y_lc(n)
	integer, intent(out) :: lb(m), j1
	real(8), intent(out) :: detl, u(n)

	integer :: i, k, kk
	real(8) :: esp, sm, b, sl

	esp = 1.d-50!12
	detl = 0.d0
	j1 = 1

	do k = 1, m

		lb(k) = 0
		sm = 0.d0

		do i = j1, n
			u(i) = z_lc(i,k)
			sm = sm + z_lc(i,k) * z_lc(i,k)
		end do

		b = dsqrt(sm)

		if(b < esp) write(*,*)k,b,sm
		if(b < esp) cycle
		if(u(j1) > 0.d0) b = -b

		u(j1) = u(j1) - b

		do kk = k+1, m
			sl = 0.d0

			do i = j1, n
				sl = sl + z_lc(i,kk) * u(i)
			end do
			sl = sl / u(j1) / b
			do i = j1, n
				z_lc(i,kk) = z_lc(i,kk) + u(i) * sl
			end do

		end do

		z_lc(j1,k) = b
		detl = detl + dlog(dabs(b))

		if(ind /= 0) then
			sl = 0.d0
			do i = j1, n
				sl = sl + y_lc(i) * u(i)
			end do
			sl = sl / u(j1) / b
			do i = j1, n
				y_lc(i) = y_lc(i) + u(i) * sl
			end do
		end if

		lb(k) = j1
		j1 = j1 + 1

!	write(*,*)j1
	end do

	j1 = j1 - 1

return
end subroutine QRDEC

!==========
subroutine RINV(x_lc, y_lc, z_lc, n1, m, ind)

	implicit none

	integer, intent(in) :: n1, m, ind
	real(8), intent(in) :: y_lc(m), z_lc(n1,m)
	real(8), intent(inout) :: x_lc(m)

	integer :: k, k2
	real(8) :: uf

	if(ind == 0) then
		do k = m, 1, -1
			uf = y_lc(k)
			do k2 = k+1, m
				uf = uf - x_lc(k2) * z_lc(k,k2)
			end do
			x_lc(k) = uf / z_lc(k,k)
		end do
		return
	end if

	do k = ind, m
		uf = y_lc(k)
		do k2 = ind, k-1
			uf = uf - x_lc(k2) * z_lc(k2,k)
		end do
		x_lc(k) = uf / z_lc(k,k)
	end do

return
end subroutine RINV

end module mod_qrdec

