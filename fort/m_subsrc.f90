!=================================================
! SUBROUTINES
!	FSNPA(we,sa,sb,xp,yp,zp,is,it)
!	SVESE(q,cm,bsv,idb,ls,nod2s)
!
!=================================================

module mod_subsrc

	implicit none

contains
!==========
subroutine FSNPA(we,sa,sb,xp,yp,zp,is,it)

	use prm_matrix
	use prm_var

	implicit none

	integer, intent(in) :: is, it
	real(8), intent(out) :: we(3), sa(3), sb(3), xp, yp, zp

	integer :: mu, mv, muu, mvv
	real(8) :: u1, u2, u3, v1, v2, v3, w1, w2, w3, uu, vv, ww

	mu = is
		if(is == ms) mu = ms - 1
	mv = it
		if(it == mt) mv = mt - 1

	muu = mu + 1
	mvv = mv + 1

	u1 = xtr(muu)     - xtr(mu)
	u2 = ytr(muu,mv ) - ytr(mu,mv)
	u3 = ztr(muu,mv ) - ztr(mu,mv)
	v2 = ytr(mu ,mvv) - ytr(mu,mv)
	v3 = ztr(mu ,mvv) - ztr(mu,mv)

	w1 = u2*v3 - u3*v2
	w2 =       - u1*v3
	w3 = u1*v2

	v1 = w2*u3 - w3*u2
	v2 = w3*u1 - w1*u3
	v3 = w1*u2 - w2*u1

	uu = dsqrt(u1*u1 + u2*u2 + u3*u3)
	vv = dsqrt(v1*v1 + v2*v2 + v3*v3)
	ww = dsqrt(w1*w1 + w2*w2 + w3*w3)

	sa(1) = u1 / uu
	sa(2) = u2 / uu
	sa(3) = u3 / uu

	sb(1) = v1 / vv
	sb(2) = v2 / vv
	sb(3) = v3 / vv

	if(vv == 0.d0) sb(1:3) = 0.d0

	we(1) = w1 / ww
	we(2) = w2 / ww
	we(3) = w3 / ww

	xp = xtr(is)
	yp = ytr(is,it)
	zp = ztr(is,it)

return
end subroutine FSNPA

!==========
subroutine SVESE(q,cm,bsv,idb,ls,nod2s)

	use prm_matrix
	use prm_var
	use prm_inv

	implicit none

	integer, intent(in) :: idb(16), ls, nod2s
	real(8), intent(in) :: bsv(16)

	real(8), intent(out) :: q(4), cm(4,4)

	integer :: i, j, ll, lp, lq, l1, l2

	do i = 1, 4

		q(i) = 0.d0

		do ll = 1, nod2s
			if(idb(ll) == 0) cycle

			j = idmp(idb(ll), i, ls)
				if(j == 0) cycle

			q(i) = q(i) + bsv(ll) * x(j)

		end do

		do j = 1, i

			cm(j,i) = 0.d0

			do lp = 1, nod2s
				if(idb(lp) == 0) cycle

				l1 = idmp(idb(lp), i, ls)
					if(l1 == 0) cycle

				do lq = 1, nod2s
					if(idb(lq) == 0) cycle

					l2 = idmp(idb(lq), i, ls)
						if(l2 == 0) cycle

					cm(j,i) = cm(j,i) + bsv(lp) * bsv(lq) * cs(l1,l2)

				end do
			end do
		end do
	end do

return
end subroutine SVESE

!==========

end module mod_subsrc

