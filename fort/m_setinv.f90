!=================================================
! SUBROUTINES
!	INPUTD(iu)
!	SETUP(flobspd,nsource,icm)
!	SETWT
!	CHGPM
!	CHGIMP
!
!=================================================

module mod_setinv

	implicit none

contains
!=========
subroutine INPUTD(iu)

	use prm_matrix
	use prm_var
	use prm_inv
	use mod_trans

	implicit none

	integer, intent(in) :: iu

	integer :: i, i1, i2, ii, j, j1, j2, k
	integer :: nr(kis), nst(15,kis), k3(kis), nl(2,kis)
	real(8) :: da(kis), dy(kis), dl(kis), aa, al, ea, dd
	real(8) :: sx, sy, sl, rx, ry, ua, ub

	read(iu,*) isn, ih, iv, il

	itotal = 0
	iinit = ih * 2 + iv

	if(iinit > kis) stop 'INSUFFICINT ARRAY DIM(KIS)'
	if(iinit > kio) stop 'INSUFFICINT ARRAY DIM(KIO)'
	if(il > kis) stop 'INSUFFICINT ARRAY DIM(KIS)'

	ea = 0.d0

	!--------------------------------------------

	do i = 1, ih

		read(iu,*) ii, aa, al, da(i), dy(i), nr(i), (nst(k,i), k=1,nr(i))
		if(ii /= i) stop ' INPUT IH DATA ERROR'

		CALL PLTXY(aa,al,st(1,i),st(2,i),0)

	end do

	do i = ih+1, isn

		read(iu,*) ii, aa, al, da(i), nr(i), (nst(k,i), k=1,nr(i))
		if(ii /= i) stop ' INPUT IV DATA ERROR'

		CALL PLTXY(aa,al,st(1,i),st(2,i),0)

		da(i) = da(i) + vrp

	end do

	do i =1, il
		read(iu,*) ii, nl(1,i), nl(2,i), dl(i), k3(i)
		if(ii /= i) stop ' INPUT DATA ERROR'
	end do

		print *, ' DATA INPUT OK. (ISN,IH,IV,IL)= ', isn, ih, iv, il
!		write(77,*), ' DATA INPUT OK. (ISN,IH,IV,IL)= ', isn, ih, iv, il


	!--- INPUT END -------------------------------

	do i = 1, ih

		if(nst(1,i) == 0) cycle

		i1 = itotal + 1
		i2 = itotal + 2

		itotal = itotal + 2

		ia(i1) = 1
		ia(i2) = 2
		ib(i1) = i
		ib(i2) = i
		ic(i1) = 0
		ic(i2) = 0

		yy(i1) = da(i)
		yy(i2) = dy(i)

		dd = da(i) * da(i) + dy(i) * dy(i)

		ea = dble(nst(1,i))
		ee(i1) = fact(1) * dsqrt(ea*ea + fabs(1)*dd)
		ee(i2) = ee(i1)

	end do

	do i = ih+1, isn

		if(nst(1,i) == 0) cycle

		itotal = itotal + 1

		ia(itotal) = 3
		ib(itotal) = i
		ic(itotal) = 0

		yy(itotal) = da(i)
		ea = dble(nst(1,i))
		ee(itotal) = fact(2) * dsqrt(ea*ea + fabs(2)*da(i)*da(i))

	end do

	do i = 1, il

		if(k3(i) == 0) cycle

		itotal = itotal + 1

		ia(itotal) = 4
		ib(itotal) = nl(1,i)
		ic(itotal) = nl(2,i)

		yy(itotal) = dl(i)
		ee(itotal) = fact(3) * dsqrt(dble(k3(i)**2) + fabs(3)*dl(i)*dl(i))
	end do

	if(itotal > kio) stop ' ITOTAL,INSUFFICIENT ARRAY DIM (KIO)'

	print *, ' TOTAL NUMBER OF DATA =', itotal
	write(77,*) ' TOTAL NUMBER OF DATA =', itotal

!==== CALL SETAG(kio,kis,itolal,ih,ia,ib,ic,nr,nst,ag,st) ===

	ag = 0.d0 !ag(1:kio,1:kio) = 0.d0

	do j = 1, itotal

		if(ia(j) == 1) cycle

		if(ia(j) == 2) then
			i = ib(j)
			ag(j-1, 2*i-1) = 1.d0
			ag(j, 2*i) = 1.d0

			if(nr(i) == 1) cycle

			ic(j-1) = 1
			ic(j) = 1

			if(nr(i) == 2) then
				j1 = nst(2,i)
				ag(j-1, 2*j1-1) = ag(j-1, 2*j1-1) - 1.d0
				ag(j, 2*j1) = ag(j, 2*j1) - 1.d0

			else if(nr(i) >= 3) then

				do k = 2, nr(i)-1

					j1 = nst(k,i)
					j2 = nst(k+1,i)

					sx = st(1,j2) - st(1,j1)
					sy = st(2,j2) - st(2,j1)
					sl = dsqrt(sx*sx + sy*sy)
					sx = sx / sl
					sy = sy / sl

					rx = (st(1,i) - st(1,j1)) / sl
					ry = (st(2,i) - st(2,j1)) / sl

					ua = rx * sx + ry * sy
					ub = rx * sy - ry * sx

					ag(j-1, 2*j1-1) = ag(j-1, 2*j1-1) + ua - 1.d0
					ag(j-1, 2*j1  ) = ag(j-1, 2*j1  ) + ub
					ag(j-1, 2*j2-1) = ag(j-1, 2*j2-1) - ua
					ag(j-1, 2*j2  ) = ag(j-1, 2*j2  ) - ub

					ag(j, 2*j1-1) = ag(j, 2*j1-1) - ub
					ag(j, 2*j1  ) = ag(j, 2*j1  ) + ua - 1.d0
					ag(j, 2*j2-1) = ag(j, 2*j2-1) + ub
					ag(j, 2*j2  ) = ag(j, 2*j2  ) - ua

				end do

			end if

		else if(ia(j) == 3) then
			i = ib(j)
			ag(j, i+ih) = 1.d0

			if(nr(i) /= 2) cycle

			ag(j, nst(2,i)+ih) = -1.d0
			ic(j) = 1

		else if(ia(j) == 4) then

			sx = st(1,ib(j)) - st(1,ic(j))
			sy = st(2,ib(j)) - st(2,ic(j))
			sl = dsqrt(sx*sx + sy*sy)
			sx = sx / sl
			sy = sy / sl

			ag(j, 2*ib(j)-1) =  sx
			ag(j, 2*ib(j)  ) =  sy
			ag(j, 2*ic(j)-1) = -sx
			ag(j, 2*ic(j)  ) = -sy

		end if
	end do

return
end subroutine INPUTD

!=========
subroutine SETUP(flobspd,nsource,icm,icheck)

	use prm_matrix
	use prm_var
	use prm_inv
	use mod_subsetup

	implicit none

	integer, intent(in) :: nsource, icheck
	character(len=24), intent(in) :: flobspd
	integer, intent(out) :: icm

	integer :: k, lsource, ls

	icm = 0

	!===== CALL INITIAL =====

	jsp = 0
	jcp = 0
	inncon = 0
	idmp = 0 !idmp(1:kfg,1:4,1:5) = 0

	detca = 0.d0
	gi = 0.d0 !gi(1:kmp,1:kmp) = 0.d0
	cs = 0.d0 !cs(1:kmp,1:kmp) = 0.d0

	do k = 1, kmp
		gi(k,k) = 1.d0
	end do

	!===== =====

	do lsource = 1, nsource

		ls = lsource

		if(ipl(ls) == 0) open(31,file=fcof(ls),status='old')
		if(ipl(ls) == 0) close(31)

		nbi = (js0(ls) + ndgs(ls)) * (jt0(ls) + ndgs(ls))

		print *, '------------------------------------------'
		print *, ' SETUP JACOBI MTRX. LSOURCE= ', ls, ' NBI=', nbi
		write(77,*) ' SETUP JACOBI MTRX. LSOURCE= ', ls, ' NBI=', nbi

		if(icheck == 0) then
			open(32,file=invjac(ls),form='unformatted')
!!			!!,status='old',form='unformatted')
			CALL INPUTG(hd_jac(ls), flobspd, ls, 32)
			close(32)
		end if

		CALL SETMP(js0(ls), jt0(ls), ndgs(ls), ls)
		CALL SETBM(ls, js0(ls), jt0(ls), ndgs(ls))

		if(nnls(1,ls) == 1) CALL INSL(rss(ls), ls)
		if(nnls(3,ls) == 1) CALL INNG(rop(ls), ls, 3)
		if(nnls(4,ls) == 1) CALL INNG(rex(ls), ls, 4)

		icm = icm + nnls(1,ls) + nnls(3,ls) + nnls(4,ls)
		print *, '------------------------------------------'

	end do

return
end subroutine SETUP

!=========
subroutine SETWT

	use prm_matrix
	use prm_var
	use prm_inv

	implicit none

	integer :: i, k, l
	real(8) :: wk(kio), ec

	do k = 1, jsp

		do i = 1, itotal
			wk(i) = 0.d0
			do l = 1, iinit
				wk(i) = wk(i) + ag(i,l) * zz(l,k)
			end do
		end do

		do i = 1, itotal
			zz(i,k) = wk(i)
		end do

	end do

	do i = 1, itotal

		ec = dabs(ee(i))
		yy(i) = yy(i) / ec

		do k = 1, jsp
			zz(i,k) = zz(i,k) / ec
		end do

	end do

	print *, " SET WEIGHT OK"

return
end subroutine SETWT

!==========
subroutine CHGPM

	use prm_matrix
	use prm_var
	use prm_inv

	implicit none

	integer :: i, j, k
	real(8) :: wk(kmp)

!  ++++  ZZ = GI * ZZ ++++

	do i = 1, itotal
		do k = 1, jsp
			wk(k) = zz(i,k)
			zz(i,k) = 0.d0
		end do

		do j = 1, jsp
			do k = 1, jsp
				zz(i,j) = zz(i,j) + wk(k) * gi(k,j)
			end do
		end do
	end do

!   ++++ CS = CS * GI ++++

	do i = 1, jcp
		do k = 1, jsp
			wk(k) = cs(i,k)
			cs(i,k) = 0.d0
		end do

		do j = 1, jsp
			do k = 1, jsp
				cs(i,j) = cs(i,j) + wk(k) * gi(k,j)
			end do
		end do
	end do

	jnp = jsp

	print *, ' TRANSFORM MODEL PARAMETERS USING INEQLITY EQS.'

return
end subroutine CHGPM
!==========
subroutine CHGIMP

	use prm_matrix
	use prm_var
	use prm_inv

	implicit none

	integer :: i, j, k, l2, l3
	real(8) :: u(kmp)

!  ++++  X = GI * X ++++

	do i = 1, jnp
		u(i) = x(i)
		x(i) = 0.d0
	end do

	do i = 1, jnp
		do k = 1, jnp
			x(i) = x(i) + gi(i,k) * u(k)
		end do
	end do

!   ++++ CS = GI * CS GI(T) ++++

	do j = 1, jmt

		do i = 1, jnp
			u(i) = cs(j,i)
			cs(j,i) = 0.d0
		end do

		do i = 1, jnp
			do k = 1, jmt
				l2 = lu(k)
				cs(j,i) = cs(j,i) + gi(i,l2) * u(k)
			end do
		end do

	end do

	do j = 1, jnp

		do i = 1, jnp
			u(i) = cs(i,j)
			cs(i,j) = 0.d0
		end do

		do i = 1, jnp
			do k = 1, jmt
				l3 = lu(k)
				cs(i,j) = cs(i,j) + gi(i,l3) * u(k)
			end do
		end do

	end do

	print *, ' MODEL PARAMETERS ARE INVERSE-TRANSFORMED.'

return
end subroutine CHGIMP
!==========

end module mod_setinv


