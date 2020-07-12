!=================================================
! SUBROUTINES
!	INPUTG(subhd_jac,flobspd,LS,IU)
!	SETMP(KS00,KT00,NDG,LS)
!	SETBM(LS,KS00,KT00,NDG)
!
!	INSL(rs,ls)
!	INNG(rs,ls,ktyp)
!
!=================================================

module mod_subsetup

	implicit none

contains
!==========
subroutine INPUTG(subhd_jac,flobspd,ls,iu)

	!===== INPUT zz(kio,kmp), slms(3,kfg,5), slmd(3,kfg,5) =====

	use prm_matrix
	use prm_var
	use prm_inv

	implicit none

	integer, intent(in) :: ls, iu
	character(len=24), intent(in) :: subhd_jac, flobspd

	integer :: i, ild, j, j1, k, ja
	real(8) :: ztem(kfg)
	character(len=24) :: ca, c1, c2, c3

	ca = '**** JACOBI MATRIX **** '

	if(nbi > kfg) stop ' INSUFFICIENT ARRAY DIM(KFG)'

	ja = 0

	read(iu) c1, c2, c3
		if(c1 /= ca .or. c2 /= subhd_jac .or. c3 /= flobspd) then
			print *, "c1 ca     ", c1, ca
			print *, "c2 hd_jac ", c2, subhd_jac
			print *, "c3 flobspd", c3, flobspd
			print *, ' INPUTG ERROR'
			stop
		end if

	print *, ' JUST =', (just(k,ls), k=1,4)
	write(77,*) ' JUST =', (just(k,ls), k=1,4)

	do i = 1, ih
		do ild = 1, 2
			do j = 1, 4

				read(iu) (ztem(j1), j1=1,nbi)
				if(just(j,ls) == 0) cycle

				ja = nbi * (just(j,ls) - 1) + jsp

				do j1 = 1, nbi
					zz(2*i - 2 + ild, ja + j1) = ztem(j1)
				end do

			end do
		end do
	end do

	do i = 1, iv
		do j = 1, 4

			read(iu) (ztem(j1), j1=1,nbi)
			if(just(j,ls) == 0) cycle

			ja = nbi * (just(j,ls) - 1) + jsp

			do j1 = 1, nbi
				zz(i + 2*ih, ja + j1) = ztem(j1)
			end do

		end do
	end do

	if(ja+nbi > kmp) stop ' INSUFFICIENT ARRAY DIM(KMP)'

	print *, ' READ JACOBI MATRIX. OK'

	read(iu) ( (slms(i,j,ls), i=1,3), (slmd(i,j,ls), i=1,3), j=1,nbi )

return
end subroutine INPUTG

!==========
subroutine SETMP(ks00,kt00,ndg,ls)

	use prm_matrix
	use prm_var
	use prm_inv

	implicit none

	integer, intent(in) :: ks00, kt00, ndg, ls

	integer :: i, j, k, ja, jb, kss, ktt, is, it

	print *, ' USED MODEL FUNCTION/IUS= ', (ius(k,ls), k=1,4)
	write(77,*) ' USED MODEL FUNCTION/IUS= ', (ius(k,ls), k=1,4)

	kss = ks00 + ndg
	ktt = kt00 + ndg

	ja = 0
	jb = 0

	do k = 1, 4

		if(just(k,ls) == 0) cycle

		j = 0

		do it = 1, ktt
			do is = 1, kss

				j  = j  + 1
				ja = ja + 1

				if (is < ius(1,ls) .or. is > ius(2,ls)) cycle
				if (it < ius(3,ls) .or. it > ius(4,ls)) cycle

				jb = jb + 1

				do i = 1, iinit
					zz(i,jb+jsp) = zz(i,ja+jsp)
				end do

				idmp(j,k,ls) = jb + jsp

			end do
		end do

	end do

	jsp = jsp + jb

	print *, ' SET MODEL PARAMETERS. JB= ', jb, ' JSP= ', jsp
	write(77,*) ' SET MODEL PARAMETERS. JB= ', jb, ' JSP= ', jsp

return
end subroutine SETMP

!==========
subroutine SETBM(ls,ks00,kt00,ndg)

	use prm_matrix
	use prm_var
	use prm_inv

	implicit none

	integer, intent(in) :: ls, ks00, kt00, ndg

	integer :: k, jc, ix, iy, j0, linv, m1, m2
	integer :: iba, ibb, ibc, ibd, kss, ktt

	print *, ' IBOC(1/CONSTRAIN TO OUTSIDE,0/NO COSTRAIN)=', (iboc(k,ls), k=1,4)
	write(77,*) ' IBOC(1/CONSTRAIN TO OUTSIDE,0/NO COSTRAIN)=', (iboc(k,ls), k=1,4)


	iba = ius(1,ls)
	ibb = ius(2,ls)
	ibc = ius(3,ls)
	ibd = ius(4,ls)

	kss = ks00 + ndg
	ktt = kt00 + ndg

	if(iboc(1,ls) == 1) iba = iba - 1
	if(iboc(2,ls) == 1) ibb = ibb + 1
	if(iboc(3,ls) == 1) ibc = ibc - 1
	if(iboc(4,ls) == 1) ibd = ibd + 1

	jc = jcp

	do k = 1, 4
		if(just(k,ls) == 0) cycle

loopy:	do iy = 1, ktt
loopx:		do ix = 1, kss

				if(ix < ius(1,ls) .or. ix > ius(2,ls)) cycle loopx
				if(iy < ius(3,ls) .or. iy > ius(4,ls)) cycle loopy

				j0 = kss * (iy-1) + ix
				linv = 0

				if(ix-1 >= iba .and. ix+1 <= ibb) then

					if(ix /= 1) m1 = idmp(j0-1, k, ls)
					if(ix == 1) m1 = 0

					if(ix /= kss) m2 = idmp(j0+1, k, ls)
					if(ix == kss) m2 = 0

					if(m1 /= 0) cs(jc+1, m1) = 1.d0
					if(m2 /= 0) cs(jc+1, m2) = 1.d0

					linv = linv + 2

				end if

				if(iy-1 >= ibc .and. iy+1 <= ibd) then

					if(iy /= 1) m1 = idmp(j0-kss, k, ls)
					if(iy == 1) m1 = 0

					if(iy /= ktt) m2 = idmp(j0+kss, k, ls)
					if(iy == ktt) m2 = 0

					if(m1 /= 0) cs(jc+1, m1) = 1.d0
					if(m2 /= 0) cs(jc+1, m2) = 1.d0

					linv = linv + 2

				end if

				if(linv == 0) cycle loopx

				cs(jc+1, idmp(j0,k,ls)) = -dble(linv)
				jc = jc + 1

			end do loopx
		end do loopy

	end do

	jcp = jc

	print *, ' CONSTRAIN MATRIX JCP= ', jcp
	write(77,*) ' CONSTRAIN MATRIX JCP= ', jcp

return
end subroutine SETBM

!==========
subroutine INSL(rs,ls)

	use prm_matrix
	use prm_var
	use prm_inv

	implicit none

	integer, intent(in) :: ls
	real(8), intent(in) :: rs

	integer :: j, j1, j2
	real(8) :: rsa, rsb, spx, spy, sqx, sqy

!++  RS=RSS IS THE PLATE CONVERGENCE DIRECTEIN (SLIP DIRECTION) ++

	rsa = (rs - 45.d0) * pi / 180.d0
	rsb = (rs + 45.d0) * pi / 180.d0

	spx = dcos(rsa)
	spy = dsin(rsa)
	sqx = dcos(rsb)
	sqy = dsin(rsb)

	do j = 1, nbi

		j1 = idmp(j,1,ls)
		j2 = idmp(j,2,ls)
		if(j1 == 0 .or. j2 == 0) cycle

		ha(j1) = 0.d0
		ha(j2) = 0.d0

		gi(j1,j1) = spx
		gi(j1,j2) = sqx
		gi(j2,j1) = spy
		gi(j2,j2) = sqy

		inncon = inncon + 1
		incc(inncon) = j1

	end do

	do j = 1, nbi

		j1 = idmp(j,1,ls)
		j2 = idmp(j,2,ls)
		if(j1 == 0 .or. j2 == 0) cycle

		inncon = inncon + 1
		incc(inncon) = j2

	end do

	print '(" SLIP DIRECTION NNLS,  RSS= ",F8.2," INCON= ",I6)', rs, inncon
	write(77,'(" SLIP DIRECTION NNLS,  RSS= ",F8.2," INCON= ",I6)') rs, inncon

return
end subroutine INSL

!==========
subroutine INNG(rs,ls,ktyp)

	use prm_matrix
	use prm_var
	use prm_inv

	implicit none

	integer, intent(in) :: ls, ktyp
	real(8), intent(in) :: rs

	integer :: j, j1

	do j = 1, nbi

		j1 = idmp(j,ktyp,ls)
		if(j1 == 0) cycle

		ha(j1) = 0.d0
		gi(j1,j1) = rs
		inncon = inncon + 1
		incc(inncon) = j1

	end do

	print *, ' NON-NEGATIVE CONSTRAINT SOURCE= ', ktyp
	write(77,*) ' NON-NEGATIVE CONSTRAINT SOURCE= ', ktyp

return
end subroutine INNG
!==========

end module mod_subsetup

