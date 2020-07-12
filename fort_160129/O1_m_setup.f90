!=================================================
! SUBROUTINES
!	SETMP(KS00,KT00,NDG,LS)
!	SETBM(LS,KS00,KT00,NDG)
!
!=================================================

module mod_setup

	implicit none

contains

!==========
subroutine SETMP(ls, ks00, kt00, ndg)

	use prm_matrix
	use prm_inv

	implicit none

	integer, intent(in) :: ls, ks00, kt00, ndg

	integer :: i, j, k, ja, jb, kss, ktt, is, it

	print *,    ' USED MODEL FUNCTION/IUS= ', (Ius(k,ls), k=1,4)
	write(77,*) ' USED MODEL FUNCTION/IUS= ', (Ius(k,ls), k=1,4)

	kss = ks00 + ndg
	ktt = kt00 + ndg

	ja = 0
	jb = 0

	do k = 1, 4

		if(Just(k,ls) == 0) cycle

		j = 0

		do it = 1, ktt
			do is = 1, kss

				j  = j  + 1
				ja = ja + 1

				if (is < Ius(1,ls) .or. is > Ius(2,ls)) cycle
				if (it < Ius(3,ls) .or. it > Ius(4,ls)) cycle

				jb = jb + 1

				do i = 1, Ih*2 + Iv
					Z_Green(i,jb+Jsp) = Z_Green(i,ja+Jsp) ! Z_Green(Itotal,Jsp)
					! ここで，Z_Green のうち，Str-slip, Dip-slip, Open, Ext の不要部分を除いている
				end do

				Idmp(j,k,ls) = jb + Jsp ! Z_Green の index

			end do
		end do

	end do

	Jsp = Jsp + jb

	print *,    ' SET MODEL PARAMETERS. JB= ', jb, ' JSP= ', Jsp
	write(77,*) ' SET MODEL PARAMETERS. JB= ', jb, ' JSP= ', Jsp

return
end subroutine SETMP

!==========
subroutine SETBM(ls, ks00, kt00, ndg)

	use prm_matrix
	use prm_inv

	implicit none

	integer, intent(in) :: ls, ks00, kt00, ndg

	integer :: k, jc, ix, iy, j0, linv, m1, m2
	integer :: iba, ibb, ibc, ibd, kss, ktt

	print *,    ' IBOC(1/CONSTRAIN TO OUTSIDE,0/NO COSTRAIN)=', (Iboc(k,ls), k=1,4)
	write(77,*) ' IBOC(1/CONSTRAIN TO OUTSIDE,0/NO COSTRAIN)=', (Iboc(k,ls), k=1,4)

	iba = Ius(1,ls)
	ibb = Ius(2,ls)
	ibc = Ius(3,ls)
	ibd = Ius(4,ls)

	kss = ks00 + ndg
	ktt = kt00 + ndg

	if(Iboc(1,ls) == 1) iba = iba - 1
	if(Iboc(2,ls) == 1) ibb = ibb + 1
	if(Iboc(3,ls) == 1) ibc = ibc - 1
	if(Iboc(4,ls) == 1) ibd = ibd + 1

	jc = Jcp

	do k = 1, 4
		if(Just(k,ls) == 0) cycle

loopy:	do iy = 1, ktt
loopx:		do ix = 1, kss

				if(ix < Ius(1,ls) .or. ix > Ius(2,ls)) cycle loopx
				if(iy < Ius(3,ls) .or. iy > Ius(4,ls)) cycle loopy

				j0 = kss * (iy-1) + ix
				linv = 0

				if(ix-1 >= iba .and. ix+1 <= ibb) then

					if(ix /= 1  ) m1 = Idmp(j0-1, k, ls)
					if(ix == 1  ) m1 = 0

					if(ix /= kss) m2 = Idmp(j0+1, k, ls)
					if(ix == kss) m2 = 0

					if(m1 /= 0  ) Cs(jc+1, m1) = 1.d0
					if(m2 /= 0  ) Cs(jc+1, m2) = 1.d0

					linv = linv + 2

				end if

				if(iy-1 >= ibc .and. iy+1 <= ibd) then

					if(iy /= 1  ) m1 = Idmp(j0-kss, k, ls)
					if(iy == 1  ) m1 = 0

					if(iy /= ktt) m2 = Idmp(j0+kss, k, ls)
					if(iy == ktt) m2 = 0

					if(m1 /= 0  ) Cs(jc+1, m1) = 1.d0
					if(m2 /= 0  ) Cs(jc+1, m2) = 1.d0

					linv = linv + 2

				end if

				if(linv == 0) print *, "linv=0 ", ix, iy
				if(linv == 0) linv = 0
				!if(linv == 0) cycle loopx
				! オリジナルから変更して，4隅のスムージングマトリクスにゼロが入るようにした。

				Cs(jc+1, Idmp(j0,k,ls)) = -dble(linv) ! smoothness matrix Cs(Jcp,Jsp)
				jc = jc + 1

			end do loopx
		end do loopy

	end do

	Jcp = jc

	print *,    ' CONSTRAIN MATRIX JCP= ', Jcp
	write(77,*) ' CONSTRAIN MATRIX JCP= ', Jcp

return
end subroutine SETBM

!==========

end module mod_setup
