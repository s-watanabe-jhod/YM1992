!=================================================
! SUBROUTINES
!	BIBSP
!	JACOBI(iu)
!	JACOBIFW(iu)
!
!=================================================

module mod_gdjac

	implicit none

contains

!==========
subroutine JACOBIFW(iu)

	use prm_matrix
	use prm_var
	use mod_jaclib
	use mod_trans

	implicit none

	integer, intent(in) :: iu
	integer :: ios
	real(8) :: xs, ys

	print *, ' HORIZONTAL DISPLACEMENTS OBS.POINTS= ', ih

	do ios = 1, ih
		CALL ROTATE(st(1,ios),st(2,ios),xs,ys,1,1, Aphi,Adlt)
		CALL HORDS(xs,ys)
		CALL DSITRNS
		write(iu) (Ds1(k,1), k=1,Nbi)
		write(iu) (Ds1(k,2), k=1,Nbi)
		write(iu) (Ds1(k,3), k=1,Nbi)
		write(iu) (Ds1(k,4), k=1,Nbi)
		write(iu) (Ds2(k,1), k=1,Nbi)
		write(iu) (Ds2(k,2), k=1,Nbi)
		write(iu) (Ds2(k,3), k=1,Nbi)
		write(iu) (Ds2(k,4), k=1,Nbi)
	end do

	print *, ' VERTICAL DISPLACEMENTS OBS.POINTS= ', iv

	do ios = ih+1, isn
		CALL ROTATE(st(1,ios),st(2,ios),xs,ys,1,1, Aphi,Adlt)
		CALL VERDS(xs,ys)
		write(iu) (Ds1(k,1), k=1,Nbi)
		write(iu) (Ds1(k,2), k=1,Nbi)
		write(iu) (Ds1(k,3), k=1,Nbi)
		write(iu) (Ds1(k,4), k=1,Nbi)
	end do

return
end subroutine JACOBIFW

!==========


end module mod_gdjac

