!=================================================
! SUBROUTINES
!	JACOBI(iu)
!
!=================================================

module mod_gdjac

	implicit none

contains
!==========
subroutine JACOBI(iu, ls)

	use prm_matrix
	use prm_inv
	use mod_jaclib
	use mod_trans
	
	implicit none
	
	integer, intent(in) :: iu, ls
	integer :: ios, i, j, k, ja, ild
	real(8) :: xs, ys, xx, yy
	
	ja = 0
	
	
	
	
	
	
	
	
	
	
	
	
	print *, ' HORIZONTAL DISPLACEMENTS OBS.POINTS= ', Ih
	
	do ios = 1, Ih
		
		CALL ROTATE(St(1,ios),St(2,ios),xs,ys,1,1, Bphi(ls),Bdlt(ls))
		
		CALL HORDS(xs,ys,Ndgs(ls),Jd2(ls),Da2(ls))
		
		!==== DSITRNS ====
		do i = 1, Nbi
			CALL ROTATE(Ds1(i,1),Ds2(i,1), xx,yy, 1,-1, Bphi(ls),Bdlt(ls))
				Ds1(i,1) = xx
				Ds2(i,1) = yy
			CALL ROTATE(Ds1(i,2),Ds2(i,2), xx,yy, 1,-1, Bphi(ls),Bdlt(ls))
				Ds1(i,2) = xx
				Ds2(i,2) = yy
			CALL ROTATE(Ds1(i,3),Ds2(i,3), xx,yy, 1,-1, Bphi(ls),Bdlt(ls))
				Ds1(i,3) = xx
				Ds2(i,3) = yy
			CALL ROTATE(Ds1(i,4),Ds2(i,4), xx,yy, 1,-1, Bphi(ls),Bdlt(ls))
				Ds1(i,4) = xx
				Ds2(i,4) = yy
		end do
		!==== /DSITRNS ====
		
		do ild = 1, 2
			do j = 1, 4
				
				if(Just(j,ls) == 0) cycle
				
				ja = Nbi * (Just(j,ls)-1) + Jsp
				
				do k = 1, Nbi
					if(ild == 1) Z_Green(2*(ios-1)+ild, ja+k) = Ds1(k,j)
					if(ild == 2) Z_Green(2*(ios-1)+ild, ja+k) = Ds2(k,j)
				end do
			end do
		end do
		
		!write(iu) (Ds1(k,1), k=1,Nbi)
		!write(iu) (Ds1(k,2), k=1,Nbi)
		!write(iu) (Ds1(k,3), k=1,Nbi)
		!write(iu) (Ds1(k,4), k=1,Nbi)
		!write(iu) (Ds2(k,1), k=1,Nbi)
		!write(iu) (Ds2(k,2), k=1,Nbi)
		!write(iu) (Ds2(k,3), k=1,Nbi)
		!write(iu) (Ds2(k,4), k=1,Nbi)
		
	end do
	
	print *, ' VERTICAL DISPLACEMENTS OBS.POINTS= ', Iv
	
	do ios = Ih+1, Isn
		
		CALL ROTATE(St(1,ios),St(2,ios),xs,ys,1,1, Bphi(ls),Bdlt(ls))
		
		CALL VERDS(xs,ys,Ndgs(ls),Jd2(ls),Da2(ls))
		
		do j = 1, 4
			
			if(Just(j,ls) == 0) cycle
			
			ja = Nbi * (Just(j,ls) - 1) + Jsp
			
			do k = 1, Nbi
				Z_Green(ios+Ih, ja+k) = Ds1(k,j)
			end do
		
		end do
		
		!write(iu) (Ds1(k,1), k=1,Nbi)
		!write(iu) (Ds1(k,2), k=1,Nbi)
		!write(iu) (Ds1(k,3), k=1,Nbi)
		!write(iu) (Ds1(k,4), k=1,Nbi)
		
	end do
	
	print *, ' AVERAGED DIRECTION OF SLIP ON THE FAULT'
	
	if(ja+Nbi > Kmp) stop ' INSUFFICIENT ARRAY DIM(KMP)'
	
	CALL SLIPM(Ndgs(ls),Jd2(ls))
	
	!==== SMITRNS ====
	do k = 1, Nbi
		CALL ROTATE(Ds1(k,1),Ds1(k,2), xx,yy, 1,-1, Bphi(ls),Bdlt(ls))
			Ds1(k,1) = xx
			Ds1(k,2) = yy
			Slms(1,k,ls) = Ds1(k,1)
			Slms(2,k,ls) = Ds1(k,2)
			Slms(3,k,ls) = Ds1(k,3)
		CALL ROTATE(Ds2(k,1),Ds2(k,2), xx,yy, 1,-1, Bphi(ls),Bdlt(ls))
			Ds2(k,1) = xx
			Ds2(k,2) = yy
			Slmd(1,k,ls) = Ds1(k,1)
			Slmd(2,k,ls) = Ds1(k,2)
			Slmd(3,k,ls) = Ds1(k,3)
	end do
	!==== /SMITRNS ====
	
	!write(iu) (Ds1(j,1),Ds1(j,2),Ds1(j,3), Ds2(j,1),Ds2(j,2),Ds2(j,3), j=1,Nbi)
	
return
end subroutine JACOBI

!==========

end module mod_gdjac

