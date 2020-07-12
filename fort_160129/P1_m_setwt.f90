!=================================================
! SUBROUTINES
!	SETWT
!
!=================================================

module mod_setwt

	implicit none

contains
!=========
subroutine SETWT

	use prm_matrix
	use prm_inv
	
	implicit none
	
	integer :: i, k, l
	real(8) :: wk(Kio)
	
	! 最初のループは Sta_corr = I のとき無意味（≒通常無意味）
	! Sta_corr は観測点間の相関
	
	do k = 1, Jsp
		
		do i = 1, Itotal
			
			wk(i) = 0.d0
			
			do l = 1, Ih*2 + Iv
				wk(i) = wk(i) + Sta_corr(i,l) * Z_Green(l,k)
			end do
			
		end do
		
		do i = 1, Itotal
			Z_Green(i,k) = wk(i)
		end do
		
	end do
	
	! Y_data と Z_Green のいずれもを E_prior で重み付けする
	! Y_data については メインプログラムで元の量に戻している。
	
	do i = 1, Itotal
		
		Y_data(i) = Y_data(i) / dabs(E_prior(i))
		
		do k = 1, Jsp
			Z_Green(i,k) = Z_Green(i,k) / dabs(E_prior(i))
		end do
		
	end do
	
	print *, " SET WEIGHT OK"
	
return
end subroutine SETWT

!==========

end module mod_setwt


