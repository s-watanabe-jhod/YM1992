!=================================================
! SUBROUTINES
!	INP2DP(iu)
!		intent(in) :: iu
!
!=================================================

module mod_inp2dp

	implicit none

contains
!==========
subroutine INP2DP(iu)

	use prm_matrix
	use prm_var

	implicit none

	integer, intent(in) :: iu
	integer :: i, kmri, kui, kvi

	read(iu,'(//1X,3(I5,1X)//)') kmri, kui, kvi
	write(*,*)kmri,kui,kvi
	read(iu,'(5(E12.6,1X))') (xi(i), i = 1, kmri)

return
end subroutine INP2DP

end module mod_inp2dp

