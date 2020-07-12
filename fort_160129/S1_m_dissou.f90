!=================================================
! SUBROUTINES
!	DISSOU(hd,nsource)
!
!=================================================

module mod_dissou

	implicit none

contains
!==========
subroutine DISSOU(hd,nsource)

	use prm_matrix
	use prm_inv
	use mod_blms
	use mod_source
	
	implicit none
	
	character(len=24), intent(in) :: hd
	integer, intent(in) :: nsource
	
	character(len=24) :: fxyz = "                        ", fls
	character(len=24) :: felp = "                        "
	character(len=24) :: fcnt = "                        "
	integer :: lsource, ls, lname, i, kmri, kui, kvi
	
	do lsource = 1, nsource
	
		ls = lsource
		
		if(Mlg(ls) == 0) cycle
		
		print *, ' -------------------------------------------'
		print *, ' SOURCE DISTRIBUTION.LSOURCE= ', ls
		
		fls = fscd(ls)
		lname = len_trim(fls)
		fxyz(1:lname-4) = fls(1:lname-4)
		felp(1:lname-4) = fls(1:lname-4)
		fcnt(1:lname-4) = fls(1:lname-4)
		fxyz(lname-3:lname+4) = "_scd.xyz"
		felp(lname-3:lname+4) = "_scd.elp"
		fcnt(lname-3:lname+4) = "_scd.cnt"
		
		open(ls+44,file=fxyz)
		open(ls+50,file=felp)
		open(ls+56,file=fcnt)
		open(ls+22,file=fls) !, status='unknown')
		
		if(Ipl(ls) == 0) open(31,file=Fcof(ls),status='old')
		
		!==== BLMS 断層面の基底関数の設定 ====
		
		if(Ipl(ls) == 0) then ! 断層面が曲面
			
			CALL STAREA(ls, Xa(ls), Ya(ls), Xb(ls), Yb(ls))     ! set Xa(ls), Ya(ls), Xb(ls), Yb(ls)
			CALL STMSAR(ls, Xfa(ls), Xfb(ls), Yfa(ls), Yfb(ls)) ! set Xfa(ls), Xfb(ls), Yfa(ls), Yfb(ls)
			
			open(31, file=Fcof(ls), status='old')
				Xcof = 0.d0
				read(31,'(//1X,3(I5,1X)//)') kmri, kui, kvi     ! 断層形状の読み込み from COF file
				read(31,'(5(E12.6,1X))') (Xcof(i), i = 1, kmri) ! 
			close(31)
			
			CALL MESH(ls) ! set Xtr(Ms), Ytr(Ms,Mt), Ztr(Ms,Mt), Da2(ls)
			
		else if(Ipl(ls) /= 0) then ! 断層面が平面
			STOP "!!!Not in service for FLAT FAULT!!!"
			! CALL STMSAR2
			! CALL MESHSS
		end if
		
		!==== /BLMS/ ====
		
		CALL SOURCE(hd,ls,ls+22,ls+44,ls+50,ls+56)
		
		close(ls+22)
		close(ls+44)
		close(ls+50)
		close(ls+56)
		
		if(Ipl(ls) == 0) close(31)
		
	end do

return
end subroutine DISSOU

!==========

end module mod_dissou


