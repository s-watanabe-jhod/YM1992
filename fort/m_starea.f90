!=================================================
! SUBROUTINES
!	STAREA
!
!=================================================

module mod_starea

	implicit none

contains

!==========
subroutine STAREA

	use prm_var
	use mod_trans

	implicit none

	real(8) :: lc_x, lc_y

	if(adlt /= 0.d0) then
		print *, " SORRY. ADLT MUST BE ZERO FOR CURVED SURFACE"
		stop
	end if

	CALL PLTXY(alat,alon,lc_x,lc_y,0)
	CALL TRANS(lc_x,lc_y,xa,ya,1,1)

	xb = xa + alen
	yb = ya + awid

	print 600, xa, xb, ya, yb
	 600 format(' ',' SET 2-D SPACE FOR FITTING'/' XMIN-XMAX ',2X,1F10.4,2X,F10.4/' YMIN-YMAX',2X,F10.4,2X,F10.4)

return
end subroutine STAREA

end module mod_starea
