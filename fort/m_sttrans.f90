!=================================================
! SUBROUTINES
!	STTRANS
!
!=================================================

module mod_sttrans

	implicit none

contains
!==========
subroutine STTRANS

	use prm_var

	implicit none

	real(8) :: rh, rd

	print 600, aphi, adlt
	 600 FORMAT(' ','SOURCE GEOMETRY (PHI,DELT)=',2(F6.2,1X))

	rh = aphi * pi/180.d0
	rd = adlt * pi/180.d0

	t11 = dcos(rh)
	t12 = dsin(rh)
	u11 = dcos(rd)
	u12 = dsin(rd)

return
end subroutine STTRANS

end module mod_sttrans
