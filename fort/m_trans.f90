!=================================================
! SUBROUTINES
!	PLTXY(alat,along,x,y,ind)
!	TRANS(u,v,ut,vt,iab,ind)
!	
!=================================================

module mod_trans

	implicit none

contains
!==========
subroutine PLTXY(alati,along,lc_x,lc_y,ind)

	use prm_var

	implicit none
!		PLTXY TRANSFORMS (X,Y) TO (ALATI,ALONG) IF IND.EQ.1
!		PLTXY TRANSFORMS (ALATI,ALONG) TO (X,Y) IF IND.EQ.0
!			WHEN ICORD.NE.0  PLTXY MAKES NO CHANGE IN TRANSFORMATION  BETWEEN
!				(X,Y) AND (ALATI,ALONG).
	real(8), intent(inout) :: alati, along, lc_x, lc_y
	integer, intent(in) :: ind

	real(8) :: a, e2, e12, d, rd
	real(8) :: rlat0, slat0, clat0, rlat, slat, clat
	real(8) :: r, an, al, v2, c1, c2, ph1, rph1, rph2, tphi1, cphi1, sphi1, sphi2

! IAU-64 (理科年表2001年に掲載) の場合?
!	a   = 6.378160d3	! 地球の赤道半径
!	e2  = 6.6944541d-3	! 離心率の2乗
!	! 扁平率 1/f = 298.25 とした場合の e^2 = 0.00669454185 のタイプミス？
!	e12 = 6.7395719d-3	! 

! WGS-84 の場合
	a   = 6.37813700000000d3	! 地球の赤道半径 km
	e2  = 6.69438002301275d-3	! 離心率の2乗
	e12 = e2 / (1.d0 - e2)

	d   = 180.d0/pi		! radian <--> degree 変換係数
	rd  = 1.d0 / d		! 

	if(ind == 1) then

		if(icord /= 0) then
			alati = lc_x
			along = lc_y
			return
		end if

		rlat0 = rd * alat0
		slat0 = dsin(rlat0)
		clat0 = dcos(rlat0)

		r = a * (1.d0 - e2)/dsqrt((1.d0 - e2*slat0**2)**3)
		an = a / dsqrt(1.d0 - e2 * slat0**2)
!		an = a / dsqrt(1.d0 - e12 * clat0**2)

		v2 = 1.d0 + e12 * clat0**2
		c1 = d/r
		c2 = d/an

		ph1 = alat0 + c1 * lc_y
		rph1 = ph1 * rd
		tphi1 = dtan(rph1)
		cphi1 = dcos(rph1)

		alati = ph1 - (c2*lc_x)**2 * v2 * tphi1 / (2.d0*d)
		along = alng0 + (c2*lc_x) / cphi1 - (c2*lc_x)**3 * (1.d0 + 2.d0 * tphi1**2) / (6.d0 * d**12 * cphi1)

		return

	else
		if(icord /= 0) then
			lc_x = alati
			lc_y = along
			return
		end if

		rlat = rd * alati
		slat = dsin(rlat)
		clat = dcos(rlat)

		v2 = 1.d0 + e12 * clat**2
		al = along - alng0

		ph1 = alati + (v2 * al**2 * slat * clat) / (2.d0 * d)
		rph1 = ph1 * rd
		rph2 = (ph1 + alat0) * 0.5d0 * rd
		sphi1 = dsin(rph1)
		sphi2 = dsin(rph2)

		r = a * (1.d0 - e2) / dsqrt((1.d0 - e2 * sphi2**2)**3)
		an = a / dsqrt(1.d0 - e2 * sphi1**2)
		c1 = d/r
		c2 = d/an

		lc_y = (ph1 - alat0) / c1
		lc_x = (al * clat) / c2 + (al**3 * clat * dcos(2.d0*rlat)) / (6.d0 * c2 * d**2)

		return

	end if

end subroutine PLTXY

!==========
subroutine TRANS(u,v,ut,vt,iab,ind)

	use prm_var

	implicit none
!--------------------------------------------------------
! (IAB,IND)
! (1, 1):  X-Y COODINATE rotation APHI deg counter-clockwisely.
!             Transform from LOCAL COOD. to FAULT COOD.
! (1,-1):  X-Y COODINATE rotation APHI deg clockwisely.
!             Transform from FAULT COOD. to LOCAL COOD.
! (2, 1):  Y-Z COODINATE rotation ADLT deg counter-clockwisely.
!             Transform from FAULT Cood. to ModelSurface COOD.
! (2,-1):  Y-Z COODINATE rotation ADLT deg clockwisely.
!             Transform from ModelSurface COOD. to FAULT COOD>
!--------------------------------------------------------
	real(8), intent(inout) :: u, v, ut, vt
	integer, intent(in) :: iab, ind

	real(8) :: p11, p12

	p11 = 0.d0
	p12 = 0.d0

	if(iab == 1) then
		p11 = t11
		p12 = t12
	else if(iab == 2) then
		p11 = u11
		p12 = u12
	else
		print *, " in TRANS, [iab] is bad value. "
		stop
	end if

	if(ind == 1) then
		ut =  p11 * u + p12 * v
		vt = -p12 * u + p11 * v
	else if(ind == -1) then
		ut =  p11 * u - p12 * v
		vt =  p12 * u + p11 * v
	end if

return
end subroutine TRANS


end module mod_trans

