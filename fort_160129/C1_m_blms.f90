!=================================================
! SUBROUTINES
!	STAREA(ls, lc_xa, lc_ya, lc_xb, lc_yb)
!	STMSAR(ls, lc_xfa, lc_yfa, lc_xfb, lc_yfb)
!	MESH
!
!=================================================

module mod_blms

	implicit none

contains

!==========
subroutine STAREA(ls, lc_xa, lc_ya, lc_xb, lc_yb)

	use prm_inv
	use mod_trans

	implicit none

	integer, intent(in)  :: ls
	real(8), intent(out) :: lc_xa, lc_ya, lc_xb, lc_yb
	real(8) :: lc_x, lc_y

	if(Bdlt(ls) /= 0.d0) stop " SORRY. ADLT MUST BE ZERO FOR CURVED SURFACE"

	CALL LL2XY(Blat(ls),Blon(ls), lc_x,lc_y, 0, Alat0,Alng0,Icord)
	CALL ROTATE(lc_x,lc_y, lc_xa, lc_ya, 1,1, Bphi(ls),Bdlt(ls))

	lc_xb = lc_xa + Blen(ls)
	lc_yb = lc_ya + Bwid(ls)

	print 600, lc_xa, lc_xb, lc_ya, lc_yb
	 600 format(' ',' SET 2-D SPACE FOR FITTING'/' XMIN-XMAX ',2X,1F10.4,2X,F10.4/' YMIN-YMAX',2X,F10.4,2X,F10.4)

return
end subroutine STAREA

!==========
subroutine STMSAR(ls, lc_xfa, lc_xfb, lc_yfa, lc_yfb)

	use prm_inv
	use mod_trans

	implicit none

	integer, intent(in)  :: ls
	real(8), intent(out) :: lc_xfa, lc_yfa, lc_xfb, lc_yfb
	real(8) :: lc_x, lc_y

	print *, ' SOURCE AREA FOR CURVED FAULT'
	print *, ' LAT,LON,DEPTH,WIDTH,LENGTH '
	print '(1X,5(F10.3,1X))', Blaf(ls),Blof(ls), Bdep(ls),Bwif(ls),Blef(ls)

	CALL LL2XY(Blaf(ls),Blof(ls), lc_x,lc_y, 0, Alat0,Alng0,Icord) ! (Blaf(ls),Blof(ls)) => (lc_x,lc_y)
	CALL ROTATE(lc_x,lc_y, lc_xfa,lc_yfa, 1,1, Bphi(ls),Bdlt(ls))

	lc_xfb = lc_xfa + Blef(ls)
	lc_yfb = lc_yfa + Bwif(ls)

return
end subroutine STMSAR

!==========
subroutine MESH(ls)

	use prm_matrix
	use prm_inv
	use mod_wbb
	use mod_mnewt

	implicit none

	integer, intent(in) :: ls

	integer :: la, lb, is, it, iss, j
	real(8) :: da
	real(8) :: x1, x2, x3, y0, y1, y2, y3, z0, z1, z2, z3
	real(8) :: cx, cy, cz, dx, dy, dz, c2, d2, cd, ar, arm
	real(8) :: xtq(0:ms), ytq(0:ms,2), ztq(0:ms,2)

	integer :: ix, iy, idb(16)
	real(8) :: px, py, bsv(16)

	da  = (Xfb(ls) - Xfa(ls)) / dble(Js0(ls)*Jd2(ls))
	Da2(ls) = da * da

	la = 1

	do is = 0, Js0(ls)*Jd2(ls)

		xtq(is)    = Xfa(ls) + da * dble(is)
		ytq(is,la) = Yfa(ls)

		! == INTDEP(xtq(is),ytq(is,la),z1) ==
		CALL WHERE1(xtq(is),    ix, px, Xa(ls), Xb(ls), Ju0(ls))
		CALL WHERE1(ytq(is,la), iy, py, Ya(ls), Yb(ls), Jv0(ls))


		CALL BSVALI(px, py, bsv, Ndgf(ls), Nod2(ls))
		CALL BSPARI(ix, iy, idb, Ndgf(ls), Ju1(ls), Jv1(ls))

		z1 = 0.d0
		do j = 1, Nod2(ls)
			z1 = z1 + bsv(j) * Xcof(idb(j))
		end do
		! ====

		ztq(is,la) = z1



	end do

	do is = 1, Js0(ls)*Jd2(ls)
		Xtr(is) = xtq(is) - da/2.d0
	end do

	do it = 1, Jt0(ls)*Jd2(ls)

		lb = 3 - la

		x1 = xtq(0)
		y1 = ytq(0,la)
		z1 = ztq(0,la)

		x2 = xtq(1)
		y2 = ytq(1,la)
		z2 = ztq(1,la)

		x3 = x1

		CALL MNEWT(ls, x1,y1,z1, x2,y2,z2, x3,y3,z3, da, Da2(ls)) ! output y3, z3

		ytq(0,lb) = y3
		ztq(0,lb) = z3

		do is = 1, Js0(ls)*Jd2(ls)

			iss = is - 1

			y0 = ytq(iss,la)
			z0 = ztq(iss,la)

			x1 = xtq(is)
			y1 = ytq(is,la)
			z1 = ztq(is,la)

			x2 = xtq(iss)
			y2 = ytq(iss,lb)
			z2 = ztq(iss,lb)

			cx = da
			cy = y1 - y0
			cz = z1 - z0

			dx = 0.d0
			dy = y2 - y0
			dz = z2 - z0


			! 二つのベクトルの面積計算 (cx,cy,cz) と (dx,dy,dz)
			! AREA2(da,cy,cz, 0.d0,dy,dz, c2,ar)
			c2 = cx*cx + cy*cy + cz*cz
			d2 = dx*dx + dy*dy + dz*dz
			cd = cx*dx + cy*dy + cz*dz
			ar = dsqrt(d2*c2 - cd*cd)
			! ====

			arm = 2.d0*Da2(ls) - ar
			x3 = x1

			CALL MNEWT(ls, x1,y1,z1, x2,y2,z2, x3,y3,z3, dy, arm) ! output y3, z3

			ytq(is,lb) = y3
			ztq(is,lb) = z3

		end do

		do is = 1, Js0(ls)*Jd2(ls)
			iss = is - 1
			Ytr(is,it) = (ytq(is,la) + ytq(iss,la) + ytq(is,lb) + ytq(iss,lb) ) / 4.d0
			Ztr(is,it) = (ztq(is,la) + ztq(iss,la) + ztq(is,lb) + ztq(iss,lb) ) / 4.d0
		end do

		la = lb

	end do

return
end subroutine MESH

!==========

end module mod_blms
