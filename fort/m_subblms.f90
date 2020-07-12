!=================================================
! SUBROUTINES
!	STMSAR
!	MESH
!	STMSAR2
!	MESHSS
!
!
!=================================================

module mod_subblms

	implicit none

contains
!==========
subroutine STMSAR

	use prm_var
	use mod_trans

	implicit none

	real(8) :: lc_x, lc_y

	print *, ' SOURCE AREA FOR CURVED FAULT'
	print *, ' LAT,LON,DEPTH,WIDTH,LENGTH '
	print '(1X,5(F10.3,1X))', alaf, alof, adep, awif, alef

	CALL PLTXY(alaf,alof,lc_x,lc_y,0) ! (alaf,alof) => (lc_x,lc_y)
	CALL TRANS(lc_x,lc_y,xfa,yfa,1,1)

	xfb = xfa + alef
	yfb = yfa + awif

return
end subroutine STMSAR

!==========
subroutine MESH

	use prm_matrix
	use prm_var
	use mod_intdep
	use mod_mnewt

	implicit none

	real(8), parameter :: zero = 0.d0

	integer :: la, lb, is, it, iss
	real(8) :: da, db, tis
	real(8) :: x1, x2, x3, y0, y1, y2, y3, z0, z1, z2, z3
	real(8) :: cy, cz, dy, dz, c2, ar, arm
	real(8) :: xtq(0:ms), ytq(0:ms,2), ztq(0:ms,2)


	da = (xfb - xfa) / dble(ms2)
	da2 = da * da
	db = 2.d0 * da2

	la = 1

	do is = 0, ms2

		tis = dble(is)
		xtq(is) = xfa + da * tis
		ytq(is,la) = yfa

		CALL INTDEP(xtq(is),ytq(is,la),z1)

		ztq(is,la) = z1

	end do

	do is = 1, ms2
		xtr(is) = xtq(is) - da/2.d0
	end do

	do it = 1, mt2

		lb = 3 - la

		x1 = xtq(0)
		y1 = ytq(0,la)
		z1 = ztq(0,la)

		x2 = xtq(1)
		y2 = ytq(1,la)
		z2 = ztq(1,la)

		x3 = x1

		CALL MNEWT(x1,y1,z1, x2,y2,z2, x3,y3,z3, da,da2)

		ytq(0,lb) = y3
		ztq(0,lb) = z3

		do is = 1, ms2

			iss = is - 1

			y0 = ytq(iss,la)
			z0 = ztq(iss,la)

			x1 = xtq(is)
			y1 = ytq(is,la)
			z1 = ztq(is,la)

			x2 = xtq(iss)
			y2 = ytq(iss,lb)
			z2 = ztq(iss,lb)

			cy = y1 - y0
			cz = z1 - z0
			dy = y2 - y0
			dz = z2 - z0

			c2 = da2 + cy*cy + cz*cz

			CALL AREA2(da,cy,cz,zero,dy,dz,c2,ar)

			arm = db - ar
			x3 = x1

			CALL MNEWT(x1,y1,z1, x2,y2,z2, x3,y3,z3, dy,arm)

			ytq(is,lb) = y3
			ztq(is,lb) = z3

		end do

		do is = 1, ms2
			iss = is - 1
			ytr(is,it) = (ytq(is,la) + ytq(iss,la) + ytq(is,lb) + ytq(iss,lb) ) / 4.d0
			ztr(is,it) = (ztq(is,la) + ztq(iss,la) + ztq(is,lb) + ztq(iss,lb) ) / 4.d0
		end do

		la = lb

	end do

return
end subroutine MESH

!==========
subroutine STMSAR2

	use prm_var
	use mod_trans

	implicit none

	real(8) :: lc_x, lc_y, y2, y3, z2

	print *, ' SOURCE AREA FOR PLANE FAULT'
	print *, ' LAT,LON,DEPTH,WIDTH,LENGTH '
	print '(1X,5(F10.3,1X))', alaf, alof, adep, awif, alef

	CALL PLTXY(alaf,alof,lc_x,lc_y,0) ! (alaf,alof) => (lc_x,lc_y)
	CALL TRANS(lc_x,lc_y,xfa,yfa,1,1)

	zda = -adep

	CALL TRANS(yfa,zda,y2,z2,2,1)

	y3 = y2 + awif

	CALL TRANS(y3,z2,yfb,zdb,2,-1)

	xfb = xfa + alef

return
end subroutine STMSAR2

!==========
subroutine MESHSS

	use prm_matrix
	use prm_var

	implicit none

	integer :: is, it
	real(8) :: dx,dy,dz, rx,ry,rz, da1, tis, tit

	dx = xfb - xfa
	dy = yfb - yfa
	dz = zdb - zda

	rx = dx / dble(ms2)
	ry = dy / dble(mt2)
	rz = dz / dble(mt2)

	do is = 1, ms2
		tis = dble(is)
		xtr(is) = (tis - 0.5d0) * rx + xfa
	end do

	do it = 1, mt2
		tit = dble(it)
		ytr(1,it) = (tit -0.5d0) * ry + yfa
		ztr(1,it) = (tit -0.5d0) * rz + zda

		do is = 2, ms2
			ytr(is,it) = ytr(1,it)
			ztr(is,it) = ztr(1,it)
		end do

	end do

	da1 = rx * dsqrt(ry*ry + rz*rz)
	da2 = dabs(da1)

return
end subroutine MESHSS
!===========

end module mod_subblms

