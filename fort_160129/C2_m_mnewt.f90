!=================================================
! SUBROUTINE
!	MNEWT(x1, y1, z1, x2, y2, z2, x3, y3, z3, di, da2l)
!
!=================================================
module mod_mnewt

	implicit none

contains
!==========
subroutine MNEWT(ls, x1,y1,z1, x2,y2,z2, x3,y3,z3, di,da2l)
	
	use prm_matrix
	use prm_inv
	use mod_wbb
	
	implicit none
	
	integer, intent(in)  :: ls
	real(8), intent(in)  :: x1, y1, z1, x2, y2, z2, x3, di, da2l
	real(8), intent(out) :: y3, z3
	
	real(8) :: ct, cx, cy, cz, c2, d2, cd
	real(8) :: dx, dl, dr, dz, du, dy
	real(8) :: arr, arl, ari
	
	integer :: ix, iy, idb(16), j
	real(8) :: px, py, bsv(16)
	
	ct = 1.d-6
	cx = x2 - x1
	cy = y2 - y1
	cz = z2 - z1
	
	dx  = x3 - x1
	dl  = di
	dr  = 0.d0
	arr = 0.d0
	
	do
	
		y3 = y1 + dl
		
		! == INTDEP(x3,y3,z3) ==
		CALL WHERE1(x3, ix, px, Xa(ls), Xb(ls), Ju0(ls))
		CALL WHERE1(y3, iy, py, Ya(ls), Yb(ls), Jv0(ls))
		
		CALL BSVALI(px, py, bsv, Ndgf(ls), Nod2(ls))
		CALL BSPARI(ix, iy, idb, Ndgf(ls), Ju1(ls), Jv1(ls))
		
		z3 = 0.d0
		do j = 1, Nod2(ls)
			z3 = z3 + bsv(j) * Xcof(idb(j))
		end do
		! ====
		
		dz = z3 - z1
		
		! 二つのベクトルの面積計算 (cx,cy,cz) と (dx,dy,dz)
		! AREA2(cx, cy, cz, dx, dl, dz, c2, arl)
		dy = dl
		c2  = cx*cx + cy*cy + cz*cz
		d2  = dx*dx + dy*dy + dz*dz
		cd  = cx*dx + cy*dy + cz*dz
		arl = dsqrt(d2*c2 - cd*cd)
		! ====
		
		if(dabs(arl-da2l) < ct) return
		if(arl < da2l) then
			dr = dl
			arr = arl
			dl = dl * da2l / arl
			cycle
		end if

		do
			du = dr + (dl - dr) * (da2l - arr) / (arl - arr)
			y3 = y1 + du
			
			! == INTDEP(x3,y3,z3) ==
			CALL WHERE1(x3, ix, px, Xa(ls), Xb(ls), Ju0(ls))
			CALL WHERE1(y3, iy, py, Ya(ls), Yb(ls), Jv0(ls))
			
			CALL BSVALI(px, py, bsv, Ndgf(ls), Nod2(ls))
			CALL BSPARI(ix, iy, idb, Ndgf(ls), Ju1(ls), Jv1(ls))
			
			z3 = 0.d0
			do j = 1, Nod2(ls)
				z3 = z3 + bsv(j) * Xcof(idb(j))
			end do
			! ====
			
			dz = z3 - z1
			
			! 二つのベクトルの面積計算 (cx,cy,cz) と (dx,dy,dz)
			! AREA2(cx, cy, cz, dx, du, dz, c2, ari)
			dy = du
			c2  = cx*cx + cy*cy + cz*cz
			d2  = dx*dx + dy*dy + dz*dz
			cd  = cx*dx + cy*dy + cz*dz
			ari = dsqrt(d2*c2 - cd*cd)
			! ====
			
			if(dabs(ari-da2l) < ct) return
			
			if(ari > da2l) then
				arl = ari
				dl = du
			else if(ari < da2l) then
				arr = ari
				dr = du
			end if

		end do

	end do

end subroutine MNEWT

end module mod_mnewt

