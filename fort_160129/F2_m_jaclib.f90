!=================================================
! SUBROUTINES
!	HORDS(xs,ys)
!	DSITRNS
!	VERDS(xs,ys)
!	SLIPM
!	SMITRNS
!
!=================================================

module mod_jaclib

	implicit none

contains
!==========
subroutine HORDS(xs,ys,nbdeg, nd2, da2_lc)

	use prm_matrix
	use prm_inv
	use mod_mkkn
	
	implicit none
	
	real(8), intent(in) :: xs, ys, da2_lc
	integer, intent(in) :: nbdeg, nd2
	
	integer :: is, it, m, mu, mv, m1, m3, ls1, lt1, j, ll, mm
	integer :: nborder2, ndd2, ms2l, mt2l
	real(8) :: u1, u2, u3, v1, v2, v3, w1, w2, w3
	real(8) :: uu, vv, ww, wl, delw
	real(8) :: q1, q2, q3, q4, q5, q6, q7, q8, q9
	real(8) :: q10, q11, q12, q13, q14, q15, q16, q17, q18
	real(8) :: dx, dy, dz, r2, rr, ggx, ggy, ggz, xg, xg2, f
	real(8) :: ggx2, ggy2, ggz2, ggxy, dex2, c1, c2, c4, cnorm
	real(8) :: p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11
	real(8) :: fde
	
	real(8), parameter :: pi = 3.14159265358979324d0
	
	gf1(1:Kfe,1:16,1:4) = 0.d0
	gf2(1:Kfe,1:16,1:4) = 0.d0
	
	ndd2 = nd2 * nd2
	ms2l = nd2 * Ks0
	mt2l = nd2 * Kt0
	
	do it =1, kt0
		
		lt1 = (it-1) * nd2
		
		do is = 1, ks0
			
			ls1 = (is-1) * nd2
			
			do m = 1, ndd2
				
				m1 = (m - int((m-1) / nd2) * nd2) + ls1
				m3 = (1 + int((m-1) / nd2)      ) + lt1
				
				mu = m1 - int((m1-1)/(ms2l-1))
				mv = m3 - int((m3-1)/(mt2l-1))
				
				u1 = xtr(mu+1) - xtr(mu)
				u2 = ytr(mu+1,mv) - ytr(mu,mv)
				u3 = ztr(mu+1,mv) - ztr(mu,mv)
				
				v1 = 0.d0
				v2 = ytr(mu,mv+1) - ytr(mu,mv)
				v3 = ztr(mu,mv+1) - ztr(mu,mv)
				
				w1 = u2 * v3 - u3 * v2
				w2 =         - u1 * v3
				w3 = u1 * v2
				
				v1 = w2 * u3 - w3 * u2
				v2 = w3 * u1 - w1 * u3
				v3 = w1 * u2 - w2 * u1
				
				uu = dsqrt(u1*u1 + u2*u2 + u3*u3)
				vv = dsqrt(v1*v1 + v2*v2 + v3*v3)
				wl = w1*w1 + w2*w2 + w3*w3
				ww = dsqrt(wl)
				delw = Delt * wl
				
				q1  = 2.d0 * u1 * w1
				q2  = 2.d0 * u2 * w2
				q3  = 2.d0 * u3 * w3
				q4  = u1 * w2 + u2 * w1
				q5  = u2 * w3 + u3 * w2
				q6  = u3 * w1 + u1 * w3
				
				q7  = 2.d0 * v1 * w1
				q8  = 2.d0 * v2 * w2
				q9  = 2.d0 * v3 * w3
				q10 = v1 * w2 + v2 * w1
				q11 = v2 * w3 + v3 * w2
				q12 = v3 * w1 + v1 * w3
				
				q13 = delw + 2.d0 * w1 * w1
				q14 = delw + 2.d0 * w2 * w2
				q15 = delw + 2.d0 * w3 * w3
				q16 = 2.d0 * w1 * w2
				q17 = 2.d0 * w2 * w3
				q18 = 2.d0 * w3 * w1
				
				dx = xs - xtr(m1)
				dy = ys - ytr(m1,m3)
				dz = -ztr(m1,m3)
				
				r2 = dx*dx + dy*dy + dz*dz
				rr = dsqrt(r2)
				
				ggx = dx / rr
				ggy = dy / rr
				ggz = dz / rr
				
				xg = ggz + 1.d0
				xg2 = xg * xg
				
				fde = 1.d0/(1.d0+Delt)
				f = 3.d0 - fde * (xg + 2.d0) / xg2 / xg
				ggx2 = ggx * ggx
				ggy2 = ggy * ggy
				ggz2 = ggz * ggz
				ggxy = ggx * ggy
				
				dex2 = fde / xg2
				c1 = f * ggx2 - 1.d0
				c2 = f * ggy2 - 1.d0
				c4 = fde - 1.d0 + 3.d0 * ggz2
				
				p1  = ggx * (c1 + 3.d0*dex2)
				p2  = ggy * (c1 + dex2)
				p3  = ggx * (c2 + dex2)
				p4  = ggy * (c2 + 3.d0*dex2)
				
				p5  = ggx * c4
				p6  = ggy * c4
				p7  = ggy * (f*ggx2 + dex2) * 2.d0
				p8  = ggx * (f*ggy2 + dex2) * 2.d0
				
				p9  = 6.d0 * ggxy * ggz
				p10 = 6.d0 * ggy2 * ggz
				p11 = 6.d0 * ggx2 * ggz
				
				gs1(m,1) = (p1*q1  + p3*q2  + p5*q3  + p7*q4  + p9 *q5  + p11*q6 ) /r2/ww/uu
				gs1(m,2) = (p1*q7  + p3*q8  + p5*q9  + p7*q10 + p9 *q11 + p11*q12) /r2/ww/vv
				gs1(m,3) = (p1*q13 + p3*q14 + p5*q15 + p7*q16 + p9 *q17 + p11*q18) /r2/wl
				gs1(m,4) = (p1 + p3 + p5) / r2
				
				gs2(m,1) = (p2*q1  + p4*q2  + p6*q3  + p8*q4  + p10*q5  + p9 *q6 ) /r2/ww/uu
				gs2(m,2) = (p2*q7  + p4*q8  + p6*q9  + p8*q10 + p10*q11 + p9 *q12) /r2/ww/vv
				gs2(m,3) = (p2*q13 + p4*q14 + p6*q15 + p8*q16 + p10*q17 + p9 *q18) /r2/wl
				gs2(m,4) = (p2 + p4 + p6) / r2
				
			end do
			
			j = (it-1) * ks0 + is
			
			nborder2 = (nbdeg+1) * (nbdeg+1)
			
			do ll = 1, nborder2
				do mm = 1, ndd2
					gf1(j,ll,1) = gf1(j,ll,1) + bi(mm,ll) * gs1(mm,1)
					gf1(j,ll,2) = gf1(j,ll,2) + bi(mm,ll) * gs1(mm,2)
					gf1(j,ll,3) = gf1(j,ll,3) + bi(mm,ll) * gs1(mm,3)
					gf1(j,ll,4) = gf1(j,ll,4) + bi(mm,ll) * gs1(mm,4)
					
					gf2(j,ll,1) = gf2(j,ll,1) + bi(mm,ll) * gs2(mm,1)
					gf2(j,ll,2) = gf2(j,ll,2) + bi(mm,ll) * gs2(mm,2)
					gf2(j,ll,3) = gf2(j,ll,3) + bi(mm,ll) * gs2(mm,3)
					gf2(j,ll,4) = gf2(j,ll,4) + bi(mm,ll) * gs2(mm,4)
				end do
			end do

		end do
	end do
	
	cnorm = da2_lc * 0.25d0 / pi
	
	CALL MKKN(ds1(1,1),gf1(1,1,1),cnorm,kfe,nbi,nbdeg,ks0,kt0)
	CALL MKKN(ds1(1,2),gf1(1,1,2),cnorm,kfe,nbi,nbdeg,ks0,kt0)
	CALL MKKN(ds1(1,3),gf1(1,1,3),cnorm,kfe,nbi,nbdeg,ks0,kt0)
	CALL MKKN(ds1(1,4),gf1(1,1,4),cnorm,kfe,nbi,nbdeg,ks0,kt0)
	CALL MKKN(ds2(1,1),gf2(1,1,1),cnorm,kfe,nbi,nbdeg,ks0,kt0)
	CALL MKKN(ds2(1,2),gf2(1,1,2),cnorm,kfe,nbi,nbdeg,ks0,kt0)
	CALL MKKN(ds2(1,3),gf2(1,1,3),cnorm,kfe,nbi,nbdeg,ks0,kt0)
	CALL MKKN(ds2(1,4),gf2(1,1,4),cnorm,kfe,nbi,nbdeg,ks0,kt0)
	
return
end subroutine HORDS

!==========
subroutine VERDS(xs,ys,nbdeg, nd2, da2_lc)

	use prm_matrix
	use prm_inv
	use mod_mkkn
	
	implicit none
	
	real(8), intent(in) :: xs, ys, da2_lc
	integer, intent(in) :: nbdeg, nd2
	
	integer :: is, it, m, mu, mv, m1, m3, ls1, lt1, j, ll, mm
	integer :: nborder2, ndd2, ms2l, mt2l
	real(8) :: u1, u2, u3, v1, v2, v3, w1, w2, w3
	real(8) :: uu, vv, ww, wl, delw
	real(8) :: q1, q2, q3, q4, q5, q6, q7, q8, q9
	real(8) :: q10, q11, q12, q13, q14, q15, q16, q17, q18
	real(8) :: dx, dy, dz, r2, rr, ggx, ggy, ggz, xg, xg2, e
	real(8) :: ggx2, ggy2, ggz2, ggxy, adex, c4, cnorm
	real(8) :: p1, p2, p3, p4, p5, p6
	real(8) :: fde
	
	real(8), parameter :: pi = 3.14159265358979324d0
	
	gf1(1:kfe,1:16,1:4) = 0.d0
	
	ndd2 = nd2 * nd2
	ms2l = nd2 * Ks0
	mt2l = nd2 * Kt0
	
	do it =1, kt0
	
		lt1 = (it-1) * nd2
		
		do is = 1, ks0
			
			ls1 = (is-1) * nd2
			
			do m = 1, ndd2
				
				m1 = (m - int((m-1) / nd2) * nd2) + ls1
				m3 = (1 + int((m-1) / nd2)      ) + lt1
				
				mu = m1 - int((m1-1)/(ms2l-1))
				mv = m3 - int((m3-1)/(mt2l-1))
				
				u1 = xtr(mu+1) - xtr(mu)
				u2 = ytr(mu+1,mv) - ytr(mu,mv)
				u3 = ztr(mu+1,mv) - ztr(mu,mv)
				
				v1 = 0.d0
				v2 = ytr(mu,mv+1) - ytr(mu,mv)
				v3 = ztr(mu,mv+1) - ztr(mu,mv)
				
				w1 = u2 * v3 - u3 * v2
				w2 =         - u1 * v3
				w3 = u1 * v2
				
				v1 = w2 * u3 - w3 * u2
				v2 = w3 * u1 - w1 * u3
				v3 = w1 * u2 - w2 * u1
				
				uu = dsqrt(u1*u1 + u2*u2 + u3*u3)
				vv = dsqrt(v1*v1 + v2*v2 + v3*v3)
				wl =  w1*w1 + w2*w2 + w3*w3
				ww = dsqrt(wl)
				delw = Delt * wl
				
				q1  = 2.d0 * u1 * w1
				q2  = 2.d0 * u2 * w2
				q3  = 2.d0 * u3 * w3
				q4  = u1 * w2 + u2 * w1
				q5  = u2 * w3 + u3 * w2
				q6  = u3 * w1 + u1 * w3
				
				q7  = 2.d0 * v1 * w1
				q8  = 2.d0 * v2 * w2
				q9  = 2.d0 * v3 * w3
				q10 = v1 * w2 + v2 * w1
				q11 = v2 * w3 + v3 * w2
				q12 = v3 * w1 + v1 * w3
				
				q13 = delw + 2.d0 * w1 * w1
				q14 = delw + 2.d0 * w2 * w2
				q15 = delw + 2.d0 * w3 * w3
				q16 = 2.d0 * w1 * w2
				q17 = 2.d0 * w2 * w3
				q18 = 2.d0 * w3 * w1
				
				dx = xs - xtr(m1)
				dy = ys - ytr(m1,m3)
				dz = -ztr(m1,m3)
				
				r2 = dx*dx + dy*dy + dz*dz
				rr = dsqrt(r2)
				
				ggx = dx / rr
				ggy = dy / rr
				ggz = dz / rr
				
				xg = ggz + 1.d0
				xg2 = xg * xg
				
				fde = 1.d0/(1.d0+Delt)
				e = 3.d0 * ggz - fde * (xg + 1.d0) / xg2
				ggx2 = ggx * ggx
				ggy2 = ggy * ggy
				ggz2 = ggz * ggz
				ggxy = ggx * ggy
				
				adex = fde / xg - ggz
				c4 = fde - 1.d0 + 3.d0 * ggz2
				
				p1  = e * ggx2 + adex
				p2  = e * ggy2 + adex
				p3  = ggz * c4
				p4  = ggxy * e * 2.d0
				
				p5  = 6.d0 * ggy * ggz2
				p6  = 6.d0 * ggx * ggz2
				
				gs1(m,1) = (p1*q1  + p2*q2  + p3*q3  + p4*q4  + p5 *q5  + p6*q6 ) /r2/ww/uu
				gs1(m,2) = (p1*q7  + p2*q8  + p3*q9  + p4*q10 + p5 *q11 + p6*q12) /r2/ww/vv
				gs1(m,3) = (p1*q13 + p2*q14 + p3*q15 + p4*q16 + p5 *q17 + p6*q18) /r2/wl
				gs1(m,4) = (p1 + p2 + p3) / r2
				
			end do
			
			j = (it-1) * ks0 + is
			
			nborder2 = (nbdeg+1) * (nbdeg+1)
			
			do ll = 1, nborder2
				do mm = 1, ndd2
					gf1(j,ll,1) = gf1(j,ll,1) + bi(MM,LL) * gs1(mm,1)
					gf1(j,ll,2) = gf1(j,ll,2) + bi(MM,LL) * gs1(mm,2)
					gf1(j,ll,3) = gf1(j,ll,3) + bi(MM,LL) * gs1(mm,3)
					gf1(j,ll,4) = gf1(j,ll,4) + bi(MM,LL) * gs1(mm,4)
				end do
			end do
			
		end do
	end do

	cnorm = da2_lc * 0.25d0 / pi

	CALL MKKN(ds1(1,1),gf1(1,1,1),cnorm,kfe,nbi,nbdeg,ks0,kt0)
	CALL MKKN(ds1(1,2),gf1(1,1,2),cnorm,kfe,nbi,nbdeg,ks0,kt0)
	CALL MKKN(ds1(1,3),gf1(1,1,3),cnorm,kfe,nbi,nbdeg,ks0,kt0)
	CALL MKKN(ds1(1,4),gf1(1,1,4),cnorm,kfe,nbi,nbdeg,ks0,kt0)

return
end subroutine VERDS

!==========
subroutine SLIPM(nbdeg, nd2)

	use prm_matrix
	use prm_inv
	use mod_mkkn

	implicit none

	integer, intent(in) :: nbdeg, nd2
	
	integer :: is, it, m, mu, mv, m1, m3, ls1, lt1, j, ll, mm
	integer :: nborder2, ndd2, ms2l, mt2l
	real(8) :: u1, u2, u3, v1, v2, v3, w1, w2, w3
	real(8) :: uu, vv, cnorm, d1, d2
	
	gf1(1:kfe,1:16,1:4) = 0.d0
	gf2(1:kfe,1:16,1:4) = 0.d0
	
	ndd2 = nd2 * nd2
	ms2l = nd2 * Ks0
	mt2l = nd2 * Kt0
	
	do it =1, kt0

		lt1 = (it-1) * nd2

		do is = 1, ks0

			ls1 = (is-1) * nd2

			do m = 1, ndd2

				m1 = (m - int((m-1) / nd2) * nd2) + ls1
				m3 = (1 + int((m-1) / nd2)      ) + lt1

				mu = m1 - int((m1-1)/(ms2l-1))
				mv = m3 - int((m3-1)/(mt2l-1))

				u1 = xtr(mu+1) - xtr(mu)
				u2 = ytr(mu+1,mv) - ytr(mu,mv)
				u3 = ztr(mu+1,mv) - ztr(mu,mv)

				v1 = 0.d0
				v2 = ytr(mu,mv+1) - ytr(mu,mv)
				v3 = ztr(mu,mv+1) - ztr(mu,mv)

				w1 = u2 * v3 - u3 * v2
				w2 =         - u1 * v3
				w3 = u1 * v2

				v1 = w2 * u3 - w3 * u2
				v2 = w3 * u1 - w1 * u3
				v3 = w1 * u2 - w2 * u1

				uu = dsqrt(u1*u1 + u2*u2 + u3*u3)
				vv = dsqrt(v1*v1 + v2*v2 + v3*v3)

				gs1(m,1) = u1 / uu
				gs1(m,2) = u2 / uu
				gs1(m,3) = u3 / uu

				gs2(m,1) = v1 / vv
				gs2(m,2) = v2 / vv
				gs2(m,3) = v3 / vv
				
			end do
			
			j = (it-1) * ks0 + is
			
			nborder2 = (nbdeg+1) * (nbdeg+1)
			
			do ll = 1, nborder2
				do mm = 1, ndd2
					gf1(j,ll,1) = gf1(j,ll,1) + bi(mm,ll) * gs1(mm,1)
					gf1(j,ll,2) = gf1(j,ll,2) + bi(mm,ll) * gs1(mm,2)
					gf1(j,ll,3) = gf1(j,ll,3) + bi(mm,ll) * gs1(mm,3)
					
					gf2(j,ll,1) = gf2(j,ll,1) + bi(mm,ll) * gs2(mm,1)
					gf2(j,ll,2) = gf2(j,ll,2) + bi(mm,ll) * gs2(mm,2)
					gf2(j,ll,3) = gf2(j,ll,3) + bi(mm,ll) * gs2(mm,3)
				end do
			end do

		end do
	end do

	cnorm = 1.d0

	CALL MKKN(ds1(1,1),gf1(1,1,1),cnorm,kfe,nbi,nbdeg,ks0,kt0)
	CALL MKKN(ds1(1,2),gf1(1,1,2),cnorm,kfe,nbi,nbdeg,ks0,kt0)
	CALL MKKN(ds1(1,3),gf1(1,1,3),cnorm,kfe,nbi,nbdeg,ks0,kt0)
	CALL MKKN(ds2(1,1),gf2(1,1,1),cnorm,kfe,nbi,nbdeg,ks0,kt0)
	CALL MKKN(ds2(1,2),gf2(1,1,2),cnorm,kfe,nbi,nbdeg,ks0,kt0)
	CALL MKKN(ds2(1,3),gf2(1,1,3),cnorm,kfe,nbi,nbdeg,ks0,kt0)

	do j = 1, nbi

		d1 = dsqrt( ds1(j,1)**2 + ds1(j,2)**2 + ds1(j,3)**2 )
		d2 = dsqrt( ds2(j,1)**2 + ds2(j,2)**2 + ds2(j,3)**2 )

		ds1(j,1) = ds1(j,1) / d1
		ds1(j,2) = ds1(j,2) / d1
		ds1(j,3) = ds1(j,3) / d1

		ds2(j,1) = ds2(j,1) / d2
		ds2(j,2) = ds2(j,2) / d2
		ds2(j,3) = ds2(j,3) / d2

	end do

return
end subroutine SLIPM

!==========

end module mod_jaclib

