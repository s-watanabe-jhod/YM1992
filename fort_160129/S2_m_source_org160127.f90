!=================================================
! SUBROUTINES
!	SOURCE(hd,ls,iu)
!
!=================================================

module mod_source

	implicit none

contains
!==========
subroutine SOURCE(hd,ls,iu,iuxyz,iuelp,iucnt)

	use prm_matrix
	use prm_inv
	use mod_wbb
	use mod_bspline
	use mod_trans
	use mod_subsrc
	use mod_qrdec
	
	implicit none
	
	character(len=24), intent(in) :: hd
	integer, intent(in) :: ls, iu, iuxyz, iuelp, iucnt
	
	real(8), parameter :: pi = 3.14159265358979324d0
	
	integer :: i, j, k, isp, node_x, node_y, is, it, k1, k2, l1, l2, n, np, nq, ml
	integer :: lx, ly, idb(16), numsrc !, ii
	integer :: js0s, js1s, jt0s, jt1s
	real(8) :: q(4), we(3), sa(3), sb(3), cm(4,4), scv(10), bsv(16), vec(3), vec_th, vec_r
	real(8) :: q_R(4), cm_R(4,4), unit(Kmp), b_sp(Kmp,Kmp), btb(Kmp,Kmp), btb_i(Kmp,Kmp), btb_sum(Kmp)
	real(8) :: b_dagger(Kmp,Kmp), bdagb(Kmp,Kmp), invtemp1, invtemp2 !, v(Kmp)
	real(8) :: elp(2), elp_th, sign, cnt_r !, celp(2): unused in this version
	real(8) :: px, py, xp1, yp1, xp, yp, zp, sa1, sa2, sb1, sb2, plat, plong
	
	scv(1:10) = 0.d0
	
	js0s = (Ius(1,ls) - Ndgs(ls) - 1) * Jd2(ls) + 1
	if(js0s <= 0) js0s = 1
	
	jt0s = (Ius(3,ls) - Ndgs(ls) - 1) * Jd2(ls) + 1
	if(jt0s <= 0) jt0s = 1
	
	js1s = Ius(2,ls) * Jd2(ls)
	if(js1s > Js0(ls)*Jd2(ls)) js1s = Js0(ls)*Jd2(ls)
	
	jt1s = Ius(4,ls) * Jd2(ls)
	if(jt1s > Jt0(ls)*Jd2(ls)) jt1s = Jt0(ls)*Jd2(ls)
	
	!=== ちょっと入れてみる 
	!  MLG: SCD 出力のために震源の大きさを計算するグリッド数（ x 方向）
	!  Jd2: Bspline のセグメントを分ける数。
	! Mlg(ls) = (js1s - js0s) / Jd2(ls)
	!===
	
	print *,    ' START OF SOURCE DISTRIBUTION COMPUTATION'
	print *,    ' MAXIMUM OUTPUT GRID POINTS(1-D)(MLG) = ', Mlg(ls)
	print *,    ' OUTPUT IN LOCAL COORDINATES '
	write(77,*) ' MAXIMUM OUTPUT GRID POINTS(1-D)(MLG) = ', Mlg(ls)
	
	!ml = (js1s - js0s) / Mlg(ls) + 1
	ml = Jd2(ls)
	lx = (js1s - js0s) / ml + 2
	ly = (jt1s - jt0s) / ml + 2
	
	node_x = Js0(ls) + Ndgs(ls) ! b-spline node の数 x 方向
	node_y = Jt0(ls) + Ndgs(ls) ! b-spline node の数 y 方向
	
	print *, lx, ly
	
	write(iu,*) "  ** SOURCE PARAMETER DISTRIBUTION **"
	write(iu,*) hd
	write(iu,*) (Just(k,ls), k=1,4), "  : USED SOURCE(STR S.,DIP S.,OPC.,EXP.)"
	write(iu,*) lx, ly, "  :NUMBER OF MESHES ALONG STRIKE- AND DIP- AXIS"
	
	write(iu,*) 
	write(iu,*) "  (I J) (XP,YP,ZP), UNIT VECTORS /STLIKE(X,Y,Z), DIP(X,Y,Z)/"
	write(iu,*) "  SOURCE MAGNITUDE (STLIKE, DIP, OPNING C., EXPLOSION)"
	write(iu,*) "  COV(S-S,S-D,D-D,S-O,D-O,O-O,S-E,D-E,O-E,E-E)"
	
	do j = 1, ly
		do i = 1, lx
			
			is = (i-1) * ml + js0s
			it = (j-1) * ml + jt0s
			
			if(is > js1s) is = js1s
			if(it > jt1s) it = jt1s
			
			k1 = (is-1) / Jd2(ls) + 1
			l1 = is - (k1-1) * Jd2(ls)
			
			px = (dble(l1) - 0.5d0) / dble(Jd2(ls))
			
			k2 = (it-1) / Jd2(ls) + 1
			l2 = it - (k2-1) * Jd2(ls)
			
			py = (dble(l2) - 0.5d0) / dble(Jd2(ls))
			
			CALL BSPARI(k1, k2, idb, Ndgs(ls), node_x, node_y)
			CALL BSVALI(px, py, bsv, Ndgs(ls), Nbod2(ls))
			
			CALL FSNPA(we, sa, sb, xp1, yp1, zp, is, it)
			
			CALL ROTATE(sa(1),sa(2),sa1,sa2,1,-1, Bphi(ls),Bdlt(ls))
			CALL ROTATE(sb(1),sb(2),sb1,sb2,1,-1, Bphi(ls),Bdlt(ls))
			CALL ROTATE(xp1,yp1,xp,yp,1,-1, Bphi(ls),Bdlt(ls))
			
			CALL SVESE(q,cm,bsv,idb,ls,Nbod2(ls))
			
			!==================
			CALL SVESE_RESOL(q_R,cm_R,bsv,idb,ls,Nbod2(ls))
			
				CALL LL2XY(plat,plong,xp,yp,1, Alat0,Alng0,Icord)
				write(98,'(6(E12.6,1X))') plat, plong, (q_R(k), k=1,4)
			
			do isp = 1, Jsp
				
				if(i >= lx-1 .or. j >= ly-1) exit ! 行列がゼロになるのを避ける。
				
				unit(1:Jsp) = 0.d0
				unit(isp)   = 1.d0
				
				CALL SVESE_I(q_R,cm_R,bsv,idb,ls,Nbod2(ls),unit)
				
				b_sp(i+(lx-2)*(j-1),              isp) = q_R(1)
				b_sp(i+(lx-2)*(j-1)+(lx-2)*(ly-2),isp) = q_R(2)
				! b_sp は 2*(lx-2)*(ly-2) × Jsp の行列（ix,iyの端を除いているのに注意）
				
			end do
			
			!==================
			
			do np = 1, 4
				n = np * (np-1) / 2
				do nq = 1, np
					scv(nq+n) = cm(nq,np)
				end do
			end do
			
			if(i == lx .or. j == ly) q(1:4) = 0.d0
			
			write(iu,'(2(I3,1X),3(F12.5,1X),2(1X,F6.4,1X,F6.4,1X,F6.4))') &
				i, j, xp, yp, zp, sa1, sa2, sa(3), sb1, sb2, sb(3)
			write(iu,'(4(E12.6,1X))' ) (q(k), k=1,4)
			write(iu,'(10(E10.4,1X))') (scv(k), k=1,10)
			
			if(i == lx .or. j == ly) then
				CALL LL2XY(plat,plong,xp,yp,1, Alat0,Alng0,Icord)
				write(iucnt,'(2X,2F10.3,1X,F10.3)') plong, plat, 0.d0
				cycle
			end if
			
			vec(1) = sa1   * q(1) + sb1   * q(2)
			vec(2) = sa2   * q(1) + sb2   * q(2)
			vec(3) = sa(3) * q(1) + sb(3) * q(2)
			
			vec_th = datan2(vec(2),vec(1)) * 180.d0 / pi
			vec_r  = dsqrt(vec(1)*vec(1) + vec(2)*vec(2))
				if(vec_r == 0.d0) then
					write(iuxyz,*)
					cycle
				end if
			
			! 誤差楕円を求める
	!		celp(1) =  scv(1) + scv(3)
	!		celp(2) = (scv(1) - scv(3))**2 + 4.d0*scv(2)**2
	!		elp(1) = celp(1) + dsqrt(celp(2))
	!		elp(2) = celp(1) - dsqrt(celp(2))
			
	!		elp(1) = dsqrt(elp(1)/2.d0)
	!		elp(2) = dsqrt(elp(2)/2.d0)
			
	!		elp_th = (elp(1)**2 - elp(2)**2) / scv(2)
	!		elp_th = datan(elp_th) * 180.d0 / pi
	!		elp_th = elp_th + bphi(ls)
			
			elp_th = bphi(ls)
			elp(1) = dsqrt(scv(1))
			elp(2) = dsqrt(scv(3))
			
			sign = dcos( (vec_th - (rss(ls)+Bphi(ls))) * pi/180.d0 )
			
			if(sign >= 0.d0) cnt_r =  vec_r
			if(sign <  0.d0) cnt_r = -vec_r
			
			CALL LL2XY(plat,plong,xp,yp,1,  Alat0,Alng0,Icord)
			write(iuxyz,'(2X,2F10.3,1X,2F15.3)') plong, plat, vec_th, vec_r
			write(iuelp,'(2X,2F10.3,1X,3F15.3)') plong, plat, elp_th, elp(1), elp(2)
			write(iucnt,'(2X,2F10.3,1X, F15.3)') plong, plat, cnt_r
		
		end do
	end do
	
	!======= 
	btb      = 0.d0
	btb_i    = 0.d0
	btb_sum  = 0.d0
	b_dagger = 0.d0
	bdagb    = 0.d0
	
	numsrc = (lx-2)*(ly-2)*2
	
	do i = 1, numsrc
		do j = 1, numsrc
			do isp = 1, Jsp
				btb(i,j)   = btb(i,j)   + b_sp(i,isp) * b_sp(j,isp) ! BTB (numsrc,numsrc)
				btb_i(i,j) = btb_i(i,j) + b_sp(i,isp) * b_sp(j,isp) ! BTB (numsrc,numsrc)
			end do
		end do
	end do
	
	do i = 1, numsrc
		if(dabs(btb_i(i,i)) < 1.d-12) print *, i, btb_i(i,i)
	end do
	
	! == 逆行列を計算？
!	do i = numsrc, 1, -1
!		ii = i
!		unit(1:ii) = 0.d0
!		unit(ii)   = 1.d0
!		
!		!CALL QRDEC(btb_i, unit, invtemp1, Kmp, Kmp, Kmp, wk1,lu, k,1)
!		
!		!print *, unit(ii)
!		
!		CALL RINV(v,unit,btb_i,Kmp,ii,0)
!		
!		do j = 1, ii
!			btb_i(i,j) = v(j) ! (BTB)^-1
!			btb_i(j,i) = v(j) ! (BTB)^-1
!		end do
!	end do
	! ==
	
	
	! == 逆行列を計算
	do k = 1, numsrc
		
		invtemp1 = btb_i(k,k)
		btb_i(k,k) = 1.d0
		
		if(dabs(invtemp1) < 1.d-12) print *, i, btb(k,k), invtemp1
		!if(dabs(invtemp1) < 1.d-12) stop "!!WARNING!! in BTB inverse "
		! invtemp1 が小さいときは，逆行列が計算できないことに注意
		
		do j = 1, numsrc
			btb_i(k,j) = btb_i(k,j)/invtemp1
		end do
		
		do i = 1, numsrc
			
			if(i /= k) then
				
				invtemp2 = btb_i(i,k)
				btb_i(i,k) = 0.d0
				
				do j = 1, numsrc
					btb_i(i,j) = btb_i(i,j) - invtemp2 * btb_i(k,j) ! (BTB)^-1
				end do
				
			end if
			
		end do
		
	end do
	! ==
	
	do i = 1, numsrc
		do j = 1, numsrc
			do k = 1, numsrc
				b_dagger(i,j) = b_dagger(i,j) + btb(j,k) * btb_i(i,k)
			end do
			write(92,*) i, j, b_dagger(i,j)
			!if(dabs(b_dagger(i,j)) > 0.1d0) write(*,*) i, j, b_dagger(i,j) ! check 用に b_dagger を借りた
		end do
	end do
	
	! ================
	b_dagger = 0.d0
	btb_sum  = 0.d0
	
	do i = 1, numsrc
		do isp = 1, Jsp
			do j = 1, numsrc
				b_dagger(i,isp) = b_dagger(i,isp) + btb_i(i,j) * b_sp(j,isp) ! BTB (numsrc,numsrc)
			end do
			write(91,*) i, isp, b_dagger(i,isp)
			!write(91,*) i, btb_sum(i)
			!btb_sum(i) = btb_sum(i) + btb(i,j)
		end do
	end do
	
	do i = 1, numsrc
		do j = 1, numsrc
			do isp = 1, Jsp
				bdagb(i,j) = bdagb(i,j) + b_sp(j,isp) * b_dagger(i,isp)
			end do
			write(92,*) i, j, bdagb(i,j)
			!if(dabs(bdagb(i,j)) > 0.01d0) write(*,*) i, j, bdagb(i,j)
		end do
	end do
	
	
	!======= 
	
return
end subroutine SOURCE

end module mod_source

