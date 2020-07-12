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
	use prm_var
	use prm_inv

	use mod_wbb
	use mod_bspline
	use mod_trans
	use mod_subsrc

	implicit none

	character(len=24), intent(in) :: hd
	integer, intent(in) :: ls, iu, iuxyz, iuelp, iucnt

	integer :: i, j, k, kss, ktt, is, it, k0, kk0, k1, k2, l1, l2, n, np, nq, ml
	integer :: lx, ly, idb(16)
	integer :: js0s, js1s, jt0s, jt1s
	real(8) :: q(4), we(3), sa(3), sb(3), cm(4,4), scv(10), bsv(16), vec(3), vec_th, vec_r
	real(8) :: elp(2), elp_th, celp(2), sign, cnt_r
	real(8) :: px, py, xp1, yp1, xp, yp, zp, sa1, sa2, sb1, sb2, plat, plong

	print *, ' START OF SOURCE DISTRIBUTION COMPUTATION'
	print *, ' MAXIMUM OUTPUT GRID POINTS(1-D)(MLG) = ', mlg(ls)
	print *, ' OUTPUT IN LOCAL COORDINATES '
	write(77,*) ' MAXIMUM OUTPUT GRID POINTS(1-D)(MLG) = ', mlg(ls)

	scv(1:10) = 0.d0

	js0s = (ius(1,ls) - nbdeg - 1) * nd2 + 1
		if(js0s <= 0) js0s = 1
	jt0s = (ius(3,ls) - nbdeg - 1) * nd2 + 1
		if(jt0s <= 0) jt0s = 1

	js1s = ius(2,ls) * nd2
		if(js1s > ms2) js1s = ms2
	jt1s = ius(4,ls) * nd2
		if(jt1s > mt2) jt1s = mt2

	ml = (js1s - js0s) / mlg(ls) + 1
	lx = (js1s - js0s) / ml + 2
	ly = (jt1s - jt0s) / ml + 2

	kss = ks0 + nbdeg
	ktt = kt0 + nbdeg

	write(iu,*) "  ** SOURCE PARAMETER DISTRIBUTION **"
	write(iu,*) hd
	write(iu,*) (just(k,ls), k=1,4), "  : USED SOURCE(STR S.,DIP S.,OPC.,EXP.)"
	write(iu,*) ls, ly, "  :NUMBER OF MESHES ALONG STRIKE- AND DIP- AXIS"

	write(iu,*)
	write(iu,*) "  (I J) (XP,YP,ZP), UNIT VECTORS /STLIKE(X,Y,Z), DIP(X,Y,Z)/"
	write(iu,*) "  SOURCE MAGNITUDE (STLIKE, DIP, OPNING C., EXPLOSION)"
	write(iu,*) "  COV(S-S,S-D,D-D,S-O,D-O,O-O,S-E,D-E,O-E,E-E)"

	do j = 1, ly
		it = (j-1) * ml + jt0s
			if(it > jt1s) it = jt1s

		k0 = (it-1) / nd2
		k2 = k0 + 1
		l2 = it - k0 * nd2
		py = (dble(l2) - 0.5d0) / dble(nd2)

		do i = 1, lx
			is = (i-1) * ml + js0s
				if(is > js1s) is = js1s

			kk0 = (is-1) / nd2
			k1 = kk0 + 1
			l1 = is - kk0 * nd2
			px = (dble(l1) - 0.5d0) / dble(nd2)

			CALL BSPARI(k1,k2,idb,nbdeg,kss,ktt)
			CALL BSVALI(px,py,bsv,nbdeg,nbod2)

			CALL FSNPA(we,sa,sb,xp1,yp1,zp,is,it)

			CALL TRANS(sa(1),sa(2),sa1,sa2,1,-1)
			CALL TRANS(sb(1),sb(2),sb1,sb2,1,-1)
			CALL TRANS(xp1,yp1,xp,yp,1,-1)

			CALL SVESE(q,cm,bsv,idb,ls,nbod2)

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
				CALL PLTXY(plat,plong,xp,yp,1)
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

			! Œë·‘È‰~‚ð‹‚ß‚é
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

			sign = dcos( (vec_th - (rss(ls)+bphi(ls))) * pi/180.d0 )

			if(sign >= 0.d0) cnt_r =  vec_r
			if(sign <  0.d0) cnt_r = -vec_r

			CALL PLTXY(plat,plong,xp,yp,1)
			write(iuxyz,'(2X,2F10.3,1X,2F15.3)') plong, plat, vec_th, vec_r
			write(iuelp,'(2X,2F10.3,1X,3F15.3)') plong, plat, elp_th, elp(1), elp(2)
			write(iucnt,'(2X,2F10.3,1X, F15.3)') plong, plat, cnt_r

		end do
	end do

return
end subroutine SOURCE

end module mod_source
