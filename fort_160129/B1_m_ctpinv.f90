!=================================================
! SUBROUTINES
!	SET_CTP(HD,flobspd,CAL,FMP,NSOURCE)
!   OBSDATAIN(iu)
!
! setup file を読み込む
!
!=================================================

module mod_ctpinv

	implicit none

contains
!===========
subroutine SET_CTP(hd,flobspd,flobspos,cal,fmp,nsource)

	use prm_inv

	implicit none

	integer, intent(out) :: nsource
	character(len=24), intent(out) :: hd, flobspd, flobspos, cal, fmp

	integer :: ls, lsin, k, kk, dummy

	! == same as SETA ==

	read(9,*) Alat0, Alng0, Icord
	read(9,'(A)') flobspd
	read(9,'(A)') flobspos
	read(9,*) Delt
	read(9,*) nsource

	print *, ' ORIGIN OF COODINATE(LAT,LON) =( ', Alat0, Alng0,')'
	print *, '      GEODETIC DATA FILE NAME =', flobspd
	print *, '     ELASTIC CONST(LAMBDA/MU) =', Delt

	print *, '  SOURCE TYPE:  /STR SLIP,DIP SLIP,OPEN C,EXP/)'
	print *, '  LIST ORDER:(STATION / DATA TYPE / SOURCE TYPE / BASE FUNCTION)'

	!----------------------------

	! == same as SETFCP but LOOP for "nsource" ==

	do ls = 1, nsource

		read(9,*)
		read(9,*) lsin

		if(lsin /= ls) stop "LS INPUT ERROR"

		read(9,*)
		read(9,'(A)') Fcof(ls)
		read(9,*) Bphi(ls), Bdlt(ls)
		read(9,*) Blat(ls), Blon(ls)
		read(9,*) Blen(ls), Bwid(ls)
		read(9,*) Ju0(ls), Jv0(ls), Ndgf(ls)
		read(9,*)
		read(9,'(A)') Hd_jac(ls)
		read(9,'(A)') Invjac(ls)
		read(9,'(A)') Fwdjac(ls)
		read(9,*) Blaf(ls), Blof(ls)
		read(9,*) Blef(ls), Bwif(ls), Bdep(ls)
		read(9,*) Js0(ls), Jt0(ls), Jd2(ls), Ndgs(ls)

		Ipl(ls) = 0

		if(Fcof(ls) == 'PLANE FAULT             ') Ipl(ls) = 1


	! +++ PARAMETER SET +++ from SETFCP in gdjac

	Nod2(ls) = (Ndgf(ls)+1) * (Ndgf(ls)+1)
	Ju1(ls)  = Ju0(ls) + Ndgf(ls)
	Jv1(ls)  = Jv0(ls) + Ndgf(ls)
	Nbod2(ls) = (Ndgs(ls)+1) * (Ndgs(ls)+1)

	! +++++++

	end do

	!---------------

	read(9,*)
	read(9,'(A)') hd

	do ls = 1, nsource
		read(9,*)
		read(9,*) lsin

		if(lsin /= ls) stop "LS INPUT ERROR"

		read(9,'(A)') Fscd(ls)
		read(9,*) Just(1,ls), Just(2,ls), Just(3,ls), Just(4,ls)
		read(9,*) dummy
		read(9,*) Rss(ls), Rop(ls), Rex(ls)
		read(9,*) Ius(1,ls), Ius(2,ls), Ius(3,ls), Ius(4,ls), Mlg(ls)
		read(9,*) Iboc(1,ls), Iboc(2,ls), Iboc(3,ls), Iboc(4,ls)

		kk = 0
		do k = 1, 4
			if(Just(k,ls) == 0) cycle
			kk = kk + 1
			Just(k,ls) = kk
		end do

	end do

	!---------------

	read(9,*)
	read(9,'(A)') cal
	read(9,'(A)') fmp

	read(9,*)
	read(9,*) Fact(1), Fact(2), Fact(3)
	read(9,*) Fabs(1), Fabs(2), Fabs(3)
	read(9,*) Vrp

	!---------------

	read(9,*)
	read(9,*) Glatmin, Glatmax, Glonmin, Glonmax

	read(9,'(A)') Fgps
	read(9,'(A)') Fpos
	read(9,'(A)') Fdum
	read(9,*) E_gsi, E_jcg
	read(9,*) Setalpha

return
end subroutine SET_CTP

!==========

subroutine OBSDATAIN(iu)

	use prm_matrix
	use prm_inv
	use mod_trans

	implicit none

	integer, intent(in) :: iu

	integer :: i, ii, i1, i2, j, j1, j2, k
	integer :: nr(Kis), nst(15,Kis)
	real(8) :: stalat, stalng, datx(Kis), daty(Kis), datz(Kis), ea
	real(8) :: sx, sy, sl, rx, ry, ua, ub

	!real(8) :: dl(Kis)             ! unused in this version
	!integer :: k3(Kis), nl(2,Kis)  ! unused in this version

	!=============

	read(iu,*) Isn, Ih, Iv, Il
	print *, Ih*2 + Iv

	if(Ih + Iv /= Isn) stop "ERROR"
	if(Isn > Kis)      stop 'STATION DIMENSION ERROR'

	if(Ih*2 + Iv > Kis)    stop 'INSUFFICINT ARRAY DIM(KIS)'
	if(Ih*2 + Iv > Kio)    stop 'INSUFFICINT ARRAY DIM(KIO)'
	if(il        > Kis)    stop 'INSUFFICINT ARRAY DIM(KIS)'

	do i = 1, Ih
		read(iu,*) ii, stalat, stalng, datx(i), daty(i), nr(i), (nst(k,i), k=1,nr(i))
		if(ii /= i) stop ' INPUT IH DATA ERROR'
		CALL LL2XY(stalat,stalng, St(1,i),St(2,i), 0, Alat0,Alng0,Icord) ! set St
	end do

	do i = Ih+1, Isn
		read(iu,*) ii, stalat, stalng, datz(i), nr(i), (nst(k,i), k=1,nr(i))
		if(ii /= i) stop ' INPUT IV DATA ERROR'
		CALL LL2XY(stalat,stalng, St(1,i),St(2,i), 0, Alat0,Alng0,Icord) ! set St
		datz(i) = datz(i) + Vrp
	end do

	if(Il /= 0) stop ' INPUT IL DATA ERROR (not support for IL data)'
	!do i =1, Il
	!	read(iu,*) ii, nl(1,i), nl(2,i), dl(i), k3(i)
	!	if(ii /= i) stop ' INPUT IL DATA ERROR'
	!end do

	print *, ' DATA INPUT OK. (ISN,IH,IV,IL)= ', Isn, Ih, Iv, Il

	!--- INPUT OBS DATA END -------------------------------

	Itotal = 0

	do i = 1, Ih

		if(nst(1,i) == 0) cycle

		i1 = Itotal + 1
		i2 = Itotal + 2

		Itotal = Itotal + 2

		Ia(i1) = 1
		Ia(i2) = 2
		Ib(i1) = i
		Ib(i2) = i
		Ic(i1) = 0
		Ic(i2) = 0

		Y_data(i1) = datx(i)
		Y_data(i2) = daty(i)

		ea = dble(nst(1,i))
		E_prior(i1) = Fact(1) * dsqrt(ea*ea + Fabs(1)*(datx(i)*datx(i) + daty(i)*daty(i)))
		E_prior(i2) = E_prior(i1)

	end do

	do i = Ih+1, Isn

		if(nst(1,i) == 0) cycle

		Itotal = Itotal + 1

		Ia(itotal) = 3
		Ib(itotal) = i
		Ic(itotal) = 0

		Y_data(Itotal) = datz(i)
		ea = dble(nst(1,i))
		E_prior(itotal) = Fact(2) * dsqrt(ea*ea + Fabs(2) * datz(i)*datz(i))

	end do

	!do i = 1, Il
	!	if(k3(i) == 0) cycle
	!	Itotal = Itotal + 1
	!	Ia(Itotal) = 4
	!	Ib(Itotal) = nl(1,i)
	!	Ic(Itotal) = nl(2,i)
	!	Y_data(itotal) = dl(i)
	!	E_prior(itotal) = Fact(3) * dsqrt(dble(k3(i)**2) + Fabs(3) * dl(i)*dl(i))
	!end do

	if(Itotal > Kio) stop ' ITOTAL, INSUFFICIENT ARRAY DIM (KIO)'

	print *,    ' TOTAL NUMBER OF DATA =', Itotal
	write(77,*) ' TOTAL NUMBER OF DATA =', Itotal

	!========================
	! OBS DATA set ここまで
	!========================

	!========================
	! SETAG Sta_corr(Kio,Kio) (old; ag) を作る
	!  nr は，固定点のステーション番号を示す nst の数を示す。
	!  nr(i) == 1 のときはどこにも固定していない系になっている。
	!  通常のプレート固定の場合は nr(*) = 1 としている（これでいいかは，要考察）。
	!  その際，Sta_corrは単位行列になる（と思う）
	!========================

	Sta_corr = 0.d0

	do j = 1, Itotal ! OBS data 分のループ（3成分独立；Ih + Ih + Iv）

		if(Ia(j) == 1) cycle ! X-ward は処理しない

		if(Ia(j) == 2) then  ! Y-ward のとき

			i = Ib(j) ! Station 番号（水平で一つ；Ih + Iv）

			Sta_corr(j-1, 2*i-1) = 1.d0
			Sta_corr(j  , 2*i  ) = 1.d0

			if(nr(i) == 1) cycle

			!====
			! 以下，nr(*) = 1 のときは不要
			!====

			Ic(j-1) = 1
			Ic(j  ) = 1

			if(nr(i) == 2) then

				j1 = nst(2,i)
				Sta_corr(j-1, 2*j1-1) = Sta_corr(j-1, 2*j1-1) - 1.d0
				Sta_corr(j  , 2*j1  ) = Sta_corr(j  , 2*j1  ) - 1.d0

			else if(nr(i) >= 3) then

				do k = 2, nr(i)-1

					j1 = nst(k  ,i)
					j2 = nst(k+1,i)

					sx = St(1,j2) - St(1,j1)
					sy = St(2,j2) - St(2,j1)
					sl = dsqrt(sx*sx + sy*sy)
					sx = sx / sl
					sy = sy / sl

					rx = (St(1,i) - St(1,j1)) / sl
					ry = (St(2,i) - St(2,j1)) / sl

					ua = rx * sx + ry * sy
					ub = rx * sy - ry * sx

					Sta_corr(j-1, 2*j1-1) = Sta_corr(j-1, 2*j1-1) + ua - 1.d0
					Sta_corr(j-1, 2*j1  ) = Sta_corr(j-1, 2*j1  ) + ub
					Sta_corr(j-1, 2*j2-1) = Sta_corr(j-1, 2*j2-1) - ua
					Sta_corr(j-1, 2*j2  ) = Sta_corr(j-1, 2*j2  ) - ub

					Sta_corr(j  , 2*j1-1) = Sta_corr(j  , 2*j1-1) - ub
					Sta_corr(j  , 2*j1  ) = Sta_corr(j  , 2*j1  ) + ua - 1.d0
					Sta_corr(j  , 2*j2-1) = Sta_corr(j  , 2*j2-1) + ub
					Sta_corr(j  , 2*j2  ) = Sta_corr(j  , 2*j2  ) - ua

				end do

			end if

		else if(Ia(j) == 3) then  ! Z-ward のとき

			i = Ib(j)
			Sta_corr(j, i+Ih) = 1.d0

			if(nr(i) /= 2) cycle

			!====
			! 以下，nr(*) = 1 のときは不要
			!====

			Sta_corr(j, nst(2,i)+Ih) = -1.d0
			Ic(j) = 1

		!Ia == 4 になることは，現状はない
		!else if(ia(j) == 4) then
		!	sx = st(1,ib(j)) - st(1,ic(j))
		!	sy = st(2,ib(j)) - st(2,ic(j))
		!	sl = dsqrt(sx*sx + sy*sy)
		!	sx = sx / sl
		!	sy = sy / sl
		!	Sta_corr(j, 2*ib(j)-1) =  sx
		!	Sta_corr(j, 2*ib(j)  ) =  sy
		!	Sta_corr(j, 2*ic(j)-1) = -sx
		!	Sta_corr(j, 2*ic(j)  ) = -sy

		end if
	end do

return
end subroutine OBSDATAIN

!==========

end module
