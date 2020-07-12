program GDFWD

!-------------------------------------------------------
!    FORWARD PROBLEM IN GEODETIC DATA
!    YOU MUST ALREADY DONE INVERSION WITH abic MINIMUM CRITERION
!    YOU MUST ALREADY PROCESSED PROGRAM GDJACFW
!-------------------------------------------------------

implicit none

character(len=24) :: flobspos, fmp, fwd, hd
character(len=24) :: fctp="                        "
character(len=24) :: fxyz="                        "
integer :: nsource, lname
character(len=6) :: cstation(3000)

call getarg(1,fctp)

open(9,file=fctp) ! fctp : ctp file

	CALL SETCPFW(hd,flobspos,fmp,fwd,nsource)

	open(33,file=flobspos,status='old')
		CALL GRDPO(33,cstation)
	close(33)

	CALL SETUPFW(flobspos,nsource, 32)

	open(34,file=fmp,status='old')
		CALL INSOL(hd, nsource, 34)
	close(34)

	CALL FORWARD(nsource)

	lname = len_trim(fwd)
	fxyz(1:lname-4) = fwd(1:lname-4)
	fxyz(lname-3:lname+4) = "_fwd.xyz"

	open(25,file=fxyz)
	open(24,file=fwd,status='unknown')
		CALL OUTFWD(hd,24,25,cstation)
	close(24)
	close(25)

close(9)

contains

!==========
subroutine SETCPFW(hd, flobspos, fmp, fwd, nsource)

	use prm_fwd

	implicit none

	integer, intent(out) :: nsource
	character(len=24), intent(out) :: hd, flobspos, fmp, fwd

	integer :: ls, lsin, kk, k
!	character(len=24) :: cal

	read(9,*) alat0, alng0, icord
	read(9,*)
	read(9,'(A)') flobspos
	read(9,*)
	read(9,*) nsource

	do ls = 1, nsource

		read(9,*)
		read(9,*) lsin
			if(lsin /= ls) stop 'LS INPUT ERROR'
		read(9,*)
		read(9,'(A)')
		read(9,*)
		read(9,*)
		read(9,*)
		read(9,*)

		read(9,*)
		read(9,'(A)') hd_jac(ls)
		read(9,*)
		read(9,'(A)') fwdjac(ls)
		read(9,*)
		read(9,*)
		read(9,*) js0(ls), jt0(ls), jd2(ls), ndgs(ls)

	end do

	read(9,*)
	read(9,'(A)') hd

	do ls = 1, nsource

		read(9,*)
		read(9,*) lsin
			if(lsin /= ls) stop 'LS INPUT ERROR'
		read(9,'(A)')
		read(9,*) just(1,ls), just(2,ls), just(3,ls), just(4,ls)
		read(9,*)
		read(9,*)
		read(9,*)
		read(9,*)

		kk = 0

		do k = 1, 4
			if(just(k,ls) == 0) cycle
			kk = kk + 1
			just(k,ls) = kk
		end do

	end do

	read(9,*)
	read(9,'(A)') !cal
	read(9,'(A)') fmp
	read(9,'(A)') fwd
	read(9,*)
	read(9,*)
	read(9,*)

		read(9,*)
	read(9,*)
	read(9,*)
	read(9,*)
	read(9,*)
	read(9,*) e_gsi, e_jcg, herr

return
end subroutine SETCPFW

!==========
subroutine GRDPO(iu,cstation)
	!----- stf にGPS基地局位置を入力 -----

	use prm_matrix
	use prm_fwd

	implicit none

	integer, intent(in) :: iu
	character(len=6), intent(out) :: cstation(3000)

	integer :: i, ii
	real(8) :: ala, alg

	read(iu,*) isf, ihf, ivf

	if(ihf+ivf /= isf) stop 'ERROR'
	if(isf > kif) stop 'STATION DIMENSION ERROR'

	do i = 1, isf
!		read(iu,'(I5,1X,F11.8,1X,F11.7,1X,A6)') ii, ala, alg, cstation(ii)
		read(iu,'(I7,1X,F12.5,1X,F12.5,4X,A6)') ii, ala, alg, cstation(ii)
		if(ii /= i) stop 'INPUT DATA ERROR'

		stf(1,i) = ala
		stf(2,i) = alg

	end do

return
end subroutine GRDPO

!==========
subroutine SETUPFW(flobspos, nsource, iu)

	use prm_matrix
	use prm_fwd

	implicit none

	integer, intent(in) :: iu, nsource
	character(len=24), intent(in) :: flobspos

	integer :: jsf, lsource, ls, nbi, i, ild, ir, j, j1, ja, jt, k
	real(8) :: ztem(kfg)
	character(len=24) :: ca, c1, c2, c3

	jsf = 0

	do lsource = 1, nsource

		ls = lsource
		nbi = (js0(ls) + ndgs(ls)) * (jt0(ls) + ndgs(ls))

		print *, '------------------------------------------'
		print *, ' SETUP JACOBI MTRX IN FWD. LSOURCE= ', ls,' NBI=', nbi

		open(iu,file=fwdjac(ls),status='old',form='unformatted')

	!======= CALL INFCJ(hd_jac(ls),flobspos,just(1,ls),nbi,jsf,ls,iu)

		ca = '**** JACOBI MATRIX **** '

		if(nbi > kfg) stop ' INSUFFICIENT ARRAY DIM(KFG)'

		read(iu) c1, c2, c3
		if(c1 /= ca     ) stop ' FILE READ ERROR (INFCJ)'
		if(c2 /= hd_jac(ls)) stop ' FILE READ ERROR (INFCJ)'
		if(c3 /= flobspos) stop ' FILE READ ERROR (INFCJ)'

		print *, ' JUST=', (just(k,ls), k=1,4)

		do i = 1, ihf
			do ild = 1, 2
				ir = 2*i -2 +ild
				do j = 1, 4
					read(iu) (ztem(j1), j1=1,nbi)
					if(just(j,ls) == 0) cycle

					ja = nbi * (just(j,ls) - 1) + jsf

					do j1 = 1, nbi
						zf(ir,ja+j1) = ztem(j1)
					end do

				end do
			end do
		end do

		do i = 1, ivf
			ir = i + 2*ihf

			do j = 1, 4
				read(iu) (ztem(j1), j1=1,nbi)
				if(just(j,ls) == 0) cycle

				ja = nbi * (just(j,ls) - 1) + jsf

				do j1 = 1, nbi
					zf(ir,ja+j1) = ztem(j1)
				end do

			end do
		end do

		do j = 1, 4
			if(just(j,ls) == 0) cycle
			jt = nbi * just(j,ls)
		end do

		jsf = jsf + jt

		if(jsf > kmp) stop ' INSUFFICIENT ARRAY DIM(KMP)'
		print *, ' READ JACOBI MATRIX. OK'

	!======= /INFCJ/

		close(iu)

	end do

return
end subroutine SETUPFW

!==========
subroutine INSOL(hd, nsource, iu) ! fmp

	use prm_matrix
	use prm_fwd

	implicit none

	integer, intent(in) :: nsource, iu
	character(len=24), intent(in) :: hd

	integer :: ls, lsa, js1, jt1, ndg1, is, it, kss, k, j
	integer :: jst1, jst2, jst3, jst4, is0, is1, it0, it1
	real(8) :: xr, error
	character(len=24) :: c1, c2
	character(len=3 ) :: chead

	er(1:kfg, 1:4, 1:5) = 0.d0
	x (1:kfg, 1:4, 1:5) = 0.d0

	read(iu,'(A3)') chead

	if(chead == ' KS') then

		print *, "chead = ' KS'"

		read(iu,*) js1, jt1, ndg1
			if(ndg1 /= ndgs(1)) stop 'INPUT ERROR NDG'
		kss = js0(1) + ndgs(1)
		read(iu,*)

		do

			k  = 0
			is = 0
			it = 0

			read(iu,*,err=12) k, is, it, xr, error

			if(k == 0 .and. is == 0 .and. it == 0) exit
			if(is == is1 .and. it == it1) exit

			j = is + (it-1) * kss
			x(j,k,1) = xr
			er(j,k,1) = error

		end do

	12	BACKSPACE(iu)

		return
	end if

! ----------

	read(iu,'(1X,A24////)') c1
	if(c1 /= hd) stop 'HD INPUT ERROR'

	do ls = 1, nsource

		read(iu,'(I5,5X,A24)') lsa, c2
		read(iu,'(3I5)'  ) js1, jt1, ndg1
		read(iu,'(4I3/)' ) jst1, jst2, jst3, jst4
		read(iu,'(4I3//)') is0, is1, it0, it1

		if(lsa /= ls     ) stop 'INPUT ERROR LS'
		if(c2  /= hd_jac(ls)) stop 'INPUT ERROR LS'

		if(js1 /= js0(ls)) stop 'INPUT ERROR JS'
		if(jt1 /= jt0(ls)) stop 'INPUT ERROR JT'

		if(ndg1 /= ndgs(ls)) stop 'INPUT ERROR NDG'

		if(jst1 /= just(1,ls)) stop 'jst1 INPUT ERROR'
		if(jst2 /= just(2,ls)) stop 'jst2 INPUT ERROR'
		if(jst3 /= just(3,ls)) stop 'jst3 INPUT ERROR'
		if(jst4 /= just(4,ls)) stop 'jst4 INPUT ERROR'

		kss = js0(ls) + ndgs(ls)

		do

			k  = 0
			is = 0
			it = 0

			read(iu,*,err=11) k, is, it, xr, error

			if(k == 0 .and. is == 0 .and. it == 0) exit
			if(is == is1 .and. it == it1) exit

			j = is + (it-1) * kss
			x (j,k,ls) = xr
			er(j,k,ls) = error

		end do

	11	BACKSPACE(iu)

	end do

return
end subroutine INSOL

!==========
subroutine FORWARD(nsource)

	use prm_matrix
	use prm_fwd

	implicit none

	integer, intent(in) :: nsource
	integer :: itot, i, jr, ls, it, is, kss, ktt, j, k

	itot = 2 * ihf + ivf

	do i = 1, itot

		jr = 0
		yf(i) = 0.d0

		do ls = 1, nsource

			kss = js0(ls) + ndgs(ls)
			ktt = jt0(ls) + ndgs(ls)

			do k = 1, 4

				if(just(k,ls) == 0) cycle

				do it = 1, ktt
					do is = 1, kss
						j = is + (it-1) * kss
						jr = jr + 1
						yf(i) = yf(i) + zf(i,jr) * x(j,k,ls)
					end do
				end do

			end do
		end do

	end do

return
end subroutine FORWARD

!==========
subroutine OUTFWD(hd,iu,iuxyz,cstation)

	use prm_matrix
	use prm_fwd

	implicit none

	integer, intent(in) :: iu, iuxyz
	character(len=24), intent(in) :: hd
	character(len=6),  intent(in) :: cstation(3000)

	integer :: n, i, ii, ier
	real(8) :: vec_th, vec_r, error

	integer, parameter :: nrnd =2400
	integer :: irnd
	real(8) :: rnd(nrnd*2)

	n = 1
	write(iu,*) isf, ihf, ivf, 0  !hd

	ier = 0
	if(e_gsi /= 0.d0) ier = ier + 1
	if(e_jcg /= 0.d0) ier = ier + 1

	if(ier == 1) stop 'CTP E_GSI or E_JCG ERROR'

	error = 0.d0

	call NRMRND(nrnd,rnd(1:nrnd),rnd(nrnd+1:nrnd+nrnd))

	irnd = 1

	do i = 1, ihf
		ii = (i-1) * 2


		if(ier == 2) then
			error = e_jcg / 1000.d0
			if(cstation(i) == "gsi") error = e_gsi / 1000.d0

			!print *, cstation(i), error


			yf(ii+1) = yf(ii+1) + error * rnd(irnd)
				irnd = irnd + 1
			yf(ii+2) = yf(ii+2) + error * rnd(irnd)
				irnd = irnd + 1

		end if

		write(iu,'(I7,4(1X,F12.5),2(1X,I3),1X,F12.5)') &
			n, stf(1,i), stf(2,i), yf(ii+1), yf(ii+2), 1, int(error*1000.d0), error*1000.d0

		vec_th = datan2(yf(ii+2),yf(ii+1)) * 180.d0 / pi
		vec_r  = dsqrt(yf(ii+1)*yf(ii+1) + yf(ii+2)*yf(ii+2))

		write(iuxyz,'(4F11.3)') stf(2,i), stf(1,i), vec_th, vec_r
		n = n + 1
	end do

	do i = ihf+1, isf
		ii = i + ihf

		if(ier == 2) then
			error = e_jcg / 1000.d0 * herr
			if(cstation(i-ihf) == "gsi") error = e_gsi / 1000.d0 * herr

			yf(ii) = yf(ii) + error * rnd(irnd)
				irnd = irnd + 1
		end if

		write(iu,'(I7,3(1X,F12.5),2(1X,I3),1X,F12.5)') &
			n, stf(1,i), stf(2,i), yf(ii), 1, int(error*1000.d0), error*1000.d0

		vec_r = abs(yf(ii))
		if(vec_r == yf(ii)) vec_th =  90.d0
		if(vec_r /= yf(ii)) vec_th = -90.d0

		write(iuxyz,'(4F11.3)') stf(2,i), stf(1,i), vec_th, vec_r
		n = n + 1
	end do

return
end subroutine OUTFWD
!==========

subroutine NRMRND(n, nrand1, nrand2)

implicit none

	integer, intent(in) :: n
	real(8), intent(out) :: nrand1(1:n), nrand2(1:n)

	real(8), parameter :: pi = acos(-1.d0)
	integer :: i

	real(8) :: mean1, sigma1, mean2, sigma2

	integer :: clock1, clock2 ! クロック取得用
	integer :: nRand ! 乱数シードのサイズ
	integer, allocatable :: seed(:) ! 乱数シード格納用
	real(8) :: rand1(n), rand2(n) ! 乱数格納用配列

	call random_seed(size=nRand)
	allocate(seed(nRand))

	call system_clock(count=clock1)
	seed = clock1
	call random_seed(put=(/seed/))
	call random_number(rand1)

	do
		call system_clock(count=clock2)
		if(clock2 == clock1) cycle
		seed = clock2
		call random_seed(put=(/seed/))
		call random_number(rand2)
		exit
	end do

!=== Box and Mueller ===

	mean1 = 0.d0
	mean2 = 0.d0
	sigma1 = 0.d0
	sigma2 = 0.d0

	do i = 1, n
		nrand1(i) = dsqrt(-2.d0*log(rand1(i))) * cos(2.d0*pi*rand2(i))
		nrand2(i) = dsqrt(-2.d0*log(rand1(i))) * sin(2.d0*pi*rand2(i))

		mean1 = mean1 + nrand1(i)/dble(n)
		mean2 = mean2 + nrand2(i)/dble(n)

	end do

	do i = 1, n
		sigma1 = sigma1 + (nrand1(i)-mean1)**2/dble(n)
		sigma2 = sigma2 + (nrand2(i)-mean2)**2/dble(n)
	end do

	sigma1 = dsqrt(sigma1)
	sigma2 = dsqrt(sigma2)

	print *, "=== MAKE RANDOM ==="
		print *, "    mean1, sigma1", real(mean1), real(sigma1)
		print *, "    mean2, sigma2", real(mean2), real(sigma2)

end subroutine NRMRND


end
