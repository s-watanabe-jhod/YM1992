program gdjac

use prm_matrix
use prm_inv

use mod_ctpinv
use mod_gdjac

use mod_blms
use mod_bibsp
!use mod_sub_flatfault
use mod_setup
use mod_setwt
use mod_sabicm
use mod_inv
use mod_dissou

!====================
implicit none


integer :: lsource, ls, iu
integer :: nsource, i, j, k, lname

character(len=24) :: flobspos, flobspd, cal, fmp, hd
character(len=24) :: fctp  = "                        "
character(len=24) :: fcalp = "                        "
character(len=24) :: fobsp = "                        "
character(len=24) :: fo_cp = "                        "

integer :: kmri, kui, kvi
real(8) :: restmp

!====================
call getarg(1,fctp)

open(77,file='../tmp/gdinv.out')

open(9,file=fctp)
	CALL SET_CTP(hd,flobspd,flobspos,cal,fmp,nsource)
close(9)

open(33,file=flobspd,status='old')
	CALL OBSDATAIN(33)
close(33)

!===== INITIAL in GDINV =====

Jsp    = 0
Jcp    = 0
Idmp   = 0 !idmp(1:kfg,1:4,1:5) = 0
Detca  = 0.d0
Cs     = 0.d0 !cs(1:kmp,1:kmp) = 0.d0
Z_Green= 0.d0
!Gi     = 0.d0 !gi(1:kmp,1:kmp) = 0.d0
!do i = 1, Kmp
!	Gi(i,i) = 1.d0
!end do

!========

do lsource = 1, nsource

	ls = lsource
	print *, ' ---------------------------------------------'
	print *, ' START OF JACOBI MATRIX CALCULATION. SOURCE=', ls

	!==== BLMS 断層面の基底関数の設定 ====

	if(Ipl(ls) == 0) then ! 断層面が曲面

		CALL STAREA(ls, Xa(ls), Ya(ls), Xb(ls), Yb(ls))     ! set Xa(ls), Ya(ls), Xb(ls), Yb(ls)
		CALL STMSAR(ls, Xfa(ls), Xfb(ls), Yfa(ls), Yfb(ls)) ! set Xfa(ls), Xfb(ls), Yfa(ls), Yfb(ls)

		open(31, file=Fcof(ls), status='old')
			Xcof = 0.d0
			read(31,'(//1X,3(I5,1X)//)') kmri, kui, kvi     ! 断層形状の読み込み from COF file
			read(31,'(5(E12.6,1X))') (Xcof(i), i = 1, kmri) !
		close(31)

		CALL MESH(ls) ! set Xtr(Ms), Ytr(Ms,Mt), Ztr(Ms,Mt), Da2(ls)

		open(29, file='xyz.dat')
		do i = 1, Js0(ls)*Jd2(ls)
		 do j = 1, Jt0(ls)*Jd2(ls)
		  write(29, *) i, j, Xtr(i), Ytr(i,j), Ztr(i,j)
		 end do
		end do
		close(29)

	else if(Ipl(ls) /= 0) then ! 断層面が平面
		STOP "!!!Not in service for FLAT FAULT!!!"
		! CALL STMSAR2
		! CALL MESHSS
	end if

	!==== /BLMS ====

	!==== BIBSP すべり分布の基底関数の設定 ====

	CALL BIBSP(Jd2(ls), Ndgs(ls), Jd2(ls)*Jd2(ls), Nbod2(ls), Bi) ! BI BSPLINE (DEGREE=NBDEG)

	!==== /BIBSP ====

	iu = 40 + lsource
	open(iu,file=Invjac(ls),form='unformatted')

	!==== OUTPAR ====
	print '(5X,"**** JACOBI MATRIX **** ",/" ",A24,/," ","SOURCE NUMBER= ",I5)', Hd_jac(ls), ls
	print '(1X,3(I5,1X),A24," (ISN,IH,IV,OBS DATA FILE)")', Isn, Ih, Iv, flobspd
	print '(1X,3(I5,1X),2X,"SOURCE (BSPL DEG ,KS ,KT)")', Ndgs(ls), Js0(ls), Jt0(ls)
	print *, '  MEAN SLIP DIRECTION OUT'
	write(iu) "**** JACOBI MATRIX **** ", Hd_jac(ls), flobspd
	!==== /OUTPAR ====

	! できればなくしたいinput変数たち (JACLIB で使用している)
	Ks0 = Js0(ls)
	Kt0 = Jt0(ls)
	! ここまで

	CALL JACOBI(iu, ls) ! 未解読（サブ）

	close(iu)

	!============================
	!  Fin   GDJAC
	!============================

	!============================
	!  Start GDINV
	!============================

	!Nbi = (Js0(ls) + Ndgs(ls)) * (Jt0(ls) + Ndgs(ls)) ! Nbi を再定義

	print *,    '------------------------------------------'
	print *,    ' SETUP JACOBI MTRX. LSOURCE= ', ls, ' NBI=', Nbi
	write(77,*) ' SETUP JACOBI MTRX. LSOURCE= ', ls, ' NBI=', Nbi

	! INPUTG は JACOBI に included

	!open(32,file=Invjac(ls),status='old',form='unformatted')
	!CALL INPUTG(Hd_jac(ls), flobspd, ls, 32)
	!close(32)

	CALL SETMP(ls, Js0(ls), Jt0(ls), Ndgs(ls)) ! 未解読
	CALL SETBM(ls, Js0(ls), Jt0(ls), Ndgs(ls)) ! 未解読

	print *, '------------------------------------------'

end do

CALL SETWT

CALL SABICM


! === CALL DRESID ===

Y_calc(1:Itotal) = 0.d0

do j = 1, Jsp
	do i = 1, Itotal
		Y_calc(i) = Y_calc(i) + Z_Green(i,j) * X_model(j)
		! モデルパラメータから推定される地殻変動データ
	end do
end do

! SETWT でつけたデータの重みを元に戻す
do i = 1, Itotal
	Y_calc(i) = Y_calc(i) * dabs(E_prior(i))
	Y_data(i) = Y_data(i) * dabs(E_prior(i))
end do

! === /DRESID ===

lname = len_trim(cal)
fcalp(1:lname-4) = cal(1:lname-4)
fobsp(1:lname-4) = cal(1:lname-4)
fo_cp(1:lname-4) = cal(1:lname-4)
fcalp(lname-3:lname+4) = "_cal.dat"
fobsp(lname-3:lname+4) = "_obs.dat"
fo_cp(lname-3:lname+4) = "_o_c.dat"

open(60,file=fcalp)
open(61,file=fobsp)
open(62,file=fo_cp)

open(21,file=cal,status='unknown')
	CALL OUTDAT(hd)
close(21)

close(60)
close(61)
close(62)

CALL COVAR(Sigma)

open(79,file="Resol.dat")

Resol  = 0.d0

do i = 1, Jsp
	do j = 1, Jsp
		restmp = 0.d0
		do k = 1, Jsp
			Resol(i,j) = Resol(i,j) + Cov(i,k) * HTH(j,k) ! Cov * HTH いずれも対称行列
		end do
		Resol(i,j) = Resol(i,j)  / Sigma / Sigma ! Resolution matrix 対称行列でないことに注意 (HT*H + alpha^2*G)^-1 * HT*H
	end do
	write(79,*) Resol(i,i) !(Resol(i,j), j=1,Jsp)
	if(Resol(i,i) < 0.d0) print *, "Negative Res(i,i) = ", Resol(i,i), " i = ", i
end do

!do i = 1, Jsp
!	do k = 1, Jdata ! Jdata まで足し合わせる
!		do j = 1, Jsp
!			tmp(i,k) = tmp(i,k) + tmpcov(i,j) * Z_Green(k,j) ! (HT*H + alpha^2*G)^-1 * HT
!		end do
!	end do
!end do
!
!do i = 1, Jsp
!	do j = 1, Jsp
!		do k = 1, Jdata
!			Resol(i,j) = Resol(i,j) + tmp(i,k) * Z_Green(k,j) ! (HT*H + alpha^2*G)^-1 * HT*H
!		end do
!	end do
!
!	if(Cov(i,i)   < 0.d0) print *, "Negative Cov(i,i) = ", Cov(i,i),   " i = ", i
!end do

close(79)

open(22,file=fmp,status='unknown')
	CALL OUTSOL(hd,flobspd,nsource)
close(22)

CALL DISSOU(hd,nsource)

close(77)

end
