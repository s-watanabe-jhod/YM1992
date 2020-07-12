!=================================================
!	2 dimensional data fitting
!=================================================

!<INPUT>
!	.ctp : control parameter file
!	FLF = *.flt :  断層面記述データファイル

!<OUTPUT>
!	fcof = *.cof : 
!

program INV2D

use prm_var
use mod_ctpin
use mod_starea
use mod_sttrans
use mod_inv2d

implicit none

integer :: i, ls, iu, ipl, nsource
character(LEN=24) :: flf, fcof, fctp="                        "

open(77,file="./inv2d.out")
call getarg(1,fctp)

open(9,file=fctp)

CALL SETA2D(nsource)	! ctpin.f
!原点と断層数をパラメータファイルから読み込む

do i = 1, nsource
	ls = i
	iu = 50 + i

	CALL SETFC2D(ls,flf,fcof,ipl)	! ctpin.f
		! (ipl=1) ==> MODEL FAULT IS FLAT 
		!各断層の位置・グリッドサイズなど情報を取得

	if(ipl == 0) then

		CALL STTRANS	! blms.f
			!断層座標への変換係数導出（u11, u12, t11, t12）
		CALL STAREA		! blms.f
			!断層座標での矩形断層領域位置の導出（xa-xb, ya-yb）

		open(11,file=flf,status="old")
			CALL INSURD(11)		! this file
		close(11)
			!gx,gy,gzに、断層座標系での断層位置（深さ）情報を代入

		CALL SETKN		! this file
		CALL SETBM2D	! this file
		CALL SABIC2D	! this file
			!断層面をb-splineで表す際のパラメータを決定

		open(iu,file=fcof)
			CALL OUT2DP(flf,iu)	! this file
		close(iu)
			! ？？jmt 個のセグメントの値を出力？？

	end if

end do

close(9)
close(77)

end program
