!===============================
! DEFINING variables
!===============================

module prm_var

	use prm_matrix

	implicit none

	!== ABICPR ==
	real(8), save :: alpha, rsq2, sigma ! 
	integer, save :: jmt, jct, jdata ! 
	!== CBSPLN ==
	integer, save :: nbdeg, nbod2, ks0, kt0, nd2, ndd2, ms2, mt2 ! 
	!== CDSYST ==
	real(8), save :: aphi, adlt, alat, alon, awid, alen ! 断層のパラメータ（ctpファイル）
	!== CONSTR2D ==
	real(8), save :: cs2d(kmr,kmr), detca ! 
	integer, save :: jcp ! 境界にないbsplineのセグメント数？？
	!== CONTID ==
	integer, save :: isn, ih, iv
	!== ESTIMA2D ==
	real(8), save :: x2d(kmr) ! 
	integer, save :: lu2d(kmr) ! 
	!== ELASTIC ==
	real(8), save :: delt, fde !ELASTIC CONSTANT delt = LAMBDA/MU
	!== ETCGRN ==
	real(8), save :: bi(ndd,16) ! bi bspline
	!== FAULTR ==
	real(8), save :: alaf, alof, awif, alef, adep ! 断層面上の基底関数を定義するもの（らしい）
	!== FLTBSP ==
	integer, save :: ndeg, nod2, ku0, kv0, ku1, kv1 ! bspline次数。nod2=(ndeg +1)**2。bsplineのセグメント数。左+次数。
	!== JACOBF ==
	real(8), save :: ds1(kfe,4), ds2(kfe,4)
	integer, save :: nbi
	!== JACOBI ==
	real(8), save :: rt(ki,kmr) ! 
	integer, save :: itotal, jmp  ! 有効（モデル断層面内に存在する）断層面座標の個数。jmp = ku1*kv1

	!== MESHPO ==
	real(8), save :: xtr(ms), ytr(ms,mt), ztr(ms,mt), da2
	!== MSAREA ==
	real(8), save :: xfa, xfb, yfa, yfb, zda, zdb ! 断層の定義領域の最大最小座標値
	!== OCOD ==
	real(8), save :: alat0, alng0 ! ローカル座標の原点
	integer, save :: icord ! ctpの座標表示が緯度経度かメートルか
	!== PARAME ==
	real(8), save :: xi(kmr)
	!== PRAREA ==
	real(8), save :: xa, xb, ya, yb ! 矩形断層領域の対角の隅
	!== STATIO ==
	real(8), save :: st(2,kis)
	!== SURFGM ==
	real(8), save :: gx(ki), gy(ki), gz(ki) ! 断層面の形状の座標値
	integer, save :: isurf ! 断層面の座標値の数
	!== USEARE2d ==
	real(8), save :: z2d(kdata2,kmr), y2d(kdata2) ! 
	!== XYTRAN ==
	real(8), save :: t11, t12, u11, u12 ! 方位座標 - 断層座標　間の座標回転の係数
	!== WORKAR ==
	real(8), save :: gf1(kfe,16,4), gf2(kfe,16,4), gs1(ndd,4), gs2(ndd,4)
	!====

	!== GPS POSITION ==
	real(8), save :: glatmin, glatmax, glonmin, glonmax, e_gsi, e_jcg, setalpha
	character(len=24), save :: fpos, fgps, fdum
	!====

	real(8), parameter :: pi = 3.14159265358979324d0
	real(8), parameter :: ground = -1.d0 ! この深さ以上で計算する

end module prm_var

