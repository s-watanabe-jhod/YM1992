!===============================
! DEFINING variables
!  prm_var で未定義のもののみ
!===============================

module prm_inv

	use prm_matrix
	use prm_var

	implicit none

	!== IDHIST ==
	character(len=24), save :: fcof(50), invjac(50), hd_jac(50), fscd(50)

	!== AAULTR ==
	real(8), save :: blaf(50), blof(50), bwif(50), blef(50), bdep(50)
	!== ABSPLN ==
	integer, save :: ndgs(50), js0(50), jt0(50), jd2(50)
	!== ADSYST ==
	real(8), save :: bphi(50), bdlt(50), blat(50), blon(50), blen(50), bwid(50)
	!== ALTBSP ==
	integer, save :: ju0(50), jv0(50), ndgf(50), ipl(50)
	!== CONSTR == + prm_var
	real(8), save :: cs(kmp,kmp) ! サイズが異なっている！ 対症療法で、cs->cs2d と使い分け
	!== DATARE ==
	integer, save :: ia(kio), ib(kio), ic(kio)
	!== ERCOVM ==
	real(8), save :: fact(3), fabs(3), vrp
	!== ESTIMA ==
	real(8), save :: x(kmp) ! 
	integer, save :: lu(kmp) ! 
	!== FREEMP ==
	real(8), save :: gi(kmp,kmp)
	integer, save :: jnp
	!== INDEXC ==
	integer, save :: just(4,50), ius(4,50), mlg(50), nnls(4,50), iboc(4,50)
	real(8), save :: rss(50), rop(50), rex(50)
	!== INEQLT ==
	real(8), save :: ha(kmp)
	integer, save :: incc(kmp), inncon
	!== MEDELR ==
	integer, save :: idmp(kfg,4,50)
	!== OBSPOS == + prm_var
	integer, save :: il !isn, ih, iv, 
		!	real(8), save :: st(2,kis)
	!== READPR == + prm_var
	real(8), save :: zz(kio,kmp), yy(kio), ee(kio)
	integer, save :: iinit, jsp
	!== RESIDU ==
	real(8), save :: ysl(kio)
	!== SLIPDR ==
	real(8), save :: slms(3,kfg,50), slmd(3,kfg,50)
	!== TRANSD ==
	real(8), save :: ag(kio,kio)
	!== USEARE ==
	real(8), save :: z(kdata,kmp), y(kdata) ! 
	!====
	
	!== SWAT 追加 ==


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

!	real(8), save ::

	!== ABICPR ==
!	real(8), save :: alpha, rsq2, sigma ! 
!	integer, save :: jmt, jct, jdata ! 
	!== CBSPLN ==
!	integer, save :: nbdeg, nbod2, ks0, kt0, nd2, ndd2, ms2, mt2 ! 
	!== CDSYST ==
!	real(8), save :: aphi, adlt, alat, alon, awid, alen ! 断層のパラメータ（ctpファイル）
	!== CONSTR2D ==
!	real(8), save :: cs2d(kmr,kmr), detca ! 
!	integer, save :: jcp ! ？淵にないbsplineのセグメント数？？
	!== CONTID ==
!	integer, save :: isn, ih, iv
	!== ESTIMA2d ==
!	real(8), save :: x2d(kmr) ! 
!	integer, save :: lu2d(kmr) ! 
	!== ELASTIC ==
!	real(8), save :: delt, fde !ELASTIC CONSTANT delt = LAMBDA/MU
	!== ETCGRN ==
!	real(8), save :: bi(ndd,16) ! bi bspline
	!== FAULTR ==
!	real(8), save :: alaf, alof, awif, alef, adep ! 断層面上の基底関数を定義するもの（らしい）
	!== FLTBSP ==
!	integer, save :: ndeg, nod2, ku0, kv0, ku1, kv1 ! bspline次数。nod2=(ndeg +1)**2。bsplineのセグメント数。左+次数。
	!== JACOBF ==
!	real(8), save :: ds1(kfe,4), ds2(kfe,4)
!	integer, save :: nbi
	!== JACOBI ==
!	real(8), save :: rt(ki,kmr) ! 
!	integer, save :: itotal, jmp  ! 有効（モデル断層面内に存在する）断層面座標の個数。jmp = ku1*kv1

	!== MESHPO ==
!	real(8), save :: xtr(ms), ytr(ms,mt), ztr(ms,mt), da2
	!== MSAREA ==
!	real(8), save :: xfa, xfb, yfa, yfb, zda, zdb ! 断層の定義領域の最大最小座標値
	!== OCOD ==
!	real(8), save :: alat0, alng0 ! ローカル座標の原点
!	integer, save :: icord ! ctpの座標表示が緯度経度かメートルか
	!== PARAME ==
!	real(8), save :: xi(kmr)
	!== PRAREA ==
!	real(8), save :: xa, xb, ya, yb ! 矩形断層領域の対角の隅
	!== STATIO ==
!	real(8), save :: st(2,kis)
	!== SURFGM ==
!	real(8), save :: gx(ki), gy(ki), gz(ki) ! 断層面の形状の座標値
!	integer, save :: isurf ! 断層面の座標値の数
	!== USEARE ==
!	real(8), save :: z2d(kdata2,kmr), y2d(kdata2) ! 
	!== XYTRAN ==
!	real(8), save :: t11, t12, u11, u12 ! 方位座標 - 断層座標　間の座標回転の係数
	!== WORKAR ==
!	real(8), save :: gf1(kfe,16,4), gf2(kfe,16,4), gs1(ndd,4), gs2(ndd,4)
	!====

end module prm_inv

