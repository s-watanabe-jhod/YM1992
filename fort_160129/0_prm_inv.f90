!===============================
! DEFINING variables
!===============================

module prm_inv

	use prm_matrix
	
	implicit none
	
	!========================
	!== ローカル座標・定数のセッティング == (SET_CTP) 
	real(8), save :: Alat0, Alng0 ! ローカル座標の原点
	integer, save :: Icord ! ctpの座標表示が緯度経度かメートルか
	real(8), save :: Delt !ELASTIC CONSTANT delt = LAMBDA/MU
	
	!== 断層のパラメータ前半（最大5枚分 b -> a とするとgdjacのパラメータ） == (SET_CTP) 
	character(len=24), save :: Fcof(5), Hd_jac(5), Invjac(5), Fwdjac(5)
	real(8), save :: Bphi(5), Bdlt(5), Blat(5), Blon(5), Blen(5), Bwid(5)
	real(8), save :: Blaf(5), Blof(5), Blef(5), Bwif(5), Bdep(5)
	integer, save :: Ju0(5), Jv0(5), Ndgf(5), Ipl(5)
	integer, save :: Js0(5), Jt0(5), Jd2(5), Ndgs(5)
	
	! 上のパラメータから計算される変数群（SET_CTP）
	integer, save :: Nod2(5), Ju1(5), Jv1(5) ! bspline次数。nod2=(ndeg +1)**2。bsplineのセグメント数。左+次数。
	integer, save :: Nbod2(5)! Ms2(5), Mt2(5) この二つは使うのやめる
	
	!== 断層のパラメータ後半（gdinv用のパラメータ） == (SET_CTP) 
	character(len=24), save :: Fscd(5)
	integer, save :: Just(4,5), Ius(4,5), Mlg(5), Iboc(4,5)
	real(8), save :: Rss(5), Rop(5), Rex(5)
	! Just: 1; Strike slip 2; Dip slip 3; Opening 4; Exp.
	! Iboc: boundary conditions
	
	!== 計算条件・出力ファイルのセッティング（gdinv用のパラメータ） == (SET_CTP) 
	! VRP: to ADD TO ALL VERTICAL DISP. DATA
	real(8), save :: Fact(3), Fabs(3), Vrp ! Vrp は，高さデータに一律で足しあげている（da(i) = da(i) + vrp）
	!== added by oldSWAT ==
	real(8), save :: Glatmin, Glatmax, Glonmin, Glonmax, E_gsi, E_jcg, Setalpha
	character(len=24), save :: Fpos, Fgps, Fdum
	
	!========================
	! OBSDATAIN
	!========================
	! station position set in OBSDATAIN
	real(8), save :: St(2,Kis)   ! station position in x, y
	integer, save :: Isn, Ih, Iv, Il ! station の数，水平，上下，レベリング，Iinit = Ih * 2 + Iv
	
	integer, save :: Itotal ! 成分ごとのデータ数（通常，Ih*2+Iv）
	integer, save :: Ia(Kio), Ib(Kio), Ic(Kio) 
	! Ia 1: x-ward, 2: y-ward, 3: z-ward) のフラグ
	! Ib Stationの通し番号（水平成分と上下成分を分けている；重みゼロがなければobsファイルの最初のカラムと同じ値のハズ）
	! Ic GPSのデータには，ゼロが入る
	
	real(8), save :: E_prior(Kio) ! a priori error    (old: ee(kio))
	real(8), save :: Y_data(Kio)  ! observation data  (old: yy(kio))
	!========================
	! SETAG
	!========================
	real(8), save :: Sta_corr(Kio,Kio) ! 通常，単位行列 (old: ag(kio,kio))
		! 固定点を置いた時に意味が出てくる（ステーション間の相関が出るということ）
		! SETWT のみで使用される
	
	
	!========================
	! BLMS
	!========================
	! set in STAREA
	real(8), save :: Xa(5), Xb(5), Ya(5), Yb(5) ! 矩形断層領域の対角の隅
	! set in STMSAR
	real(8), save :: Xfa(5), Xfb(5), Yfa(5), Yfb(5), Zda(5), Zdb(5) ! 断層の定義領域の最大最小座標値
	! input from COF file (set after STMSAR)
	real(8), save :: Xcof(Kmr) ! 断層形状の読み込み from COF file
	! set in MESH
	real(8), save :: Xtr(Ms), Ytr(Ms,Mt), Ztr(Ms,Mt), Da2(5)
	
	!========================
	! BIBSP
	!========================
	real(8), save :: Bi(ndd,16) ! bi bspline
	
	!========================
	! JACOBI
	!========================
	integer, save :: Nbi ! = (Js0(ls) + Ndgs(ls)) * (Jt0(ls) + Ndgs(ls))
	real(8), save :: Ds1(Kfe,4), Ds2(Kfe,4) ! Z_Green, Slms, Slmd を計算するためのTemp領域
	real(8), save :: Slms(3,Kfg,5), Slmd(3,Kfg,5) ! F77時代からコメントアウトされていて，使われていないようである。
	real(8), save :: Z_Green(kio,kmp) ! (old: zz(kio))
	! MKKN で使用
	real(8), save :: gf1(Kfe,16,4), gf2(Kfe,16,4), gs1(Ndd,4), gs2(Ndd,4)
	
	!========================
	! SETUP で initialize, SETMP SETBM で定義
	!========================
	integer, save :: Jsp ! = Jmt = M "MODEL PARAMETERS" in SETMP モデルパラメータの数
	integer, save :: Jcp ! = Jct = P "CONSTRAIN MATRIX" in SETBM 拘束マトリクス Cs(Jcp,Jsp) の成分の数 rank(G) 
	integer, save :: Idmp(Kfg,4,5) 
	real(8), save :: Cs(Kmp,Kmp), Cov(Kmp,Kmp), Resol(Kmp,Kmp), HTH(Kmp,Kmp) ! 最初は smoothness matrix ; COVAR の出力
	real(8), save :: Detca ! 何かの行列式の対数。プログラム中ではずっとゼロ
	!real(8), save :: Gi(Kmp,Kmp) ! 単位行列 I: 使われていない？ <- もしかするとYabuMatsuの行列 G？
	
	real(8), save :: Z_HaG(Kdata,Kmp), Y(Kdata) ! Y 使わないデータを取り除いたデータ; Z = H + alpha*G
		! QRDEC 内で Z_HaG は変換される
	
	!========================
	! SABICM
	!========================
	real(8), save :: Alpha, Rsq2, Sigma ! 
	integer, save :: Jdata ! N in YabuMatsu
	
	!== ESTIMA ==
	real(8), save :: X_model(Kmp) ! 
	real(8), save :: Y_calc(Kio)
	
	
	!=====
	! できればなくしたいinput変数たち (JACLIBで使用している)
	!=====
	integer, save :: Ks0, Kt0
	
	
	
	
	!====
	
	!real(8), parameter :: ground = -1.d0 ! この深さ以上で計算する
	
end module prm_inv

