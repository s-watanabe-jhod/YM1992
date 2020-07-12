#コンパイラ
FC = gfortran

#使用ソースファイル(*.o)
O_INV2D = prm_matrix.o prm_var.o m_ctpin.o m_trans.o m_starea.o m_sttrans.o m_bspline.o m_wbb.o m_qrdec.o m_rf2d.o m_inv2d.o inv2d.o

O_INTP2D = prm_matrix.o prm_var.o m_ctpin.o m_trans.o m_starea.o m_sttrans.o m_bspline.o m_wbb.o m_intdep.o m_inp2dp.o intp2d.o

O_GDJAC = prm_matrix.o prm_var.o m_mkkn.o m_trans.o m_ctpin.o m_starea.o m_sttrans.o m_bspline.o m_wbb.o m_intdep.o m_mnewt.o m_subblms.o m_inp2dp.o m_jaclib.o m_gdjac.o gdjac.o

O_GDINV = prm_matrix.o prm_var.o prm_inv.o m_trans.o m_ctpinv.o m_starea.o m_sttrans.o m_bspline.o m_wbb.o m_intdep.o m_mnewt.o m_subblms.o m_inp2dp.o m_qrdec.o m_nnls.o m_subsetup.o m_setinv.o m_subabic.o m_sabicm.o m_inv.o m_subsrc.o m_source.o m_dissou.o gdinv.o

O_GDJACFW = prm_matrix.o prm_var.o m_mkkn.o m_trans.o m_ctpin.o m_starea.o m_sttrans.o m_bspline.o m_wbb.o m_intdep.o m_mnewt.o m_subblms.o m_inp2dp.o m_jaclib.o m_gdjac.o gdjacfw.o

O_GDFWD = prm_matrix.o prm_fwd.o gdfwd.o

O_MKSCD = prm_matrix.o prm_var.o prm_inv.o m_trans.o m_ctpinv.o m_starea.o m_sttrans.o m_bspline.o m_wbb.o m_intdep.o m_mnewt.o m_subblms.o m_inp2dp.o m_qrdec.o m_nnls.o m_subsetup.o m_setinv.o m_subsrc.o m_source.o m_dissou.o mkscd.o

#コンパイルオプション
# -Wall	: コンパイル時に全ての警告メッセージを有効
# -O2	: 初期化されていない変数を検出
FCFLAGS = -Waliasing -Wampersand -Wconversion -Wsurprising -Wline-truncation -Wunused
# FCFLAGS = -Wall -O2
# SMEMORY = -Wl,--stack,4194304

#サフィックス
.SUFFIXES: .f90

#アウトプット
inv2d: $(O_INV2D)
		$(FC) $(FCFLAGS) $(SMEMORY) -o ../exe/inv2d.exe $(O_INV2D)

intp2d: $(O_INTP2D)
		$(FC) $(FCFLAGS) $(SMEMORY) -o ../exe/intp2d.exe $(O_INTP2D)

gdjac: $(O_GDJAC)
		$(FC) $(FCFLAGS) $(SMEMORY) -o ../exe/gdjac.exe $(O_GDJAC)

gdinv: $(O_GDINV)
		$(FC) $(FCFLAGS) $(SMEMORY) -o ../exe/gdinv.exe $(O_GDINV)

gdjacfw: $(O_GDJACFW)
		$(FC) $(FCFLAGS) $(SMEMORY) -o ../exe/gdjacfw.exe $(O_GDJACFW)

gdfwd: $(O_GDFWD)
		$(FC) $(FCFLAGS) $(SMEMORY) -o ../exe/gdfwd.exe $(O_GDFWD)

mkscd: $(O_MKSCD)
		$(FC) $(FCFLAGS) $(SMEMORY) -o ../exe/mkscd.exe $(O_MKSCD)

#メインプログラム
inv2d.o:		inv2d.f90
			$(FC) $(FCFLAGS) -c inv2d.f90

intp2d.o:		intp2d.f90
			$(FC) $(FCFLAGS) -c intp2d.f90

gdjac.o:		gdjac.f90
			$(FC) $(FCFLAGS) -c gdjac.f90

gdinv.o:		gdinv.f90
			$(FC) $(FCFLAGS) -c gdinv.f90

gdjacfw.o:		gdjacfw.f90
			$(FC) $(FCFLAGS) -c gdjacfw.f90

gdfwd.o:		gdfwd.f90
			$(FC) $(FCFLAGS) -c gdfwd.f90

mkscd.o:		mkscd.f90
			$(FC) $(FCFLAGS) -c mkscd.f90

#パラメータファイル
prm_matrix.o:	prm_matrix.f90
			$(FC) $(FCFLAGS) -c prm_matrix.f90
prm_var.o: prm_var.f90
			$(FC) $(FCFLAGS) -c prm_var.f90
prm_inv.o: prm_inv.f90
			$(FC) $(FCFLAGS) -c prm_inv.f90
prm_fwd.o: prm_fwd.f90
			$(FC) $(FCFLAGS) -c prm_fwd.f90

#モジュールのソースたち
m_ctpin.o:		m_ctpin.f90
			$(FC) $(FCFLAGS) -c m_ctpin.f90

m_ctpinv.o:		m_ctpinv.f90
			$(FC) $(FCFLAGS) -c m_ctpinv.f90

m_sttrans.o:	m_sttrans.f90
			$(FC) $(FCFLAGS) -c m_sttrans.f90

m_starea.o:		m_starea.f90
			$(FC) $(FCFLAGS) -c m_starea.f90

m_trans.o:		m_trans.f90
			$(FC) $(FCFLAGS) -c m_trans.f90

m_inv2d.o:		m_inv2d.f90
			$(FC) $(FCFLAGS) -c m_inv2d.f90

m_wbb.o:		m_wbb.f90
			$(FC) $(FCFLAGS) -c m_wbb.f90

m_bspline.o:	m_bspline.f90
			$(FC) $(FCFLAGS) -c m_bspline.f90

m_rf2d.o:		m_rf2d.f90
			$(FC) $(FCFLAGS) -c m_rf2d.f90

m_qrdec.o:		m_qrdec.f90
			$(FC) $(FCFLAGS) -c m_qrdec.f90

m_inp2dp.o:		m_inp2dp.f90
			$(FC) $(FCFLAGS) -c m_inp2dp.f90

m_intdep.o:		m_intdep.f90
			$(FC) $(FCFLAGS) -c m_intdep.f90

m_mkkn.o:		m_mkkn.f90
			$(FC) $(FCFLAGS) -c m_mkkn.f90

m_mnewt.o:		m_mnewt.f90
			$(FC) $(FCFLAGS) -c m_mnewt.f90

m_subblms.o:		m_subblms.f90
			$(FC) $(FCFLAGS) -c m_subblms.f90

m_jaclib.o:		m_jaclib.f90
			$(FC) $(FCFLAGS) -c m_jaclib.f90

m_gdjac.o:		m_gdjac.f90
			$(FC) $(FCFLAGS) -c m_gdjac.f90

m_setinv.o:		m_setinv.f90
			$(FC) $(FCFLAGS) -c m_setinv.f90

m_subsetup.o:	m_subsetup.f90
			$(FC) $(FCFLAGS) -c m_subsetup.f90

m_sabicm.o:		m_sabicm.f90
			$(FC) $(FCFLAGS) -c m_sabicm.f90

m_subabic.o:	m_subabic.f90
			$(FC) $(FCFLAGS) -c m_subabic.f90

m_nnls.o:		m_nnls.f90
			$(FC) $(FCFLAGS) -c m_nnls.f90

m_dissou.o:		m_dissou.f90
			$(FC) $(FCFLAGS) -c m_dissou.f90

m_source.o:		m_source.f90
			$(FC) $(FCFLAGS) -c m_source.f90

m_subsrc.o:		m_subsrc.f90
			$(FC) $(FCFLAGS) -c m_subsrc.f90

m_inv.o:		m_inv.f90
			$(FC) $(FCFLAGS) -c m_inv.f90

#削除ツール
clean:
			rm -v *.mod
			rm -v *.o
