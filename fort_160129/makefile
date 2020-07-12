#コンパイラ
FC = gfortran

#使用ソースファイル(*.o)

# O_GDJAC = prm_matrix.o prm_var.o m_mkkn.o Z_m_trans.o m_ctpin.o m_starea.o m_bspline.o m_wbb.o m_intdep.o m_mnewt.o m_subblms.o m_jaclib.o m_gdjac.o gdjac.o
O_GDJACFW = 0_prm_matrix.o prm_var.o m_mkkn.o Z_m_trans.o m_ctpin.o m_starea.o m_bspline.o m_wbb.o m_intdep.o m_mnewt.o m_subblms.o m_jaclib.o m_gdjac.o gdjacfw.o
O_GDFWD = 0_prm_matrix.o 0_prm_fwd.o gdfwd.o
O_MKSCD = 0_prm_matrix.o prm_var.o 0_prm_inv.o Z_m_trans.o m_ctpinv.o m_starea.o m_bspline.o m_wbb.o m_intdep.o m_mnewt.o m_subblms.o m_qrdec.o m_subsetup.o m_setinv.o m_subsrc.o m_source.o m_dissou.o mkscd.o


O_GDINV =  0_prm_matrix.o 0_prm_inv.o Z_m_trans.o Y2_m_bspline.o Y1_m_wbb.o X_m_qrdec.o S3_m_subsrc.o S2_m_source.o C2_m_mnewt.o C1_m_blms.o S1_m_dissou.o R1_m_inv.o Q2_m_subabic.o Q1_m_sabicm.o P1_m_setwt.o O1_m_setup.o F3_m_mkkn.o F2_m_jaclib.o F1_m_jacobi.o E1_m_bibsp.o B1_m_ctpinv.o A1_gdinv.o


#コンパイルオプション
# -Wall	: コンパイル時に全ての警告メッセージを有効
# -O2	: 初期化されていない変数を検出
#FCFLAGS = -Wall -O2
FCFLAGS = -Waliasing -Wampersand -Wconversion -Wsurprising -Wline-truncation -Wunused
SMEMORY =
#-Wl,--stack,4194304

#サフィックス
.SUFFIXES: .f90

#アウトプット

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

#gdjac.o:		gdjac.f90
#			$(FC) $(FCFLAGS) -c gdjac.f90

A1_gdinv.o:		A1_gdinv.f90
			$(FC) $(FCFLAGS) -c A1_gdinv.f90

gdjacfw.o:		gdjacfw.f90
			$(FC) $(FCFLAGS) -c gdjacfw.f90

gdfwd.o:		gdfwd.f90
			$(FC) $(FCFLAGS) -c gdfwd.f90

mkscd.o:		mkscd.f90
			$(FC) $(FCFLAGS) -c mkscd.f90

#パラメータファイル
0_prm_matrix.o:	0_prm_matrix.f90
			$(FC) $(FCFLAGS) -c 0_prm_matrix.f90
0_prm_inv.o: 0_prm_inv.f90
			$(FC) $(FCFLAGS) -c 0_prm_inv.f90
0_prm_fwd.o: 0_prm_fwd.f90
			$(FC) $(FCFLAGS) -c 0_prm_fwd.f90

#モジュールのソースたち

B1_m_ctpinv.o:		B1_m_ctpinv.f90
			$(FC) $(FCFLAGS) -c B1_m_ctpinv.f90

C1_m_blms.o:		C1_m_blms.f90
			$(FC) $(FCFLAGS) -c C1_m_blms.f90

C2_m_mnewt.o:		C2_m_mnewt.f90
			$(FC) $(FCFLAGS) -c C2_m_mnewt.f90

E1_m_bibsp.o:		E1_m_bibsp.f90
			$(FC) $(FCFLAGS) -c E1_m_bibsp.f90

F1_m_jacobi.o:		F1_m_jacobi.f90
			$(FC) $(FCFLAGS) -c F1_m_jacobi.f90

F2_m_jaclib.o:		F2_m_jaclib.f90
			$(FC) $(FCFLAGS) -c F2_m_jaclib.f90

F3_m_mkkn.o:		F3_m_mkkn.f90
			$(FC) $(FCFLAGS) -c F3_m_mkkn.f90

O1_m_setup.o:		O1_m_setup.f90
			$(FC) $(FCFLAGS) -c O1_m_setup.f90

P1_m_setwt.o:		P1_m_setwt.f90
			$(FC) $(FCFLAGS) -c P1_m_setwt.f90

Q1_m_sabicm.o:		Q1_m_sabicm.f90
			$(FC) $(FCFLAGS) -c Q1_m_sabicm.f90

Q2_m_subabic.o:		Q2_m_subabic.f90
			$(FC) $(FCFLAGS) -c Q2_m_subabic.f90

R1_m_inv.o:		R1_m_inv.f90
			$(FC) $(FCFLAGS) -c R1_m_inv.f90

S1_m_dissou.o:		S1_m_dissou.f90
			$(FC) $(FCFLAGS) -c S1_m_dissou.f90

S2_m_source.o:		S2_m_source.f90
			$(FC) $(FCFLAGS) -c S2_m_source.f90

S3_m_subsrc.o:	S3_m_subsrc.f90
			$(FC) $(FCFLAGS) -c S3_m_subsrc.f90

X_m_qrdec.o:		X_m_qrdec.f90
			$(FC) $(FCFLAGS) -c X_m_qrdec.f90

Y1_m_wbb.o:		Y1_m_wbb.f90
			$(FC) $(FCFLAGS) -c Y1_m_wbb.f90

Y2_m_bspline.o:		Y2_m_bspline.f90
			$(FC) $(FCFLAGS) -c Y2_m_bspline.f90

Z_m_trans.o:		Z_m_trans.f90
			$(FC) $(FCFLAGS) -c Z_m_trans.f90

#削除ツール
clean:
			rm -v *.mod
			rm -v *.o
