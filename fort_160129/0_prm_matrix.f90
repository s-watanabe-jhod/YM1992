
!===== パラメータを決めるモジュール ====

module prm_matrix

	implicit none

	integer, parameter :: Ks = 30, Kt = 30, Ku = 20, Kv = 20, Nd = 100
	integer, parameter :: Kis = 2000, Kif = 2000
	integer, parameter :: Kfg = 3000, Kmp = 3000, Kfo = 3000, Kio = 3000
	integer, parameter :: Kfe = (Ks+3) * (Kt+3), Kmr = (Ku+3) * (Kv+3)
	integer, parameter :: Kmpa = Kmp+1, Kdata = Kmp + Kio
	integer, parameter :: Ms = Ks*Nd, Mt = Kt*Nd, Ndd = Nd*Nd

end module prm_matrix

! ki 断層面を決める節点の数
