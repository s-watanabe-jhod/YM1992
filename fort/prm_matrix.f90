
!===== パラメータを決めるモジュール ====

module prm_matrix

	implicit none

	integer, parameter :: ki = 5400, ks = 30, kt = 30, ku = 20, kv = 20, nd = 100
	integer, parameter :: kis = 4000, kif = 4000
	integer, parameter :: kfg = 4000, kmp = 4000, kfo = 4000, kio = 4000
	integer, parameter :: kfe = (ks+3) * (kt+3), kmr = (ku+3) * (kv+3)
	integer, parameter :: kmpa = kmp+1, kdata = kmp + kio
	integer, parameter :: kdata2 = kmp + ki
	integer, parameter :: ms = ks*nd, mt = kt*nd, ndd = nd*nd

end module prm_matrix

! ki 断層面を決める節点の数

