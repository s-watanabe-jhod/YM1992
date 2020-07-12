!===============================
! DEFINING variables
!  except for prm_var
!===============================

module prm_fwd

	use prm_matrix

	!== ESTMPR ==
	real(8) x(kfg,4,50), er(kfg,4,50)
	!== FORWDI ==
	character(len=24), save :: fwdjac(50), hd_jac(50)
	!== FWPRM ==
	integer, save :: isf, ihf, ivf
	real(8), save :: zf(kfo,kmp), yf(kfo), stf(2,kif)
	!== JUSTIN ==
	integer :: ndgs(50), js0(50), jt0(50), jd2(50), just(4,50)
	!== OCOD ==
	real(8), save :: alat0, alng0
	integer, save :: icord
	real(8), save :: xtr(ms), ytr(ms,mt), ztr(ms,mt), da2

	real(8), save :: e_gsi, e_jcg, herr, setalpha

	real(8), parameter :: pi = 3.14159265358979324d0

end module prm_fwd
