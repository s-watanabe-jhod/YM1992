!===============================
! DEFINING variables
!  except for prm_var
!===============================

module prm_fwd

	use prm_matrix

	!== ESTMPR ==
	real(8) x(kfg,4,5), er(kfg,4,5)
	!== FORWDI ==
	character(len=24), save :: fwdjac(5), hd_jac(5)
	!== FWPRM ==
	integer, save :: isf, ihf, ivf
	real(8), save :: zf(kfo,kmp), yf(kfo), stf(2,kif)
	!== JUSTIN ==
	integer :: ndgs(5), js0(5), jt0(5), jd2(5), just(4,5)
	!== OCOD ==
	real(8), save :: alat0, alng0
	integer, save :: icord
	real(8), save :: xtr(ms), ytr(ms,mt), ztr(ms,mt), da2

	real(8), save :: e_gsi, e_jcg, herr, setalpha

end module prm_fwd