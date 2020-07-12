!=================================================
! SUBROUTINES
!	DISSOU(hd,nsource)
!
!=================================================

module mod_dissou

	implicit none

contains
!==========
subroutine DISSOU(hd,nsource)

	use prm_var
	use prm_inv

	use mod_sttrans
	use mod_starea
	use mod_inp2dp
	use mod_subblms

	use mod_source

	implicit none

	character(len=24), intent(in) :: hd
	integer, intent(in) :: nsource

	character(len=24) :: fxyz = "                        ", fls
	character(len=24) :: felp = "                        "
	character(len=24) :: fcnt = "                        "
	integer :: lsource, ls, lname

	do lsource = 1, nsource

		ls = lsource

		alaf = blaf(ls)
		alof = blof(ls)
		awif = bwif(ls)
		alef = blef(ls)
		adep = bdep(ls)

		aphi = bphi(ls)
		adlt = bdlt(ls)
		alat = blat(ls)
		alon = blon(ls)
		awid = bwid(ls)
		alen = blen(ls)

		ndeg = ndgf(ls)
		nod2 = (ndeg + 1) * (ndeg + 1)
		ku0 = ju0(ls)
		kv0 = jv0(ls)
		ku1 = ku0 + ndeg
		kv1 = kv0 + ndeg

		nbdeg = ndgs(ls)
		nbod2 = (nbdeg + 1) * (nbdeg + 1)
		ks0 = js0(ls)
		kt0 = jt0(ls)

		nd2 = jd2(ls)
		ndd2 = nd2 * nd2
		ms2 = nd2 * ks0
		mt2 = nd2 * kt0

		if(mlg(ls) == 0) cycle

		print *, ' -------------------------------------------'
		print *, ' SOURCE DISTRIBUTION.LSOURCE= ', ls

		fls = fscd(ls)
		lname = len_trim(fls)
		fxyz(1:lname-4) = fls(1:lname-4)
		felp(1:lname-4) = fls(1:lname-4)
		fcnt(1:lname-4) = fls(1:lname-4)
		fxyz(lname-3:lname+4) = "_scd.xyz"
		felp(lname-3:lname+4) = "_scd.elp"
		fcnt(lname-3:lname+4) = "_scd.cnt"

		open(ls+44,file=fxyz)
		open(ls+50,file=felp)
		open(ls+56,file=fcnt)
		open(ls+22,file=fls) !, status='unknown')
		if(ipl(ls) == 0) open(31,file=fcof(ls),status='old')

		CALL STTRANS

	!==== BLMS ífëwñ ÇÃäÓíÍä÷êîÇÃê›íËÅH ====

		if(ipl(ls) == 0) then ! ífëwñ Ç™ã»ñ 
			CALL STAREA
			CALL STMSAR
			CALL INP2DP(31)
			CALL MESH
		else if(ipl(ls) /= 0) then ! ífëwñ Ç™ïΩñ 
			CALL STMSAR2
			CALL MESHSS
		end if

	!==== /BLMS/ ====

		CALL SOURCE(hd,ls,ls+22,ls+44,ls+50,ls+56)

		close(ls+22)
		close(ls+44)
		close(ls+50)
		close(ls+56)
		if(ipl(ls) == 0) close(31)

	end do

return
end subroutine DISSOU

!==========

end module mod_dissou


