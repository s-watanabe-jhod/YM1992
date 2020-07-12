!=================================================
! SUBROUTINES
!	SETCP(HD,flobspd,CAL,FMP,NSOURCE)
!
! ctp file Çì«Ç›çûÇﬁ
!
!=================================================

module mod_ctpinv

	implicit none

contains
!===========
subroutine SETCP(hd,flobspd,cal,fmp,nsource)

	use prm_var
	use prm_inv

	implicit none

	integer, intent(out) :: nsource
	character(len=24), intent(out) :: hd, flobspd, cal, fmp

	integer :: ls, lsin, k, kk

	read(9,*) alat0, alng0, icord
	read(9,'(A)') flobspd
	read(9,*) 
	read(9,*) 
	read(9,*) nsource

	do ls = 1, nsource

		read(9,*) 
		read(9,*) lsin

		if(lsin /= ls) stop "LS INPUT ERROR"

		read(9,*) 
		read(9,'(A)') fcof(ls)
		read(9,*) bphi(ls), bdlt(ls)
		read(9,*) blat(ls), blon(ls)
		read(9,*) blen(ls), bwid(ls)
		read(9,*) ju0(ls), jv0(ls), ndgf(ls)
		read(9,*) 
		read(9,'(A)') hd_jac(ls)
		read(9,'(A)') invjac(ls)
		read(9,*) 
		read(9,*) blaf(ls), blof(ls)
		read(9,*) blef(ls), bwif(ls), bdep(ls)
		read(9,*) js0(ls), jt0(ls), jd2(ls), ndgs(ls)

		ipl(ls) = 0

		if(fcof(ls) == 'PLANE FAULT             ') ipl(ls) = 1

	end do

	!---------------

	read(9,*) 
	read(9,'(A)') hd

	write(*,*)nsource
	do ls = 1, nsource
		read(9,*) 
		read(9,*) lsin

		if(lsin /= ls) stop "LS INPUT ERROR"

		read(9,'(A)') fscd(ls)
		read(9,*) just(1,ls), just(2,ls), just(3,ls), just(4,ls)
		read(9,*) nnls(1,ls), nnls(2,ls), nnls(3,ls), nnls(4,ls)
		read(9,*) rss(ls), rop(ls), rex(ls)
		read(9,*)ius(1,ls),ius(2,ls),ius(3,ls),ius(4,ls),mlg(ls)
		read(9,*) iboc(1,ls), iboc(2,ls), iboc(3,ls), iboc(4,ls)

		kk = 0

		do k = 1, 4

			if(just(k,ls) == 0) cycle

			kk = kk + 1
			just(k,ls) = kk

		end do

	end do

	!---------------

	read(9,*) 
	read(9,'(A)') cal
	read(9,'(A)') fmp

	read(9,*) 
	read(9,*) fact(1), fact(2), fact(3)
	read(9,*) fabs(1), fabs(2), fabs(3)
	read(9,*) vrp

	!---------------

	read(9,*) 
	read(9,*) glatmin, glatmax, glonmin, glonmax
	
	read(9,'(A)') fgps
	read(9,'(A)') fpos
	read(9,'(A)') fdum
	read(9,*) e_gsi, e_jcg
	read(9,*) setalpha

return
end subroutine SETCP

!==========

end module

