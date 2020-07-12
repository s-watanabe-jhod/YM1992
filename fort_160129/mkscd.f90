!===== MAKING SCD file from FMP file =====

program mkscd

use prm_var

use mod_ctpinv
use mod_setinv

use mod_dissou

implicit none

integer :: nsource, icm, igsi, i, ll
real(8) :: glat(3000), glon(3000)
character(len=24) :: flobspd, cal, fmp, hd, fctp = "                        "
character(len=6) :: cstation(3000)

open(77,file="../tmp/mkscd.out")

call getarg(1,fctp)

open(9,file=fctp)
CALL SETCP(hd,flobspd,cal,fmp,nsource)

open(33,file=fdum,status='old')
	CALL INPUTD(33)
close(33)

CALL SETUP(fdum,nsource,icm,1)

if(nsource /= 1) stop 'nsource /= 1: this program cannnot make multi-source!'

open(22,file=fmp)
CALL MKFMP(22)
close(22)

open(22,file=fmp)
CALL RDFMP(22)
close(22)

CALL DISSOU(hd,nsource)

print *, " -------------------------------------------"
print *, " SET GSI GPS POSITION"
print *, "  LATI min-max: ", int(glatmin), int(glatmax)
print *, "  LONG min-max: ", int(glonmin), int(glonmax)

open(98,file=fgps)
open(99,file=fpos)

read(98,*) igsi

ll = 1

do i = 1, igsi
	read(98,'(F11.7,1X,F11.8,1X,A6)') glon(ll), glat(ll), cstation(ll)
	if(glon(ll) < glonmin .or. glon(ll) > glonmax) cycle
	if(glat(ll) < glatmin .or. glat(ll) > glatmax) cycle
	ll = ll + 1
end do

write(99,*) 2*(ll-1), ll-1, ll-1, 0

do i = 1, ll-1
	write(99,'(I7,1X,F12.5,1X,F12.5,4X,A6)') i, glat(i), glon(i), cstation(i)
end do

do i = 1, ll-1
	write(99,'(I7,1X,F12.5,1X,F12.5,4X,A6)') i+(ll-1), glat(i), glon(i), cstation(i)
end do

close(98)
close(99)

close(77)
close(9)

stop


contains
!==========
subroutine RDFMP(ifn)

	use prm_matrix
	use prm_inv
	use prm_var

	implicit none

	integer, intent(in) :: ifn
	integer :: n, m, is, it, n1, n2, n3, iso, is1, it1

	n = 1

	read(ifn,*) 
	read(ifn,*) n1, n2, n3
	print *, "KS KT DEG", n1, n2, n3
	read(ifn,*) 

	do m = 1, 2
		do it = 1, n2
			do is = 1, n1
				read(ifn,*) iso, is1, it1, x(n), cs(n,n)
				if(is1 /= is) stop "is"
				if(it1 /= it) stop "it"
				if(iso > 4) stop "iso"
				n = n + 1
			end do
		end do
	end do 

return
end subroutine RDFMP

!==========
subroutine MKFMP(ifn)

	use prm_inv

	implicit none

	integer, intent(in) :: ifn
	integer :: m, is, it, n1, n2, n3, iss, itt
	real(8) :: xconst1, xconst2, const

	n1 = js0(1)
	n2 = jt0(1)
	n3 = ndgs(1)
	xconst1 = 0.002d0
	xconst2 = 0.2d0

	write(ifn,*) "KS KT DEG"
	write(ifn,*) n1, n2, n3
	write(ifn,*) " SOURCE    IS    IT         M.P.         E.E. "

	m = 1

	open(97,file="./init.dat")
	read(97,*)

	do it = 1, n2
		do is = 1, n1
			read(97,*) itt, iss, const
			if(itt/=it .or. iss/=is) stop "it or is error in init.dat"
			write(ifn,'(3I4,2F10.3)') m, is, it, 0.d0, 0.d0
		end do
	end do

	close(96)

	m = 2

	open(97,file="./init.dat")
	read(97,*)

	do it = 1, n2
		do is = 1, n1
			read(97,*) itt, iss, const
			if(itt/=it .or. iss/=is) stop "it or is error in init.dat"
			write(ifn,'(3I4,2F10.3)') m, is, it, -const, 0.d0
		end do
	end do

	close(97)

	write(ifn,*) "end"

return
end subroutine MKFMP

!==========

end


