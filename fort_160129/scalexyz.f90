program scalexyz

implicit none

real(8) :: dscale

integer :: i, ios, j, n = 1000000
real(8) :: dat(1:5)
character(len=30) :: fin  = "                              "
character(len=30) :: fout = "                              "
character(len=30) :: ein  = "                              "
character(len=30) :: eout = "                              "
character(len=10) :: cdsc = "          "

call getarg(1,fin )
call getarg(2,fout)
call getarg(3,ein )
call getarg(4,eout)
call getarg(5,cdsc)

read(cdsc,*) dscale

open(10,file=fin )
open(20,file=fout)

do i = 1, n

	read(10,*,iostat=ios) (dat(j),j=1,4)
	 if(ios == -1) exit
	dat(4) = dat(4)/dscale
	write(20,'(4F20.8)') (dat(j),j=1,4)

end do

close(10)
close(20)

if(ein(1:4)=="none") goto 100

open(11,file=ein )
open(21,file=eout)

do i = 1, n

	read(11,*,iostat=ios) (dat(j),j=1,5)
	 if(ios == -1) exit
	dat(4) = dat(4)/dscale
	dat(5) = dat(5)/dscale
	write(21,'(5F20.8)') (dat(j),j=1,5)

end do

close(11)
close(21)

100 continue

end
