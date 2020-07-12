program grid_camp

implicit none

real(8), parameter :: rng(1:4) = (/148.d0,149.d0, 43.d0,46.d0/)
integer, parameter :: nlon = 10, nlat = 30
real(8), parameter :: dlon = (rng(2)-rng(1))/nlon, dlat = (rng(4)-rng(3))/nlat
!
! rng [領域指定] 経度(min,max) 緯度(min,max) [125,155],[20,50]
! (nlat+1)*(nlon+1) [grid node個数]
! 

integer :: ii, jj
real(8) :: alon, alat, dep
character(len=24) :: fflt = "                        "
character(len=24) :: fdat = "                        "

if (dlon <= 0.d0 .or. dlat <= 0.d0) stop "range error!"

call getarg(1,fflt)
call getarg(2,fdat)

open(10,file=fdat)
open(20,file=fflt)

write(10,*) "** CAMP STANDARD MODEL / TOHOKU TR. **"
write(10,*) "???? : the number of data"

do ii = 0, nlon
	alon = rng(1) + dble(ii) * dlon

	do jj = 0, nlat
		alat = rng(3) + dble(jj) * dlat

		CALL BDY(alon,alat,dep)

		if(dep ==  -10.d0) cycle
		if(dep == -100.d0) dep = 0.d0
		write(10,*) alon, alat, -dep
		write(20,*) alat, alon, -dep

	end do
	print *, ii, "/", nlon, "  Fin"

end do

close(10)
close(20)

end


