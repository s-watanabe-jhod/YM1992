program grid_camp

implicit none

real(8), parameter :: rng(1:4) = (/122.d0,125.d0, 22.d0,26.d0/)
integer, parameter :: nlon = 30, nlat = 40
real(8), parameter :: dlon = (rng(2)-rng(1))/nlon, dlat = (rng(4)-rng(3))/nlat
!
! rng [領域指定] 経度(min,max) 緯度(min,max) [125,155],[20,50]
! (nlat+1)*(nlon+1) [grid node個数]
! 

integer :: ii, jj
real(8) :: alon, alat, dep
character(len=24) :: finp = "                        "

if (dlon <= 0.d0 .or. dlat <= 0.d0) stop "range error!"

call getarg(1,finp)

open(10,file=finp)

do ii = 0, nlon
	alon = rng(1) + dble(ii) * dlon

	do jj = 0, nlat
		alat = rng(3) + dble(jj) * dlat

		CALL BDY(alon,alat,dep)

		if(dep <= -3.d0) cycle
		write(10,*) alat, alon, -dep

	end do
	print *, ii, "/", nlon, "  Fin"

end do

close(10)

end


