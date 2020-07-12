program mkgpspos

implicit none

real(8) :: amax, amin, omax, omin, alat, alon, dlat, dlon
integer :: lat, lon, i, j, ii

amin = 31.d0 	! LAT min
amax = 36.d0 	! LAT max
omin = 131.d0 	! LON min
omax = 138.d0 	! LON max

lat  = 25 		! number of grid
lon  = 35 		! number of grid

!=========

dlat = (amax-amin) / dble(lat-1)
dlon = (omax-omin) / dble(lon-1)

if(lat*lon > 1000) stop 'grid number is too large'
if(amin > amax) stop 'lat error'
if(omin > omax) stop 'lon error'

open(10,file="./check/gps_grid.dat")
open(11,file="./check/grid_dummyg.dat")

write(10,*) lat*lon
write(11,*) lat*lon, lat*lon, 0, 0

ii = 0

do i = 1, lat
	alat = amin + dlat * dble(i-1)

	do j = 1, lon
		ii = ii + 1
		alon = omin + dlon * dble(j-1)
		write(10,'(F11.7,1X,F11.8,1X,"#  gsi ")') alon, alat
		write(11,'(I5,1X,F11.7,1X,F11.8,1X,"0.0 0.0  1  1")') ii, alon, alat
	end do

end do

close(10)
close(11)
print *, "FILE : ./check/gps_grid.dat"

end

