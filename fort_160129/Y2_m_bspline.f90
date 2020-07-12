!=================================================
! SUBROUTINE
!	BSPLINE(px,bx,ndeg)
!
!	no common
!=================================================

module func_bspline

	implicit none

contains

function b21(xx)
	implicit none
	real(8), intent(in) :: xx
	real(8) :: b21
		b21 = xx
	return
end function

function b31(xx)
	implicit none
	real(8), intent(in) :: xx
	real(8) :: b31
		b31 = xx * xx /2.d0
	return
end function

function b32(xx)
	implicit none
	real(8), intent(in) :: xx
	real(8) :: b32
		b32 = ((-2.d0 * xx +2.d0) * xx +1.d0) / 2.d0
	return
end function

function b41(xx)
	implicit none
	real(8), intent(in) :: xx
	real(8) :: b41
		b41 = xx * xx * xx / 6.d0
	return
end function

function b42(xx)
	implicit none
	real(8), intent(in) :: xx
	real(8) :: b42
		b42 = (((-3.d0 * xx + 3.d0) * xx + 3.d0) * xx + 1.d0) / 6.0d0
	return
end function

end module func_bspline

!***************************************

module mod_bspline

	implicit none

contains
!==========
subroutine BSPLINE(px,bx,ndeg)

	use func_bspline
	implicit none

	integer, intent(in) :: ndeg
	real(8), intent(in) :: px
	real(8), intent(out) :: bx(4)

	real(8) :: pr

	pr = 1.d0 - px

	if(ndeg == 0) then
		bx(1) = 1.d0
	else if(ndeg == 1) then
		bx(1) = b21(px)
		bx(2) = b21(pr)
	else if(ndeg == 2) then
		bx(1) = b31(px)
		bx(2) = b32(px)
		bx(3) = b31(pr)
	else if(ndeg == 3) then
		bx(1) = b41(px)
		bx(2) = b42(px)
		bx(3) = b42(pr)
		bx(4) = b41(pr)
	end if

return
end subroutine BSPLINE

end module mod_bspline

