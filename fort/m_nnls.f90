!*************************************************
! module mod_nnls
!	NNLS(el,fs,ux,iuna,iula,iiq,eg,fb,u,p,la,lb)
!	
!
!*************************************************

module mod_nnls

	implicit none

contains
!=========
subroutine rtNNLS(el,fs,ux,iuna,iula,iiq,eg,fb,u,p,la,lb)

	use mod_qrdec

	implicit none

	integer, intent(in) :: iuna, iula, iiq
	real(8), intent(in) :: el(iuna,iiq), fs(iula)

	integer, intent(out) :: la(iiq), lb(iiq)
	real(8), intent(out) :: fb(iula), ux(iiq), u(iula), p(iiq)

	real(8), intent(inout) :: eg(iuna,iiq)

	integer :: i, j, jm, iiqa, iite, lea, j1
	real(8) :: ea, emi, q, qs, detl

!   *** PROGRAM NNLS ***
!   *** search UX WHICH NORM(EL*UX-FS) minimize subject to UX > 0 ***
!   *** UX  : SOLUTION ***
!   *** FB  : RESIDUAL FS-EL*UX  ***

!  ++++ INITIALIZE ++++

	print *, ' START OF PROGRAM NNLS, IIQ= ', iiq

	iiqa = iiq * 2
	iite = 0

	ux(1:iiq) = 0.d0
	la(1:iiq) = 1

lp1:do

		iite = iite + 1

!  +++ FB = FS-EL*UX,   P = EL/T*FB +++

		do i = 1, iula
			fb(i) = fs(i)
			do j = 1, iiq
				fb(i) = fb(i) - el(i,j) * ux(j)
			end do
		end do

		do j = 1, iiq
			p(j) = 0.d0
			do i = 1, iula
				p(j) = p(j) + el(i,j) * fb(i)
			end do
		end do

		ea = 0.d0
		lea = 0

		do j = 1, iiq
			if(la(j) == 0) cycle
			if( p(j) < ea) cycle

			lea = j
			ea = p(j)
		end do

		if(lea == 0) return

		la(lea) = 0

lp2:	do

!  +++ ( FB - EG*P ) MINIMIZE/  P: SELECTED PARAMETERS +++

			do i = 1, iula
				fb(i) = fs(i)
			end do

			jm = 0

			do j = 1, iiq
				if(la(j) == 1) cycle

				jm = jm + 1

				do i = 1, iula
					eg(i,jm) = el(i,j)
				end do
			end do

! ----- CALL LSQS(eg,fb,p,iuna,iula,jm,lb,u) -----

			CALL QRDEC(eg,fb,detl,iuna,iula,jm,u,lb,j1,1)

			if(j1 /= jm) stop ' ERROR LSQS UNDERDETERMINED '

			CALL RINV(p,fb,eg,iuna,jm,0)

! ----- /LSQS/ -----

			if(iite >= iiqa) stop ' IITE LIMIT OVER'

			emi = 1.d2

			do j = iiq, 1, -1

				if(la(j) == 0) then

					p(j) = p(jm)
					jm = jm - 1

					if(p(j) < emi) emi = p(j)

				else if(la(j) == 1) then

					p(j) = 0.d0

				end if

			end do

			if(emi <= 0.d0) then

				q = 1.d2
				lea = 0

				do j = 1, iiq
					if(la(j) == 1) cycle
					if( p(j) > 0.d0) cycle

					qs = ux(j) / (ux(j) - p(j))
					if(qs >= q) cycle

					q = qs
					lea = j
				end do

				do j = 1, iiq
					ux(j) = ux(j) + q *(p(j) - ux(j))
					if(ux(j) > 0.d0) cycle

					ux(j) = 0.d0
					la(j) = 1
				end do

				la(lea) = 1
				ux(lea) = 0.d0

				cycle lp2

			else if(emi > 0.d0) then

				do j = 1, iiq
					ux(j) = p(j)
				end do

				cycle lp1

			end if

		end do lp2

	end do lp1

end subroutine rtNNLS

!=========

end module mod_nnls
