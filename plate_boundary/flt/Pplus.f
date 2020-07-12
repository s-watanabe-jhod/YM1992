CProgram B2C
CImplicit none

C     MAIN PROGRAM
C
      PARAMETER (NP=10,NDMAX=10000)
C
	INTEGER xx,M,N,T,tt,F(NDMAX),count,MON(NDMAX),Ma(10)
	REAL(KIND=16) X(NDMAX,3),Y(NDMAX,3)
C
	CHARACTER STAN*40
	CHARACTER ctmp*1
	OPEN(UNIT=10,FILE='./gphs.flt')
	OPEN(UNIT=20,FILE='./gphs+0.flt')
	!OPEN(UNIT=10,FILE='./gphs+.flt')
	!OPEN(UNIT=20,FILE='./gphs+3.flt')
	read(10,'(a40)')STAN
	write(20,'(a40)')STAN
	read(10,*)T
	M=0
10	M=M+1
	read(10,*,END=20)(X(M,xx),xx=1,3)
	GO TO 10
20	M=M-1
	close(10)
	
	count=0
	do tt=1,M-1
	count=count+1
	Y(count,1)=X(tt,1)
	Y(count,2)=X(tt,2)
!	if(X(tt,3).gt.-10..and.X(tt,3).lt.0.)then
!	X(tt,3)=X(tt,3)-0.5*(10+X(tt,3))
!	endif
	Y(count,3)=X(tt,3)
	if(X(tt,1).gt.X(tt+1,1))then
	count=count+1
	Y(count,1)=X(tt+1,1)-1.8!0.1
	Y(count,2)=X(tt+1,2)
	Y(count,3)=-X(tt+1,3)*0.5!2.0
!	count=count+1
!	Y(count,1)=X(tt+1,1)-1.2!0.1
!	Y(count,2)=X(tt+1,2)
!	Y(count,3)=-X(tt+1,3)*0.3!2.0
!	count=count+1
!	Y(count,1)=X(tt+1,1)-0.6!0.1
!	Y(count,2)=X(tt+1,2)
!	Y(count,3)=-X(tt+1,3)*0.3!2.0
	endif
	enddo

	write(20,'(i5)')count
	do tt=1,count
!	write(20,'(2(f6.2,1x),f14.8)')(Y(tt,xx),xx=1,3)
	write(20,'(2(f6.2,1x),f14.8)')Y(tt,1),Y(tt,2),Y(tt,3)
	enddo

	END
