CProgram B2C
CImplicit none

C     MAIN PROGRAM
C
      PARAMETER (NP=10,NDMAX=5100)
C
	INTEGER xx,M,N,T,tt,F(NDMAX),tmp,MON(NDMAX),Ma(10)
	REAL(KIND=16) X(NDMAX,3),Y(NDMAX,3),A(NDMAX,3),NTMP(3)
C
	CHARACTER STA(NDMAX,10)*16,tmp6*6,tmp16*16,tmp10*10,tmp7*7
	CHARACTER ctmp*1
      OPEN(UNIT=10,FILE='./gphs.dat')
      OPEN(UNIT=20,FILE='./gphs2.dat')
	do M=1,5029
	read(10,*)(NTMP(N),N=1,3)
	write (20,'(f6.2,1x,f6.2,1x,f10.6)')NTMP(2),NTMP(1),NTMP(3)
	enddo
	stop
	end
