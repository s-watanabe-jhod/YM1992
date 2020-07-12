      program camp
      implicit NONE   

      logical FOREVER
      parameter (FOREVER=.true.)
      real*8 lon,lat
      real*8 x,y
      integer p,pmax,ptmp
      parameter (pmax=50000)
      real*8 a(pmax)
      real*8 k(pmax),l(pmax)
      real*8 dk,dl
      parameter (dk=8.0d0)
      parameter (dl=8.0d0)
      integer flag(pmax)
      real*8 Nxk,Nyl
      real*8 z

      do while(FOREVER)

         read(*,*)lon,lat
         if(lon.eq.0.0d0)then
            goto999
         endif
         if((lon.lt.125.0d0).or.(lon.gt.155.0d0)
     &        .or.(lat.lt.20.0d0).or.(lat.gt.50.0d0))then
            write(*,*)'OUT OF THE MODEL REGION'
            goto 992
         endif
         call Lambert(lon,lat,x,y)

         open(unit=90,file='camp_pac.bdy',status='old')
         do p=1,pmax
            read(90,*,END=990)k(p),l(p),a(p),flag(p)
            ptmp=p
         enddo
 990     continue
         close (unit=90)
         
         z=0.0d0
         do p=1,ptmp
            if(((x-2.0d0*dk.lt.k(p)).and.(x+2.0d0*dk.gt.k(p))).and.
     &           ((y-2.0d0*dl.lt.l(p)).and.(y+2.0d0*dl.gt.l(p))))then
               call spline(x,k(p),dk,Nxk)
               call spline(y,l(p),dl,Nyl)
               z=z+a(p)*Nxk*Nyl
               if(flag(p).eq.0)then
                  z=0.0d0
                  goto 991
               endif
            endif
         enddo
 991     continue
         if(z.gt.100.0d0)then
            z=0.0d0
         endif
         if(z.le.0.0d0)then
            write(*,*)'OUT OF THE MODEL REGION'
            goto 992
         endif
         write(*,*)z
 992     continue
      enddo
 999  continue

      stop
      end


      subroutine spline(s,sj,dsj,N)
      implicit NONE   

      real*8 s
      real*8 N
      real*8 sj,dsj
      real*8 sj2
      real*8 M

      sj2=sj+2.0D+00*dsj 

      if((s.ge.sj2-4.0D+00*dsj).and.(s.lt.sj2-3.0D+00*dsj))then
         M=(1.0D+00/(24.0D+00*(dsj**4.0d0)))
     &        *((s-sj2+4.0d0*dsj)**3.0d0)
      elseif((s.ge.sj2-3.0d0*dsj).and.(s.lt.sj2-2.0d0*dsj))then
         M=(1.0d0/(24.0d0*(dsj**4.0d0)))
     &        *(((s-sj2+4.0d0*dsj)**2.0d0)*(sj2-2.0d0*dsj-s)
     &        +(s-sj2+4.0d0*dsj)*(s-sj2+3.0d0*dsj)*(sj2-dsj-s)
     &        +(sj2-s)*((s-sj2+3.0d0*dsj)**2.0d0))
      elseif((s.ge.sj2-2.0d0*dsj).and.(s.lt.sj2-dsj))then
         M=(1.0d0/(24.0d0*(dsj**4.0d0)))
     &        *((s-sj2+4.0d0*dsj)*((sj2-dsj-s)**2.0d0)
     &        +(sj2-s)*(s-sj2+3.0d0*dsj)*(sj2-dsj-s)
     &        +((sj2-s)**2.0d0)*(s-sj2+2.0d0*dsj))
      elseif((s.ge.sj2-dsj).and.(s.lt.sj2))then
         M=(1.0d0/(24.0d0*(dsj**4.0d0)))*((sj2-s)**3.0d0)
      else
         M=0.0d0
      endif

      N=4.0d0*dsj*M

      return
      end


      subroutine Lambert(lambda,phi,x,y)
      implicit NONE

      real*8 PI 
      parameter (PI=3.1415926535897932d0)
      real*8 phi,lambda
      real*8 phi_0,lambda_0
      real*8 phi_1,phi_2
      real*8 theta
      real*8 n,rho
      real*8 rho_0
      real*8 F
      real*8 R
      real*8 x,y,xtmp,ytmp

      phi_0=35.0d0
      lambda_0=140.0d0
      phi_1=30.0d0
      phi_2=40.0d0
      R=6371.0d0

      phi_0=phi_0*((2.0d0*PI)/360.0d0)
      phi_1=phi_1*((2.0d0*PI)/360.0d0)
      phi_2=phi_2*((2.0d0*PI)/360.0d0)
      lambda_0=lambda_0*((2.0d0*PI)/360.0d0)
      n=log(cos(phi_1)/cos(phi_2))
     &     /log(tan((PI/4.0d0)+(phi_2/2.0d0))
     &     /tan((PI/4.0d0)+(phi_1/2.0d0)))
      F=(cos(phi_1)*(tan((PI/4.0d0)+(phi_1/2.0d0)))**n)/n
      rho_0=R*F/(tan((PI/4.0d0)+(phi_0/2.0d0))**n)
      phi=phi*((2.0d0*PI)/360.0d0)
      lambda=lambda*((2.0d0*PI)/360.0d0)
      rho=R*F/(tan((PI/4.0d0)+(phi/2.0d0))**n)
      theta=n*(lambda-lambda_0)
      xtmp=rho*sin(theta)
      ytmp=rho_0-rho*cos(theta)
      x=-ytmp
      y=-xtmp

      return
      end
