c==============================================================================
      real function ykarea(r,amax,a,b,r1)
c==============================================================================
c
c     if amax=0.0, this function returns the integral  x*rho(x) dx 
c     for the limits 0<x<r, where rho(x) is the function derived by 
c     Yusifov & Kukuk (2004)
c
c     if amax>0.0, this function returns the value of r when the 
c     integral is equal to amax.
c
      implicit none
      real r,a,b,r1,x,dx,r0,rsun,amax
      parameter(dx=0.01,rsun=8.5)

      r0=rsun+r1
      ykarea=0.0
      do x=0,r,dx
         ykarea=ykarea+(((x+r1)/r0)**a)*exp(-1.0*b*(x-rsun)/r0)*x*dx
         if (amax.gt.0.0.and.ykarea.gt.amax) then
            ykarea=x
            return
         endif
      enddo
      end
c==============================================================================
      real function ykr(seed)
c==============================================================================
c
c     uses ykarea to draw a random deviate from Yusifov & Kukuk's 
c     radial distribution. a=1.64,b=4.01,r1=0.55
c
      integer seed
      real ykr0
      ykr=ykr0(seed,1.64,4.01,0.55)
c      ykr= 5.0
      end function
c==============================================================================
      real function llfr(seed)
c==============================================================================
c
c     uses ykarea to draw a random deviate from Yusifov & Kukuk's 
c     radial distribution. a=1.64,b=4.01,r1=0.0
c
      integer seed
      real ykr0
      llfr=ykr0(seed,3.51,7.89,0.0)
      end
c==============================================================================
      real function ykr0(seed,a,b,r1)
c==============================================================================
c
c     uses ykarea to draw a random deviate from Yusifov & Kukuk's 
c     radial distribution. a=1.64,b=4.01,r1=0.55
c
      implicit none
      real amax,area,ykarea,a,b,r1,psrran
      integer seed
      logical first
      data first/.true./
      save


      if (first) then
         amax=ykarea(500.0,0.0,a,b,r1)
c         write(*,*) "first!!!"
         first=.false.
      endif
      area = psrran(seed) * amax
      ykr0 = ykarea(500.0,area,a,b,r1)
      end
c==============================================================================
