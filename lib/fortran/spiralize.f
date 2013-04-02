c==============================================================================
      subroutine spiralize(r,seed,x,y)
c==============================================================================
c
c     Creates a spiral arm structure following the procedure given in fk06
c
      implicit none
      integer seed,arm
      real r,x,y,k(4),r0(4),theta0(4),theta,angle,twopi,dr,dx,dy
      real psr_normal,psrran
      parameter(twopi=6.283185307)
      data k/4.25,4.25,4.89,4.89/
      data r0/3.48,3.48,4.9,4.9/
      data theta0/1.57,4.71,4.09,0.95/
      save
c
c     first select one of 4 modeled spiral arms and calculate theta given r
c
      arm = int(psrran(seed)*4.0)+1
      theta=k(arm)*log(r/r0(arm))+theta0(arm)
c
c     apply blurring scheme 
c
      angle=twopi*psrran(seed)*exp(-0.35*r)
      if (psrran(seed).lt.0.5) angle=angle*(-1.0)
      theta=theta+angle
c
c     finally, blur around centroid with normal distribution
c
      dr=abs(psr_normal(seed,0.0,0.5*r))
      angle=psrran(seed)*twopi
      dx=dr*cos(angle)
      dy=dr*sin(angle)
c
c     combined results give x and y
c
      x=r*cos(theta)+dx
      y=r*sin(theta)+dy
      end
