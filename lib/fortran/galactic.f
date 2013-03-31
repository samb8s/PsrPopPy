c==============================================================================
      subroutine galactic(x, y, z, l, b)
c==============================================================================
c
c     Converts from x,y,z in kpc to l and b in degrees assuming
c     that the GC is at (0,0,0) and the Sun is at (0,rsun,0)
c
      implicit none
      real x, y, z, l, b, d, dcb, rsun, rad2deg
      parameter(rsun=8.5,rad2deg=57.2957795131)
      if (x.eq.0.and.y.eq.rsun.and.z.eq.0) stop 'silly pos to galactic'
      d=sqrt(x*x+(rsun-y)*(rsun-y)+z*z)
      b=asin(z/d)
      dcb=d*cos(b)
      if (y.le.rsun) then
         if (abs(x/dcb).ge.1.0) then 
            l=1.57079632679
         else
            l=asin(x/dcb)
         endif
      else 
         if (abs(x/dcb).ge.1.0) then 
            l=0.0
         else
            l=acos(x/dcb)
         endif
         l=l+1.57079632679
         if (x.lt.0) l=l-6.28318530718
      endif
      l=l*rad2deg
      b=b*rad2deg
      end
c==============================================================================
