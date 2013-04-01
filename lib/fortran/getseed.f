ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getseed(seed)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     gets a starting seed for the random number generator based on 
c     the ship's clock
c
      implicit none
      integer yy,mm,dd,hh,mi,ss,seed,i,org,nano
      character*14 curtime
      real psrran,r
      call clock(yy,mm,dd,hh,mi,ss,nano)
      write(curtime,'(i4)') nano
      read(curtime,*) seed
      org=seed
c      do i=1,ss
c      r=psrran(seed)
c      enddo
      seed=org
      end
