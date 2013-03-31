ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine getseed(seed)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     gets a starting seed for the random number generator based on 
c     the ship's clock
c
      implicit none
      integer yy,mm,dd,hh,mi,ss,seed,i,org
      character*14 curtime
      real psrran,r
      call clock(yy,mm,dd,hh,mi,ss)
      write(curtime,'(2i2.2)') ss,mm
      read(curtime,*) seed
      org=seed
      do i=1,ss
      r=psrran(seed)
      enddo
      seed=org
      end
