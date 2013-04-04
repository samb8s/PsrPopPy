ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine clock(yy,mm,dd,hh,mi,ss,nano)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Returns the present time from the ship's clock
c      
c     yy - i4 - year   
c     mm - i4 - month
c     dd - i4 - day  
c     hh - i4 - hour
c     mi - i4 - minute
c     ss - i4 - second
c      
c     Creation date: 98/06/08 (dlorimer@naic.edu)
c      
c==============================================================================
c      
      implicit none
      integer yy,mm,dd,hh,mi,ss,nano
      integer sysdate
      integer system

      sysdate = system("date +'%N' | awk '{print substr($1,6,4)}' >
     & tmp.date")
      open(unit=2,file="tmp.date")
      read(2,*) nano
      sysdate = system('rm -f tmp.date')
c      integer iarray(3)
      
c      call idate(iarray)
c      mm=iarray(2)
c      dd=iarray(1)
c      yy=iarray(3)
c      call itime(iarray)
c      hh=iarray(1)
c      mi=iarray(2)
c      ss=iarray(3)
      end
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
