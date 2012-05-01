c==============================================================================
      subroutine readfn(filename,pmin,pmax,pdist,nbins)
c==============================================================================
      implicit none
      character*(*) filename
      real pdist(*),pmin,pmax
      integer nbins

      integer lun,i,istat

      i=0
      call glun(lun)
c      write(*,*) 'Opening ',filename
      open(unit=lun,file=filename,status='old',iostat=istat)
      if (istat.ne.0) stop 'Error opening file!'
      read(lun,*,iostat=istat) pmin,pmax
c      write(*,*) pmin,pmax
      do while(istat.eq.0)
         i=i+1
         read(lun,*,iostat=istat) pdist(i) 
c         write(*,*) pdist(i)
      enddo
      close(unit=lun)
      nbins=i-1
      if (i.lt.1) stop 'Error reading file!'
      end
