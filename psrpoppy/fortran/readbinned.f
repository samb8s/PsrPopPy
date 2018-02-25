ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine readbinned(filename,binned,xmin,xmax,nbins,maxbin)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     reads in a binned data file created by bindat and returns 
c     the binned data, ranges, number of bins and computs the
c     maximum binvalue. Data are then passed to the deviate function
c
      implicit none
      character filename*(*),line*80
      integer binned(*),nbins,maxbin,i,lbin
      real xmin,xmax
      logical filex

      inquire(file=filename,exist=filex)

      if (.not.filex) then
         write(*,*) "could not open file: ",filename
         stop "readbinned..."
      endif

      maxbin=0
      call glun(lbin)
      open(unit=lbin,file=filename,status="old")
      read(lbin,'(a)') line
      read(lbin,*) i,xmin,xmax,nbins
      do i=1,nbins
         read(lbin,*) binned(i)
         maxbin=max(maxbin,binned(i))
      enddo
      close(unit=lbin)
      end
