c---------------------------------------
      real function psr_tsky(l, b, freq, inpath, linpath)
c---------------------------------------
c
c     Returns tsky for l,b.
c     If the first entry, reads the data from file.
c
      implicit none 
      logical first
      real l, b, nsky(90, 180),freq
      integer j,i,nl,lun,lpth
      character*120 path,tskyfile
      character inpath*(*)
      integer linpath
c
c     Check for first entry..
c
      data first / .true. /
      save first,nsky

c      write(*,*) linpath
c      write(*,*) inpath(1:linpath)
      if (first) then
c
c       Read in catalogue
c
        call glun(lun)
c        tskyfile=path(1:lpth)//'/lookuptables/tsky.ascii'
        tskyfile = inpath(1:linpath)//'/lookuptables/tsky.ascii'
        open(unit=lun, status='old', file=tskyfile, err=999)
        read(unit=lun, fmt=1000, end=998) ((nsky(i,j),j = 1, 180)
     &  ,i = 1, 90)
 1000   format(16f5.1)
        close(unit=lun) 
        first = .false.
      end if
c
c     Convert to standard l,b
c
      j = b + 91.5
      if (j .gt. 180) j = 180
      nl = l - 0.5
      if (l .lt. 0.5) nl = 359
      i = (nl / 4) + 1
c
c     Read off tsky from array converting from
c     408MHz n.b assume spectral index of -2.6
c
      psr_tsky = nsky(i,j) * (freq/408.0)**(-2.6)
      return 
c
c     Error messages, prog terminaes
c
  998 write(unit=*, fmt=1010) 
 1010 format(/40h ***** Unexpected end to TSKY file *****)
  999 write(unit=*, fmt=1020) 
 1020 format(/37h ***** Unable to open TSKY file *****)
      stop 
      end
