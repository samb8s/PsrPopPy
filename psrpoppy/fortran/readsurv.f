c=============================================================================
      subroutine readsurv(filename,nsurv)
c=============================================================================
      implicit none
      include 'survey.inc'
      integer lun,nsurv,i
      character*80 string,filename*(*),dmfile
      logical filex
      call glun(lun)
      open(lun,file=filename,status="old")
      string='#'
      do while(string(1:1).eq.'#')
         read(lun,'(a)') string
      enddo
      read(string,*) beta(nsurv)
      read(lun,*) gain(nsurv)
      read(lun,*) tint(nsurv)
      read(lun,*) tsamp(nsurv)
      read(lun,*) tsys(nsurv)
      read(lun,*) nu(nsurv)
      read(lun,*) bw(nsurv)
      read(lun,*) dnu(nsurv)
      read(lun,*) npol(nsurv)
      read(lun,*) fwhm(nsurv)
      read(lun,*) ramin(nsurv)
      read(lun,*) ramax(nsurv)
      read(lun,*) decmin(nsurv)
      read(lun,*) decmax(nsurv)
      read(lun,*) lomin(nsurv)
      read(lun,*) lomax(nsurv)
      read(lun,*) abmin(nsurv)
      read(lun,*) abmax(nsurv)
      read(lun,*) sfrac(nsurv)
      read(lun,*) snmin(nsurv)
      close(lun)
      dmfile=filename(1:index(filename,' ')-1)//'.DMs'
      inquire(file=dmfile,exist=filex)
      if (filex) then
         open(lun,file=dmfile)
         string='#'
         do while(string(1:1).eq.'#')
            read(lun,'(a)') string
         enddo
         read(string,*)       dmdd(nsurv,1),acts(nsurv,1)
         do i=2,maxdms
            read(lun,*,end=1) dmdd(nsurv,i),acts(nsurv,i)
         enddo
 1       ndms(nsurv)=i-1
      else
         do i=1,maxdms
            dmdd(nsurv,i)=0.0
            acts(nsurv,i)=0.0
         enddo
         ndms(nsurv)=0
      endif
      end
