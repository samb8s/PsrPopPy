ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      program survey
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     observe model populations produced by "populate" for an arbitrary survey
c
c     There are three main usages:
c
c     survey name_of_popfile filename
c
c     reads survey parameters from single line of ASCII file
c     beta,gain,tint,tsamp,tsys,bw,nu,dnu,fwhm,par
c
c     survey name_of_popfile -parname1 parvalue1 -parname2 parvalue2 ....
c
c     specify the parameters on the command line, where the -parnames may be
c     any one of
c
c     beta,gain,tint,tsamp,tsys,bw,nu,dnu,fwhm,par
c
c     survey name_of_popfile filename -parname1 parvalue1 ....
c
c     reads parameters from a file then overrides given parameters
c
c     beta - survey degradation factor (1.0=perfect; 1.25=1-bit sampling...)
c     gain - telescope gain (K/Jy)
c     tint - integration time (s)
c     tsamp- sampling time (ms)
c     bw   - total bandwidth (MHz)
c     nu   - centre frequency (MHz)
c     nch  - number of channels across bw
c     npol - number of polarizations
c     fwhm - full-width half max of (assumed Gaussian) telescope beam (arcmin)
c
c
c     Apr 1, 2006 --- three further command-line switches are now available
c
c     -pms --- period in ms
c     -lum --- luminosity in mJy kpc2
c     -wms --- intrinsic with in ms
c     these will override the values in the parameter file to allow a scale
c     factor calculation
c
c     2011/10/12 - sbates adding scindx variable
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include 'survey.inc'
      logical filex,modified(maxsurv),sfan
      integer i,j,nmax,ngen,istat,maxdet,lpop,seed,ndet,ndets(maxsurv),
     &        smeared(maxsurv),tooweak(maxsurv),ldet,lsrc,lpopf,ntop,
     &        notvis(maxsurv),iargc,narg,nsurv,discovered(maxsurv),lpca
      parameter(nmax=1000000)
      integer detected(nmax),lout,lcline,lpth, lcarg
      integer gpsArray(nmax)
      real freq,ldeg(nmax),bdeg(nmax),ppsr(nmax),width(nmax),dmp(nmax),
     &     tau(nmax),tsky(nmax),dtrue(nmax),lpsr(nmax),s(nmax),si(nmax),
     &     sn(maxsurv),x,y,z,surveysn,par,flux,snmax,psr_dist,pdot,sfac,
     &     weff(maxsurv),wmax,dderi(nmax),wms,pms,lum,dcp,psr_normal,bf,
     &     gpsA(nmax), gpsB(nmax), gpsC, lgnu1, lgnu2, scindx
      character*80 string,value,source,popfile,cline*240,temp,path,outf
      character*80 carg
      common /outputtxt/ lout
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      outf="survey.det"
      ntop=0
      pdot=0.0
      lpca=0
      lout=6
      sfan=.false.
      par=0.0
      wms=0.0
      pms=0.0
      lum=0.0
      dcp=0.0
      scindx = -3.86
      do i=1,maxsurv
         modified(i)=.false.
         surveyfile(i)=''
         beta(i)=1.0
         gain(i)=1000.0
         tint(i)=10000.0
         tsys(i)=1.0
         tsamp(i)=0.001
         bw(i)=3000.0
         nu(i)=10000.0
         fwhm(i)=1.5
         ramin(i)=0.0
         ramax(i)=360.0
         decmin(i)=-90.0
         decmax(i)=90.0
         lomin(i)=-180.0
         lomax(i)=180.0
         abmin(i)=0.0
         abmax(i)=90.0
         sfrac(i)=1.0
      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      narg=iargc()
      if (narg.lt.2) then
         write(*,*)
         write(*,*) 'SURVEY: A program to simulate radio pulsar surveys'
         write(*,*)
         write(*,*) 'usage: survey name_of_popfile survey1 survey2...'
         write(*,*)
         stop
      endif
      call getarg(1,popfile)
      lpopf=index(popfile,' ')-1
      i=2
      nsurv=0
      do while (i.le.narg)
         call getarg(i,string)
         inquire(file=string,exist=filex)
         if (.not.filex) then
                call getpath(path,lpth)
                temp=path(1:lpth)//'/surveys/'//string
                inquire(file=temp,exist=filex)
	        path=path(1:lpth)//'/surveys/'
		lpth=lpth+9
	 else 
		path = './'
		lpth = 2
	 endif
         if (string.eq."ppdot.dat".and.filex) then
            call glun(lpca)
            open(lpca,file="ppdot.dat")
            sfan=.true.
         else if (string.eq."-o") then
	    i=i+1
	    call getarg(i,outf)
         else if (string.eq."-scindx") then
        i=i+1
        call getarg(i,carg)
        read(carg,*) scindx
         else if (string.eq."-n") then
	    i=i+1
	    call getarg(i,value)
            read(value,*) ntop
            nsurv=1
            surveyfile(1)='NOSURV'
         else if (filex) then
            nsurv=nsurv+1
            surveyfile(nsurv)=string
            call readsurv(path(1:lpth)//surveyfile(nsurv),nsurv)
         else
            i=i+1
            call getarg(i,value)
            read(value,*) par
            if (string.eq."-pms") then
               pms=par
               sfan=.true.
            else if (string.eq."-wms") then
               wms=par
            else if (string.eq."-lum") then
               lum=par
            else if (string.eq."-dcp") then
               dcp=par
            else
               modified(nsurv)=.true.
               if (string.eq."-beta") then
                  beta(nsurv)=par
               else if (string.eq."-gain") then
                  gain(nsurv)=par
               else if (string.eq."-tint") then
                  tint(nsurv)=par
               else if (string.eq."-tsys") then
                  tsys(nsurv)=par
               else if (string.eq."-tsamp") then
                  tsamp(nsurv)=par
               else if (string.eq."-bw") then
                  bw(nsurv)=par
               else if (string.eq."-nu") then
                  nu(nsurv)=par
               else if (string.eq."-dnu") then
                  dnu(nsurv)=par
               else if (string.eq."-npol") then
                  npol(nsurv)=par
               else if (string.eq."-fwhm") then
                  fwhm(nsurv)=par
               else if (string.eq."-ramin") then
                  ramin(nsurv)=par
               else if (string.eq."-ramax") then
                  ramax(nsurv)=par
               else if (string.eq."-decmin") then
                  decmin(nsurv)=par
               else if (string.eq."-decmax") then
                  decmax(nsurv)=par
               else if (string.eq."-lomin") then
                  lomin(nsurv)=par
               else if (string.eq."-lomax") then
                  lomax(nsurv)=par
               else if (string.eq."-abmin") then
                  abmin(nsurv)=par
               else if (string.eq."-abmax") then
                  abmax(nsurv)=par
               else if (string.eq."-sfrac") then
                  sfrac(nsurv)=par
               else 
                  stop 'unknown parameter'
               endif
               if (nsurv.eq.0) nsurv=1
            endif
         endif
         i=i+1
      enddo
      do i=1,nsurv
         if (modified(i)) then
         call glun(lpop)
         open(lpop,file="modified."//surveyfile(i),status="unknown")
         write(lpop,*) beta(i), "  ! survey degradation factor"
         write(lpop,*) gain(i), "  ! antenna gain (K/Jy)"
         write(lpop,*) tint(i), "  ! integration time (s)"
         write(lpop,*) tsamp(i),"  ! sampling time (ms)"
         write(lpop,*) tsys(i), "  ! system temperature (K)"
         write(lpop,*) nu(i),   "  ! centre frequency (MHz)"
         write(lpop,*) bw(i),"     ! bandwidth (MHz)"
         write(lpop,*) dnu(i),"    ! channel bandwidth (MHz)"
         write(lpop,*) npol(i),"    ! number of polarizations"
         write(lpop,*) fwhm(i),"   ! full-width half maximum (arcmin)"
         write(lpop,*) ramin(i),"   ! minimum RA (deg)"
         write(lpop,*) ramax(i),"   ! maximum RA (deg)"
         write(lpop,*) decmin(i)," ! minimum DEC (deg)"
         write(lpop,*) decmax(i)," ! maximum DEC (deg)"
         write(lpop,*) lomin(i),"   ! minimum Galactic longitude (deg)"
         write(lpop,*) lomax(i),"   ! maximum Galactic longitude (deg)"
         write(lpop,*) abmin(i),"! minimum abs(Galactic latitude) (deg)"
         write(lpop,*) abmax(i),"! maximum abs(Galactic latitude) (deg)"
         write(lpop,*) sfrac(i),"  ! fractional sky coverage (0-1)"
         write(lpop,*) snmin(i),"  ! minimum signal-to-noise ratio"
         inquire(file=dmlist(i),exist=filex)
         if (filex) write(lpop,*) dmlist(i)(1:index(dmlist(i),' ')-1),
     &                        "  ! list of trial DMs and tsamps"
         close(lpop)
         endif
      enddo
c
c     get a seed for the random number generator based on current time
c
      if (.not.sfan) write(*,*) "Fast sampled radio survey simulator"
      call getseed(seed)
      if (.not.sfan) write(*,*) "Starting seed for random numbers:",seed
c
c     Read the sources generated by populate
c
      call glun(lpop)
      open(lpop,file=popfile,form="unformatted")
      read(lpop) lcline,cline
      read(lpop) ngen,freq
      bf=0.0
      do i=1,ngen
         read(lpop) ldeg(i),bdeg(i),ppsr(i),width(i),dmp(i),tau(i),
     &                tsky(i),dtrue(i),dderi(i),lpsr(i),s(i),si(i),
     &                gpsArray(i), gpsA(i), gpsB(i)
         bf=bf+0.03+0.09*(log10(ppsr(i)/10000.0))**2
      enddo
      close(lpop)
      bf=bf/real(ngen)
      if (ppsr(1).eq.0.0) source = 'RRATs'
      if (ppsr(1).gt.0.0) source = 'pulsars'
      lsrc=index(source,' ')-1
      if (i-1.ne.ngen) stop "error reading psrpop.all"
      if (.not.sfan) write(*,*) "Read in",ngen," ",source(1:lsrc),
     & " from "//popfile(1:lpopf)
      if (.not.sfan) write(*,'(a)') " Generated via: "//cline(1:lcline)
 1    format(3(f9.2,1x),2(f6.1,1x),f11.1,1x,
     &       f5.2,1x,f7.1,1x,5(f7.2,1x),f10.2,1x,i3)
      if (.not.sfan) then
         call glun(ldet)
         open(ldet,file=outf,status="unknown")
         write(ldet,'(''# '',a)') source(1:index(source,' '))
         write(ldet,'(''# P (ms) * 0 1000'')')
         write(ldet,'(''# DM (cm\\\\u-3\\\\d pc) * 0 1000'')')
         write(ldet,'(''# W (ms) * 0 100'')')
         write(ldet,'(''# l (deg) * -180 180'')')
         write(ldet,'(''# b (deg) * -90 90'')')
         write(ldet,'(''# S (\\\\gmJy) * 0 1e6'')')
         write(ldet,'(''# SI * -3 0'')')
         write(ldet,'(''# S/N * 0 100'')')
         write(ldet,'(''# X (kpc) * -15 15'')')
         write(ldet,'(''# Y (kpc) * -15 15'')')
         write(ldet,'(''# Z (kpc) * -1 1'')')
         write(ldet,'(''# D (kpc) * 0 15'')')
         write(ldet,'(''# R (kpc) * 0 15'')')
         write(ldet,'(''# L (mJy kpc\\\\u2\\\\d) * 0 1000'')')
         write(ldet,'(''# code'')')
         do i=1,nsurv
            write(ldet,'(''# SURVEY: '',a)') 
     &           surveyfile(i)(1:index(surveyfile(i),' ')-1)
         enddo
      endif
      do i=1,nsurv
         discovered(i)=0
         ndets(i)=0
         smeared(i)=0
         tooweak(i)=0
         notvis(i)=0
      enddo
 444  ndet=0
      if (lpca.gt.0) then
         read(lpca,*,err=555,end=555) pms,pdot,lum
         if (lum.le.0.0) goto 444
      endif
c
c      LOOP OVER THE GENERATED PULSARS AND SEE IF DETECTED
c
      do i=1,ngen
         detected(i)=0
         snmax=0.0
         do j=1,nsurv
            if (lum.gt.0.0) s(i)=lum/dtrue(i)/dtrue(i) 
c       EDIT HERE TO ADD SUPPORT FOR GPS SPECTRA - FREQ IN GHZ
c       SBATES OCTOBER 2011
            if (gpsArray(i) < 1) then
              flux=s(i)*(nu(j)/freq)**si(i) ! scale flux to survey freq
            else
              lgnu1 = log10(freq/1000.0)
              lgnu2 = log10(nu(j)/1000.0)
c       CALCULATE C, USING REFERENCE FREQUENCY
              gpsC = log10(s(i))-(gpsA(i)*lgnu1**2)-gpsB(i)*lgnu1
c       CALCULATE FLUX, USING SURVEY FREQUENCY
              flux = 10.0**(gpsA(i)*(lgnu2**2)+gpsB(i)*lgnu2+gpsC)

            endif
            if (pms.gt.0.0) ppsr(i)=pms
            if (dcp.gt.0.0) then
               weff(j)=(dcp/100)*ppsr(i)**0.9
               weff(j)=log10(width(i))
               weff(j)= 10.0**psr_normal(seed,weff(j),0.27)
            else if (wms.gt.0.0) then
               weff(j) = wms
            else
               weff(j)=width(i)  
            endif
            sn(j)=surveysn(seed,ppsr(i),weff(j),dmp(i),tau(i),
     &      flux,ldeg(i),bdeg(i),tsky(i),freq,j, scindx)
            if (sn(j).gt.snmax) then
		snmax=sn(j)
		wmax=weff(j)
	    endif
            if (sn(j).gt.snmin(j)) then
		write(44,*) flux,s(i)
               if (detected(i).eq.0) discovered(j)=discovered(j)+1
               detected(i)=detected(i)+2**(j-1) ! assign binary detection code
               ndets(j)=ndets(j)+1
            else if (sn(j).eq.0.0) then
               notvis(j)=notvis(j)+1
            else if (sn(j).eq.-1.0) then
               smeared(j)=smeared(j)+1
            else 
               tooweak(j)=tooweak(j)+1
            endif
         enddo
         if (ntop.gt.0) detected(i)=1
         if (detected(i).gt.0) then
            ndet=ndet+1
            call calc_xyz(ldeg(i), bdeg(i), dderi(i), x, y, z)
c            call calc_xyz(ldeg(i), bdeg(i), dtrue(i), x, y, z)
            if (ldeg(i).gt.180.0) ldeg(i)=ldeg(i)-360.0
            if (.not.sfan)
     &      write(ldet,1) ppsr(i),dmp(i),wmax,ldeg(i),bdeg(i),
     &      1000*s(i),si(i),snmax,x,y,z,dderi(i),sqrt(x*x+y*y),
     &      s(i)*dderi(i)*dderi(i),gpsArray(i)
         endif
         if (ndet.eq.ntop.and.ntop.gt.0) stop 
      enddo
      if (sfan) then
         if (ndet.eq.0) ndet=-1 ! signifies a pulsar that was not found at all!
         sfac=real(ngen)/real(ndet)
         if (pdot.eq.0) then
            write(*,*) pms,sfac
         else
            write(*,*) pms,sfac,pdot,lum
         endif
         goto 444
      endif
      if (sfan) stop
      close(ldet)
      open(ldet,file="survey.info",status='unknown')
      write(ldet,*) ngen,ndet
      write(*,*) ngen," ",source(1:lsrc)," in model Galaxy"
      write(*,*) ndet," ",source(1:lsrc)," found"
      write(*,*) " Survey      Dis    Det    Out  Smear   Weak"
      do i=1,nsurv
         write(*,2) surveyfile(i)(1:index(surveyfile(i),' ')-1),
     &   discovered(i),ndets(i)-discovered(i),notvis(i),
     &   smeared(i),tooweak(i)
         write(ldet,2) surveyfile(i)(1:index(surveyfile(i),' ')-1),
     &   discovered(i),ndets(i)-discovered(i),notvis(i),
     &   smeared(i),tooweak(i)
      end do
      close(ldet)
 2    format(a8,3x,5(i6,1x))
 555  end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
