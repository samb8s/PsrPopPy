c==============================================================================
                            program populate
c==============================================================================
c
c   This program is trying to update the Lyne, Manchester & Taylor (1985)
c   analysis to determine the underlying Galactic distribution of pulsars.
c
c   A model UNDERLYING population is drawn from up to 5 distribution functions:
c
c     R : Galactocentric Radius
c     z : Height above the galactic plane
c     L : Radio "Luminosity" at 400 MHz 
c     P : Period
c     W : Intrinsic pulse width
c
c   These functions are supplied by the user, the program then applies a
c   model of the Parkes multibeam survey on the data to derive the model
c   OBSERVED sample which can then be directly compared to the real sample.
c
c==============================================================================
      implicit none
      include 'survey.inc'
      real pi,rsun,rad,psr_exp
      parameter(pi=3.1415927,rsun=8.5,rad=57.29577951)

      integer seed, maxpsrs, i, j, ngen, unbeamed, dmod, lpth
      parameter (maxpsrs=1000000)

      real sm,smtau,smth,smiso,tauiss,fghz,surveysn,logtau
      character lim*1,path*80,fixedfile*80
      integer maxf,nf,fn, gpsArray(maxpsrs), gpsCount, gpsSwitch
      parameter (maxf=64)
      real fl(maxf),fb(maxf),fd(maxf)
      real*4 r0,lrad,brad,dmax,scindx,
     &     x(maxpsrs), y(maxpsrs), rd(maxpsrs), x2, y2, z2,
     &     z(maxpsrs), ppsr(maxpsrs), lpsr(maxpsrs), lumchoose,
     &     tsky(maxpsrs), width(maxpsrs,10), spindex(maxpsrs),
     &     alpha(maxpsrs), rho(maxpsrs,10), tau(maxpsrs),
     &     dmp(maxpsrs), dderi(maxpsrs), sn(maxpsrs,10),
     &     ldeg(maxpsrs),bdeg(maxpsrs),freq, dc, drawgamma,
     &     psr_tsky, psrran, fbea, psr_normal, pmean, psigma,
     &     psr_dist, theta0, derr, lderi, dm, wms,gamma_a,gamma_m,
     &     dtrue(maxpsrs), bfact, s(maxpsrs),simean,sistdv,
     &     gpsA, gpsB, gpsFrac, flux, minfreq
      logical detected(maxpsrs),yk04,isotropic,spiral,slab,pgamma,
     &      llnorm,beaming,filexists,cc97,lumfile,scattering,
     &      pnorm,plnorm,fixed,llf06,disk,uniformr
      logical beamresults(maxpsrs,10)
      integer minfreqint
c
c     Model Observed population
c
      integer ndet,n2det
      real rdet(maxpsrs), zdet(maxpsrs), pdet(maxpsrs),duty
c
c     Initial values...
c
      real draw1d
      integer maxbins
      parameter(maxbins=256)
c
c     Radial distribution
c
      integer nrbins
      real rprob(maxbins),rmin,rmax,rscale,ykr,llfr
c
c     Z distribution
c
      integer nzbins
      real zprob(maxbins),zmin,zmax,zscale
c
c     Period distribution
c
      integer npbins
      real pprob(maxbins),pmin,pmax,pstatic
c
c     Luminosity distribution
c
      integer nlbins
      real lprob(maxbins),lmin,lmax,loglmax,lslope,ntot,n,lmean,lsigma
      real lintercept,lfunc
c
c     pulse width distribution
c
      integer nwbins
      real wprob(maxbins),wmin,wmax
c
c     sys stuff
c
      real slim
      integer narg,iargc,istat,system,lun,nsurv,lout,lcarg,lcline
      character*80 carg,cline*240,string,outfile,temp
      character*3 uprow
      common /outputtxt/ lout
      character (len=80) :: outfiles(3)
      include 'version.inc'
c      character (len=80) :: outfiles(3) = (/"psrpop.model1",
c     &                                     'psrpop.model2',
c     &                                     'psrpop.model3'/)
      lout=6
      do i=1,maxsurv
         surveyfile(i)=''
      enddo
      outfile='psrpop.model'
      outfiles(1)='psrpop.model1'
      outfiles(2)='psrpop.model2'
      outfiles(3)='psrpop.model3'
      uprow = char(27)//'[A'
      narg=iargc()
      zscale=0.0
      pstatic=-1.0
      dmax=0.0
      duty=-1.0
      scindx=-3.86
      simean=-1.6
      sistdv=0.35
      freq=1374.0
      fixed=.false.
      llnorm=.false.
      yk04=.false.
      llf06=.false.
      uniformr=.true.
      cc97=.false.
      pgamma=.false.
      pnorm=.false.
      plnorm=.false.
      scattering=.true.
      lumfile=.true.
      isotropic=.false.
      slab=.false.
      disk=.false.
      spiral=.false.
      seed=-1
      rscale=0.0
      gpsSwitch = 0
      gpsA = 0.
      gpsB = 0.
      call getarg(1,carg)
      if (carg.eq.'-version') then
         write(*,'(''POPULATE part of PSRPOP version '',f4.1)') 
     &   psrpop_version
         stop
      endif
      if (narg.lt.1.or.carg.eq.'-help') then
        write(*,*)
        write(*,*) "POPULATE: A program to generate pulsar populations"
        write(*,*)
        write(*,'('' usage: populate {options} {survey file(s)}'')')
        write(*,*)
        write(*,'('' -o file - output file (def=psrpop.model)'')')
        write(*,'('' -n ngen - number of pulsars to generate/detect'')')
        write(*,'('' -d dmod - distance model (lmt=0;bwhv=1;cwfsr=2;'',
     &            ''tc93=3;ne2001=4=default)'')')
        write(*,'('' -e derr - distance error (fraction; def=0.0)'')')
        write(*,*)"-z zkpc - use exponential scale height (def=model)"
        write(*,'('' -s seed - seed for random number generator '',
     &            ''(def=based on time of day)'')')
        write(*,*)"-f file - fix sky positions according to file"
        write(*,*)"-p dist - period distribution (def=from file)"
        write(*,*)"-w duty - width from duty cycle in % (def=from file)"
        write(*,*)"-r dist - radial distribution (def=from file)"
        write(*,*)"-si m s - Mean and sigma of spectral",
     &            " indices (def=",simean,sistdv,")"
        write(*,*)"-sc ind - Spectral index for scattering",
     &            " (def=",scindx,")"
        write(*,*)"-rf mhz - Reference frequency for luminosities",
     &            " (def=",freq," MHz)"
        write(*,*)"-rrats  - generate RRATs with no spin period info"
        write(*,*)"-quiet  - suppress output text to file populate.out"
        write(*,*)"-spiral - generates spiral arms a la fk06"
        write(*,*)"-scatt0 - set all scattering times to zero"
        write(*,*) 'add "-gps N a b" to add N% GPS pulsars'
        write(*,*) 'with ax^2 + bx + c spectrum shape'
        write(*,*)
c        write(*,*) "Example 1. Generate 1000 pulsars with LMT distmod"
c        write(*,*) "> populate -n 1000 -d 0 -e 0.1"
c        write(*,*)
c        write(*,*) "Example 2. Generate until 1000 pulsars are detected"
c        write(*,*) "> populate -n 1000 -d 0 -e 0.1 PMSURV"
c        write(*,*) "N.B. PMSURV is a separate (pregenerated) ASCII file"
c        write(*,*)
        stop
      else
         dmod=4
         derr=0.0
         nsurv=0
         lcline=0
c
c        write out command line options to an ASCII file
c
         call glun(lun)
         open(lun,file="populate.cmd",status="unknown",access="append")
         do i=1,narg
            call getarg(i,carg)
            lcarg=index(carg,' ')
            cline=cline(1:lcline)//carg(1:lcarg)
            lcline=lcline+lcarg
         enddo
         write(lun,'("populate ",a)') cline(1:lcline)
         close(lun)
c
c        Now parse the command line
c
         i=1
         do while (i.le.narg)
            call getarg(i,carg)
            inquire(file=carg,exist=filexists)
            if (.not.filexists) then
		call getpath(path,lpth)
		temp=path(1:lpth)//'/surveys/'//carg
		inquire(file=temp,exist=filexists)
	    else
		temp=carg
	    endif
            if (filexists) then
               nsurv=nsurv+1
               surveyfile(nsurv)=carg
               call readsurv(temp,nsurv)
            else if (carg.eq.'-quiet') then
               call glun(lout)
               open(lout,file='populate.out',status='unknown')
               write(lout,'("populate ",a)') cline(1:lcline)
            else if (carg.eq.'-rrats') then
               pstatic=0.0
            else if (carg.eq.'-n') then
               i=i+1
               call getarg(i,carg)
               read(carg,*) ngen
               if (ngen.gt.maxpsrs) stop 'max number of psrs reached!'
            else if (carg.eq.'-o') then
               i=i+1
               call getarg(i,outfile)
            else if (carg.eq.'-d') then
               i=i+1
               call getarg(i,carg)
               read(carg,*) dmod
            else if (carg.eq.'-dmax') then
               i=i+1
               call getarg(i,carg)
               read(carg,*) dmax
            else if (carg.eq.'-f') then
               i=i+1
               call getarg(i,fixedfile)
               inquire(file=fixedfile,exist=filexists)
               if (.not.filexists) then
                  write(*,*) 'fixed file does not exist'
                  stop
               else
                  call glun(lun)
                  open(lun,file=fixedfile,status='old')
                  j=0
                  do while (.true.)
                     j=j+1
                     if (j.gt.maxf) stop 'too many fixed file entries'
                     read(lun,*,end=222,err=111) fl(j),fb(j),fd(j)
                  enddo
 111              stop 'error reading fixed file'
 222              close(lun)
                  nf=j-1
                  fn=0
                  if (nf.lt.1) stop 'nothing in fixed file'
                  fixed=.true.
               endif
            else if (carg.eq.'-p') then
               i=i+1
               call getarg(i,carg)
               if (carg.eq.'cc97') then
                  cc97=.true.
               else if (carg.eq.'gamma') then
                  pgamma=.true.
                  call getarg(i+1,carg)
                  read(carg,*) gamma_a
                  call getarg(i+2,carg)
                  read(carg,*) gamma_m
                  i=i+2
               else if (carg.eq.'norm') then
                  pnorm=.true.
                  call getarg(i+1,carg)
                  read(carg,*) pmean
                  call getarg(i+2,carg)
                  read(carg,*) psigma
                  i=i+2
               else if (carg.eq.'lnorm') then
                  plnorm=.true.
                  call getarg(i+1,carg)
                  read(carg,*) pmean
                  call getarg(i+2,carg)
                  read(carg,*) psigma
                  i=i+2
               else
                  read(carg,*) pstatic
               endif
            else if (carg.eq.'-l') then
               lumfile=.false.
               call getarg(i+1,carg)
               if (carg.eq.'lnorm') then
                  llnorm=.true.
                  call getarg(i+2,carg)
                  read(carg,*) lmean
                  call getarg(i+3,carg)
                  read(carg,*) lsigma
               else
                  read(carg,*) lmin
                  call getarg(i+2,carg)
                  read(carg,*) lmax
                  call getarg(i+3,carg)
                  read(carg,*) lslope
               endif
               i=i+3
            else if (carg.eq.'-w') then
               i=i+1
               call getarg(i,carg)
               read(carg,*) duty
            else if (carg.eq.'-isotropic') then
               isotropic=.true.
            else if (carg.eq.'-slab') then
               slab=.true.
            else if (carg.eq.'-disk') then
               disk=.true.
            else if (carg.eq.'-spiral') then
               spiral=.true.
            else if (carg.eq.'-r') then
               i=i+1
               call getarg(i,carg)
               if (carg.eq.'yk04') then
                  yk04=.true.
               else if (carg.eq.'lfl+06') then
                  llf06=.true.
               else if (carg.eq.'uniform') then
                  uniformr=.true.
               else
                  read(carg,*) rscale
               endif
            else if (carg.eq.'-e') then
               i=i+1
               call getarg(i,carg)
               read(carg,*) derr
            else if (carg.eq.'-z') then
               i=i+1
               call getarg(i,carg)
               read(carg,*) zscale
            else if (carg.eq.'-rf') then
               i=i+1
               call getarg(i,carg)
               read(carg,*) freq
            else if (carg.eq.'-s') then
               i=i+1
               call getarg(i,carg)
               read(carg,*) seed
            else if(carg.eq.'-sc') then
                i=i+1
                call getarg(i, carg)
                read(carg,*) scindx
            else if (carg.eq.'-si') then
               i=i+1
               call getarg(i,carg)
               read(carg,*) simean
               i=i+1
               call getarg(i,carg)
               read(carg,*) sistdv
            else if (carg.eq.'-gps') then
                gpsSwitch = 1
                i = i+1
                call getarg(i, carg)
                read(carg,*) gpsFrac
                i = i + 1
                call getarg(i, carg)
                read(carg,*) gpsA
                i = i+1
                call getarg(i, carg)
                read(carg,*) gpsB
            else if (carg.eq.'-scatt0') then
               scattering=.false.
            endif
            i=i+1
         enddo
      endif
c
c     decide and notify user of various population options
c
      if (cc97) then
         write(lout,*) 'Using cc97 MSP period distribution'
      else if (pgamma) then
         write(lout,*) 'Using gamma period distribution'
         write(lout,*) 'A = ',gamma_a,' M = ',gamma_m
      else if (pnorm) then
         write(lout,*) 'Using normal period distribution'
         write(lout,*) 'Pmean = ',pmean,' ms Psigma = ',psigma,' ms'
      else if (plnorm) then
         write(lout,*) 'Using log-normal period distribution'
         write(lout,*) 'Pmean = ',10.0**pmean,' ms Psigma = ',
     &                  10.0**psigma,' ms'
      else if (pstatic.eq.0.0) then
         write(lout,*) 'Generating RRATs instead of pulsars'
      else if (pstatic.gt.0.0) then
         write(lout,*) 'Using fixed period of ',pstatic,' ms'
      else
         write(lout,*) 'Using period distribution from pdist.dat'
         call readfn("pdist.dat",pmin,pmax,pprob,npbins)
      endif
      if (duty.eq.-1.0) then
         write(lout,*) 'Using pulse-width distribution from wdist.dat'
         call readfn("wdist.dat",wmin,wmax,wprob,nwbins)
      else if (duty.gt.0.0) then
         write(lout,*) 'Using fixed duty cycle of ',duty,'%'
      else
         write(lout,*) 'Using beaming model to calculate pulse widths'
      endif
      if (isotropic) then
         write(lout,*) 'Generating isotropic population'
      else if (fixed) then
         write(lout,*) 'Generating fixed position population'
      else if (slab) then
         write(lout,*) 'Generating slab population'
      else if (disk) then
         write(lout,*) 'Generating disk population'
      else
         if (rscale.ne.0.0) then
            write(lout,*) 'Using Gaussian radial distribution',
     &           ' with scale length:',rscale,' kpc'
         else if (yk04) then
            write(lout,*) 'Using Yusifov & Kucuk radial distribution'
         else if (llf06) then
            write(lout,*) 'Using Lorimer et al. radial distribution'
         else if (uniformr) then
            write(lout,*) 'Using uniform spatial distribution'
         else
            write(lout,*) 'Using radial distribution from rdist.dat'
            call readfn("rdist.dat",rmin,rmax,rprob,nrbins)
         endif
         if (spiral) 
     &      write(lout,*) 'Spiral-arm structure included'
         if (zscale.ne.0.0) then
            write(lout,*) 'Using exponential z-scale height of',
     &           zscale,' kpc'
         else
            write(lout,*) 'Using z distribution from zdist.dat'
            call readfn("zdist.dat",zmin,zmax,zprob,nzbins)
         endif
      endif
      write(lout,*) 'Spectral index mean and standard deviation',
     &          simean,sistdv
      write(lout,*) 'Spectral index for scattering',
     &          scindx
      if (lumfile) then
         write(lout,*) 'Using luminosity distribution from ldist.dat'
         call readfn("ldist.dat",lmin,lmax,lprob,nlbins)
      else if (llnorm) then
         write(lout,*) 'Log-normal luminosity distribution',lmean,lsigma
      else if (lmin.eq.lmax) then
         write(lout,*) 'Assigning fixed luminosity value:',lmin,
     &       ' mJy kpc^2'
      else
         write(lout,*) 'Luminosity function min, max and slope:',
     &   lmin,lmax,lslope
      endif
      write(lout,*) 'Reference frequency',freq,' MHz.'
      fghz=freq/1000.0
      if (.not.scattering) write(lout,*) 'No scattering in this galaxy!'
      if (seed.eq.-1) call getseed(seed)
      write(lout,*) 'Initial seed for random numbers:',seed
c
c     Generate the pulsars
c
      if (nsurv.eq.0) then
        if (pstatic.eq.0.0) then
           write(lout,*) 'About to generate',ngen,' RRATs...'
        else
           write(lout,*) 'About to generate',ngen,' pulsars...'
        endif
        ndet=0
        n2det=0
      else
        if (pstatic.eq.0.0) then
           write(lout,*) 'Attempting to detect',ngen,' RRATs...'
        else
           write(lout,*) 'Attempting to detect',ngen,' pulsars...'
        endif
        n2det=ngen
        ngen=maxpsrs
      endif
      unbeamed=0
c      write(lout,*)

c     GET MINIMUM FREQUENCY AND CORRESPONDING SURVEYINT

      minfreq = nu(1)
      minfreqint = 1
      do j=1, nsurv
        if (nu(j) .lt. minfreq) then
          minfreq = nu(j) ! this is the lowest frequency
          minfreqint = j
        endif
      enddo

      do i = 1, ngen
c
c       intial period distribution
c
 1      continue
        if (cc97) then
           ppsr(i)=0.65/(1.0-psrran(seed)) ! cc97 MSP period distribution
           if (ppsr(i).gt.30.0) goto 1 ! impose hard cut-off at 30 ms
           if (ppsr(i).lt.1.0) goto 1  ! and another one at 1 ms
        else if (pgamma) then
           ppsr(i)=1000.0*drawgamma(0.0,10.0,gamma_a,gamma_m,seed)
        else if (pnorm) then
           ppsr(i)=psr_normal(seed,pmean,psigma)
           if (ppsr(i).le.0.0) goto 1
        else if (plnorm) then
           ppsr(i)=10.0**psr_normal(seed,pmean,psigma)
        else if (pstatic.gt.0.0) then
           ppsr(i)=pstatic
        else
c          period distribution from something like the observations
           ppsr(i)=10.0**draw1d(pmin,pmax,pprob,npbins,seed)
        endif
	write(55,*) ppsr(i)
c
c       inclination angle in range 0<alpha<90 degrees in uniform omega
c
        alpha(i)=abs(asin(2.0*psrran(seed)-1.0)*rad)
        if (duty.eq.-1.0) then
           width(i,1)=draw1d(wmin,wmax,wprob,nwbins,seed)
        else if (duty.gt.0.0) then
 777       width(i,1)=(duty/100)*ppsr(i)**0.9
           width(i,1)=log10(width(i,1))
           width(i,1)= 10.0**psr_normal(seed,width(i,1),0.3)
           if (width(i,1)/ppsr(i).lt.0.005) goto 777
c           width(i)=(duty/100)*ppsr(i)
        else
c          Hacking the code here to include a frequency-dependant rho
c          S BATES 2012

           do j=1, nsurv
             call RFMrhowidth(ppsr(i), alpha(i), minfreq, rho(i,j),
     &                        width(i,j), seed)
           enddo
c           call getrhowidth(ppsr(i),alpha(i),seed,rho(i),width(i))
c           if (width(i,minfreqint).eq.0.0.and.rho(i,minfreqint).eq.0.0) 
c     &         goto 1
c
c          Figure out if the pulsar is beaming towards us first!
c          AT EACH FREQUENCY - AND STORE IN ARRAY
           do j=1, nsurv
             beamresults(i,j) = beaming(seed,alpha(i),rho(i,j))
           enddo
c           if (.not.beaming(seed,alpha(i),rho(i,minfreqint))) then

c          CHECK IF BEAMING AT LOWEST FREQUENCY (WIDEST BEAM)
c          IF NOT, WE CAN SKIP EVERYTHING ELSE
           if (.not.beamresults(i,minfreqint)) then
              unbeamed=unbeamed+1
              goto 1
           endif
           do j=1, nsurv
             width(i,j)=ppsr(i)*(width(i,j)/360.0) ! width in ms
           enddo
        endif
        if (pstatic.eq.0.0) ppsr(i)=0.0
c
c       Get initial positions of pulsars 
c
 888    if (isotropic) then
           bdeg(i)=asin(psrran(seed))*180.0/pi
           if (psrran(seed).lt.0.5) bdeg(i)=bdeg(i)*(-1.0)
           ldeg(i)=360.0*psrran(seed)
           call calc_xyz(ldeg(i),bdeg(i),1.0,x(i),y(i),z(i))
        else if (fixed) then
           fn=fn+1
           if (fn.gt.nf) fn=1
           ldeg(i)=fl(fn)
           bdeg(i)=fb(fn)
           dtrue(i)=fd(fn)
           call calc_xyz(ldeg(i),bdeg(i),dtrue(i),x(i),y(i),z(i))
        else if (slab.or.disk) then
           x(i)=-15.0+psrran(seed)*30.0
           y(i)=-15.0+psrran(seed)*30.0
           z(i)=-5.0+psrran(seed)*10.0
           if (disk) z(i)=0.0
           call galactic(x(i),y(i),z(i),ldeg(i),bdeg(i))
        else
           if (rscale.gt.0.0) then
              r0 = rscale*sqrt(-1.0*log(1.0-psrran(seed))) 
           else if (yk04) then
              r0 = ykr(seed)
           else if (llf06) then
              r0 = llfr(seed)
           else if (uniformr) then
              x(i) = -3.0 + 6.0*psrran(seed)
              y(i) =  5.5 + 6.0*psrran(seed)
              goto 63
           else
              r0 = draw1d(rmin,rmax,rprob,nrbins,seed)
           endif
           if (spiral) then
              call spiralize(r0,seed,x(i),y(i))
           else
              theta0= 2.0*pi*psrran(seed)
              x(i) = r0 * cos(theta0)
              y(i) = r0 * sin(theta0)
           endif
 63        if (zscale.eq.0.0) then
              z(i) = draw1d(zmin,zmax,zprob,nzbins,seed)
           else
              z(i) = psr_exp(seed,0.0,zscale)
           endif
c           write(*,*) x(i),y(i),z(i)
c fudge for galactic centre
c          x(i) = psr_normal(seed,0.0,0.07)
c          y(i) = psr_normal(seed,0.0,0.07)
c          z(i) = psr_normal(seed,0.0,0.07)
        endif
c
c       Calculate true/derived distances from sun and Galactic coordinates
c
        dtrue(i) = sqrt(x(i)**2.0+(y(i)-rsun)**2.0+z(i)**2.0)
        if (dmax.gt.0.and.dtrue(i).gt.dmax) goto 888
        dderi(i) = psr_normal(seed,dtrue(i),derr*dtrue(i))
        call galactic(x(i),y(i),z(i),ldeg(i),bdeg(i))
        if (ldeg(i).gt.180.0) ldeg(i)=ldeg(i)-360.0
c
c       fudge for galactic centre
c
c	if (abs(ldeg(i)).gt.0.3.or.abs(bdeg(i)).gt.0.3) goto 888
c
c       calculate DM from derived distance
c
        dmp(i)=dm(dderi(i),ldeg(i),bdeg(i),dmod,smtau)
        smtau=0.0
        if (scattering) then
           if (smtau.eq.0.0) then
c
c          model scattering time in ms uses best-fit relation 
c          from Bhat et al. 2004 with a scatter about the scatter
c
           logtau=-6.46+0.154*log10(dmp(i))+
     &          1.07*(log10(dmp(i)))**2.0+scindx*log10(fghz)
           tau(i)=10.0**psr_normal(seed,logtau,0.8)
           else
c
c          use NE2001 scattering time
c
           tau(i)=tauiss(dtrue(i),smtau,fghz)
           endif
        else
           tau(i)=0.0 ! -noscatter option has been selected
        endif
c
c fudge for gc pulsars
c
c        tau(i)=3e5*(fghz**-4.4) 
c
c       Sky Temperature
c
	tsky(i)=psr_tsky(ldeg(i),bdeg(i),freq)
c
c       Pulsar luminosity drawn from a power-law luminosity function...
c
        if (lumfile) then
           lpsr(i) = 10.0**(draw1d(lmin,lmax,lprob,nlbins,seed))
       else if (llnorm) then
           lpsr(i) = 10.0**psr_normal(seed,lmean,lsigma)
        else if (lmin.eq.lmax) then
           lpsr(i) = lmin
        else
c           n=0.1
c           do while (n.lt.1.0) 
c              n=psrran(seed)*ntot
c              lpsr(i) = lmax*10.0**(log10(n)/lslope)
c           enddo
           lpsr(i)=lfunc(lmin,lslope,lmax,seed)
        endif
        if (isotropic.or.slab) lpsr(i)=1.0e4
c        if (disk) lpsr(i)=1.0
        s(i)=lpsr(i)/dtrue(i)/dtrue(i)
        spindex(i) = psr_normal(seed,simean,sistdv)
        if (nsurv.gt.0) then
           detected(i)=.false.
           do j=1,nsurv
              wms=width(i,j) ! do not pass down width(i) as these get modified
              flux=s(i)*(nu(j)/freq)**spindex(i) ! scale flux to survey freq
c             IS PULSAR BEAMING? IF NOT, SET S/N = 0
              if (.not.beamresults(i,j)) then
                sn(i,j) = 0.
              else
                sn(i,j)=surveysn(seed,ppsr(i),wms,dmp(i),tau(i),
     &             flux,ldeg(i),bdeg(i),tsky(i),freq,j, scindx)
              endif
c             IS S/N OVER THRESHOLD?
              if (sn(i,j).ge.snmin(j)) then
                 if (.not.detected(i)) ndet=ndet+1
                 write(lout,*) uprow,
     &           nint(100*real(ndet)/real(n2det)),
     &           '% of population (i.e.',
     &           i,' pulsars) generated'
                 detected(i)=.true.
                 if (ndet.eq.n2det) goto 555
              endif
           enddo
        else
           detected(i)=.true.
        endif
      enddo
c
c     dump output file as normal
c
 555  continue
      if (nsurv.gt.0) then
         write(lout,*) uprow,"Detected",ndet," out of",i," pulsars...",
     &                  "                       "
c         write(lout,*) "A further",unbeamed," pulsars were beaming away"
c         write(lout,*) "Mean beaming fraction:",
c     &               real(ndet)/real(i+unbeamed)
         ngen=i
      else
         write(lout,*) "Generated",ngen," pulsars..."
      endif

c   Before writing out the file, randomly
c   pick gpsFrac % pulsars to behave
c   like GPS sources SDBATES 2011

c     INTIALISE GPSARRAY TO ALL ZEROS
      do i=1, ngen
        gpsArray(i) = 0
      enddo
        
      if (gpsSwitch.gt.0) then
          gpsCount = 0
          write(*,*) "Generating GPS pulsars"
          do while (real(gpsCount)/real(ngen) .lt. gpsFrac)
            i = int(psrran(seed) * (ngen)) + 1
            gpsArray(i) = 1
            gpsCount = 0
c       SUM THE GPSARRAY TO SEE IF WE REACHED GPSFRAC
            do j=1, ngen
              gpsCount = gpsCount + gpsArray(j)
            enddo 
c           write(*,*) real(gpsCount), real(ngen), gpsFrac
          enddo
          write(*,*) gpsCount, " GPS sources produced"
      endif

c     TO CHECK HOW RANDOMLY DISTRIBUTED THE GPS SOURCES ARE
c      open(97, file="randomnums.txt", access="append", status="new")
c      do i=1,ngen
c        write(97,*) gpsArray(i)
c      enddo
c      close(97)

      do j=1, nsurv
        write(lout,'('' Writing binary population file: '',a)') 
     &       outfile(1:index(outfile,' ')-1)
        call glun(lun)
        open(lun,file=outfiles(j),form="unformatted",status="unknown")
        cline='populate '//cline
        write(lun) lcline+9,cline
        write(lun) ngen,freq
        do i=1,ngen
	     lpsr(i)=s(i)*dderi(i)*dderi(i)
         write(lun) ldeg(i),bdeg(i),ppsr(i),width(i,minfreqint),dmp(i),
     &              tau(i),
     &              tsky(i),dtrue(i),dderi(i),lpsr(i),s(i),spindex(i),
     &              gpsArray(i), gpsA, gpsB
        enddo
        close(lun)
      enddo
      open(lun,file='populate.info',status='unknown')
      write(lun,*) ngen,ndet
      close(lun)
      do j=1, nsurv
        write(*,*) nu(j)
      enddo
      end
c==============================================================================
