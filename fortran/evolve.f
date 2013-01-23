cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      program evolve

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c Evolve creates a population of observable pulsars by following various
c spatial and period evolution models.

c DECLARATIONS AND INITIAL VALUES

      implicit none
      include 'survey.inc'
      real psrran,ykr,psr_exp,psr_tsky,psr_normal,dm,tauiss,surveysn
      real draw1d
      external poft
      integer maxpsrs
      parameter (maxpsrs=8000000)

c Variable Declarations

      real*4 x(maxpsrs),y(maxpsrs),z(maxpsrs),d(maxpsrs)
      real*4 ldeg(maxpsrs), bdeg(maxpsrs)  
      real*4 Vx(maxpsrs), Vy(maxpsrs), Vz(maxpsrs)
      real*4 pnot(maxpsrs), bfield(maxpsrs), alpha(maxpsrs)
      real*4 lcorr(maxpsrs), lum(maxpsrs), s(maxpsrs), sn(maxpsrs)
      real*4 simean, sistdv, spindex(maxpsrs), pwidth(maxpsrs), wms
      real*4 freq, fghz, tsky(maxpsrs), dmp(maxpsrs), beta1(maxpsrs)
      real*4 tau(maxpsrs),smtau, r0, alpha0(maxpsrs), age(maxpsrs)
      real*4 rho1(maxpsrs),r
      doubleprecision pt(maxpsrs), pdot(maxpsrs)
      logical bound, detected(maxpsrs), psrgood
      integer ngen,dmod,numnotbeaming,numdead,ndet,n2det,nint
      integer ngood
      integer npulses(maxpsrs)
      real fprand,n,braking,chi,numgen,coschi,randomx,randomy
      real sinchi,td, sinchinot, birthrate

c Temp Variable Declarations

      real tempage,temppnot,tempbfield,tempfp,temppdeath,rho
      real result,tempconstant,a,m,min,minp,tempmin, v, tempbeta
      doubleprecision temppt, fpavg, templum, tempnum
      doubleprecision temppdot, temppdotnot, lamda, rholamda
      integer temppulse
 
c Constant Declarations

      integer lun, i, j, lcline, narg, lcarg, lpth, nsurv
      integer initial,final,time1,time2,hours,minutes,seconds
      character*80 path,cline*240,carg,temp
      character*4 spindown, beaming
      character*8 vbirth
      character*6 align
      character*52 outfile
      character*3 uprow
      logical deathline, spiralarms, filexists
      real seed, pi, tmax, c, rsun, psrrad
      real vmu, vsigma, bmu, bsigma, zscale
      real pmu, psigma, pdeathline,ep,epdot
      real lnot, lcorrmu, lcorrsigma,lumprob
      doubleprecision k, kprime, loginertia

c Constant Assignments

      call getseed(seed)
      pi = 3.141592
      c = 3e10
      loginertia = 45
      rsun = 8.5
      psrrad = 1e6
      numdead = 0
      numnotbeaming = 0
      ngood = 0
      ndet = 0
      nsurv = 0
      lcline = 0

      k = 8.0*pi**2.0*psrrad**6.0/3.0/c**3.0/(10**loginertia)

c Variable Assignments (User can change some of these)

      tmax = 1e9
      vmu = 0.0
      vsigma = 180.0
      pmu = 300.0
      psigma = 150.0
      bmu = 12.65
      bsigma = 0.55
      zscale = .050
      lnot = 0.18
      lcorrmu = 0.0
      lcorrsigma = 0.8
      simean = -1.6
      sistdv = 0.35
      dmod = 4 
      pdeathline = 0.17e12
      ep = -1.5
      epdot = 0.5 
      freq = 1374.0
      fghz = freq/1000.0
      ngen = 1000
      spindown = 'fk06'
      braking = 3.0
      n = 3.0
      align = 'orthog'
      chi = 90.0
      coschi = 0.0
      sinchi = 0.0
      td = 0.0
      vbirth = 'exp'
      deathline = .true.
      beaming = 'tm98'
      spiralarms = .true.
      psrgood = .false.
      narg = iargc()
      uprow = char(27)//'[A'
      fpavg = 0.0

c GET PARAMETERS FROM USER
      
      time1 = time()

      call glun(lun)
      open(lun,file="evolve.cmd",status="unknown",access="append")
      do j=1,narg
         call getarg(j,carg)
         lcarg=index(carg,' ')
         cline = cline(1:lcline) // carg(1:lcarg)
         lcline = lcline+lcarg
      enddo
      write(lun,'("evolve ",a)') cline(1:lcline)
      close(lun)

      do j=1,maxsurv
         surveyfile(j) = ''
      enddo

c Show various options and switches for the program

      call getarg(1,carg)
      if (narg.lt.1.or.carg.eq.'-help') then
         write(*,*)
         write(*,*) "EVOLVE: A program to model pulsar populations."
         write(*,*)
         write(*,'('' usage: evolve {options} '')')
         write(*,*)
         write(*,*) "-npsrs ngen     - number of psrs to generate",
     &        " or detect (def = 1000)"
         write(*,*) "-pnot pmu psigma - initial period distribution",
     &        " in ms (def = 300.0, 150.0)"
         write(*,*) "-spindown model - fk06 or cs06",
     &        " (def = fk06)"
         write(*,*) "-braking n      - braking index (def = 3.0, 0 for",
     &        " random)"
         write(*,*) "-align chi      - orthog, random, rand45, rand07,",
     &        " rand08,rand09,rand10 (def = orthog)"
         write(*,*) "-bmu bmu bsigma - log of average magnetic field",
     &        " in Gauss (def = 12.65, 0.60)"
         write(*,*) "-vbirth dist    - velocity distribution function",
     &         " as normal or exponential (def = exp)"
         write(*,*) "-vmu vmu vsigma - average velocity in km/s (def ",
     &        "= 0.0, 180.0)"
         write(*,*) "-lum alpha beta - luminosity power law",
     &        " (def -> alpha = -1.5, beta = 0.5; 0 0 for random)"
         write(*,*) "-nodeathline    - turn off the deathline"
         write(*,*) "-beaming model  - tm98, bi90, lm88, nv83, cons,",
     &        " wj08, or none (def = tm98)" 
         write(*,*) "-nospiralarms   - turn off spiralarms"
         write(*,*)
         stop
      endif

c Assign values to certain variables depending on the user's inputs

      j=1
      write (*,*)
      do while (j.le.narg)
         call getarg(j,carg)
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
            write (*,*) "Using survey: ", surveyfile(nsurv)
         elseif (carg.eq.'-npsrs') then
            j=j+1
            call getarg(j,carg)
            read (carg,*) ngen
         elseif (carg.eq.'-spindown') then
            j=j+1
            call getarg(j,carg)
            read (carg,*) spindown
         elseif (carg.eq.'-braking') then
            j=j+1
            call getarg(j,carg)
            read (carg,*) braking
         elseif (carg.eq.'-align') then
            j=j+1
            call getarg(j,carg)
            read (carg,*) align
            if(align.eq.'rand07') td = 1e7
            if(align.eq.'rand08') td = 1e8
            if(align.eq.'rand09') td = 1e9
            if(align.eq.'rand10') td = 1e10
         elseif (carg.eq.'-vbirth') then
            j=j+1
            call getarg(j,carg)
            read (carg,*) vbirth
         elseif (carg.eq.'-lum') then
            j=j+1
            call getarg(j,carg)
            read (carg,*) ep
            j=j+1
            call getarg(j,carg)
            read (carg,*) epdot
         elseif (carg.eq.'-pnot') then
            j=j+1
            call getarg(j,carg)
            read (carg,*) pmu
            j=j+1
            call getarg(j,carg)
            read (carg,*) psigma
         elseif (carg.eq.'-bmu') then
            j=j+1
            call getarg(j,carg)
            read (carg,*) bmu
            j=j+1
            call getarg(j,carg)
            read (carg,*) bsigma
         elseif (carg.eq.'-vmu') then
            j=j+1
            call getarg(j,carg)
            read (carg,*) vmu
            j=j+1
            call getarg(j,carg)
            read (carg,*) vsigma
         elseif (carg.eq.'-nodeathline') then
            deathline = .false.
         elseif (carg.eq.'-beaming') then
            j=j+1
            call getarg(j,carg)
            read (carg,*) beaming
         elseif (carg.eq.'-nospiralarms') then
            spiralarms = .false.
         endif   
         j=j+1
      enddo

c Display the user's population model choices

      if(nsurv.gt.0) then
         n2det = ngen
         ngen = 1e9
         write (*,*) "Number of Puslars to Detect: ", n2det
      else
         write (*,*) "Number of Pulsars to Generate: ", ngen
      endif
      write (*,*) "Spindown Model: ", spindown
      write (*,*) "Initial Velocity Distribution: ", vbirth
      write (*,*) "Beaming Model: ", beaming
      write (*,*) "Braking Index: ", braking
      write (*,*) "Alignment Angle: ", align
      if(.not.deathline) write (*,*) "Not using a death line."
      if(.not.spiralarms) write (*,*) "Not using spiral arms."
      write (*,*) 

c GENERATE PULSAR ARRAY

c Here we will generate pulsar data for ngen different pulsars.

 111  i = 1
      write (*,*)
      do i = 1,ngen
         if(i.eq.58e6) GOTO 111
      if(i.eq.1000) then
         write (*,*)
         write (*,*)
      elseif(i.gt.1000.and.mod(i,1000).eq.0) then
         write (*,*) uprow, i, " pulsars generated, with ", ngood,
     &        " good ones."
      endif

c PICK INITIAL AGE (age in years)

      tempage = psrran(seed)*tmax

c PICK INITIAL PERIOD (pnot in ms)

      temppnot = -1
      do while (temppnot.le.0)
         if(spindown.eq.'fk06') then
            temppnot = psr_normal(seed, pmu, psigma)
         elseif(spindown.eq.'cs06') then
c            temppnot = psrran(seed)*190+10
            temppnot = psr_normal(seed, pmu, psigma)
         else 
            stop ' Not a valid Spindown option!!!! '
         endif
      enddo 

c PICK INITIAL MAGNETIC FIELD (bfield in Gauss)

      tempbfield = 10**(psr_normal(seed, bmu, bsigma))

c PICK ALIGNMENT ANGLE (in degrees)

      if(align.eq.'orthog') then
         chi = 90.0
         sinchinot = 1.0
         sinchi = 1.0
      elseif(align.eq.'random') then
         coschi = psrran(seed)
         chi = acos(coschi)*180.0/pi
         sinchinot = sin(chi*pi/180.0)
         sinchi = sinchinot
      elseif(align.eq.'rand45') then
         coschi = psrran(seed)*(2.0-sqrt(2.0))/2.0+(sqrt(2.0)/2.0)
         chi = acos(coschi)*180.0/pi
         sinchinot = sin(chi*pi/180.0)
      else
         coschi = psrran(seed)
         chi = acos(coschi)*180.0/pi
         sinchinot = sin(chi*pi/180.0)
         sinchi = sinchinot*exp(-tempage/td)
         chi = asin(sinchi)*180.0/pi
      endif

c CALCULATE BRAKING INDEX

      if(braking.eq.0) then
         n = 2.5 + psrran(seed)*0.5
      else
         n = braking
      endif

c DETERMINE IF THE PULSAR IS DEAD     

c Faucher-Giguere and Kaspi spindown model

      if(spindown.eq.'fk06') then
	 kprime = k*tempbfield**2.0
	temppt=((temppnot/1e3)**(n-1.0)+(n-1.0)*kprime*tempage*
     &	      365.25*24.0*3.6e3*(sinchinot**2.0))**(1.0/(n-1.0))		
	temppt=temppt*1000.0
	temppdot = sinchi**2.0*kprime*(temppt/1e3)**(2.0-n)
	temppdeath = 3.2*10**19.0*sqrt(temppdot*temppt/1e3)/
     &      ((temppt/1e3)**2.0)
         if(deathline) then
            if(temppdeath.lt.pdeathline) then
               numdead = numdead + 1
               psrgood = .false.
            else
               psrgood = .true.
            endif
         else 
            psrgood = .true.   
         endif

c Contopolous and Spitkovsky spindown model 

      elseif(spindown.eq.'cs06') then
         temppdeath=(0.81*tempbfield/1e9/temppnot)**(2.0/(n+1.0))*1e3
         initial = temppnot
         final = 1e5   
         a = coschi**2.0/temppdeath
         min = 1e24
         minp = temppnot
         tempconstant=3.3e-40*tempbfield**2*(temppnot)**(n-3.0)*
     &	    tempage*365.25*24*3.6e9
         call qsimp(poft,initial,temppdeath,result,a,n)
         if(result.lt.tempconstant.or.result.gt.1e14) then
            temppt = 1e6
            goto 123
         endif
         do m = initial,final+1
            call qsimp(poft,initial,real(m),result,a,n)
            if(result.gt.1e14) then
               temppt = 1e6
               goto 123
            endif
            tempmin = abs(result-tempconstant)
            if(tempmin.le.min) then
               min = tempmin
               minp = m
            else
               temppt = minp
               goto 123
            endif
            if(m.eq.(final+1)) temppt = 1e6
         enddo
 123     if(temppt.gt.temppdeath.and.deathline) then
            numdead = numdead + 1
            psrgood = .false.
         else
            psrgood = .true.
            temppdot = 3.3e-40*tempbfield**2.0*(1e3/temppnot)*
     &           ((temppt/temppnot)**(2.0-n))*(1.0-a*temppt)
         endif         
      endif

c DETERMINE IF THE PULSAR IS BEAMING TOWARDS US

      if(psrgood) then
         rho = 5.4*(temppt/1e3)**(-1.0/2.0)
         temppulse = 0
         if(beaming.eq.'bi90') then
            tempfp=(4.0/pi)*sin((6.2*(temppt/1e3)**(-1.0/2.0))*pi/180.0)
         elseif(beaming.eq.'lm88') then
            tempfp=(4.0/pi)*sin(6.5*((temppt/1e3)**(-1.0/3.0))*pi/180.0)
         elseif(beaming.eq.'tm98') then
            tempfp = 0.09*((log10(temppt/1e3)-1.0)**2.0)+0.03
         elseif(beaming.eq.'nv83') then
            tempfp = (temppt/1e3)**(-0.65)/4.9
         elseif(beaming.eq.'cons') then
            tempfp = 0.2
         elseif(beaming.eq.'none') then
            tempfp = 1
         elseif(beaming.eq.'wj08') then
            if(chi.lt.1e-5) chi = 0
	    tempbeta = psrran(seed)
	    tempbeta = acos(tempbeta)*180.0/pi
            if(tempbeta.le.rho) then
               tempfp=1		
            else
               tempfp=0
            endif
            if(chi.gt.rho.and.(chi+rho).lt.90.0) then
               tempfp=2.0*sin(chi*pi/180.0)*sin(rho*pi/180.0)
            elseif(chi.gt.rho.and.(chi+rho).gt.90.0) then
               tempfp=cos((chi-rho)*pi/180.0)
            elseif(chi.le.rho.and.(chi+rho).lt.90.0) then
               tempfp=1-cos((chi+rho)*pi/180.0)
            else
               tempfp=1
            endif
         else
            stop "Not a valid beaming model!!!! "
         endif
   
         fprand = psrran(seed)

         if(tempfp.lt.fprand) then
            numnotbeaming = numnotbeaming + 1
            psrgood = .false.
         else
            ngood = ngood + 1
            psrgood = .true.
            temppulse = temppulse + 1
         endif

         fpavg = fpavg + tempfp

      endif

c IF IT IS A GOOD PULSAR, THEN DO THE REST OF THE CALCULATIONS

      if(psrgood) then
         detected(ngood) = .false.
         bfield(ngood) = tempbfield
         npulses(ngood) = temppulse
         alpha0(ngood) = asin(sinchinot)*180.0/pi
         alpha(ngood) = chi
         age(ngood) = tempage
         if(beaming.ne.'wj08') then
            beta1(ngood) = psrran(seed)*rho
         else
            beta1(ngood) = tempbeta
         endif
         rho1(ngood) = rho

c	write(62,*) temppt/1e3, tempfp

c PULSAR ROTATIONAL EVOLUTION (pt in ms, pdot in s/s, pwidth in ms)

         pnot(ngood) = temppnot
         pt(ngood) = temppt
         pdot(ngood) = temppdot
         pwidth(ngood) = 0.05*pt(ngood)

c PICK INITIAL POSITION (x,y,z in kpc)
         
         if(spiralarms) then
            r0 = ykr(seed)
            call spiralize(r0,seed,x(ngood),y(ngood))
         else
            x(ngood) = -20.0 + psrran(seed)*40.0
            y(ngood) = -20.0 + psrran(seed)*40.0
         endif
         z(ngood) = psr_exp(seed,0.0,zscale)

c PICK INITIAL VELOCITY (Vx,Vy,Vz in km/s)

         if(vbirth.eq.'gaussian') then
            Vx(ngood) = psr_normal(seed, vmu, vsigma)
            Vy(ngood) = psr_normal(seed, vmu, vsigma)
            Vz(ngood) = psr_normal(seed, vmu, vsigma)
         elseif(vbirth.eq.'exp') then
            Vx(ngood) = psr_exp(seed, vmu, vsigma)
            Vy(ngood) = psr_exp(seed, vmu, vsigma)
            Vz(ngood) = psr_exp(seed, vmu, vsigma)

c         elseif(vbirth.eq.'other') then

c            v=draw1d(0,1000.0,vel0609,50,seed)
c            v=draw1d(0,1000.0,vel1913,50,seed)

c            randomx = psrran(seed)*2*pi
c            randomy = psrran(seed)*2 - 1
c            Vx(ngood) = v*cos(randomx)*randomy
c            Vy(ngood) = v*sin(randomx)*randomy
c            Vz(ngood) = v*cos(asin(randomy))

         else
            stop ' Not a valid velocity distribution '
         endif

c PULSAR SPATIAL EVOLUTION
         
         call vxyz(0.005,x(ngood),y(ngood),z(ngood),tempage/(1e6),
     &        Vx(ngood),Vy(ngood),Vz(ngood),
     &        x(ngood),y(ngood),z(ngood),Vx(ngood),Vy(ngood),
     &        Vz(ngood),bound)

c CALCULATE LUMINOSITY (lum in mJy kpc^2)

         if(ep.eq.0.0.and.epdot.eq.0.0) then
            randomy = 1e5
            lumprob = -1
            do while (randomy.gt.lumprob)
               randomx = psrran(seed)*1e4
               randomy = psrran(seed)*18.5
               if (randomx.lt.0.1) then
                  lumprob = -1
               elseif (randomx.gt.2) then
                  lumprob = 1.662*randomx**(-2.0)
               else
                  lumprob = randomx**(-19.0/15.0)
               endif
            enddo
            lum(ngood) = randomx
c         elseif(spindown.eq.'fk06') then
         else
            lcorr(ngood)=psr_normal(seed, lcorrmu, lcorrsigma)
            lum(ngood)=log10(lnot*(pdot(ngood)*1e15)**epdot*
     &           (pt(ngood)/1e3)**ep)+lcorr(ngood)
            lum(ngood)=10**lum(ngood)
c         elseif(spindown.eq.'cs06') then
c            lamda=psrran(seed)*20
c            rholamda=500*lamda**2*exp(-1*lamda)
c            tempnum=psrran(seed)*300
c            do while(tempnum.gt.rholamda)
c               lamda=psrran(seed)*20
c               rholamda=500*lamda**2*exp(-1*lamda)
c               tempnum=psrran(seed)*300
c            enddo
c            templum=(1.0/3.0)*log10(pdot(ngood)/(pt(ngood)/1e3)**3.0)
c     &            + 6.64
c            lum(ngood)=10**(lamda/3.6+templum-1.8)
c            lum(ngood)=lum(ngood)*((freq/400.0)**spindex(ngood))
         endif

c CALCULATE OBSERVED FLUX (s in mJy)

         call galactic(x(ngood),y(ngood),z(ngood),ldeg(ngood),
     &        bdeg(ngood))
         if(ldeg(ngood).gt.180.0) ldeg(ngood)=ldeg(ngood)-360.0
         d(ngood)=sqrt(x(ngood)**2.0+(y(ngood)-rsun)**2.0+
     &        z(ngood)**2.0)
         s(ngood) = lum(ngood)/(d(ngood)**2.0)


c    I don't understand this bit at all!! :oS
         
c	 if(d(ngood).lt.40.or.x(ngood)**2.0+(y(ngood)-rsun)**2.0.lt.50) then
c	    s(ngood) = 100000
c	    lum(ngood) = 100000
c	 else
c	    s(ngood) = 0.0
c	    lum(ngood) = 0.0
c	 endif

c CALCULATE TSKY AND SPECTRAL INDEX

         tsky(ngood) = psr_tsky(ldeg(ngood),bdeg(ngood),freq)
         spindex(ngood) = psr_normal(seed,simean,sistdv)

c CALCULATE DM AND TAU

         dmp(ngood)=dm(d(ngood),ldeg(ngood),bdeg(ngood),dmod,smtau)
         tau(ngood) = tauiss(d(ngood),smtau,fghz)

c ATTEMPT TO DETECT THE PULSARS (IF DESIRED)

         if (nsurv.gt.0) then
            detected(ngood)=.false.
            do j=1,nsurv
               wms=pwidth(ngood)  ! do not pass down pwidth(ngood) as it gets modified
               sn(ngood)=surveysn(seed,real(pt(ngood)),wms,dmp(ngood),
     &         tau(ngood),s(ngood),ldeg(ngood),bdeg(ngood),tsky(ngood),
     &         freq,j)
               if (sn(ngood).ge.snmin(j)) then
                  if (.not.detected(ngood)) ndet=ndet+1
                  if (i.gt.1000) then
                     write(*,*) nint(100*real(ndet)/real(n2det)),
     &                    '% of population (i.e.',
     &                    ndet,' pulsars) detected',uprow
                  endif
                  detected(ngood)=.true.
                  if (ndet.eq.n2det) then
                     write (*,*)
                     goto 555
                  endif
               endif
            enddo
c         else
c            detected(ngood)=.true.
	 else
	    r=sqrt(x(ngood)**2.0+(y(ngood)-rsun)**2.0)
	    write(*,*) d(ngood), r
	    if (d(ngood).lt.40.0.and.r.lt.50.0) then
	       detected(ngood)=.true.
	       s(ngood) = 999.9
               if (i.gt.1000) then
                  write(*,*) nint(100*real(ndet)/real(n2det)),
     &                 '% of population (i.e.',
     &                 ndet,' pulsars) detected',uprow
               endif
               if (.not.detected(ngood)) ndet=ndet+1
               if (ndet.eq.n2det) then
                  write (*,*)
                  goto 555
               endif
	    else
	       detected(ngood)=.false.
	       s(ngood) = 0.0
	    endif
         endif

c END THE "if psr good" LOOP AND THE "do while" LOOP

      endif 

      enddo

c SEND THE DATA TO AN OUTPUT FILE
     
 555  call glun(lun)

c      write(outfile,10) "psrpop","_",spindown,"_",beaming,"_",
c     &     nint(10.0*braking),"_",nint(-10.0*ep),nint(10.0*epdot),
c     &     "_",align,"_",nint(pmu),nint(psigma),"_",nint(100.0*bmu),
c     &     nint(100.0*bsigma),".model"

      write(outfile,10) "evolve.model"

 10   format(A12)

c 10   format(A6,A1,A4,A1,A4,A1,I2.2,A1,I2.2,I2.2,A1,A6,A1,I3.3,I3.3,A1,
c     &     I4.4,I3.3,A6)

      write(*,*) outfile
      open(lun,file=outfile,form="unformatted",status="unknown")
      cline = 'evolve '//cline
      write (lun) lcline+7, cline, 1
      numgen = numdead+numnotbeaming+ngood
      write (lun) ngood, freq, numgen
      do i=1,ngood
         write (lun) ldeg(i),bdeg(i),real(pt(i)),real(pdot(i)),
     &        pwidth(i),alpha(i),beta1(i),dmp(i),tau(i),tsky(i),d(i),
     &        d(i),lum(i),s(i),spindex(i),rho1(i),age(i)  
      enddo
      close(lun)

      write (*,*)
      write (*,*) "Total number of pulsars generated = ", numgen
      write (*,*) "Number of dead pulsars = ", numdead
      write (*,*) "Number of non-beaming pulsars = ", numnotbeaming
      if(nsurv.gt.0) then 
         write (*,*) "Number of good but not detected pulsars = ", 
     &        ngood-ndet
         write (*,*) "Number of detected pulsars = ", ndet
      else
         write (*,*) "Number of good pulsars = ", ngood
      endif
      write (*,*) "Average beaming fraction = ",
     &     fpavg/real(numnotbeaming+ngood)
      write (*,*) 

      time2 = time()
      hours = (time2-time1)/3600
      minutes = ((time2-time1)-hours*3600)/60
      seconds = (time2-time1)-hours*3600-minutes*60
      write(*,*) "Time elapsed: ",hours," hours,",minutes,
     &     " minutes, and",seconds," seconds."
      write(*,*) "Total number of seconds: ",time2-time1 
      write(*,*) "Birthrate (psrs/century) is: ", numgen*100.0/tmax

c THIS IS THE END OF THE PROGRAM

      stop
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
