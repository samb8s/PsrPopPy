      program namedpsrs
      implicit none
      integer maxpsrs,lun,i,seed,ngen,dmod,narg,istat
      parameter(maxpsrs=10000)
      real ldeg(maxpsrs),bdeg(maxpsrs),p0(maxpsrs),width(maxpsrs),
     &     dmp(maxpsrs),tau(maxpsrs),tsky(maxpsrs),dtrue(maxpsrs),
     &     lpsr(maxpsrs),s(maxpsrs),spindex(maxpsrs),freq,fghz,
     &     psr_normal,psr_tsky,deg2rad,sm,smtau,smtheta,smiso,
     &     dist(maxpsrs),x(maxpsrs),y(maxpsrs),z(maxpsrs),pb,p1(maxpsrs)
      integer pksmb,swmb,phmb,swhi
      integer scod(maxpsrs)
      character*1 limit,cline*240,psr*12
      pksmb=1
      swmb=2
      phmb=4
      swhi=8
      call getarg(1,cline)
      narg=iargc()
      if (cline.eq.''.or.narg.lt.2) stop 'usage: namedpsrs dmodel file'
      read(cline,*) dmod
      call glun(lun)
      deg2rad=3.1415927/180.0
      open(lun,file='psrcat.in',status='new')
      close(lun)
      call getarg(2,cline)
      open(lun,file=cline,status='old')
      istat=0
      do while (istat.eq.0)
         read(lun,*,iostat=istat) psr
         seed=system('psrcat -c 
     &      "gl gb p0 p1 w50 dm dist_dm s1400 pb survey name" -nonumber
     &       -o short -nohead '//psr//' |sed s/"*"/0/g>>psrcat.in')
      enddo
      close(lun)

      seed=system('sort -o psrcat.in -u psrcat.in')
      freq=1374.0
      fghz=freq/1000.0
      open(lun,file="psrcat.in",status="unknown")
	i=1
      do while (.true.)
         read(lun,'(a)',end=999) cline
         read(cline,*) ldeg(i),bdeg(i),p0(i),p1(i),width(i),dmp(i),
     &                dtrue(i),s(i),pb
           if (dmp(i).gt.0.0) then
              call dmdist(dmp(i),ldeg(i),bdeg(i),dmod,dist(i))
           else
              dist(i)=0.0
           endif
           if (dist(i).gt.0.0) then
              call calc_xyz(ldeg(i),bdeg(i),dist(i),x(i),y(i),z(i))
              p0(i)=p0(i)*1000.0
              tsky(i)=psr_tsky(ldeg(i),bdeg(i),freq)
              tau(i)=10.0**(-6.46+0.154*log10(dmp(i))+
     &             1.07*(log10(dmp(i)))**2.0-3.86*log10(fghz)) ! tau in ms
              spindex(i) = psr_normal(seed,-1.6,0.35)
	      scod(i)=0
	      if (pksmb.gt.0.and.index(cline,'pksmb').gt.0) 
     &        scod(i)=scod(i)+pksmb
	      if (swmb.gt.0.and.index(cline,'swmb').gt.0)  
     &        scod(i)=scod(i)+swmb
	      if (phmb.gt.0.and.index(cline,'phmb').gt.0)  
     &        scod(i)=scod(i)+phmb
	      if (swhi.gt.0.and.index(cline,'swhi').gt.0)  
     &        scod(i)=scod(i)+swmb
              scod(i)=1
              i=i+1
        else
c           write(55,*) cline
	endif
      enddo
 999  ngen=i-1
      close(lun)
      write(*,*) ngen," pulsars from PSRCAT -> pksmb.real"
      open(lun,file="pksmb.real",status="unknown")
      write(lun,'(''# pulsars'')')
      write(lun,'(''# P (ms) * 0 1000'')')
      write(lun,'(''# DM (cm\\\\u-3\\\\d pc) * 0 1000'')')
      write(lun,'(''# W (ms) * 0 100'')')
      write(lun,'(''# l (deg) * -180 180'')')
      write(lun,'(''# b (deg) * -90 90'')')
      write(lun,'(''# S (\\\\gmJy) * 0 1e6'')')
      write(lun,'(''# SI * -3 0'')')
      write(lun,'(''# S/N * 0 100'')')
      write(lun,'(''# X (kpc) * -15 15'')')
      write(lun,'(''# Y (kpc) * -15 15'')')
      write(lun,'(''# Z (kpc) * -1 1'')')
      write(lun,'(''# D (kpc) * 0 15'')')
      write(lun,'(''# R (kpc) * 0 15'')')
      write(lun,'(''# L (mJy kpc\\\\u2\\\\d) * 0 1000'')')
      write(lun,'(''# code'')')
      do i=2,narg
         call getarg(i,cline)
         write(lun,'(''# SURVEY: '',a)') cline(1:index(cline,' ')-1)
      enddo
 1    format(3(f9.2,1x),2(f6.1,1x),f11.1,1x,
     &       f5.2,1x,f7.1,1x,5(f7.2,1x),f10.2,1x,i3)
      do i=1,ngen
         if (ldeg(i).gt.180.0) ldeg(i)=ldeg(i)-360.0
         write(lun,1) p0(i),dmp(i),width(i)*1.06,ldeg(i),bdeg(i),
     &        1000*s(i),0,0,x(i),y(i),z(i),dist(i),
     &        sqrt(x(i)*x(i)+y(i)*y(i)),s(i)*dist(i)*dist(i),scod(i)
      enddo
      close(lun)
      open(lun,file="ppdot.dat",status="unknown")
      do i=1,ngen
         if (p1(i).gt.0) write(lun,*) p0(i),p1(i),s(i)*dist(i)*dist(i)
      enddo
      close(lun)
      seed=system('rm -f psrcat.in psrcat.number')
      end
