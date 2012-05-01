      program realpsrs
      implicit none
      character*240 string
      integer maxpsrs,lun,i,seed,ngen,ncat
      parameter(maxpsrs=10000)
      real ldeg(maxpsrs),bdeg(maxpsrs),ppsr(maxpsrs),width(maxpsrs),
     &     dmp(maxpsrs),tau(maxpsrs),tsky(maxpsrs),dtrue(maxpsrs),
     &     lpsr(maxpsrs),s(maxpsrs),spindex(maxpsrs),freq,fghz,
     &     psr_normal,psr_tsky,dderi(maxpsrs)
      call glun(lun)
      seed=system('psrcat  -nonumber -o short -nohead | wc -l >
     &             psrcat.number')
      open(lun,file="psrcat.number")
      read(lun,*) ncat
      close(lun)
      seed=system('psrcat -c "gl gb p0 w50 dm dist r_lum14 s1400" 
     & -nonumber -o short -nohead |grep -v "*">psrcat.in')
      freq=1374.0
      fghz=freq/1000.0
      open(lun,file="psrcat.in",status="unknown")
      do i=1,maxpsrs
         read(lun,*,end=999) ldeg(i),bdeg(i),ppsr(i),width(i),dmp(i),
     &                dtrue(i),lpsr(i),s(i)
         dderi(i)=dtrue(i)
         ppsr(i)=ppsr(i)*1000.0
         width(i)=ppsr(i)*0.04
         tsky(i)=psr_tsky(ldeg(i),bdeg(i),freq)
         tau(i)=10.0**(-6.46+0.154*log10(dmp(i))+
     &        1.07*(log10(dmp(i)))**2.0-3.86*log10(fghz)) ! tau in ms
        spindex(i) = -1.5 ! psr_normal(seed,-1.6,0.35)
      enddo
 999  ngen=i-1
      close(lun)
      write(*,*) ngen," from",ncat," pulsars from PSRCAT -> psrpop.real"
      open(lun,file="psrpop.real",form="unformatted",status="unknown")
      string="the pulsar catalogue"
	i=20
      write(lun) i,string
      write(lun) ngen,freq
      do i=1,ngen
         write(lun) ldeg(i),bdeg(i),ppsr(i),width(i),dmp(i),tau(i),
     &            tsky(i),dtrue(i),dderi(i),lpsr(i),s(i),spindex(i)
      enddo
      close(lun)
      seed=system('rm -f psrcat.in psrcat.number')
      end
